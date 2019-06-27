from __future__ import print_function
import cython
from astropy.io import fits
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

######################################################################
# Class to hold the parameters of the optimization, visible to both C
# for fast access and publicly to python.
######################################################################

cdef class Params:
    cdef public double sau, esino, ecoso, inc, asc, lam, mpri, msec, jit
    cdef public double ecc, per, arg
    
    def __init__(self, par):

        cdef extern from "math.h":
            double atan2(double x, double y)
            double sqrt(double x)
        
        self.sau = par[0]
        self.esino = par[1]
        self.inc = par[2]
        self.asc = par[3]
        self.ecoso = par[4]
        self.lam = par[5]
        self.mpri = par[6]
        self.msec = par[7]
        self.jit = par[8]
        self.ecc = self.ecoso**2 + self.esino**2
        self.per = sqrt(self.sau*self.sau*self.sau/(self.mpri + self.msec))*365.25
        self.arg = atan2(self.esino, self.ecoso)

######################################################################
# A small structure to hold the important components of Mirek's HTOF
# routine to avoid python overheads.
######################################################################

cdef class AstrometricFitter:
    cdef public int nepochs, npar
    cdef double [:, :] chi2_matrix
    cdef double [:, :] ra_solution_vecs
    cdef double [:, :] dec_solution_vecs

    def __init__(self, fitter):
        self.chi2_matrix = fitter.fitter._chi2_matrix
        self.ra_solution_vecs = fitter.fitter.astrometric_solution_vector_components['ra'].T.copy(order='C')
        self.dec_solution_vecs = fitter.fitter.astrometric_solution_vector_components['dec'].T.copy(order='C')
        self.nepochs = fitter.fitter.astrometric_solution_vector_components['ra'].shape[0]
        self.npar = fitter.fitter.astrometric_solution_vector_components['ra'].shape[1]
        
######################################################################
# Structure to hold all of the data, also includes a method for
# loading in data.  Perhaps add a covariance multiplier below, for
# bright Gaia stars?
######################################################################

cdef class Data:
    cdef double [:] epochs
    cdef double [:] RV
    cdef double [:] RV_err
    cdef int [:] RVinst
    cdef double [:] relsep
    cdef double [:] PA
    cdef double [:] relsep_err
    cdef double [:] PA_err
    cdef public int nRV, nAst, nHip1, nHip2, nGaia, nTot, nInst
    cdef public double pmra_H, pmdec_H, pmra_HG, pmdec_HG, pmra_G, pmdec_G
    cdef public double plx, plx_err
    cdef double [:, :] Cinv_H
    cdef double [:, :] Cinv_HG
    cdef double [:, :] Cinv_G
    cdef public double refep
    cdef public double epRA_H, epDec_H, epRA_G, epDec_G, dt_H, dt_G
    cdef public int use_abs_ast

    def __init__(self, Hip, RVfile, relAstfile,
                 use_epoch_astrometry=False,
                 epochs_Hip1=None, epochs_Hip2=None, epochs_Gaia=None,
                 refep=None): #2455197.5000):
        try:
            rvdat = np.genfromtxt(RVfile)
            rvep = rvdat[:, 0]
            self.RV = rvdat[:, 1]
            self.RV_err = rvdat[:, 2]
            self.nRV = rvdat.shape[0]
            print("Loading RV data from file " + RVfile)
        except:
            print("Unable to load RV data from file " + RVfile)
            self.nRV = 0
            self.nInst = 1
            rvep = []
        try:
            self.RVinst = (rvdat[:, 3]).astype(np.int32)
            # Check to see that the column we loaded was an integer
            assert np.all(self.RVinst == rvdat[:, 3])
            self.nInst = np.amax(rvdat[:, 3]) + 1
            print("Loaded data from %d RV instruments." % (self.nInst))
        except:
            print("Unable to read RV instruments from fourth column.")
            print("Assuming all data are from one instrument.")
            self.RVinst = (rvdat[:, 2]*0).astype(np.int32)
            self.nInst = 1

        try:
            reldat = np.genfromtxt(relAstfile,
                                   usecols=(1,2,3,4,5), skip_header=1)
            relep = (reldat[:, 0] - 2000)*365.25 + 2451544.5
            self.relsep = reldat[:, 1]
            self.relsep_err = reldat[:, 2]
            self.PA = reldat[:, 3]*np.pi/180
            self.PA_err = reldat[:, 4]*np.pi/180
            self.nAst = reldat.shape[0]
            print("Loaded %d relative astrometric data points from file " % (self.nAst) + relAstfile)
        except:
            print("Unable to load relative astrometry data from file " + relAstfile)
            self.nAst = 0
            relep = []

        try:
            t = fits.open('HGCA_vDR2_corrected.fits')[1].data
            t = t[np.where(t['hip_id'] == Hip)]
            assert len(t) > 0
            print("Loading absolute astrometry data for Hip %d" % (Hip))
            self.use_abs_ast = 1
        except:
            print("Unable to load absolute astrometry data for Hip %d" % (Hip))
            self.use_abs_ast = 0
            self.epochs = np.asarray(list(rvep) + list(relep))
            self.nTot = len(self.epochs)

            #########################################################
            # Note: we still need parallax to run this calculation
            # without the full HGCA.  With this exception, the lines
            # below should be sufficient to make the fit run without
            # the HGCA.  The main program can access and modify the
            # parallax and parallax error.
            #########################################################
            
            self.pmra_H = self.pmra_G = self.pmra_HG = 0
            self.pmdec_H = self.pmdec_G = self.pmdec_HG = 0
            self.epRA_H = self.epDec_H = 1991.25
            self.epRA_G = self.epDec_G = 2015.5
            if refep is not None:
                self.refep = refep
            else:
                self.refep = self.epochs[0]
            self.Cinv_H = np.zeros((2, 2))
            self.Cinv_G = np.zeros((2, 2))
            self.Cinv_HG = np.identity(2)
            return
        
        self.plx = 1e-3*t['parallax_gaia']
        self.plx_err = 1e-3*t['parallax_gaia_error']
        self.pmra_H = 1e-3*t['pmra_hip']
        self.pmdec_H = 1e-3*t['pmdec_hip']
        self.pmra_HG = 1e-3*t['pmra_hg']
        self.pmdec_HG = 1e-3*t['pmdec_hg']
        self.pmra_G = 1e-3*t['pmra_gaia']
        self.pmdec_G = 1e-3*t['pmdec_gaia']
        self.epRA_H = t['epoch_ra_hip']
        self.epDec_H = t['epoch_dec_hip']
        self.epRA_G = t['epoch_ra_gaia']
        self.epDec_G = t['epoch_dec_gaia']
        if refep is not None:
            self.refep = refep
        else:
            self.refep = self.epochs[0]
        
        if not use_epoch_astrometry:
            self.nHip1 = self.nHip2 = self.nGaia = 6
            self.dt_H = 3.36
            self.dt_G = 1.83
            ep_2010 = 2455197.5000
            
            dmurH_epc = (self.epRA_H - 2010.0)*365.25 + ep_2010
            dmudH_epc = (self.epDec_H - 2010.0)*365.25 + ep_2010
            dmurG_epc = (self.epRA_G - 2010.0)*365.25 + ep_2010
            dmudG_epc = (self.epDec_G - 2010.0)*365.25 + ep_2010

            absasteps = np.asarray([dmurH_epc - self.dt_H*0.5*365.25, dmurH_epc,
                                    dmurH_epc + self.dt_H*0.5*365.25,
                                    dmudH_epc - self.dt_H*0.5*365.25, dmudH_epc,
                                    dmudH_epc + self.dt_H*0.5*365.25,
                                    dmurH_epc - self.dt_H*0.5*365.25, dmurH_epc,
                                    dmurH_epc + self.dt_H*0.5*365.25,
                                    dmudH_epc - self.dt_H*0.5*365.25, dmudH_epc,
                                    dmudH_epc + self.dt_H*0.5*365.25,
                                    dmurG_epc - self.dt_G*0.5*365.25, dmurG_epc,
                                    dmurG_epc + self.dt_G*0.5*365.25,
                                    dmudG_epc - self.dt_G*0.5*365.25, dmudG_epc,
                                    dmudG_epc + self.dt_G*0.5*365.25])
        else:
            self.nHip1 = epochs_Hip1.shape[0]
            self.nHip2 = epochs_Hip2.shape[0]
            self.nGaia = epochs_Gaia.shape[0]
            absasteps = np.asarray(list(epochs_Hip1) + list(epochs_Hip2) + list(epochs_Gaia))
            
        self.epochs = np.asarray(list(rvep) + list(relep) + list(absasteps))
        self.nTot = len(self.epochs)
        
        eRA, eDec, corr = [1e-3*t['pmra_hip_error'], 1e-3*t['pmdec_hip_error'], t['pmra_pmdec_hip']]
        C_H = np.asarray([[eRA**2, eRA*eDec*corr], [eRA*eDec*corr, eDec**2]])
        eRA, eDec, corr = [1e-3*t['pmra_hg_error'], 1e-3*t['pmdec_hg_error'], t['pmra_pmdec_hg']]
        C_HG = np.asarray([[eRA**2, eRA*eDec*corr], [eRA*eDec*corr, eDec**2]])
        eRA, eDec, corr = [1e-3*t['pmra_gaia_error'], 1e-3*t['pmdec_gaia_error'], t['pmra_pmdec_gaia']]
        C_G = np.asarray([[eRA**2, eRA*eDec*corr], [eRA*eDec*corr, eDec**2]])
        
        self.Cinv_H = np.linalg.inv(C_H.reshape(2, 2)).astype(float)
        self.Cinv_HG = np.linalg.inv(C_HG.reshape(2, 2)).astype(float)
        self.Cinv_G = np.linalg.inv(C_G.reshape(2, 2)).astype(float)
        

cdef class Model:

    cdef public int nEA, nRV, nAst, nHip1, nHip2, nGaia
    cdef public double pmra_H, pmra_HG, pmra_G, pmdec_H, pmdec_HG, pmdec_G
    cdef double *EA
    cdef double *sinEA
    cdef double *cosEA
    cdef double *RV
    cdef double *relsep
    cdef double *PA
    cdef double *dRA_H1
    cdef double *dDec_H1
    cdef double *dRA_H2
    cdef double *dDec_H2
    cdef double *dRA_G
    cdef double *dDec_G
    
    def __init__(self, Data data):
        self.nEA = data.nTot
        self.nRV = data.nRV
        self.nAst = data.nAst
        self.nHip1 = data.nHip1
        self.nHip2 = data.nHip2
        self.nGaia = data.nGaia
        
        self.EA = <double *> PyMem_Malloc(self.nEA * sizeof(double))
        self.sinEA = <double *> PyMem_Malloc(self.nEA * sizeof(double))
        self.cosEA = <double *> PyMem_Malloc(self.nEA * sizeof(double))

        self.RV = <double *> PyMem_Malloc(self.nRV * sizeof(double))
        
        self.relsep = <double *> PyMem_Malloc(self.nAst * sizeof(double))
        self.PA = <double *> PyMem_Malloc(self.nAst * sizeof(double))

        self.dRA_H1 = <double *> PyMem_Malloc(self.nHip1 * sizeof(double))
        self.dDec_H1 = <double *> PyMem_Malloc(self.nHip1 * sizeof(double))
        self.dRA_H2 = <double *> PyMem_Malloc(self.nHip2 * sizeof(double))
        self.dDec_H2 = <double *> PyMem_Malloc(self.nHip2 * sizeof(double))
        self.dRA_G = <double *> PyMem_Malloc(self.nGaia * sizeof(double))      
        self.dDec_G = <double *> PyMem_Malloc(self.nGaia * sizeof(double))

        
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

def calc_EA_RPP(Data data, Params par, Model model):

    cdef extern from "calcEA.h":
        double shortsin(double x)
        double MAmod(double x)
        double EAstart(double M, double ecc)
        void getbounds(double x[], double y[], double e)
    cdef extern from "math.h":
        double fabs(double x)
        double sqrt(double x)

    cdef double *EA_tab = <double *> PyMem_Malloc(6*13 * sizeof(double))
    cdef double *bounds = <double *> PyMem_Malloc(13 * sizeof(double))

    getbounds(bounds, EA_tab, par.ecc)
    
    cdef double pi = 3.14159265358979323846264338327950288
    cdef double pi_d_4 = 0.25*pi
    cdef double pi_d_2 = 0.5*pi
    cdef double threepi_d_4 = 0.75*pi
    cdef double twopi = 2*pi
    cdef int i, j, k, n

    cdef double one_over_ecc = 1/max(1e-17, par.ecc)

    n = data.nTot

    #################################################################
    # These are for computing the mean anomaly.  Cheaper to do this
    # within the calculation of eccentric anomaly and not write to its
    # own array, since we never need it again.
    #################################################################
    
    cdef double _MA, _EA, sinEA, cosEA, dEA, dEAsq_d6
    cdef int EAsign
    cdef double num, denom
    cdef double twopi_d_per = twopi/par.per
    cdef double MA0 = par.lam - par.arg
    cdef double dMA, MA_last
    
    cdef double one_sixth = 1./6

    if par.ecc < 0.78:
        for i in range(n):

            #############################################################
            # Ensure _MA is between 0 and pi. This enables the use of
            # shorter Taylor series for sine, simplifies conversion
            # between sine and cosine, and facilities the use of a lookup
            # table. Note that pi:2pi maps trivially onto 0:pi (boolean
            # wrapped).
            #############################################################

            if i > 0:
                _MA = MA_last + twopi_d_per*(data.epochs[i] - data.epochs[i-1])
                _MA = MAmod(_MA)
            else: 
                _MA = MAmod(twopi_d_per*(data.epochs[i] - data.refep) + MA0)
            MA_last = _MA
        
            if _MA > pi:
                EAsign = -1
                _MA = twopi - _MA
            else:
                EAsign = 1

            #############################################################
            # Use the lookup table for the initial guess.
            #############################################################

            for j in range(11, -1, -1):
                if _MA > bounds[j]:
                    break
            
            k = 6*j
            dx = _MA - bounds[j]
            _EA = EA_tab[k] + dx*(EA_tab[k + 1] + dx*(EA_tab[k + 2] + dx*(EA_tab[k + 3] + dx*(EA_tab[k + 4] + dx*EA_tab[k + 5]))))
            
            #############################################################
            # For sinEA, since _EA in [0,pi], sinEA should always be >=0
            # (no sign ambiguity).  sqrt is much cheaper than sin.  If
            # |cos|>|sin|, compute them in reverse order, again using sqrt
            # to avoid a trig call.  Also, use trig identities, sin with a
            # low argument and the series expansion to minimize
            # computational cost.
            #############################################################

            if not _EA > pi_d_4:
                sinEA = shortsin(_EA)
                cosEA = sqrt(1 - sinEA**2)
            elif _EA < threepi_d_4:
                cosEA = shortsin(pi_d_2 - _EA)
                sinEA = sqrt(1 - cosEA**2)
            else:
                sinEA = shortsin(pi - _EA)
                cosEA = -sqrt(1 - sinEA**2)

            num = (_MA - _EA)*one_over_ecc + sinEA
            denom = one_over_ecc - cosEA
            
            # Second order approximation.
            
            dEA = num*denom/(denom**2 + 0.5*sinEA*num)
        
            #############################################################
            # Apply our correction to EA, sinEA, and cosEA using
            # series.  Go to second order, since that was our level of
            # approximation above and will get us to basically machine
            # precision for eccentricities below 0.78.
            #############################################################

            model.EA[i] = EAsign*(_EA + dEA)
            model.sinEA[i] = EAsign*(sinEA*(1 - 0.5*dEA*dEA) + dEA*cosEA)
            model.cosEA[i] = cosEA*(1 - 0.5*dEA*dEA) - dEA*sinEA

    #####################################################################
    # Higher eccentricities will require a third-order correction to
    # achieve machine precision for all values of the eccentric
    # anomaly.  In the singular corner, they also use a series
    # expansion rather than the piecewise polynomial fit.
    #####################################################################
            
    else:
        for i in range(n):
            if i > 0:
                _MA = MA_last + twopi_d_per*(data.epochs[i] - data.epochs[i-1])
                _MA = MAmod(_MA)
            else: 
                _MA = MAmod(twopi_d_per*(data.epochs[i] - data.refep) + MA0)
            MA_last = _MA
               
            if _MA > pi:
                EAsign = -1
                _MA = twopi - _MA
            else:
                EAsign = 1

            #############################################################
            # Use the lookup table for the initial guess as long as we
            # are not in the singular corner.
            #############################################################

            if _MA + (1 - par.ecc) > 0.25:
                for j in range(11, -1, -1):
                    if _MA > bounds[j]:
                        break

                k = 6*j
                dx = _MA - bounds[j]
                _EA = EA_tab[k] + dx*(EA_tab[k + 1] + dx*(EA_tab[k + 2] + dx*(EA_tab[k + 3] + dx*(EA_tab[k + 4] + dx*EA_tab[k + 5]))))
                
            #############################################################
            # Use the series expansions in the singular corner.
            #############################################################

            else:
                _EA = EAstart(_MA, par.ecc)

            if not _EA > pi_d_4:
                sinEA = shortsin(_EA)
                cosEA = sqrt(1 - sinEA**2)
            elif _EA < threepi_d_4:
                cosEA = shortsin(pi_d_2 - _EA)
                sinEA = sqrt(1 - cosEA**2)
            else:
                sinEA = shortsin(pi - _EA)
                cosEA = -sqrt(1 - sinEA**2)

            num = (_MA - _EA)*one_over_ecc + sinEA
            denom = one_over_ecc - cosEA
            
            # Second order approximation.
            
            dEA = num*denom/(denom**2 + 0.5*sinEA*num)

            # This brings the scheme to third order.
        
            dEA = num/(denom + dEA*(0.5*sinEA + one_sixth*cosEA*dEA))
                    
            dEAsq_d6 = dEA*dEA*one_sixth
        
            #############################################################
            # Apply our correction to EA, sinEA, and cosEA using
            # series.  Go to third order, since that was our level of
            # approximation above and will get us to basically machine
            # precision at the higher eccentricities.
            #############################################################

            model.EA[i] = EAsign*(_EA + dEA)
            model.sinEA[i] = EAsign*(sinEA*(1 - 3*dEAsq_d6) + dEA*(1 - dEAsq_d6)*cosEA)
            model.cosEA[i] = cosEA*(1 - 3*dEAsq_d6) - dEA*(1 - dEAsq_d6)*sinEA
        
    PyMem_Free(bounds)
    PyMem_Free(EA_tab)
    return

    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

######################################################################
# Use the Thieles-Innes equations to compute the relative positions of
# the two objects (to be scaled by the parallax).  This will modify
# the relative separations, position angles, and RA and Dec offsets
# (for absolute astrometry) directly in the model structure.
######################################################################

def calc_offsets(Data data, Params par, Model model):

    cdef extern from "math.h" nogil:
        double sin(double _x)
        double cos(double _x)
        double sqrt(double _x)
        double atan2(double _y, double _x)

    cdef double a_1 = -par.sau/(1. + par.mpri/par.msec)
    cdef double cosarg = cos(par.arg)
    cdef double sinarg = sin(par.arg)
    cdef double cosasc = cos(par.asc)
    cdef double sinasc = sin(par.asc)    
    cdef double cosinc = cos(par.inc)

    cdef double A = a_1*(cosarg*cosasc - sinarg*sinasc*cosinc)
    cdef double B = a_1*(cosarg*sinasc + sinarg*cosasc*cosinc)
    cdef double F = a_1*(-sinarg*cosasc - cosarg*sinasc*cosinc)
    cdef double G = a_1*(-sinarg*sinasc + cosarg*cosasc*cosinc)
    
    cdef int n = data.nTot - data.nRV
    cdef int i, i1, i2
    cdef double X, Y, dRA, dDec, sqonemeccsqr = sqrt(1 - par.ecc**2)

    for i in range(data.nAst):
        X = model.cosEA[i + data.nRV] - par.ecc
        Y = model.sinEA[i + data.nRV]*sqonemeccsqr

        dRA = B*X + G*Y
        dDec = A*X + F*Y
        
        model.relsep[i] = -sqrt(dRA**2 + dDec**2)*par.sau/a_1
        model.PA[i] = atan2(-dRA, -dDec)

    i1 = data.nRV + data.nAst
    i2 = i1 + data.nHip1
    for i in range(i1, i2):
        X = model.cosEA[i] - par.ecc
        Y = model.sinEA[i]*sqonemeccsqr

        model.dRA_H1[i - i1] = B*X + G*Y
        model.dDec_H1[i - i1] = A*X + F*Y

    i1 = i2
    i2 = i1 + data.nHip2
    for i in range(i1, i2):
        X = model.cosEA[i] - par.ecc
        Y = model.sinEA[i]*sqonemeccsqr

        model.dRA_H2[i - i1] = B*X + G*Y
        model.dDec_H2[i - i1] = A*X + F*Y
    
    i1 = i2
    i2 = i1 + data.nGaia
    for i in range(i1, i2):
        X = model.cosEA[i] - par.ecc
        Y = model.sinEA[i]*sqonemeccsqr

        model.dRA_G[i - i1] = B*X + G*Y
        model.dDec_G[i - i1] = A*X + F*Y
    
    return 

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

######################################################################
# HTOF fit to the epoch astrometry.  This is effectively a translation
# of the actual HTOF package to avoid python overheads.
######################################################################

def calc_PMs_epoch_astrometry(Data data, Model model, AstrometricFitter Hip1,
                              AstrometricFitter Hip2, AstrometricFitter Gaia):

    cdef extern from "lstsq.h":
        void lstsq_C(double A_in[], double b[], int m, int n, double coef[])

    cdef double *b_Hip1 = <double *> PyMem_Malloc(Hip1.npar * sizeof(double))
    cdef double *b_Hip2 = <double *> PyMem_Malloc(Hip2.npar * sizeof(double))
    cdef double *b_Gaia = <double *> PyMem_Malloc(Gaia.npar * sizeof(double))

    cdef double *res_Hip1 = <double *> PyMem_Malloc(Hip1.npar * sizeof(double))
    cdef double *res_Hip2 = <double *> PyMem_Malloc(Hip2.npar * sizeof(double))
    cdef double *res_Gaia = <double *> PyMem_Malloc(Gaia.npar * sizeof(double))

    cdef double *chi2mat_Hip1 = <double *> PyMem_Malloc(Hip1.npar*Hip1.npar * sizeof(double))
    cdef double *chi2mat_Hip2 = <double *> PyMem_Malloc(Hip2.npar*Hip2.npar * sizeof(double))
    cdef double *chi2mat_Gaia = <double *> PyMem_Malloc(Gaia.npar*Gaia.npar * sizeof(double))
    
    cdef int i, j
    cdef double x
    cdef double RA_H1, Dec_H1, pmra_H1, pmdec_H1
    cdef double RA_H2, Dec_H2, pmra_H2, pmdec_H2
    cdef double RA_G, Dec_G, pmra_G, pmdec_G
    
    for i in range(Hip1.npar):
        x = 0
        for j in range(Hip1.nepochs):
            x += model.dRA_H1[j]*Hip1.ra_solution_vecs[i, j]
            x += model.dDec_H1[j]*Hip1.dec_solution_vecs[i, j]
        b_Hip1[i] = x
    
    for i in range(Hip2.npar):
        x = 0
        for j in range(Hip2.nepochs):
            x += model.dRA_H2[j]*Hip2.ra_solution_vecs[i, j]
            x += model.dDec_H2[j]*Hip2.dec_solution_vecs[i, j]
        b_Hip2[i] = x
    
    for i in range(Gaia.npar):
        x = 0
        for j in range(Gaia.nepochs):
            x += model.dRA_G[j]*Gaia.ra_solution_vecs[i, j]
            x += model.dDec_G[j]*Gaia.dec_solution_vecs[i, j]
        b_Gaia[i] = x

    for i in range(Hip1.npar):
        for j in range(Hip1.npar):
            chi2mat_Hip1[i*Hip1.npar + j] = Hip1.chi2_matrix[i, j]
    for i in range(Hip2.npar):
        for j in range(Hip2.npar):
            chi2mat_Hip2[i*Hip2.npar + j] = Hip2.chi2_matrix[i, j]
    for i in range(Gaia.npar):
        for j in range(Gaia.npar):
            chi2mat_Gaia[i*Gaia.npar + j] = Gaia.chi2_matrix[i, j]
            

    lstsq_C(chi2mat_Hip1, b_Hip1, Hip1.npar, Hip1.npar, res_Hip1)
    RA_H1 = res_Hip1[0]
    Dec_H1 = res_Hip1[1]
    pmra_H1 = res_Hip1[2]*365.25
    pmdec_H1 = res_Hip1[3]*365.25
    
    lstsq_C(chi2mat_Hip2, b_Hip2, Hip2.npar, Hip2.npar, res_Hip2)
    RA_H2 = res_Hip2[0]
    Dec_H2 = res_Hip2[1]
    pmra_H2 = res_Hip2[2]*365.25
    pmdec_H2 = res_Hip2[3]*365.25
    
    lstsq_C(chi2mat_Gaia, b_Gaia, Gaia.npar, Gaia.npar, res_Gaia)
    RA_G = res_Gaia[0]
    Dec_G = res_Gaia[1]
    pmra_G = res_Gaia[2]*365.25
    pmdec_G = res_Gaia[3]*365.25
    
    model.pmra_H = 0.4*pmra_H1 + 0.6*pmra_H2
    model.pmdec_H = 0.4*pmdec_H1 + 0.6*pmdec_H2
    model.pmra_G = pmra_G
    model.pmdec_G = pmdec_G
    model.pmra_HG = (RA_G - (0.4*RA_H1 + 0.6*RA_H2))/(data.epRA_G - data.epRA_H)
    model.pmdec_HG = (Dec_G - (0.4*Dec_H1 + 0.6*Dec_H2))/(data.epDec_G - data.epDec_H)
    
    PyMem_Free(b_Hip1)
    PyMem_Free(b_Hip2)
    PyMem_Free(b_Gaia)
    
    PyMem_Free(res_Hip1)
    PyMem_Free(res_Hip2)
    PyMem_Free(res_Gaia)
    
    PyMem_Free(chi2mat_Hip1)
    PyMem_Free(chi2mat_Hip2)
    PyMem_Free(chi2mat_Gaia)

    return

    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

######################################################################
# This piece does the Simpson's rule-type calculation currently
# implemented in Trent's version.  It should generally be replaced by
# the proper calculation done by Mirek's HTOF.
######################################################################

def calc_PMs_no_epoch_astrometry(Data data, Model model):

    cdef double RA_H, Dec_H, RA_G, Dec_G
    
    RA_H = (4*model.dRA_H1[1] + model.dRA_H1[0] + model.dRA_H1[2])/6
    Dec_H = (4*model.dDec_H1[4] + model.dDec_H1[3] + model.dDec_H1[5])/6
    
    RA_G = (4*model.dRA_G[1] + model.dRA_G[0] + model.dRA_G[2])/6
    Dec_G = (4*model.dDec_G[4] + model.dDec_G[3] + model.dDec_G[5])/6

    model.pmra_HG = (RA_G - RA_H)/(data.epRA_G - data.epRA_H)
    model.pmdec_HG = (Dec_G - Dec_H)/(data.epDec_G - data.epDec_H)

    model.pmra_H = (model.dRA_H1[2] - model.dRA_H1[0])/data.dt_H
    model.pmdec_H = (model.dDec_H1[5] - model.dDec_H1[3])/data.dt_H

    model.pmra_G = (model.dRA_G[2] - model.dRA_G[0])/data.dt_G
    model.pmdec_G = (model.dDec_G[5] - model.dDec_G[3])/data.dt_G

    return
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

#######################################################################
# Compute the RVs.  Store these RVs directly in the model object.
# Note that we never explicitly compute the true anomaly; we use trig
# identities to take a shortcut and save a factor of about 20 in
# runtime.
#######################################################################

def calc_RV(Data data, Params par, Model model):

    cdef extern from "math.h" nogil:
        double sin(double _x)
        double cos(double _x)
        double fabs(double _x)
        double sqrt(double _x)

    cdef double pi = 3.14159265358979323846264338327950288
    cdef double pi_d_2 = pi/2.
    
    cdef double ecccosarg = par.ecc*cos(par.arg)
    cdef double sqrt1pe = sqrt(1 + par.ecc)
    cdef double sqrt1me = sqrt(1 - par.ecc)
    cdef double RVampl = 2*pi*par.sau*sin(par.inc)/par.per/sqrt(1 - par.ecc**2)*1731458.33*par.msec/(par.mpri + par.msec)

    cdef double cosarg = cos(par.arg)
    cdef double sinarg = sin(par.arg)
    cdef double sqrt1pe_div_sqrt1me = sqrt((1 + par.ecc)/(1 - par.ecc))
    cdef double TA, ratio, fac, tanEAd2
    
    cdef int i

    ##################################################################
    # Trickery with trig identities.  The code below is mathematically
    # identical to the use of the true anomaly.  If sin(EA) is small
    # and cos(EA) is close to -1, no problem as long as sin(EA) is not
    # precisely zero (set tan(EA/2)=1e100 in this case).  If sin(EA)
    # is small and EA is close to zero, use the fifth-order Taylor
    # expansion for tangent.  This is good to ~1e-15 for EA within
    # ~0.015 of 0.  Assume eccentricity is not precisely unity (this
    # should be forbidden by the priors).  Very, very high
    # eccentricities (significantly above 0.9999) may be problematic.
    # This routine assumes range reduction of the eccentric anomaly to
    # (-pi, pi] and will throw an error if this is violated.
    ##################################################################

    cdef double one_d_24 = 1./24
    cdef double one_d_240 = 1./240
    
    for i in range(data.nRV):

        if fabs(model.sinEA[i]) > 1.5e-2:
            tanEAd2 = (1 - model.cosEA[i])/model.sinEA[i]
        elif model.EA[i] < -pi or model.EA[i] > pi:
            raise ValueError("EA input to calc_RV must be betwen -pi and pi.")
        elif fabs(model.EA[i]) < pi_d_2:
            EA = model.EA[i]
            tanEAd2 = EA*(0.5 + EA**2*(one_d_24 + one_d_240*EA**2))
        elif model.sinEA[i] != 0:
            tanEAd2 = (1 - model.cosEA[i])/model.sinEA[i]
        else:
            tanEAd2 = 1e100
            
        ratio = sqrt1pe_div_sqrt1me*tanEAd2
        fac = 2/(1 + ratio**2)
        model.RV[i] = RVampl*(cosarg*(fac - 1) - sinarg*ratio*fac + ecccosarg)

    # Don't use the following: we do about 20 times better above.
    #for i in range(data.nRV):
    #    TA = 2*atan2(sqrt1pe*sin(model.EA[i]/2), sqrt1me*cos(model.EA[i]/2))
    #    model.RV[i] = RVampl*(cos(TA + par.arg) + par.ecc*cos(par.arg))

    return
        

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

######################################################################
# Compute the log likelihood from the RVs, relative separation,
# position angle, and absolute astrometry, all of which should be in
# model.  The absolute astrometry should now be six proper motions (RA
# and Dec for Hip, HG, and Gaia).  The likelihood is marginalized over
# the RV zero points (one per instrument), the parallax, and the
# center-of-mass proper motion of the system.  The same jitter is
# assumed to apply to all RV instruments.  These marginalizations,
# especially over parallax and center of mass proper motion, are
# responsible for most of the code complexity (but add little to the
# computational cost).
######################################################################

def calcL(Data data, Params par, Model model, freemodel=True):

    cdef int i, j
    cdef double lnL, ivar, dRV
    lnL = 0
    
    cdef extern from "math.h" nogil:
        double sin(double _x)
        double cos(double _x)
        double log(double _x)
        double pow(double _x, double _y)
        double atan2(double _y, double _x)

    cdef extern from "lstsq.h":
        void lstsq_C(double A_in[], double b[], int m, int n, double coef[])
        
    cdef double jitsq = pow(10., par.jit)
    cdef double rv_ivar = 1
    
    cdef double *A = <double *> PyMem_Malloc(data.nInst * sizeof(double))
    if not A:
        raise MemoryError()
    cdef double *B = <double *> PyMem_Malloc(data.nInst * sizeof(double))
    if not B:
        raise MemoryError()
    cdef double *C = <double *> PyMem_Malloc(data.nInst * sizeof(double))
    if not C:
        raise MemoryError()
    
    for i in range(data.nInst):
        A[i] = B[i] = C[i] = 0
    
    for i in range(data.nRV):
        ivar = 1./(data.RV_err[i]**2 + jitsq)
        dRV = data.RV[i] - model.RV[i]
        rv_ivar *= ivar
        
        # prevent underflow
        if rv_ivar < 1e-200:
            lnL += log(rv_ivar)
            rv_ivar = 1

        C[data.RVinst[i]] += dRV**2*ivar
        B[data.RVinst[i]] -= 2*dRV*ivar
        A[data.RVinst[i]] += ivar

    ##################################################################
    # Marginalize out RV, and add the jitter component of the variance
    # to the log likelihood.
    ##################################################################

    for i in range(data.nInst):
        if A[i] == 0:
            continue
        lnL -= -B[i]**2/(4*A[i]) + C[i] + log(A[i])
        
    lnL += log(rv_ivar)

    PyMem_Free(A)
    PyMem_Free(B)
    PyMem_Free(C)

    ##################################################################
    # Ok, tricky part below.  We will take care of the mean proper
    # motion and parallax.  Once we have the matrix equation M*x=b,
    # with x=(plx,pmra,pmdec), we just solve for the best-fit
    # parameters and account for the normalization with the
    # determinant of M.
    ##################################################################

    cdef double *M = <double *> PyMem_Malloc(9 * sizeof(double))
    if not M:
        raise MemoryError()
    cdef double *b = <double *> PyMem_Malloc(3 * sizeof(double))
    if not b:
        raise MemoryError()
    cdef double *res = <double *> PyMem_Malloc(3 * sizeof(double))
    if not res:
        raise MemoryError()
        
    for i in range(3):
        b[i] = res[i] = 0
        for j in range(3):
            M[i*3 + j] = 0
    
    for i in range(data.nAst):
        M[0] += model.relsep[i]**2/data.relsep_err[i]**2
        b[0] += model.relsep[i]*data.relsep[i]/data.relsep_err[i]**2
        lnL -= (atan2(sin(model.PA[i] - data.PA[i]),
                      cos(model.PA[i] - data.PA[i])))**2/data.PA_err[i]**2

    M[0] += model.pmra_H**2*data.Cinv_H[0, 0]
    M[0] += 2*model.pmra_H*model.pmdec_H*data.Cinv_H[1, 0]
    M[0] += model.pmdec_H**2*data.Cinv_H[1, 1]

    M[0] += model.pmra_HG**2*data.Cinv_HG[0, 0]
    M[0] += 2*model.pmra_HG*model.pmdec_HG*data.Cinv_HG[1, 0]
    M[0] += model.pmdec_HG**2*data.Cinv_HG[1, 1]

    M[0] += model.pmra_G**2*data.Cinv_G[0, 0]
    M[0] += 2*model.pmra_G*model.pmdec_G*data.Cinv_G[1, 0]
    M[0] += model.pmdec_G**2*data.Cinv_G[1, 1]

    M[0] += 1/data.plx_err**2
    
    M[1] = model.pmra_H*data.Cinv_H[0, 0] + model.pmdec_H*data.Cinv_H[0, 1]
    M[1] += model.pmra_HG*data.Cinv_HG[0, 0] + model.pmdec_HG*data.Cinv_HG[0, 1]
    M[1] += model.pmra_G*data.Cinv_G[0, 0] + model.pmdec_G*data.Cinv_G[0, 1]
    
    M[2] = model.pmdec_H*data.Cinv_H[1, 1] + model.pmra_H*data.Cinv_H[0, 1]
    M[2] += model.pmdec_HG*data.Cinv_HG[1, 1] + model.pmra_HG*data.Cinv_HG[0, 1]
    M[2] += model.pmdec_G*data.Cinv_G[1, 1] + model.pmra_G*data.Cinv_G[0, 1]
    
    b[0] += model.pmra_H*data.pmra_H*data.Cinv_H[0, 0]
    b[0] += model.pmdec_H*data.pmdec_H*data.Cinv_H[1, 1]
    b[0] += (model.pmdec_H*data.pmra_H + model.pmra_H*data.pmdec_H)*data.Cinv_H[0, 1]

    b[0] += model.pmra_HG*data.pmra_HG*data.Cinv_HG[0, 0]
    b[0] += model.pmdec_HG*data.pmdec_HG*data.Cinv_HG[1, 1]
    b[0] += (model.pmdec_HG*data.pmra_HG + model.pmra_HG*data.pmdec_HG)*data.Cinv_HG[0, 1]

    b[0] += model.pmra_G*data.pmra_G*data.Cinv_G[0, 0]
    b[0] += model.pmdec_G*data.pmdec_G*data.Cinv_G[1, 1]
    b[0] += (model.pmdec_G*data.pmra_G + model.pmra_G*data.pmdec_G)*data.Cinv_G[0, 1]

    b[0] += data.plx/data.plx_err**2
    
    M[1*3] = M[1]
    M[2*3] = M[2]

    M[1*3 + 1] = data.Cinv_H[0, 0] + data.Cinv_HG[0, 0] + data.Cinv_G[0, 0]
    M[1*3 + 2] = data.Cinv_H[0, 1] + data.Cinv_HG[0, 1] + data.Cinv_G[0, 1]
    M[2*3 + 1] = M[1*3 + 2]
    M[2*3 + 2] = data.Cinv_H[1, 1] + data.Cinv_HG[1, 1] + data.Cinv_G[1, 1]

    b[1] = data.pmra_H*data.Cinv_H[0, 0] + data.pmdec_H*data.Cinv_H[0, 1]
    b[1] += data.pmra_HG*data.Cinv_HG[0, 0] + data.pmdec_HG*data.Cinv_HG[0, 1]
    b[1] += data.pmra_G*data.Cinv_G[0, 0] + data.pmdec_G*data.Cinv_G[0, 1]

    b[2] = data.pmdec_H*data.Cinv_H[1, 1] + data.pmra_H*data.Cinv_H[0, 1]
    b[2] += data.pmdec_HG*data.Cinv_HG[1, 1] + data.pmra_HG*data.Cinv_HG[0, 1]
    b[2] += data.pmdec_G*data.Cinv_G[1, 1] + data.pmra_G*data.Cinv_G[0, 1]

    cdef double plx_best, pmra_best, pmdec_best, deltaRA, deltaDec, detM
    cdef double plx_best2, pmra_best2, pmdec_best2

    ##################################################################
    # We will need the determinant of M later; M might be destroyed by
    # the SVD.  Compute it now without numpy.
    ##################################################################
    
    detM = M[0*3 + 0]*(M[1*3 + 1]*M[2*3 + 2] - M[1*3 + 2]*M[2*3 + 1])
    detM -= M[0*3 + 1]*(M[1*3 + 0]*M[2*3 + 2] - M[1*3 + 2]*M[2*3 + 0])
    detM += M[0*3 + 2]*(M[1*3 + 0]*M[2*3 + 1] - M[1*3 + 1]*M[2*3 + 0])

    ##################################################################
    # Solve the matrix equation for the parallax, pmra, and pmdec.
    ##################################################################
    
    lstsq_C(M, b, 3, 3, res)
    plx_best = res[0]
    pmra_best = res[1]
    pmdec_best = res[2]

    PyMem_Free(M)
    PyMem_Free(b)
    PyMem_Free(res)
    
    ##################################################################
    # Now take care of the rest of the log likelihood.
    ##################################################################
    
    cdef double dlnL = 0
    
    for i in range(data.nAst):
        dlnL -= (model.relsep[i]*plx_best - data.relsep[i])**2/data.relsep_err[i]**2
        
    deltaRA = plx_best*model.pmra_H - data.pmra_H + pmra_best
    deltaDec = plx_best*model.pmdec_H - data.pmdec_H + pmdec_best
    dlnL -= deltaRA**2*data.Cinv_H[0, 0]
    dlnL -= deltaDec**2*data.Cinv_H[1, 1]
    dlnL -= 2*deltaRA*deltaDec*data.Cinv_H[0, 1]
    
    deltaRA = plx_best*model.pmra_HG - data.pmra_HG + pmra_best
    deltaDec = plx_best*model.pmdec_HG - data.pmdec_HG + pmdec_best
    dlnL -= deltaRA**2*data.Cinv_HG[0, 0]
    dlnL -= deltaDec**2*data.Cinv_HG[1, 1]
    dlnL -= 2*deltaRA*deltaDec*data.Cinv_HG[0, 1]
    
    deltaRA = plx_best*model.pmra_G - data.pmra_G + pmra_best
    deltaDec = plx_best*model.pmdec_G - data.pmdec_G + pmdec_best
    dlnL -= deltaRA**2*data.Cinv_G[0, 0]
    dlnL -= deltaDec**2*data.Cinv_G[1, 1]
    dlnL -= 2*deltaRA*deltaDec*data.Cinv_G[0, 1]
    
    dlnL -= (data.plx - plx_best)**2/data.plx_err**2
    
    dlnL -= log(detM)
    lnL += dlnL

    if freemodel:
        PyMem_Free(model.EA)
        PyMem_Free(model.sinEA)
        PyMem_Free(model.cosEA)
        PyMem_Free(model.RV)
        PyMem_Free(model.relsep)
        PyMem_Free(model.PA)
        PyMem_Free(model.dRA_H1)
        PyMem_Free(model.dDec_H1)
        PyMem_Free(model.dRA_H2)
        PyMem_Free(model.dDec_H2)
        PyMem_Free(model.dRA_G)
        PyMem_Free(model.dDec_G)
        
    return lnL/2

    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

######################################################################
# Priors for eccentricity, masses, inclination, semimajor axis.  Prior
# for parallax is included in likelihood so that we can marginalize
# the quantity out.
######################################################################

def lnprior(Params par):

    cdef extern from "math.h" nogil:
        double sin(double _x)
        double log(double _x)

    cdef double pi = 3.14159265358979323846264338327950288 # np.pi
    cdef double zeroprior = -np.inf
    
    if par.sau <= 0 or par.mpri <= 0 or par.msec < 1e-4 or par.ecc >= 1:
        return zeroprior
    if par.inc <= 0 or par.inc >= pi or par.asc < -pi or par.asc >= 3*pi:
        return zeroprior
    if par.lam < -pi or par.lam >= 3*pi or par.jit < -20 or par.jit > 20:
        return zeroprior

    return log(sin(par.inc)*1./(par.sau*par.msec*par.mpri))
