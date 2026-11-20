
############### Original unchanged orbit_plot ##################


import os
import numpy as np
from astropy.time import Time
import astropy
import warnings
from random import randrange
from orvara import corner_modified
from scipy.interpolate import interp1d
from scipy import stats, ndimage
from scipy.stats import kde
from orvara import orbit
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from orvara.format_fits import pull_chain_params, burnin_chain

try:
    # astropy re-exported ErfaWarning until 6.1; erfa is the real home.
    from erfa import ErfaWarning
except ImportError:                                    # pragma: no cover
    from astropy.utils.exceptions import ErfaWarning

rc('font',**{'size':11.5})

#: Proper-motion errors at or above this (mas/yr) mark a placeholder row
#: written by orvara.gaia_accel in gaia_only mode, not a measurement.
SYNTHETIC_PM_ERROR = 1e6


class Orbit:

    def __init__(self, OP, step='best', epochs='custom'):

        data = orbit.Data(OP.Hip, OP.HGCAFile, OP.RVfile, OP.relAstfile, verbose=False)

        if isinstance(epochs, list) or isinstance(epochs, np.ndarray):
            data.custom_epochs(epochs, iplanet=OP.iplanet)
        elif epochs == 'custom':
            data.custom_epochs(OP.epoch, iplanet=OP.iplanet)
        elif epochs == 'observed':
            pass
        else:
            raise ValueError("epochs should be 'custom', 'observed', or array-like")
        model = orbit.Model(data)

        if step == 'best':
            step = np.where(OP.chain['lnp'] == np.amax(OP.chain['lnp']))[0][0]
            
        self.plx = OP.chain['plx_ML'][step]
        jit_per_inst = 'jitter0' in [col.name for col in OP.chain.columns]
        if jit_per_inst:
            nj = data.nInst
        else:
            nj = 1
        for i in range(OP.nplanets):
            par = orbit.Params(pull_chain_params(OP.chain.columns, step), 
                               iplanet=i, nplanets=OP.nplanets, 
                               ninst_RV=data.nInst, ninst_jit=nj)
            if i == OP.iplanet:
                self.par = par
            orbit.calc_EA_RPP(data, par, model)
            orbit.calc_offsets(data, par, model, i)
            orbit.calc_RV(data, par, model)
            mu_RA, mu_Dec = model.return_proper_motions(par)
            if i == 0:
                self.mu_RA = mu_RA
                self.mu_Dec = mu_Dec
            else:
                self.mu_RA += mu_RA
                self.mu_Dec += mu_Dec
            
        # most likely orbit for delRA and delDec
        self.dRAs_G, self.dDecs_G, self.dRAs_H1, self.dDecs_H1, self.dRAs_H2, self.dDecs_H2 = model.return_dRA_dDec()
        self.dras, self.ddecs = self.dRAs_G*self.plx*1000, self.dDecs_G*self.plx*1000
        self.RV = model.return_RVs()
        self.relsep = model.return_relsep()*self.plx
        self.PA = (model.return_PAs()*180/np.pi) % 360

        self.mu_RA_CM = 1e3*OP.chain['pmra_ML'][step]
        self.mu_Dec_CM = 1e3*OP.chain['pmdec_ML'][step]
        
        self.mu_RA = 1e3*self.mu_RA*self.plx + self.mu_RA_CM
        self.mu_Dec = 1e3*self.mu_Dec*self.plx + self.mu_Dec_CM
        try:
            self.offset = OP.calc_RV_offset(step)
        except Exception:
            self.offset = 0
        if OP.cmref == 'msec_solar':
            self.colorpar = self.par.msec
        elif OP.cmref == 'msec_jup':
            self.colorpar = self.par.msec*1989/1.898
        elif OP.cmref == 'ecc':
            self.colorpar = self.par.ecc

        model.free()


class OrbitPlots:

###################################### Initialize Class ############################################

    def __init__(self):
        pass        

    def start(self):
        self.cmlabel_dic = {'msec_jup': r'$\mathrm{M_{comp} (M_{Jup})}$','msec_solar': r'$\mathrm{M_{comp} (M_{\odot})}$', 'ecc': 'Eccentricity'}
        self.color_list = ['r', 'g', 'b', 'y', 'c', 'b']
        
        ############################### load in data #######################
        # define epochs
        self.epoch, self.epoch_calendar = self.define_epochs()
        # load mcmc data
        self.chain = self.load_mcmc_data()
        self.rand_idx = [randrange(len(self.chain['lnp'])) for i in range(self.num_orbits)]

        # load observed RV data
        try:
            self.epoch_obs, self.RV_obs, self.RV_obs_err, self.nInst, self.epoch_obs_dic, self.RV_obs_dic, self.RV_obs_err_dic = self.load_obsRV_data()
            self.have_RVdat = True
        except Exception:
            self.have_RVdat = False
        # load relative astrometry data:
        try: #if os.access(self.relAstfile,os.R_OK):
            self.have_reldat = True
            self.ep_relAst_obs, self.relsep_obs, self.relsep_obs_err, self.PA_obs, self.PA_obs_err, self.ast_indx = self.load_relAst_data(self.iplanet)
        except Exception: #else:
            self.have_reldat = False
        # load HGCA data:
        self.ep_mualp_obs, self.ep_mudec_obs, self.mualp_obs, self.mudec_obs, self.mualp_obs_err, self.mudec_obs_err = self.load_HGCA_data()        
        
        ################################ set colorbar ######################
        # setup the normalization and the colormap

        msec_name = 'msec%d' % (self.iplanet)
        if self.cmref == 'msec_jup':
            vmin = 1989/1.898*stats.scoreatpercentile(self.chain[msec_name], 1)
            vmax = 1989/1.898*stats.scoreatpercentile(self.chain[msec_name], 99)
        elif self.cmref == 'msec_solar':
            vmin = stats.scoreatpercentile(self.chain[msec_name], 1)
            vmax = stats.scoreatpercentile(self.chain[msec_name], 99)
        elif self.cmref == 'ecc':
            ip = '%d' % (self.iplanet)
            ecc = self.chain['esino' + ip]**2 + self.chain['ecoso' + ip]**2
            vmin = stats.scoreatpercentile(ecc, 1)
            vmax = stats.scoreatpercentile(ecc, 99)
        else:
            raise ValueError("Reference parameter for color should be 'msec_jup', 'msec_solar', or 'ecc'")
            
        self.normalize = mcolors.Normalize(vmin=vmin, vmax=vmax)
        self.colormap = getattr(cm, self.color_map)
        self.sm = cm.ScalarMappable(norm=self.normalize, cmap=self.colormap)
        self.sm.set_array(np.linspace(vmin, vmax, 10))
        
        print("Generating plots for target " + self.title)

###################################### Define Functions ############################################
    def JD_to_calendar(self, JD):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=ErfaWarning)
            return Time(JD, format='jd').decimalyear

    def calendar_to_JD(self, year):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=ErfaWarning)
            return Time(year, format='decimalyear').jd

    def define_epochs(self):
        """
            Function to define a custom range of epochs
        """
        start_epoch = self.calendar_to_JD(self.start_epoch)
        end_epoch = self.calendar_to_JD(self.end_epoch)
        range_epoch = end_epoch - start_epoch
        epoch = np.linspace(start_epoch, end_epoch + 0.5*range_epoch, self.num_steps)
        epoch_calendar = np.zeros(len(epoch))
        for i in range(len(epoch_calendar)):
            epoch_calendar[i] = self.JD_to_calendar(epoch[i])
        return epoch, epoch_calendar

    def load_mcmc_data(self):
        """
            Function to load in the MCMC chain from fit_orbit
        """
        chain = fits.open(self.MCMCfile)[1].data
        chain = burnin_chain(chain.columns, self.burnin, reshape=True)
        self.lnp = chain['lnp']
        self.nplanets = int(fits.open(self.MCMCfile)[0].header['nplanets'])
        return chain
        
    def load_obsRV_data(self):
        """
            Function to load in the observed RV data
        """
        rvdat = np.genfromtxt(self.RVfile)
        epoch_obs = rvdat[:, 0]
        RV_obs = rvdat[:, 1]
        RV_obs_err = rvdat[:, 2]
        try:
            RVinst = (rvdat[:, 3]).astype(np.int32)
            self.RVinst = RVinst
            # Check to see that the column we loaded was an integer
            assert np.all(RVinst == rvdat[:, 3])
            nInst = int(np.amax(rvdat[:, 3]) + 1)
            
            self.multi_instr = True
        except Exception:
            self.multi_instr = False
            nInst = 1
            self.RVinst = (RV_obs*0).astype(int)
    
        idx_dic = {}
        epoch_obs_dic = {}
        RV_obs_dic = {}
        RV_obs_err_dic = {}
        
        for i in range(nInst):
            idx_dic[i] = (np.where(self.RVinst == i)[0])
            epoch_obs_dic[i] = epoch_obs[idx_dic[i]]
            RV_obs_dic[i] = RV_obs[idx_dic[i]]
            RV_obs_err_dic[i] = RV_obs_err[idx_dic[i]]
        return epoch_obs, RV_obs, RV_obs_err, nInst, epoch_obs_dic, RV_obs_dic, RV_obs_err_dic

    def load_relAst_data(self, iplanet=None):
        """
            Function to load in the relative astrometry data
        """
        reldat = np.genfromtxt(self.relAstfile, usecols=(0,1,2,3,4))
        if len(reldat.shape) == 1:
            reldat = np.reshape(reldat, (1, -1))

        try:
            icompanion = np.genfromtxt(self.relAstfile, usecols=(6)).astype(int)
        except Exception:
            icompanion = np.zeros(reldat.shape[0]).astype(int)

        indx = np.where(icompanion == iplanet)
        if iplanet is not None:
            reldat = reldat[indx]
        if len(reldat) == 0:
            return None

        # Try to guess whether we should assume the epochs of the
        # relative astrometry file to be decimal years or JD.
        if np.median(reldat[:, 0]) < 3000:
            relep = (reldat[:, 0] - 2000)*365.25 + 2451544.5
        else:
            relep = reldat[:, 0]

        ep_relAst_obs = relep #reldat[:, 0]

        relsep_obs = reldat[:, 1]
        relsep_obs_err = reldat[:, 2]
        PA_obs = reldat[:, 3]
        PA_obs_err = reldat[:, 4]
        
        return ep_relAst_obs, relsep_obs, relsep_obs_err, PA_obs, PA_obs_err, indx
    
    def load_HGCA_data(self):
        """
        Load proper motions, NSS astrometric acceleration (and jerk for Acc9)
        from the one-star HGCA FITS file built by gaia_accel_gen_fit.

        Sets on self:
          self.gaia_npar          : 5, 7, or 9
          self.ep_gaia_acc        : decimal year of Gaia central epoch (RA)
          self.nss_accra,  _err   : NSS accel in RA  (mas/yr²)
          self.nss_accdec, _err   : NSS accel in Dec (mas/yr²)
          self.nss_jerkra, _err   : NSS jerk  in RA  (mas/yr³), None if npar<9
          self.nss_jerkdec,_err   : NSS jerk  in Dec (mas/yr³), None if npar<9
        """
        t = fits.open(self.HGCAFile)[1].data
        try:
            self.have_pmdat = True
            # np.where returns an array of indices; NumPy >= 2.0 refuses to
            # convert a length-1 array to a Python int, and this whole block
            # is inside a try/except, so getting it wrong silently drops the
            # proper-motion panels instead of raising.
            matches = np.where(t['hip_id'] == self.Hip)[0]
            if len(matches) == 0:
                raise KeyError(f'HIP {self.Hip} is not in {self.HGCAFile}')
            i = int(matches[0])
            ep_mualp_obs = np.array([t['epoch_ra_hip'][i], t['epoch_ra_gaia'][i]])
            ep_mudec_obs = np.array([t['epoch_dec_hip'][i], t['epoch_dec_gaia'][i]])
            mualp_obs = np.array([t['pmra_hip'][i], t['pmra_gaia'][i]])
            mualp_obs_err = np.array([t['pmra_hip_error'][i], t['pmra_gaia_error'][i]])
            mudec_obs = np.array([t['pmdec_hip'][i], t['pmdec_gaia'][i]])
            mudec_obs_err = np.array([t['pmdec_hip_error'][i], t['pmdec_gaia_error'][i]])
            for k in range(len(ep_mualp_obs)):
                ep_mualp_obs[k] = self.calendar_to_JD(ep_mualp_obs[k])
                ep_mudec_obs[k] = self.calendar_to_JD(ep_mudec_obs[k])

            # A Gaia-only NSS file (orvara.gaia_accel, mode = gaia_only) carries
            # a placeholder Hipparcos row with 1e7 mas/yr errors so that it
            # drops out of the likelihood.  Drop it from the plot too, or the
            # error bar sets the y-range to +/- 1e7.
            real = mualp_obs_err < SYNTHETIC_PM_ERROR
            ep_mualp_obs, ep_mudec_obs = ep_mualp_obs[real], ep_mudec_obs[real]
            mualp_obs, mualp_obs_err = mualp_obs[real], mualp_obs_err[real]
            mudec_obs, mudec_obs_err = mudec_obs[real], mudec_obs_err[real]
            # With a single proper-motion epoch (a Gaia-only NSS source) the two
            # free barycentric proper-motion parameters reproduce the Gaia value
            # exactly, so a proper-motion O-C carries no information.
            self.n_pm_epochs = int(np.count_nonzero(real))

            # ── NSS acceleration / jerk (Acc7 or Acc9) ──────────────────────
            col_names = t.names if hasattr(t, 'names') else t.dtype.names
            self.gaia_npar = int(t['gaia_npar'][i]) if 'gaia_npar' in col_names else 5
            # Central epoch of the Gaia window, per coordinate.  For an NSS
            # file this is 2016.0 in both.
            self.ep_gaia_acc = float(t['epoch_ra_gaia'][i])
            self.ep_gaia_acc_dec = float(t['epoch_dec_gaia'][i])
            # Window length: the NSS model's 2 x 517.5 d for a one-star NSS
            # file, the HGCA's effective baseline otherwise.
            try:
                from orvara.main import identify_hgca_version
                self.dt_gaia = identify_hgca_version(self.HGCAFile)[1]
            except Exception:
                self.dt_gaia = 2.76

            if self.gaia_npar >= 7 and 'accra_gaia' in col_names:
                self.nss_accra = float(t['accra_gaia'][i])
                self.nss_accra_err = float(t['accra_gaia_error'][i])
                self.nss_accdec = float(t['accdec_gaia'][i])
                self.nss_accdec_err = float(t['accdec_gaia_error'][i])
                self.have_nss_acc = True
            else:
                self.have_nss_acc = False
                self.nss_accra = self.nss_accra_err = None
                self.nss_accdec = self.nss_accdec_err = None

            if self.gaia_npar == 9 and 'jerkra_gaia' in col_names:
                self.nss_jerkra = float(t['jerkra_gaia'][i])
                self.nss_jerkra_err = float(t['jerkra_gaia_error'][i])
                self.nss_jerkdec = float(t['jerkdec_gaia'][i])
                self.nss_jerkdec_err = float(t['jerkdec_gaia_error'][i])
                self.have_nss_jerk = True
            else:
                self.have_nss_jerk = False
                self.nss_jerkra = self.nss_jerkra_err = None
                self.nss_jerkdec = self.nss_jerkdec_err = None

        except Exception as e:
            self.have_pmdat = False
            self.n_pm_epochs = 0
            self.have_nss_acc = False
            self.have_nss_jerk = False
            return [None, None, None, None, None, None]
        return ep_mualp_obs, ep_mudec_obs, mualp_obs, mudec_obs, mualp_obs_err, mudec_obs_err
    
    def calc_RV_offset(self, step):
        """
            Function to calculate the offset of the observed RV data
        """
        try:
            # calculate the offsets of the RV curves
            assert self.multi_instr
            offset_dic = {}
            for i in range(self.nInst):
                offset = self.chain['RV_ZP_%d_ML' % (i)][step]
                offset_dic[i] = offset
        except Exception:
            offset = offset_dic = [self.chain['RV_ZP_0_ML'][step]]
        return offset_dic


    def _year_ticks(self, ax, nbins=6):
        """
        Readable epoch ticks.  Over a span of 4 yr or more, label whole years
        only and at most ~nbins of them; matplotlib's default otherwise puts a
        label every 2.5 yr ('2012.5 2015.0 ...'), which overlap.
        """
        lo, hi = ax.get_xlim()
        whole = (hi - lo) >= 4
        ax.xaxis.set_major_locator(
            MaxNLocator(nbins=nbins, integer=whole, steps=[1, 2, 2.5, 5, 10]))
        if whole:
            ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    def align_ylabels(self, fig, axes, pad=0.02):
        xmin = np.inf
        for axis in axes:
            xmin = min(xmin, axis.yaxis.get_ticklabel_extents(fig.canvas.get_renderer())[0].get_points()[0, 0])
        for axis in axes:
            bb = axis.get_window_extent(fig.canvas.get_renderer()).get_points()
            dx = bb[1, 0] - bb[0, 0]
            x0 = bb[0, 0]
            axis.yaxis.set_label_coords((xmin - x0)/dx - pad, 0.5)
        return


########################################## Plotting ################################################

# 1. Astrometric plots

    ###############################################################################################
    ############################## plot astrometric orbits ######################
    
    def closed_orbit(self, par, plx, nodes=False, n=1000):
        """
        Compute a closed orbit (EA: 0 -> 2pi) to trace out.
        If nodes is specified, get the nodes, periastron position, and 
        a point just off to compute the direction of motion.
        """
        A = np.cos(par.arg)*np.cos(par.asc) 
        A -= np.sin(par.arg)*np.sin(par.asc)*np.cos(par.inc)
        B = np.cos(par.arg)*np.sin(par.asc)
        B += np.sin(par.arg)*np.cos(par.asc)*np.cos(par.inc)
        F = -np.sin(par.arg)*np.cos(par.asc) 
        F -= np.cos(par.arg)*np.sin(par.asc)*np.cos(par.inc)
        G = -np.sin(par.arg)*np.sin(par.asc)
        G += np.cos(par.arg)*np.cos(par.asc)*np.cos(par.inc)

        if nodes:
            eccterm = np.sqrt((1 - par.ecc)/(1 + par.ecc))
            EA = [-1e-3, 0, 2*np.arctan(eccterm*np.tan((np.pi - par.arg)/2.)),
                  2*np.arctan(eccterm*np.tan(-par.arg/2.))]
        else:
            EA = np.linspace(0, 2*np.pi, n)
        
        X = np.cos(EA) - par.ecc
        Y = np.sin(EA)*np.sqrt(1 - par.ecc**2)
        
        dra = (B*X + G*Y)*(par.sau)*plx #changed
        ddec = (A*X + F*Y)*(par.sau)*plx #changed
        
        return dra, ddec

    def astrometry(self):

        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)

        # plot the num_orbits randomly selected curves
        for i in range(self.num_orbits):

            orb = Orbit(self, step=self.rand_idx[i])
            dra, ddec = self.closed_orbit(orb.par, orb.plx)

            ax.plot(dra, ddec, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.4, linewidth=0.8)

        #plot the most likely orbit
        
        orb_ml = Orbit(self, step='best')
        dra, ddec = self.closed_orbit(orb_ml.par, orb_ml.plx)
        ax.plot(dra, ddec, color='black')

        # plot the relAst data points

        if self.have_reldat:
            ra_obs = self.relsep_obs * np.sin(self.PA_obs*np.pi /180.)
            dec_obs = self.relsep_obs * np.cos(self.PA_obs*np.pi /180.)
            ax.scatter(ra_obs, dec_obs, s=45, facecolors=self.marker_color, edgecolors='none', zorder=99)
            ax.scatter(ra_obs, dec_obs, s=45, facecolors='none', edgecolors='k', zorder=100)

        # plot the predicted positions (set in config.ini)
        epoch_int = []
        for year in self.epoch_calendar:
            epoch_int.append(int(year))

        epochs = [self.calendar_to_JD(float(year)) for year in self.predicted_ep]
        orb_ml = Orbit(self, 'best', epochs=epochs)

        # define some functions to use later
        def calc_linear(x,y):
            x1, x2 = x[0], x[1]
            y1, y2 = y[0], y[1]
            m = (y1- y2)/(x1 - x2)
            b = y1 - m*x1
            return m,b
                    
        # new method to rotate the labels according to angle of normal to the curve tangent
        
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xlim=[x0 - 0.15*(x1 - x0), x1 + 0.15*(x1 - x0)]
        ylim=[y0 - 0.15*(y1 - y0), y1 + 0.15*(y1 - y0)]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        i = 0
        for year in self.predicted_ep:
            year = int(year)

            # rotate the labels
            # create a list containing the data for the most likely orbits for all the predicted epochs because the predicted epochs lie on the most likely orbit only

            t = np.asarray(self.epoch_calendar)
            
            y, x = [orb_ml.relsep[i]*np.cos(orb_ml.PA[i]*np.pi/180),
                    orb_ml.relsep[i]*np.sin(orb_ml.PA[i]*np.pi/180)]
            _dy, _dx = [-orb_ml.mu_Dec[i] + orb_ml.mu_Dec_CM,
                        -orb_ml.mu_RA[i] + orb_ml.mu_RA_CM]
            i += 1

            m, b = calc_linear([x, x + _dx], [y, y + _dy])
            if not np.isfinite(m) or not np.isfinite(b):
                # some ill constrained orbits will have issues plotting predicted epochs.
                # skip if we find such an epoch.
                print(f'Failed to plot the predicted epoch {year} in the "Astrometry orbits plot." ')
                continue
            # plot the predicted epochs data points on the most likely orbit
            ax.scatter(x, y, s=55, facecolors='none', edgecolors='k', zorder=100)
            # calculate the angle between the tangent line and the x-axis
            y_intercept = m*xlim[0]+b
            x_intercept = (ylim[0]-b)/m
            angle = np.arctan(np.abs(y_intercept - ylim[0])/np.abs(x_intercept - xlim[0]))*180./np.pi
            angle -= 90

            # calculate the line of semi-major and semi-minor axes of the ecplise
            max_y = ddec[(np.where(dra == max(dra)))[0][0]]
            min_y = ddec[(np.where(dra == min(dra)))[0][0]]
            max_x = dra[(np.where(ddec == max(ddec)))[0][0]]
            min_x = dra[(np.where(ddec == min(ddec)))[0][0]]
            # calculate the slopes and y-intercepts of the semi-major and semi-minor axes
            m_semimajor, b_semimajor = calc_linear([min_x,max_x],[min(ddec),max(ddec)])
            m_semiminor, b_semiminor = calc_linear([min(dra),max(dra)],[min_y,max_y])

            # if the point you are labeling has x value larger than the corresponding x value on the semimajor axis (same y), and has y value larger than the corresponding y value on the semiminor axis (same x), then it's in the first quardret and the labels should be aligned to the left and rotated by -angle (here the angle is already corrected by 90 degrees in line 486)
            if x < (y-b_semimajor)/m_semimajor and y > m_semiminor*x + b_semiminor:
                ax.annotate('  ' + str(year), xy=(x, y), verticalalignment='bottom', horizontalalignment='left',rotation =-angle, size=9.5)
            # if the point you are labeling has x value larger than the corresponding x value on the semimajor axis (same y), and has y value smaller than the corresponding y value on the semiminor axis (same x), then it's in the fourth quardret and the labels should be aligned to the left rotated by +angle
            elif x < (y-b_semimajor)/m_semimajor and y < m_semiminor*x + b_semiminor:
                ax.annotate('  ' + str(year), xy=(x, y), verticalalignment='top', horizontalalignment='left',rotation =angle, size=9.5)
            # if the point you are labeling has x value smaller than the corresponding x value on the semimajor axis (same y), and has y value larger than the corresponding y value on the semiminor axis (same x), then it's in the second quardret and the labels should be aligned to the right and rotated by -angle
            elif x > (y-b_semimajor)/m_semimajor and y < m_semiminor*x + b_semiminor:
                ax.annotate(str(year)+'  ', xy=(x, y), verticalalignment='top', horizontalalignment='right',rotation= -angle, size=9.5)
            # if the point you are labeling has x value smaller than the corresponding x value on the semimajor axis (same y), and has y value smaller than the corresponding y value on the semiminor axis (same x), then it's in the third quardret and the labels should be aligned to the right and rotated by +angle
            elif x > (y-b_semimajor)/m_semimajor and y > m_semiminor*x + b_semiminor:
                ax.annotate(str(year)+'  ', xy=(x, y), verticalalignment='bottom', horizontalalignment='right', rotation=angle, size=9.5)

        # plot line of nodes, periastron and the direction of motion of the companion
        dra_nd, ddec_nd = self.closed_orbit(orb_ml.par, orb_ml.plx, nodes=True)

        ax.plot(dra_nd[2:], ddec_nd[2:], 'k--',linewidth = 1)
        ax.plot([0, dra_nd[1]], [0, ddec_nd[1]], 'k:')
        # add arrow
        arrow = mpatches.FancyArrowPatch((dra_nd[0], ddec_nd[0]), (dra_nd[1], ddec_nd[1]),  arrowstyle='->', mutation_scale=25, zorder=100)

        ax.add_patch(arrow)
        ax.plot(0, 0, marker='*',  color='black', markersize=10)
        
        if self.show_title:
            ax.set_title('Astrometric Orbits')
        if self.add_text:
            ax.text(self.x_text,self.y_text, self.text_name, horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,fontsize=15)
        if self.usecolorbar:
            cbar = fig.colorbar(self.sm, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
            cbar.ax.get_yaxis().labelpad=20
        
        ax.set_aspect(np.abs((x0-x1)/(y0-y1)))
        # invert axis
        ax.invert_xaxis()
        # set ticks
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        # set labels and title
        ax.set_xlabel(r'$\mathrm{\Delta \alpha}$ (arcsec)', fontsize=14)
        ax.set_ylabel(r'$\mathrm{\Delta \delta}$ [arcsec]', fontsize=14)

        print("Plotting Astrometry orbits, your plot is generated at " + self.outputdir)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            plt.tight_layout()
            plt.savefig(os.path.join(self.outputdir,'astrometric_orbit_' + self.title)+'.pdf', bbox_inches='tight',transparent=True) # or +'.png'

# 2. RV orbits plot

    ################################################################################################
    ########################### plot the RV full orbits #######################
    
    def RV_fullorbit(self):
        fig = plt.figure(figsize=(5.7, 6.5))
        ax = fig.add_axes((0.15, 0.3, 0.8, 0.6))
        orb_ml = Orbit(self, 'best')
                
        # plot the 50 randomly selected curves
        for i in range(self.num_orbits):
            orb = Orbit(self, step=self.rand_idx[i])
            orb.RV += orb_ml.offset[0] - orb.offset[0]
                
            ax.plot(self.epoch_calendar, orb.RV, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)

        # plot the most likely one
        ax.plot(self.epoch_calendar, orb_ml.RV, color='black')

        # plot the observed data points (RV & relAst)
        if self.have_RVdat:
            rv_epoch_list = []
            for i in range(self.nInst):
                epoch_obs_Inst = np.zeros(len(self.epoch_obs_dic[i]))
                for j in range(len(self.epoch_obs_dic[i])):
                    epoch_obs_Inst[j] = self.JD_to_calendar(self.epoch_obs_dic[i][j])
                rv_epoch_list.append(epoch_obs_Inst)

            jit_ml = orb_ml.par.return_jitters()
        
            for i in range(self.nInst):
                ax.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], yerr=np.sqrt(self.RV_obs_err_dic[i]**2 + jit_ml[i]**2),
                            fmt=self.color_list[i]+'o', ecolor='black', alpha = 0.8, zorder = 299)
                ax.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], facecolors='none', edgecolors='k', alpha = 0.8, zorder=300)
           
        if self.set_limit:
            ax.set_xlim(float(self.user_xlim[0]), float(self.user_xlim[1]))
            ax.set_ylim(float(self.user_ylim[0]),float(self.user_ylim[1]))
            
        ax.set_xlim(self.start_epoch, self.end_epoch)
        self._year_ticks(ax)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))

        if self.show_title:
            ax.set_title('RV vs. Epoch')
        if self.add_text:
            ax.text(self.x_text,self.y_text, self.text_name, horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,fontsize=15)
        if self.usecolorbar:
            cbar_ax = fig.add_axes([1.3, 0.325, 0.03, 0.55])
            cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
            cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
            cbar.ax.get_yaxis().labelpad=20
            fig.delaxes(cbar_ax)
            
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax.set_xlabel('Epoch (year)', labelpad = 10, fontsize=13)
        ax.set_ylabel('RV (m/s)', fontsize=13)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            plt.tight_layout()
            print("Plotting RV orbits, your plot is generated at " + self.outputdir)
            plt.savefig(os.path.join(self.outputdir, 'RV_orbit_' + self.title)+'.pdf', transparent=True, bbox_inches='tight', dpi=200) # or +'.pdf'


# 3. relative RV plot

    ################################################################################################
    ############### plot the RV and O-C #####################

    def RV(self):
        if not self.have_RVdat:
            print("Cannot plot RV with O-C: no RV observations")
            return
        fig = plt.figure(figsize=(4.7, 5.5))
        ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
        ax2 = fig.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height

        orb_ml = Orbit(self, 'best')

        for i in range(self.num_orbits):
            orb = Orbit(self, step=self.rand_idx[i])
            orb.RV -= orb.offset[0] - orb_ml.offset[0]
            ax1.plot(self.epoch_calendar, orb.RV, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
            ax2.plot(self.epoch_calendar, orb.RV - orb_ml.RV, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)

        ax1.plot(self.epoch_calendar, orb_ml.RV, color='black')
        ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))
        
        rv_epoch_list = []
        all_RV_eps = []
        for i in range(self.nInst):
            epoch_obs_Inst = np.zeros(len(self.epoch_obs_dic[i]))
            for j in range(len(self.epoch_obs_dic[i])):
                epoch_obs_Inst[j] = self.JD_to_calendar(self.epoch_obs_dic[i][j])
            rv_epoch_list.append(epoch_obs_Inst)
            all_RV_eps += list(epoch_obs_Inst)
            
        # plot which instrument, self.whichInst=1 means 1st instrument, etc.

        orb_ml_obs = Orbit(self, 'best', epochs='observed')
        
        #if self.whichInst == str('All'):
        #    print('You have chosen to plot RV for all the Instruments')
        all_OC = []
        all_OC_err = []
        
        for i in range(self.nInst):
            plot_this = True
            if not self.whichInst == str('All'):
                plot_this = False
                whichInst = int(self.whichInst)
                if i + 1 == whichInst and i < self.nInst:
                    plot_this = True
            if not plot_this:
                continue
            
            jit_ml = orb_ml.par.return_jitters()
            ax1.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], yerr=np.sqrt(self.RV_obs_err_dic[i]**2 + jit_ml[i]**2), fmt=self.color_list[i]+'o', ecolor='black', capsize=3, alpha = 0.8, zorder=199+i)#, ecolor='black', markersize = 1, elinewidth = 0.3, capsize=1, capthick = 0.3, zorder = 200+i, alpha = 0.8)
            ax1.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], s=45, facecolors='none', edgecolors='k', zorder=200+i, alpha = 0.8)
            
            OC = self.RV_obs_dic[i] + orb_ml_obs.offset[i] - orb_ml_obs.RV[self.RVinst == i]
            y_err = self.RV_obs_err_dic[i]
            all_OC += list(OC)
            all_OC_err += list(y_err)
            
            ax2.errorbar(rv_epoch_list[i], OC, yerr=np.sqrt(y_err**2 + jit_ml[i]**2), fmt=self.color_list[i]+'o', ecolor='black', capsize=3)
            ax2.scatter(rv_epoch_list[i], OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=0.5)

            #else:
            #    print('ValueError: Please enter a valid instrument number between 1 and '+str(self.nInst)+ ' or enter All to plot the observed data points from all Instruments')
            #    raise SystemExit
        
        # axes settings
        # ax1
        ax2.sharex(ax1)
        range_ep_obs = max(all_RV_eps) - min(all_RV_eps)
        RV_obs_max = max([max(self.RV_obs_dic[i] + orb_ml.offset[i]) for i in range(self.nInst)])
        RV_obs_min = min([min(self.RV_obs_dic[i] + orb_ml.offset[i]) for i in range(self.nInst)])
        range_RV_obs = RV_obs_max - RV_obs_min
        
        ax1.set_xlim(min(all_RV_eps) - range_ep_obs/20., max(all_RV_eps) + range_ep_obs/20.)
        ax1.set_ylim(RV_obs_min - range_RV_obs/10., RV_obs_max + range_RV_obs/10.)
        ax1.tick_params(labelbottom=False)   # not NullFormatter: sharex shares the formatter and would blank the O-C axis too
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax1.set_ylabel('RV (m/s)', fontsize=13)
        # ax2

        min_OC = np.amin(np.asarray(all_OC) - np.asarray(all_OC_err))
        max_OC = np.amax(np.asarray(all_OC) + np.asarray(all_OC_err))
        min_OC = min(min_OC, -max_OC)
        max_OC = max(max_OC, -min_OC)
        range_OC = max_OC - min_OC
        
        ax2.set_ylim(min_OC - range_OC/7., max_OC + range_OC/7.)

        if self.set_limit:
            ax2.set_xlim(float(self.user_xlim[0]), float(self.user_xlim[1]))
            ax1.set_ylim(float(self.user_ylim[0]),float(self.user_ylim[1]))
            
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())

        ax2.set_xticks(ax2.get_xticks()[::2])
            
        ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax2.set_xlabel('Epoch (yr)', fontsize=13)
        ax2.set_ylabel('O-C', fontsize=13)
        
        # from advanced plotting settings in config.ini
        if self.show_title:
            ax1.set_title('RV vs. Epoch')
        if self.add_text:
            ax1.text(self.x_text,self.y_text, self.text_name, horizontalalignment='left',verticalalignment='center', transform=ax1.transAxes,fontsize=15)
        if self.usecolorbar:
            cbar_ax = fig.add_axes([1.55, 0.16, 0.05, 0.7])
            cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
            cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
            cbar.ax.get_yaxis().labelpad=20
            fig.delaxes(cbar_ax)

        print("Plotting RV, your plot is generated at " + self.outputdir)
        try:
            self.align_ylabels(fig, [ax1, ax2])
        except Exception:
            pass                            

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            plt.savefig(os.path.join(self.outputdir, 'RV_OC_' + self.title + '_Inst' + str(self.whichInst) +'.pdf'), transparent=True, bbox_inches='tight', dpi=200)
################################################################################################



# 4. relative separation plot

    ################################################################################################
    ##################### plot Relative Separation vs. Epoch and O-C ############
    
    def relsep(self):
    
        if self.have_reldat:
            fig = plt.figure(figsize=(4.7, 5.5))
            ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
            ax2 = fig.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height

            orb_ml = Orbit(self, step='best')
            # plot the randomly selected curves
                        
            for i in range(self.num_orbits):
                orb = Orbit(self, step=self.rand_idx[i]) 
                
                ax1.plot(self.epoch_calendar, orb.relsep, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax2.plot(self.epoch_calendar, orb.relsep - orb_ml.relsep, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)

            ax1.plot(self.epoch_calendar, orb_ml.relsep, color='black')
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

            # plot the observed data points
            orb_ml = Orbit(self, step='best', epochs='observed')

            ep_relAst_obs_calendar = []
            for i in range(len(self.ep_relAst_obs)):
                ep_relAst_obs_calendar.append(self.JD_to_calendar(self.ep_relAst_obs[i]))
            ax1.errorbar(ep_relAst_obs_calendar, self.relsep_obs, yerr=self.relsep_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=4, markersize=5, zorder=299)
            ax1.scatter(ep_relAst_obs_calendar, self.relsep_obs, s=60, facecolors=self.marker_color, edgecolors='k', alpha=1, zorder=300)

            dat_OC = self.relsep_obs - orb_ml.relsep[self.ast_indx]
            
            ax2.errorbar(ep_relAst_obs_calendar, dat_OC, yerr=self.relsep_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=4, markersize=5, zorder=299)
            ax2.scatter(ep_relAst_obs_calendar, dat_OC, s=60, facecolors=self.marker_color, edgecolors='k', alpha=1, zorder=300)
                
            # axes settings
            # ax1
            ax2.sharex(ax1)
            range_eprA_obs = max(ep_relAst_obs_calendar) - min(ep_relAst_obs_calendar)
            range_relsep_obs = max(self.relsep_obs)  - min(self.relsep_obs)
            ax1.set_xlim(min(ep_relAst_obs_calendar) - range_eprA_obs/8., max(ep_relAst_obs_calendar) + range_eprA_obs/8.)
            ax1.set_ylim(min(self.relsep_obs) - range_relsep_obs/2., max(self.relsep_obs) + range_relsep_obs/2.)
            ax1.tick_params(labelbottom=False)   # not NullFormatter: sharex shares the formatter and would blank the O-C axis too
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            lab1 = ax1.set_ylabel('Separation (arcsec)', fontsize=13)
            # ax2
            range_datOC = max(dat_OC) - min(dat_OC)
            if np.abs(min(dat_OC)) <= np.abs(max(dat_OC)):
                ax2.set_ylim(-np.abs(max(dat_OC)) - range_datOC, max(dat_OC) + range_datOC)
            elif np.abs(min(dat_OC)) > np.abs(max(dat_OC)):
                ax2.set_ylim(min(dat_OC) - range_datOC, np.abs(min(dat_OC)) + range_datOC)
                
            if self.set_limit:
                ax2.set_xlim(float(self.user_xlim[0]), float(self.user_xlim[1]))
                ax1.set_ylim(float(self.user_ylim[0]),float(self.user_ylim[1]))
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            self._year_ticks(ax2)
            ax2.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
            lab2 = ax2.set_ylabel('O-C', fontsize=13)
            #ax2.get_yaxis().set_label_coords()

            # from advanced plotting settings in config.ini
            if self.show_title:
                ax1.set_title('Relative Separation vs. Epoch')
            if self.add_text:
                ax1.text(self.x_text,self.y_text, self.text_name, horizontalalignment='left',verticalalignment='center', transform=ax1.transAxes,fontsize=15)
            if self.usecolorbar:
                cbar_ax = fig.add_axes([1.55, 0.16, 0.05, 0.7])
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                cbar.ax.get_yaxis().labelpad=20
                fig.delaxes(cbar_ax)

        else:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_axes((0.15, 0.1, 0.8, 0.8))
            
            # plot the 50 randomly selected curves
            
            for i in range(self.num_orbits):
                orb = Orbit(self, step=self.rand_idx[i]) 

                ax.plot(self.epoch_calendar, orb.relsep, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                
            # plot the most likely one
            orb_ml = Orbit(self, step='best') 
            ax.plot(self.epoch_calendar, orb_ml.relsep, color='black')
                
            # axes settings
            ax.set_xlim(self.start_epoch, self.end_epoch)
            self._year_ticks(ax)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax.set_ylabel('Separation (arcsec)', fontsize=13)
            ax.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)    
            plt.tight_layout()
            print("Plotting Separation, your plot is generated at " + self.outputdir)
            try:
                self.align_ylabels(fig, [ax1, ax2])
            except Exception:
                pass
            plt.savefig(os.path.join(self.outputdir, 'relsep_OC_' + self.title)+'.pdf', transparent=True, bbox_inches='tight', dpi=200)
################################################################################################



# 5. position angle plot

    ################################################################################################
    ############### plot the position angle and O-C ###############

    def PA(self):
    
        if self.have_reldat:
            fig = plt.figure(figsize=(4.7, 5.5))
            ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
            ax2 = fig.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height

            orb_ml = Orbit(self, step='best')

            for i in range(self.num_orbits):

                orb = Orbit(self, step=self.rand_idx[i])

                ax1.plot(self.epoch_calendar, orb.PA, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax2.plot(self.epoch_calendar, orb.PA - orb_ml.PA, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)

            # plot the highest likelihood orbit
            ax1.plot(self.epoch_calendar, orb_ml.PA, color='black')
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

            # plot the observed data points
            orb_ml = Orbit(self, 'best', epochs='observed')

            ep_relAst_obs_calendar = []
            for i in range(len(self.ep_relAst_obs)):
                ep_relAst_obs_calendar.append(self.JD_to_calendar(self.ep_relAst_obs[i]))
            ax1.errorbar(ep_relAst_obs_calendar, self.PA_obs, yerr=self.PA_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=4, markersize=5, zorder=299)
            ax1.scatter(ep_relAst_obs_calendar, self.PA_obs, s=60, facecolors=self.marker_color, edgecolors='k', alpha=1, zorder=300)
            dat_OC = self.PA_obs - orb_ml.PA[self.ast_indx]

            ax2.errorbar(ep_relAst_obs_calendar, dat_OC, yerr=self.PA_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=4, markersize=5, zorder=299)
            ax2.scatter(ep_relAst_obs_calendar, dat_OC, s=60, facecolors=self.marker_color, edgecolors='k', alpha=1, zorder=300)

            # axes settings
            # ax1
            ax2.sharex(ax1)
            range_eprA_obs = max(ep_relAst_obs_calendar) - min(ep_relAst_obs_calendar)
            range_PA_obs = max(self.PA_obs)  - min(self.PA_obs)
            ax1.set_xlim(min(ep_relAst_obs_calendar) - range_eprA_obs/8., max(ep_relAst_obs_calendar) + range_eprA_obs/8.)
            ax1.set_ylim(min(self.PA_obs) - range_PA_obs/5., max(self.PA_obs) + range_PA_obs/5.)
            ax1.tick_params(labelbottom=False)   # not NullFormatter: sharex shares the formatter and would blank the O-C axis too
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_ylabel(r'Position Angle ($^{\circ}$)', labelpad=5, fontsize=13)
            # ax2
            range_datOC = max(dat_OC) - min(dat_OC)
            if np.abs(min(dat_OC)) <= np.abs(max(dat_OC)):
                ax2.set_ylim(-np.abs(max(dat_OC)) - range_datOC, max(dat_OC) + range_datOC)
            elif np.abs(min(dat_OC)) > np.abs(max(dat_OC)):
                ax2.set_ylim(min(dat_OC) - range_datOC, np.abs(min(dat_OC)) + range_datOC)
                
            if self.set_limit:
                ax2.set_xlim(float(self.user_xlim[0]), float(self.user_xlim[1]))
                ax1.set_ylim(float(self.user_ylim[0]),float(self.user_ylim[1]))
                
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            self._year_ticks(ax2)
            ax2.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
            ax2.set_ylabel('O-C', labelpad=6, fontsize=13)

            # from advanced plotting settings in config.ini
            if self.show_title:
                ax1.set_title('Position angle vs. Epoch')
            if self.add_text:
                ax1.text(self.x_text,self.y_text, self.text_name, horizontalalignment='left',verticalalignment='center', transform=ax1.transAxes,fontsize=15)
            if self.usecolorbar:
                cbar_ax = fig.add_axes([1.55, 0.16, 0.05, 0.7])
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                cbar.ax.get_yaxis().labelpad=20
                fig.delaxes(cbar_ax)

        else:
            fig = plt.figure(figsize=(5, 5.5))
            ax = fig.add_axes((0.15, 0.1, 0.8, 0.8))

            # plot the 50 randomly selected curves
            for i in range(self.num_orbits):
                orb = Orbit(self, step=self.rand_idx[i]) 

                ax.plot(self.epoch_calendar, orb.PA, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)

            orb_ml = Orbit(self, 'best')
            # plot the most likely one
            ax.plot(self.epoch_calendar, orb_ml.PA, color='black')

            ax.set_xlim(self.start_epoch, self.end_epoch)
            self._year_ticks(ax)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax.set_title(self.title)
            ax.set_ylabel(r'Position Angle ($^{\circ}$)', fontsize=13)
            ax.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            plt.tight_layout()
            print("Plotting Position Angle, your plot is generated at " + self.outputdir)
            try:
                self.align_ylabels(fig, [ax1, ax2])
            except Exception:
                pass
            plt.savefig(os.path.join(self.outputdir,'PA_OC_' + self.title)+'.pdf',bbox_inches='tight', dpi=200, transparent=True)
#############################################################################with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            ###################




# 6. Proper motion plots

    ################################################################################################
    ############### plot the proper motions, together or separately ###############
    
    def _sky_path_offsets(self, step, epochs_decimalyear):
        """
        Reflex-motion offsets (dRA*, dDec) of the primary at the given epochs,
        for one step of the chain, in orvara's internal units (multiply by the
        parallax in arcsec x 1000 for mas).
        """
        epochs = np.array([self.calendar_to_JD(ep) for ep in epochs_decimalyear])
        data = orbit.Data(self.Hip, self.HGCAFile, self.RVfile, self.relAstfile,
                          verbose=False)
        data.custom_epochs(epochs, iplanet=self.iplanet)
        model = orbit.Model(data)
        jit_per_inst = 'jitter0' in [col.name for col in self.chain.columns]
        nj = data.nInst if jit_per_inst else 1
        for i in range(self.nplanets):
            par = orbit.Params(pull_chain_params(self.chain.columns, step),
                               iplanet=i, nplanets=self.nplanets,
                               ninst_RV=data.nInst, ninst_jit=nj)
            orbit.calc_EA_RPP(data, par, model)
            # calc_offsets is what actually writes dRA_G / dDec_G.  Calling
            # calc_EA_RPP alone leaves them at the zeros set in Model.__init__,
            # which made the model acceleration identically zero and the O-C
            # equal to the full measured acceleration.
            orbit.calc_offsets(data, par, model, i)
            par.free()
        dRA_G, dDec_G, _, _, _, _ = model.return_dRA_dDec()
        model.free()
        return np.asarray(dRA_G, dtype=float), np.asarray(dDec_G, dtype=float)

    def model_accel_jerk(self, step, dt_gaia=None):
        """
        Model astrometric acceleration -- and, for a nine-parameter source, jerk
        -- at the Gaia central epoch, for one step of the chain.

        The sky path is sampled at the same epochs the likelihood uses in
        ``orvara.orbit.calc_PMs_no_epoch_astrometry``, and reduced with the same
        finite-difference coefficients, so this panel and chi^2_G describe the
        same quantity:

        * five- or seven-parameter source: three points at -T/2, 0, +T/2,
          acceleration by the quadratic difference, no jerk;
        * nine-parameter source: four points at -T/2, -T/6, +T/6, +T/2,
          acceleration and jerk by the cubic differences.

        :param step: index into the chain.
        :param dt_gaia: Gaia window length in years; defaults to the value
            identified from the HGCA file (the NSS model's 2 x 517.5 d for a
            one-star NSS file).
        :returns: (accra, accdec, jerkra, jerkdec); accelerations in mas/yr^2,
            jerks in mas/yr^3, zero when the source has no jerk term.
        """
        if dt_gaia is None:
            dt_gaia = getattr(self, 'dt_gaia', 2.76)
        h = dt_gaia / 2.0
        npar = getattr(self, 'gaia_npar', 5)
        if npar == 9:
            sample_offsets = np.array([-h, -h / 3.0, h / 3.0, h])
        else:
            sample_offsets = np.array([-h, 0.0, h])
        ep_dec = getattr(self, 'ep_gaia_acc_dec', self.ep_gaia_acc)

        # RA and Dec have their own central epochs in the HGCA convention, so
        # the path is evaluated once per coordinate, as the likelihood does.
        dRA_G = self._sky_path_offsets(step, self.ep_gaia_acc + sample_offsets)[0]
        dDec_G = self._sky_path_offsets(step, ep_dec + sample_offsets)[1]
        scale = 1e3 * self.chain['plx_ML'][step]        # AU/plx units -> mas
        dRA_G = scale * np.asarray(dRA_G, dtype=float)
        dDec_G = scale * np.asarray(dDec_G, dtype=float)

        if npar == 9:
            f0, f1, f2, f3 = dRA_G[:4]
            g0, g1, g2, g3 = dDec_G[:4]
            accra = (f0 - f1 - f2 + f3) * 9.0 / (8.0 * h * h)
            accdec = (g0 - g1 - g2 + g3) * 9.0 / (8.0 * h * h)
            jerkra = (-f0 + 3.0 * f1 - 3.0 * f2 + f3) * 27.0 / (8.0 * h ** 3)
            jerkdec = (-g0 + 3.0 * g1 - 3.0 * g2 + g3) * 27.0 / (8.0 * h ** 3)
        else:
            f0, f1, f2 = dRA_G[:3]
            g0, g1, g2 = dDec_G[:3]
            accra = (f0 + f2 - 2.0 * f1) / (h * h)
            accdec = (g0 + g2 - 2.0 * g1) / (h * h)
            jerkra = jerkdec = 0.0
        return accra, accdec, jerkra, jerkdec

    def model_window_mean_pm(self, step='best', n=401):
        """
        Model proper motion averaged over the Gaia window, per coordinate.

        A Gaia DR3 NSS proper motion is the mean over 2016.0 +/- 517.5 d, and the
        likelihood compares it with (f(+T/2) - f(-T/2)) / T of the model sky
        path (see ``orvara.orbit.calc_PMs_no_epoch_astrometry``), i.e. with the
        time-average of mu(t), not its value at 2016.0.  For an Acceleration9
        source the two differ by jerk (T/2)^2 / 6 (4.3 mas/yr for J1717 in Dec).

        :returns: (mean mu_RA*, mean mu_Dec) in mas/yr.
        """
        h = 0.5 * getattr(self, 'dt_gaia', 2.76)
        out = []
        for ep_c, use_ra in ((self.ep_gaia_acc, True),
                             (getattr(self, 'ep_gaia_acc_dec', self.ep_gaia_acc), False)):
            t_cal = ep_c + np.linspace(-h, h, n)
            jd = np.asarray(self.calendar_to_JD(t_cal), dtype=float)
            orb = Orbit(self, step, epochs=jd)
            mu = np.asarray(orb.mu_RA if use_ra else orb.mu_Dec, dtype=float)
            # trapezoid mean on the uniform grid
            out.append((0.5 * mu[0] + mu[1:-1].sum() + 0.5 * mu[-1]) / (n - 1))
        return out[0], out[1]

    def proper_motions(self):
        if self.have_pmdat:
            if self.pm_separate:
                fig = plt.figure(figsize=(4.7, 8.0))
                ax1 = fig.add_axes((0.15, 0.55, 0.80, 0.38))
                ax2 = fig.add_axes((0.15, 0.39, 0.80, 0.14))
                ax_a1 = fig.add_axes((0.15, 0.08, 0.80, 0.27))
                fig1 = plt.figure(figsize=(4.7, 8.0))
                ax3 = fig1.add_axes((0.15, 0.55, 0.80, 0.38))
                ax4 = fig1.add_axes((0.15, 0.39, 0.80, 0.14))
                ax_a3 = fig1.add_axes((0.15, 0.08, 0.80, 0.27))
            else:
                fig = plt.figure(figsize=(13, 8.0))
                ax1 = fig.add_axes((0.08, 0.55, 0.37, 0.38))
                ax2 = fig.add_axes((0.08, 0.39, 0.37, 0.14))
                ax_a1 = fig.add_axes((0.08, 0.08, 0.37, 0.27))
                ax3 = fig.add_axes((0.55, 0.55, 0.37, 0.38))
                ax4 = fig.add_axes((0.55, 0.39, 0.37, 0.14))
                ax_a3 = fig.add_axes((0.55, 0.08, 0.37, 0.27))

            # ── epoch arrays for observed data points ─────────────────────
            ep_mualp_obs_calendar = []
            ep_mudec_obs_calendar = []
            t_RA = np.asarray(self.ep_mualp_obs).flatten()
            t_Dec = np.asarray(self.ep_mudec_obs).flatten()
            orb_ml_RA = Orbit(self, 'best', epochs=t_RA)
            orb_ml_Dec = Orbit(self, 'best', epochs=t_Dec)
            for i in range(len(self.ep_mualp_obs)):
                ep_mualp_obs_calendar.append(self.JD_to_calendar(self.ep_mualp_obs[i]))
            for i in range(len(self.ep_mudec_obs)):
                ep_mudec_obs_calendar.append(self.JD_to_calendar(self.ep_mudec_obs[i]))

            t1_RA = min(ep_mualp_obs_calendar); t2_RA = max(ep_mualp_obs_calendar)
            dt_RA = t2_RA - t1_RA
            t1_Dec = min(ep_mudec_obs_calendar); t2_Dec = max(ep_mudec_obs_calendar)
            dt_Dec = t2_Dec - t1_Dec
            # A Gaia-only source has a single proper-motion epoch, so span the
            # Gaia window (plus margin) rather than a zero-width range.
            if dt_RA == 0:
                half = getattr(self, 'dt_gaia', 2.76)
                t1_RA, t2_RA, dt_RA = t1_RA - half, t2_RA + half, 2 * half
            if dt_Dec == 0:
                half = getattr(self, 'dt_gaia', 2.76)
                t1_Dec, t2_Dec, dt_Dec = t1_Dec - half, t2_Dec + half, 2 * half
            xl_RA = [t1_RA - dt_RA /8., t2_RA + dt_RA /8.]
            xl_Dec = [t1_Dec - dt_Dec/8., t2_Dec + dt_Dec/8.]
            for ax in [ax1, ax2, ax_a1]:
                ax.set_xlim(xl_RA)
            for ax in [ax3, ax4, ax_a3]:
                ax.set_xlim(xl_Dec)

            orb_ml = Orbit(self, 'best')
            minT = min(t1_RA - dt_RA/5., t1_Dec - dt_Dec/5.)
            maxT = max(t2_RA + dt_RA/5., t2_Dec + dt_Dec/5.)
            T = np.asarray(self.epoch_calendar)
            indx = (T > minT) * (T < maxT) > 0

            # ── random orbit samples ──────────────────────────────────────
            sampled_orbits = []
            for i in range(self.num_orbits):
                orb = Orbit(self, step=self.rand_idx[i])
                sampled_orbits.append(orb)
                c = self.colormap(self.normalize(orb.colorpar))
                ax1.plot(T[indx], orb.mu_RA[indx], color=c, alpha=0.3)
                ax3.plot(T[indx], orb.mu_Dec[indx], color=c, alpha=0.3)
                ax2.plot(T[indx], orb.mu_RA[indx] - orb_ml.mu_RA[indx], color=c, alpha=0.3)
                ax4.plot(T[indx], orb.mu_Dec[indx] - orb_ml.mu_Dec[indx], color=c, alpha=0.3)

            # ── best-fit curves ───────────────────────────────────────────
            ax1.plot(T[indx], orb_ml.mu_RA[indx], color='black', label=r'Model $\mu(t)$')
            ax3.plot(T[indx], orb_ml.mu_Dec[indx], color='black', label=r'Model $\mu(t)$')

            # Name the proper-motion points after what they are: a Gaia-only NSS
            # source has just the Gaia window-mean proper motion.
            if getattr(self, 'have_nss_acc', False) and getattr(self, 'n_pm_epochs', 2) < 2:
                pm_data_label = 'Gaia NSS (Acc%d)' % self.gaia_npar
            elif getattr(self, 'have_nss_acc', False):
                pm_data_label = 'Hipparcos / Gaia NSS (Acc%d)' % self.gaia_npar
            else:
                pm_data_label = 'Hipparcos / Gaia'
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch_calendar)), 'k--', dashes=(5, 5))
            ax4.plot(self.epoch_calendar, np.zeros(len(self.epoch_calendar)), 'k--', dashes=(5, 5))

            # ── HGCA proper motion data points ────────────────────────────
            ax1.errorbar(ep_mualp_obs_calendar, self.mualp_obs,
                         yerr=self.mualp_obs_err, color=self.marker_color,
                         fmt='o', ecolor='black', capsize=5, markersize=1, zorder=299)
            ax1.scatter(ep_mualp_obs_calendar, self.mualp_obs, s=80,
                        facecolors=self.marker_color, edgecolors='k', zorder=300,
                        label=pm_data_label)
            ax3.errorbar(ep_mudec_obs_calendar, self.mudec_obs,
                         yerr=self.mudec_obs_err, color=self.marker_color,
                         fmt='o', ecolor='black', capsize=5, markersize=1, zorder=299)
            ax3.scatter(ep_mudec_obs_calendar, self.mudec_obs, s=80,
                        facecolors=self.marker_color, edgecolors='k', zorder=300,
                        label=pm_data_label)

            # ── PM O-C ────────────────────────────────────────────────────
            dat_OC_RA = np.array(self.mualp_obs - orb_ml_RA.mu_RA, dtype=float)
            dat_OC_Dec = np.array(self.mudec_obs - orb_ml_Dec.mu_Dec, dtype=float)
            if self.have_nss_acc:
                # The Gaia NSS proper motion is a window mean, so compare it with
                # the model's window mean (what chi^2_G uses), not with mu(2016.0).
                mean_RA, mean_Dec = self.model_window_mean_pm('best')
                h_win = 0.5 * getattr(self, 'dt_gaia', 2.76)
                ep_dec_c = getattr(self, 'ep_gaia_acc_dec', self.ep_gaia_acc)
                iG_RA = int(np.argmin(np.abs(np.asarray(ep_mualp_obs_calendar) - self.ep_gaia_acc)))
                iG_Dec = int(np.argmin(np.abs(np.asarray(ep_mudec_obs_calendar) - ep_dec_c)))
                dat_OC_RA[iG_RA] = self.mualp_obs[iG_RA] - mean_RA
                dat_OC_Dec[iG_Dec] = self.mudec_obs[iG_Dec] - mean_Dec
                ax1.hlines(mean_RA, self.ep_gaia_acc - h_win, self.ep_gaia_acc + h_win,
                           colors='k', linestyles=':', lw=1.8, zorder=298,
                           label='Model, mean over Gaia window')
                ax3.hlines(mean_Dec, ep_dec_c - h_win, ep_dec_c + h_win,
                           colors='k', linestyles=':', lw=1.8, zorder=298,
                           label='Model, mean over Gaia window')
                # model entries first, data last, as in the acceleration panel
                for ax_pm in (ax1, ax3):
                    h, l = ax_pm.get_legend_handles_labels()
                    order = sorted(range(len(l)), key=lambda k: l[k] == pm_data_label)
                    ax_pm.legend([h[k] for k in order], [l[k] for k in order],
                                 fontsize=8, loc='best', framealpha=0.8)
            ax2.errorbar(ep_mualp_obs_calendar, dat_OC_RA,
                         yerr=self.mualp_obs_err, color=self.marker_color,
                         fmt='o', ecolor='black', capsize=5, markersize=1, zorder=299)
            ax2.scatter(ep_mualp_obs_calendar, dat_OC_RA, s=80,
                        facecolors=self.marker_color, edgecolors='k', zorder=300)
            ax4.errorbar(ep_mudec_obs_calendar, dat_OC_Dec,
                         yerr=self.mudec_obs_err, color=self.marker_color,
                         fmt='o', ecolor='black', capsize=5, markersize=1, zorder=299)
            ax4.scatter(ep_mudec_obs_calendar, dat_OC_Dec, s=80,
                        facecolors=self.marker_color, edgecolors='k', zorder=300)

            # ── NSS acceleration panel ────────────────────────────────────
            # The model acceleration is drawn as a function of time: the time
            # derivative of the proper-motion curves in the panels above
            # (black = maximum-likelihood orbit, coloured = posterior draws).
            # The O-C below compares the DR3 *coefficient* with the model
            # coefficient obtained by the same finite-difference reduction the
            # likelihood uses (OrbitPlots.model_accel_jerk), so it is the
            # quantity that enters chi^2_G.
            acc_ep_cal = self.ep_gaia_acc   # decimal year
            Tc = T[indx]

            def _accel_curve(mu):
                """d(mu)/dt on the plotted epoch grid, mas/yr -> mas/yr^2."""
                return np.gradient(np.asarray(mu, dtype=float)[indx], Tc)

            for orb in sampled_orbits:
                c = self.colormap(self.normalize(orb.colorpar))
                ax_a1.plot(Tc, _accel_curve(orb.mu_RA), color=c, alpha=0.3)
                ax_a3.plot(Tc, _accel_curve(orb.mu_Dec), color=c, alpha=0.3)
            ml_curve_RA, ml_curve_Dec = _accel_curve(orb_ml.mu_RA), _accel_curve(orb_ml.mu_Dec)
            ax_a1.plot(Tc, ml_curve_RA, color='black',
                       label=r'Model $\mathrm{d}\mu/\mathrm{d}t$')
            ax_a3.plot(Tc, ml_curve_Dec, color='black',
                       label=r'Model $\mathrm{d}\mu/\mathrm{d}t$')

            # ML coefficient comparison (OrbitPlots.model_accel_jerk) is not drawn:
            # the red segment's slope is the measured jerk, and the numbers belong in
            # the beststep / posterior output rather than on the figure.

            nss_color = '#d62728'   # red: distinct from the viridis posterior draws
            ENVELOPE_LABEL = r'Gaia NSS 1$\sigma$ envelope'
            envelope_handles = {}
            if self.have_nss_acc:
                half_window = 0.5 * getattr(self, 'dt_gaia', 2.76)   # yr
                t_seg = acc_ep_cal + np.linspace(-half_window, half_window, 50)
                dt_seg = t_seg - acc_ep_cal
                for ax_a, val, err, jval, jerr in [
                    (ax_a1, self.nss_accra,  self.nss_accra_err,
                     self.nss_jerkra,  self.nss_jerkra_err),
                    (ax_a3, self.nss_accdec, self.nss_accdec_err,
                     self.nss_jerkdec, self.nss_jerkdec_err)]:
                    # central value at the DR3 reference epoch
                    ax_a.errorbar([acc_ep_cal], [val], yerr=[err],
                                  color=nss_color, fmt='D', ecolor='black',
                                  markeredgecolor='black',
                                  capsize=5, markersize=7, zorder=301,
                                  label='Gaia NSS (Acc%d)' % self.gaia_npar)
                    if self.have_nss_jerk and jval is not None:
                        # nine-parameter solution: the DR3 model acceleration is
                        # a + j (t - t0) across the window; draw that segment
                        # with its 1-sigma band.
                        seg = val + jval * dt_seg
                        seg_err = np.sqrt(err**2 + (jerr * dt_seg)**2)
                        ax_a.plot(t_seg, seg, color=nss_color, lw=2.2, zorder=300,
                                  label=r'Gaia NSS $a + j\,(t-t_0)$')
                        env_fill = ax_a.fill_between(t_seg, seg - seg_err, seg + seg_err,
                                                     color=nss_color, alpha=0.18, lw=0, zorder=298,
                                                     label=ENVELOPE_LABEL)
                        env_edge, = ax_a.plot(t_seg, seg - seg_err, color=nss_color, lw=1.3,
                                              ls='--', zorder=299)
                        ax_a.plot(t_seg, seg + seg_err, color=nss_color, lw=1.3, ls='--', zorder=299)
                        # legend entry = band + dashed edge, so the dashed lines are labelled
                        envelope_handles[ax_a] = (env_fill, env_edge)
                    else:
                        # seven-parameter solution: constant over the window
                        ax_a.hlines(val, t_seg[0], t_seg[-1], color=nss_color,
                                    lw=2.2, zorder=300,
                                    label=r'Gaia NSS $a$ (Acc7, constant)')

                # y-range from the maximum-likelihood curve and the measurement,
                # so a few wide posterior draws cannot flatten the panel
                for ax_a, curve, val, err in [
                    (ax_a1, ml_curve_RA,  self.nss_accra,  self.nss_accra_err),
                    (ax_a3, ml_curve_Dec, self.nss_accdec, self.nss_accdec_err)]:
                    lo = min(np.min(curve), val - 2 * err)
                    hi = max(np.max(curve), val + 2 * err)
                    pad = 0.15 * (hi - lo) if hi > lo else 1.0
                    ax_a.set_ylim(lo - pad, hi + pad)

            # ── axis labels and formatting ─────────────────────────────────
            if self.set_limit:
                for ax in [ax1, ax3]:
                    ax.set_xlim(float(self.user_xlim[0]), float(self.user_xlim[1]))
                    ax.set_ylim(float(self.user_ylim[0]), float(self.user_ylim[1]))

            ax2.sharex(ax1)
            ax_a1.sharex(ax1)
            ax4.sharex(ax3)
            ax_a3.sharex(ax3)
            for ax in [ax2, ax_a1]:
                ax.set_xlim(xl_RA)
            for ax in [ax4, ax_a3]:
                ax.set_xlim(xl_Dec)

            ax1.set_ylabel(r'$\mu_{\alpha*}$ (mas yr$^{-1}$)', labelpad=6, fontsize=13)
            ax3.set_ylabel(r'$\mu_{\delta}$ (mas yr$^{-1}$)', labelpad=6, fontsize=13)
            ax2.set_ylabel('O-C (PM)', labelpad=6, fontsize=10)
            ax4.set_ylabel('O-C (PM)', labelpad=6, fontsize=10)
            ax_a1.set_ylabel(r'$\alpha_{\alpha*}$ (mas yr$^{-2}$)', labelpad=6, fontsize=11)
            ax_a3.set_ylabel(r'$\alpha_{\delta}$ (mas yr$^{-2}$)', labelpad=6, fontsize=11)

            for ax in [ax1, ax3, ax2, ax4]:
                ax.tick_params(labelbottom=False)   # see RV(): a shared NullFormatter would also blank the bottom panel
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both',
                               left=True, right=True, bottom=True, top=True)
            for ax in [ax_a1, ax_a3]:
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both',
                               left=True, right=True, bottom=True, top=True)
                ax.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)

            y1_RA = max(np.abs(dat_OC_RA) + max(self.mualp_obs_err))
            y1_Dec = max(np.abs(dat_OC_Dec) + max(self.mudec_obs_err))
            ax2.set_ylim(-1.2*y1_RA,  1.2*y1_RA)
            ax4.set_ylim(-1.2*y1_Dec, 1.2*y1_Dec)

            # Gaia-only NSS source (Acc7 or Acc9, no Hipparcos): the only proper
            # motion is the Gaia one, which the free barycentric proper motion
            # absorbs, so its O-C is zero by construction.  Drop the O-C panel
            # and let the proper-motion panel take its space; the orbit is
            # tested by the acceleration (and jerk) panel below.
            if getattr(self, 'n_pm_epochs', 2) < 2:
                for top, oc in ((ax1, ax2), (ax3, ax4)):
                    b_top, b_oc = top.get_position(), oc.get_position()
                    top.set_position([b_top.x0, b_oc.y0, b_top.width, b_top.y1 - b_oc.y0])
                    oc.remove()

            if self.have_nss_acc:
                # 'best' keeps the legend off the NSS segment, which a fixed
                # upper-left legend covers when the acceleration is near its maximum
                for ax_a in (ax_a1, ax_a3):
                    handles, labels = ax_a.get_legend_handles_labels()
                    if ax_a in envelope_handles:
                        handles = [envelope_handles[ax_a] if lab == ENVELOPE_LABEL else h
                                   for h, lab in zip(handles, labels)]
                    ax_a.legend(handles, labels, fontsize=8, loc='best', framealpha=0.8)
            if self.show_title:
                ax1.set_title('Right Ascension vs. Epoch')
                ax3.set_title('Declination vs. Epoch')
            if self.add_text:
                for ax in [ax1, ax3]:
                    ax.text(self.x_text, self.y_text, self.text_name,
                            horizontalalignment='left', verticalalignment='center',
                            transform=ax.transAxes, fontsize=15)

            if self.usecolorbar:
                if self.pm_separate:
                    for figure, axes_list, name in [
                        (fig,  [ax1, ax2, ax_a1], 'RA_'),
                        (fig1, [ax3, ax4, ax_a3], 'Dec_')]:
                        cbar_ax = figure.add_axes([1.55, 0.16, 0.05, 0.7])
                        cbar = figure.colorbar(self.sm, ax=cbar_ax, fraction=12)
                        cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref],
                                           rotation=270, fontsize=13)
                        cbar.ax.get_yaxis().labelpad = 20
                        figure.delaxes(cbar_ax)
                        try:
                            self.align_ylabels(figure, [a for a in axes_list if a in figure.axes], pad=0.02)
                        except Exception:
                            pass
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=UserWarning)
                            figure.savefig(
                                os.path.join(self.outputdir,
                                             'ProperMotions_' + name + self.title) + '.pdf',
                                transparent=True, bbox_inches='tight', dpi=200)
                else:
                    cbar_ax = fig.add_axes([1.46, 0.16, 0.05, 0.7])
                    cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                    cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref],
                                       rotation=270, fontsize=13)
                    cbar.ax.get_yaxis().labelpad = 20
                    fig.delaxes(cbar_ax)
                    try:
                        self.align_ylabels(fig, [a for a in (ax1, ax2, ax_a1) if a in fig.axes], pad=0.02)
                        self.align_ylabels(fig, [a for a in (ax3, ax4, ax_a3) if a in fig.axes], pad=0.02)
                    except Exception:
                        pass
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        fig.savefig(
                            os.path.join(self.outputdir,
                                         'Proper_Motions_' + self.title) + '.pdf',
                            transparent=True, bbox_inches='tight', dpi=200)
        else:
            fig = plt.figure(figsize=(11, 5.5))
            ax1 = fig.add_axes((0.10, 0.1, 0.35, 0.77))
            ax2 = fig.add_axes((0.60, 0.1, 0.35, 0.77))
            orb_ml = Orbit(self, 'best')
            for i in range(self.num_orbits):
                orb = Orbit(self, step=self.rand_idx[i])
                ax1.plot(self.epoch_calendar, orb.mu_RA,
                         color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax2.plot(self.epoch_calendar, orb.mu_Dec,
                         color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
            ax1.plot(self.epoch_calendar, orb_ml.mu_RA, color='black')
            ax2.plot(self.epoch_calendar, orb_ml.mu_Dec, color='black')
            for ax, ylabel in [(ax1, r'$\mu_{\alpha*}$ (mas yr$^{-1}$)'),
                               (ax2, r'$\mu_{\delta}$ (mas yr$^{-1}$)')]:
                ax.set_xlim(self.start_epoch, self.end_epoch)
                self._year_ticks(ax)
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both',
                               left=True, right=True, bottom=True, top=True)
                ax.set_xlabel('date (yr)')
                ax.set_ylabel(ylabel, labelpad=6)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            plt.tight_layout()
            print("Plotting Proper Motions, your plot is generated at " + self.outputdir)
            if self.pm_separate and not self.usecolorbar:
                fig.savefig(os.path.join(self.outputdir,
                                         'ProperMotions_RA_'  + self.title) + '.pdf',
                            transparent=True, bbox_inches='tight', dpi=200)
                fig1.savefig(os.path.join(self.outputdir,
                                          'ProperMotions_Dec_' + self.title) + '.pdf',
                             transparent=True, bbox_inches='tight', dpi=200)
            elif not self.usecolorbar:
                fig.savefig(os.path.join(self.outputdir,
                                         'Proper_Motions_' + self.title) + '.pdf',
                            transparent=True, bbox_inches='tight', dpi=200)

################################################################################################

#astrometric prediction plot


    def astrometric_prediction_plot(self, nbins=300, smooth=True):
        # NOTE ONLY WORKS IN DECIMAL YEAR
        if self.num_orbits < 500:
            print('WARNING: You will want to use more than {self.num_orbits} orbits to predict the '
                  'position. Ideally you should use 500 or more. Set num_orbits=500 or more in the config file.')
        fig, ax = plt.subplots(figsize=(5, 5))
        JDepochs = self.calendar_to_JD(np.array([self.predicted_ep_ast]).flatten())
        ra, dec = self.astrometric_prediction(JDepochs)
        ra, dec = ra.flatten(), dec.flatten()
        print(f'The mean RA and Dec position of the companion at {self.predicted_ep_ast} is:')
        ra_err, dec_err = np.std(ra), np.std(dec)
        mean_ra, mean_dec = np.mean(ra), np.mean(dec)
        ra_range, dec_range = (mean_ra-5*ra_err, mean_ra+5*ra_err), (mean_dec-5*dec_err,mean_dec+5*dec_err)
        print(f'RA: {round(mean_ra, 4)} +- {round(ra_err, 4)} mas, 'f'Dec: {round(mean_dec, 4)} +- {round(dec_err, 4)} mas')
        k = kde.gaussian_kde([ra, dec])
        xi, yi = np.mgrid[min(ra_range):max(ra_range):nbins * 1j, min(dec_range):max(dec_range):nbins * 1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()])).reshape(xi.shape)
        if smooth:
            sigma = ra_err/(np.ptp(ra_range)/nbins) # size of the ra error on the position, in terms of bins in the KDE.
            sigma /= 3  # we want the kernal of the smoothing window to be a little smaller than our "beam" size.
            zi = ndimage.gaussian_filter(zi, sigma=sigma)
        # solving for the values that correspond to 1, 2 and 3 sigma contours.
        dens_sorted = np.sort(zi.flatten())
        cdf = np.cumsum(dens_sorted)/np.sum(dens_sorted)
        cdf_func = interp1d(cdf, dens_sorted)
        levels = [cdf_func(p) for p in [1 - 0.9973, 1 - 0.954, 1 - 0.683]]
        # plotting the contours and density plot.
        ax.contour(xi, yi, zi, levels=levels, colors=['k', 'C0', 'b'])
        # antialiased=False removes the faint white seams between filled levels.
        # (ContourSet.collections was removed in Matplotlib 3.10; before that
        # the seams were hidden by setting each collection's edgecolor.)
        cnt = ax.contourf(xi, yi, zi, levels=100, cmap=cm.hot_r, antialiased=False)
        if hasattr(cnt, 'collections'):          # Matplotlib < 3.10
            for c in cnt.collections:
                c.set_edgecolor("face")
        #

        # Mark the central star if (0, 0) is within the axis limits
        if ra.min() * ra.max() < 0 and dec.min() * dec.max() < 0:
            ax.plot(0, 0, marker='*', markersize=15, color='c')

        ax.set_aspect('equal')  # change me to 'auto' if the plot is too squeezed.
        ax.set_xlabel(r'$\mathrm{\Delta \alpha}$ (mas)', fontsize=14)
        ax.set_ylabel(r'$\mathrm{\Delta \delta}$ (mas)', fontsize=14)
        ax.annotate('Location of %s, %.2f' % (self.text_name, self.predicted_ep_ast), fontsize=14, xy=(0.1, 0.9), xycoords='axes fraction')
        ax.invert_xaxis()
        plt.tight_layout()
        plt.savefig(os.path.join(self.outputdir, 'astrometric_prediction_' + self.title)+'.pdf', transparent=True, pad_inches=0)

    def astrometric_prediction(self, JDepochs):
        """
        Version of astrometric_prediction that uses orvara's orbit internals so
        that the 3-body approximation is actually used.

        Note that this function can be evaluating at arbitrary many epochs (JDepochs just has to be an array of
        Julian Dates). The evaluation time of this only scales like how many orbits you use (self.num_orbits),
        not how many epochs you want to predict. You can basically predict arbitrarily many epochs for free.

        returns: array of RA and array of DEC values. Each array is shape Nxm where N is the number of
        sampled orbits (e.g. 1000) and m is the number of times in JDepochs.

        Computation time of this function is independent of m, but scales O(N).
        """
        ra, dec = [], []
        #sep, pa = [], []  # we dont calculate sep, pa anymore because pa is ill behaved near the star.
        # I.e. just binning pa, we end up binning negative and positive position angles together
        # and averaging to 0 which is nonsensical.
        print(f'simulating {self.num_orbits} orbits for astrometric prediction. Will take roughly 4 seconds per '
              f'100 orbits.')
        for i in range(self.num_orbits):
            orb = Orbit(self, step=self.rand_idx[i], epochs=JDepochs)
            dec.append(orb.relsep * np.cos(orb.PA * np.pi / 180))
            ra.append(orb.relsep * np.sin(orb.PA * np.pi / 180))
            #sep.append(orb.relsep)
            #pa.append(orb.PA)  # we do not save these because PA becomes ill defined when we are close to the star and
            # have a large scatter in separation.
        """
        ra is an array like:
        [[ra_at_epoch0_and_orbit0, ra_at_epoch1_and_orbit0, ...],
         [ra_at_epoch0_and_orbit1, ra_at_epoch1_and_orbit1, ...],
         ...
         [ra_at_epoch0_and_orbitN, ra_at_epoch1_and_orbitN, ...]]
        """
        # convert to mas.
        ra, dec = np.array(ra)*1000, np.array(dec)*1000
        return ra, dec

    def astrometric_prediction_dict(self, JDepochs):
        """
        :param: JDepochs: array of julian dates for which you want the predicted position of the companion.

        returns: dict of mean ra and mean dec, with Guassian 1 standard deviation on those values. Returns
        also, the correlation between ra and dec.
        """
        ra, dec = self.astrometric_prediction(JDepochs)
        # calculating the ra, dec, pa and separation mean values and errors
        mean_ra, mean_dec = np.mean(ra, axis=0), np.mean(dec, axis=0)
        ra_err, dec_err = np.std(ra, axis=0), np.std(dec, axis=0)
        ra_dec_corr = np.mean((ra - mean_ra)*(dec - mean_dec)/(ra_err * dec_err), axis=0)
        return {'ra': mean_ra, 'dec': mean_dec, 'ra_err': ra_err, 'dec_err': dec_err,
                'ra_dec_correlation_coefficient': ra_dec_corr}


# 7. Corner plot
   ################################################################################################
    ############### plot a nicer corner plot###############

    def plot_corner(self, **kwargs):
        rcParams["lines.linewidth"] = 1.0
        rcParams["axes.labelpad"] = 20.0
        rcParams["xtick.labelsize"] = 10.0
        rcParams["ytick.labelsize"] = 10.0
        
        chain = self.chain

        Mpri = chain['mpri']
        npl = '%d' % (self.iplanet)
        if self.cmref == 'msec_solar':
            Msec = chain['msec' + npl]         # in M_{\odot}
            labels=[r'$\mathrm{M_{pri}\, (M_{\odot})}$', r'$\mathrm{M_{sec}\, (M_{\odot})}$', 'a (AU)', r'e', r'$\mathrm{i\, (^{\circ})}$']
        else:
            Msec = chain['msec' + npl]*1989/1.898
            labels=[r'$\mathrm{M_{pri}\, (M_{\odot})}$', r'$\mathrm{M_{sec}\, (M_{Jup})}$', 'a (AU)', r'e', r'$\mathrm{i\, (^{\circ})}$']
        Semimajor = chain['sau' + npl]
        Ecc = chain['esino' + npl]**2 + chain['ecoso' + npl]**2
        Inc = chain['inc' + npl]
        
        chain = np.vstack([Mpri, Msec, Semimajor, Ecc, Inc * 180/np.pi]).T
        
        # in corner_modified, the error is modified to keep 2 significant figures
        figure = corner_modified.corner(chain, labels=labels, quantiles=[0.16, 0.5, 0.84], range=[0.999 for l in labels], verbose=False, show_titles=True,
                                        title_kwargs={"fontsize": 12}, hist_kwargs={"lw":1.}, label_kwargs={"fontsize":15}, xlabcord=(0.5,-0.45), ylabcord=(-0.45,0.5), **kwargs)

        print("Plotting Corner plot, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'Corner_' + self.title)+'.pdf', transparent=True)

###################################################################################################
###################################################################################################

#diagnostic plots

    def plot_chains(self,labels=None,thin=1,alpha=0.1):
        #labels=['Mpri','Msec','a',r'$\mathrm{\sqrt{e}\, sin\, \omega}$',r'$\mathrm{\sqrt{e}\, cos\, \omega}$','inc','asc','lam']
        print("Generating diagnostic plots to check convergence.")

        chain = fits.open(self.MCMCfile)[1].data
        chain = burnin_chain(chain.columns, self.burnin, reshape=False)
        nwalkers, nsteps = chain['lnp'].shape
        ndim = len(pull_chain_params(chain.columns, tuple([[0],[0]]))) + 1
        fig, ax = plt.subplots(nrows=ndim,sharex=True, figsize=(10,10))
        ijit = 0
        for i in range(len(chain.columns)):
            for walker in range(nwalkers):
                ax[i].plot(chain.columns[i].array[walker],
                           color="black", alpha=alpha, lw=0.5);
            ax[i].set_ylabel(chain.columns[i].name,fontsize=15,labelpad = 10)
            ax[i].margins(y=0.1)
            for label in ax[i].get_yticklabels():
                label.set_fontsize(15)
                
            ax[i].yaxis.set_label_coords(-0.1, 0.5)
            if chain.columns[i].name == 'lnp':
                break

        ax[i].set_xlabel("sample",fontsize=15)
        ax[i].minorticks_on()
        ax[0].set_title("Overview of chains",y=1.03,fontsize=15)
        for label in ax[i].get_xticklabels():
                label.set_fontsize(15)
        plt.savefig(os.path.join(self.outputdir, 'Diagnostic_' + self.title)+'.png')


#save data

    def safe_get(self, chain, key):
        """
        A wrapper for self.chain[key] that won't crash if key does not exist
        in the chain. This is needed because astropy.fits.BinTableHDU objects do
        not have a .get() method like dictionaries, that allow you to return e.g. None
        """
        try:
            values = chain[key]
        except KeyError:
            values = None
        return values

    def save_data(self):
    
        def print_best_chisq():
            par_label = ['plx_ML', 'pmra_ML', 'pmdec_ML', 'chisq_sep', 'chisq_PA', 'chisq_H',
                         'chisq_HG', 'chisq_G', 'chisq_relRV']
            print("Saving beststep parameters to " + self.outputdir)
            text_file = open(os.path.join(self.outputdir, 'beststep_params_' + self.title) +'.txt', "w")
            indx = np.where(self.lnp == np.amax(self.lnp))[0][0]
            for i in range(len(par_label)):
                chisq_val = self.safe_get(self.chain, par_label[i])
                str_to_write = " None "
                if chisq_val is not None:
                    str_to_write = "  %.10g" % chisq_val[indx]
                text_file.write(par_label[i])
                text_file.write(str_to_write)
                text_file.write("\n")

            for i in range(len(self.chain.columns)):
                try:
                    offset = self.chain['RV_ZP_%d_ML' % (i)][indx]
                    text_file.write('RV_ZP_%d_ML' % (i))
                    text_file.write("  %.10g" % (offset))
                    text_file.write("\n")
                except Exception:
                    break
            text_file.close()
        
        def print_par_values(x,m):
            q_16, q_50, q_84 = corner_modified.quantile(x, m)
            q_m, q_p = q_50-q_16, q_84-q_50

            # modified to keep 2 significant figures in the errors
            idecimal_m = np.floor(np.log10(float('%.1g'%(q_m))))
            idecimal_p = np.floor(np.log10(float('%.1g'%(q_p))))

            if idecimal_m < 2:
                fmt_m_e = "{{0:{0}}}".format(".%df"%(-idecimal_m + 1)).format
            else:
                fmt_m_e = "{{0:{0}}}".format(".0f").format
            if idecimal_p < 2:
                fmt_p_e = "{{0:{0}}}".format(".%df"%(-idecimal_p + 1)).format
            else:
                fmt_p_e = "{{0:{0}}}".format(".0f").format
            min_decimals = min(idecimal_m, idecimal_p)

            if min_decimals < 2:
                fmt = "{{0:{0}}}".format(".%df"%(-min_decimals + 1)).format
            else:
                fmt = "{{0:{0}}}".format(".0f").format

            title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
            title = title.format(fmt(q_50), fmt_m_e(q_m), fmt_p_e(q_p))
            if m == [0.16, 0.5, 0.84]:
                return(title)
            else:
                return(round(q_16,3), round(q_84,3))
        
        def print_posterior(perc_sigmas):
        
            chain = self.chain
            npl = '%d' % (self.iplanet)
            
            #save posterior and derived parameters
            if 'jitter' in [col.name for col in chain.columns]:
                RV_Jitter = [print_par_values(chain['jitter'],perc_sigmas)]
                jitter_labels = ['Jit (m/s)']
            else:
                i = 0
                RV_Jitter = []
                jitter_labels = []
                while 'jitter%d' % (i) in [col.name for col in chain.columns]:
                    RV_Jitter += [print_par_values(chain['jitter%d' % (i)],perc_sigmas)]
                    jitter_labels += ['Jit, Inst #%d (m/s)' % (i)]
                    i += 1

            Mpri = print_par_values(chain['mpri'],perc_sigmas)
            if self.cmref == 'msec_solar':
                Msec = print_par_values(chain['msec' + npl],perc_sigmas)
                unit = '(solar)'
            else:
                Msec = print_par_values(chain['msec' + npl]*1989/1.898,perc_sigmas)
                unit = '(jup)'
            a = print_par_values(chain['sau' + npl],perc_sigmas)
            sqesino = print_par_values(chain['esino' + npl],perc_sigmas)
            sqecoso = print_par_values(chain['ecoso' + npl],perc_sigmas)
            inc = print_par_values((chain['inc' + npl]*180/np.pi)%(180),perc_sigmas)
            asc = print_par_values((chain['asc' + npl]*180/np.pi)%(360),perc_sigmas)
            lam = print_par_values((chain['lam' + npl]*180/np.pi)%(360),perc_sigmas)
            # Multiply by 1000 for units of mas
            plx = print_par_values(1e3*chain['plx_ML'],perc_sigmas)
            period_data = np.sqrt(chain['sau' + npl]**3/(chain['mpri'] + chain['msec' + npl]))
            period = print_par_values(period_data,perc_sigmas)
            omega_data=(np.arctan2(chain['esino' + npl],chain['ecoso' + npl])%(2*np.pi))
            omega = print_par_values(omega_data*180/np.pi, perc_sigmas)
            e = print_par_values(chain['esino' + npl]**2 + chain['ecoso' + npl]**2,perc_sigmas)
            sma = print_par_values(1e3*chain['plx_ML']*chain['sau' + npl],perc_sigmas)
            t0_data = 2455197.5 + 365.25*period_data/(2*np.pi)*((omega_data -chain['lam' + npl])%(2*np.pi)) #reference epoch 2455197.5
            t0 = print_par_values(t0_data,perc_sigmas)
            q = print_par_values(chain['msec' + npl]/chain['mpri'],perc_sigmas)
            
            label = jitter_labels
            label += ['Mpri (solar)', 'Msec '+unit, 'a (AU)', 'sqesinw','sqecosw', 'inclination (deg)','ascending node (deg)', 'mean longitude (deg)','parallax (mas)', 'period (yrs)', 'argument of periastron (deg)', 'eccentricity', 'semimajor axis (mas)', 'T0 (JD)', 'mass ratio' ]
            result = RV_Jitter
            result += [Mpri, Msec, a, sqesino, sqecoso, inc, asc, lam, plx, period, omega, e, sma, t0, q]

            print("Saving posterior parameters to " + self.outputdir)
            text_file = open(os.path.join(self.outputdir, 'posterior_params_' + self.title) +'.txt', "w")

            for i in range(len(label)):
                text_file.write(label[i])
                text_file.write("         ")
                text_file.write(str(result[i]))
                text_file.write("\n")
            text_file.close()
        
        if self.save_params:
            print_best_chisq()
            print_posterior(list([float(err) for err in self.err_margin]))
            
#######
# end of code
