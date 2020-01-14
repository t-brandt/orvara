import os
import numpy as np
from random import randrange
import emcee, corner
from scipy.interpolate import interp1d
import scipy.optimize as op
from orbit3d import orbit
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
from geomdl import BSpline
from geomdl import utilities
from geomdl import operations

"""
    Example:
    OPs = orbit_plots.OrbitPlots(title, Hip, start_ep, end_ep, predicted_ep, cmref, num_lines, cm_name, usecolorbar, colorbar_size, colorbar_pad, burnin, user_xlim, user_ylim, user_xtext, user_ytext, steps, mcmcfile, RVfile, AstrometryFile, HGCAFile, outputdir)
    OPs.astrometry()
    OPs.RV()
    OPs.relRV()
    OPs.relsep()
    OPs.PA()
    OPs.proper_motions()
    OPs.plot_corner()
"""

# Progress:
# Finalized all of the preliminary plots: OPs.astrometry(), OPs.RV(), OP.relRV(), OPs.relsep() , OPs.PA(), OPs.proper_motions(), OPS.plot_corner()
# deleted the 'manually changing the axes ticks'

# TO DO:
# Fix the scipy method to rotate the astrometric orbits plots predicted epoch labels
# Fix the excepts, or delete them

class OrbitPlots:

###################################### Initialize Class ############################################
    def __init__(self, title, Hip, start_ep, end_ep, predicted_ep, cmref, num_lines, cm_name, usecolorbar, colorbar_size, colorbar_pad, marker_color, burnin, set_limit, user_xlim, user_ylim, show_title, add_text, user_text_name, user_xtext, user_ytext, steps, mcmcfile, RVfile, AstrometryFile, HGCAFile, outputdir, whichInstrument, pm_separate):
        
        self.title = title
        self.Hip = Hip
        self.start_epoch = start_ep
        self.end_epoch = end_ep
        self.predicted_ep = predicted_ep
        self.cmref = cmref
        self.num_lines = num_lines
        self.cm_name = cm_name
        self.marker_color = marker_color
        self.burnin = burnin
        self.MCMCfile = mcmcfile
        self.RVfile = RVfile
        self.relAstfile = AstrometryFile
        self.HGCAFile = HGCAFile
        self.outputdir = outputdir
        self.set_limit = set_limit
        self.user_xlim = user_xlim
        self.user_ylim = user_ylim
        self.steps = steps
        self.usecolorbar = usecolorbar
        self.colorbar_size = colorbar_size
        self.colorbar_pad = colorbar_pad
        self.show_title = show_title
        self.add_text = add_text
        self.text_name = user_text_name
        self.x_text = user_xtext
        self.y_text = user_ytext
        self.pm_separate = pm_separate
        self.whichInst = whichInstrument
        
        self.cmlabel_dic = {'msec': r'$\mathrm{M_{comp} (M_{Jup})}$', 'ecc': 'Eccentricity'}
        self.color_list = ['r', 'g', 'b', 'y', 'c', 'b']
        
        ############################### load in data #######################
        # define epochs
        self.epoch, self.epoch_calendar = self.define_epochs()
        # load mcmc data
        self.chain, self.beststep, self.extras = self.load_mcmc_data()
        # load observed RV data
        self.epoch_obs, self.RV_obs, self.RV_obs_err, self.nInst, self.epoch_obs_dic, self.RV_obs_dic, self.RV_obs_err_dic = self.load_obsRV_data()
        # load relative astrometry data:
        self.ep_relAst_obs, self.relsep_obs, self.relsep_obs_err, self.PA_obs, self.PA_obs_err = self.load_relAst_data()
        # load HGCA data:
        self.ep_mualp_obs, self.ep_mudec_obs, self.mualp_obs, self.mudec_obs, self.mualp_obs_err, self.mudec_obs_err = self.load_HGCA_data()
        
        ############################## calculate orbits ####################
        # calcualte the best fit orbit
        self.data, self.dras_ml, self.ddecs_ml, self.RV_ml, self.mu_RA_ml, self.mu_Dec_ml, self.relsep_ml, self.mualp_ml, self.mudec_ml, self.f_RVml, self.f_relsepml, self.f_mualpml, self.f_mudecml, self.PA_ml, self.TA_ml, self.node0, self.node1, self.idx_pras, self.pras = self.bestfit_orbit()  # plx in units of arcsecs
        # calculate more random orbits drawn from the mcmc chian
        self.RV_dic, self.dras_dic, self.ddecs_dic, self.relsep_dic, self.PA_dic, self.mualp_dic, self.mudec_dic, self.dic_keys, self.RV_dic_vals, self.dras_dic_vals, self.ddecs_dic_vals, self.relsep_dic_vals, self.PA_dic_vals, self.mualp_dic_vals, self.mudec_dic_vals = self.random_orbits()
        
        ################################ set colorbar ######################
        # setup the normalization and the colormap
        self.nValues = np.array(self.dic_keys)*1989/1.898
        self.normalize = mcolors.Normalize(vmin=self.nValues.min(), vmax=self.nValues.max())
        self.colormap = getattr(cm, self.cm_name)
        self.sm = cm.ScalarMappable(norm=self.normalize, cmap=self.colormap)
        self.sm.set_array(self.nValues)
        
        print("Generating plots for target " + self.title)

###################################### Define Functions ############################################
    def JD_to_calendar(self, JD):
        """
            Function to convert Julian Date to Calendar Date
        """
        a = int(JD + 0.5)
        if a < 2299161:
            c = a + 1524
        else:
            b = int((a - 1867216.25)/36524.25)
            c = a + b - int(b/4) + 1525
        d = int((c - 122.1)/365.25)
        e = int(365.25*d)
        f = int((c - e)/30.6001)

        D = c - e - int(30.6001*f) + ((JD + 0.5) - int(JD + 0.5))
        M = f - 1 - 12*int(f/14.)
        Y = d - 4715 - int((7 + M)/10.)
        year = Y + M/12. + D/365.25
        return year

    def calendar_to_JD(self, year, M=1, D=1):
        """
            Function to convert Calendar Date to Julian Date
        """
        if M <= 2:
            y = year - 1
            m = M + 12
        else:
            y = year
            m = M
        if year <= 1583: # more precisely should be less than or equal to 10/4/1582
            B = -2
        else:
            B = int(y/400.) - int(y/100.)
        UT = 0
        JD = int(365.25*y) + int(30.6001*(m+1)) + B + 1720996.5  + D + UT/24.
        return JD

    def define_epochs(self):
        """
            Function to define a custom range of epochs
        """
        start_epoch = self.calendar_to_JD(self.start_epoch)
        end_epoch = self.calendar_to_JD(self.end_epoch)
        range_epoch = end_epoch - start_epoch
        epoch = np.linspace(start_epoch - 0.1*range_epoch, end_epoch + 0.5*range_epoch, self.steps)
        epoch_calendar = np.zeros(len(epoch))
        for i in range(len(epoch_calendar)):
            epoch_calendar[i] = self.JD_to_calendar(epoch[i])
        return epoch, epoch_calendar

    def load_mcmc_data(self):
        """
            Function to load in the MCMC chain from fit_orbit
        """
        source = self.MCMCfile.split('_')[0]
        chain, lnp, extras = [fits.open(self.MCMCfile)[i].data for i in range(3)]
        beststep = np.where(lnp == lnp.max())
        return chain, beststep, extras
        
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
            # Check to see that the column we loaded was an integer
            assert np.all(RVinst == rvdat[:, 3])
            nInst = int(np.amax(rvdat[:, 3]) + 1)
            
            self.multi_instr = True
            idx_dic = {}
            epoch_obs_dic = {}
            RV_obs_dic = {}
            RV_obs_err_dic = {}
            
            for i in range(nInst):
                idx_dic[i] = (np.where(RVinst == i)[0][0], np.where(RVinst == i)[0][-1])
                epoch_obs_dic[i] = epoch_obs[idx_dic[i][0]: idx_dic[i][-1] + 1]
                RV_obs_dic[i] = RV_obs[idx_dic[i][0]: idx_dic[i][-1] + 1]
                RV_obs_err_dic[i] = RV_obs_err[idx_dic[i][0]: idx_dic[i][-1] + 1]
        except:
            self.multi_instr = False
        return epoch_obs, RV_obs, RV_obs_err, nInst, epoch_obs_dic, RV_obs_dic, RV_obs_err_dic

    def load_relAst_data(self):
        """
            Function to load in the relative astrometry data
        """
        try:
            reldat = np.genfromtxt(self.relAstfile, usecols=(1,2,3,4,5), skip_header=1)
            self.have_reldat = True
            
            ep_relAst_obs = reldat[:, 0]
            for i in range(len(ep_relAst_obs)):
                ep_relAst_obs[i] = self.calendar_to_JD(ep_relAst_obs[i])
            relsep_obs = reldat[:, 1]
            relsep_obs_err = reldat[:, 2]
            PA_obs = reldat[:, 3]
            PA_obs_err = reldat[:, 4]
        except:
            self.have_reldat = False
        return ep_relAst_obs, relsep_obs, relsep_obs_err, PA_obs, PA_obs_err
    
    def load_HGCA_data(self):
        """
            Function to load in the epoch astrometry data
        """
        t = fits.open(self.HGCAFile)[1].data
        try:
            self.have_pmdat = True
            i = int(np.where(t['hip_id'] == self.Hip)[0])
            ep_mualp_obs = np.array([t['epoch_ra_hip'][i], t['epoch_ra_gaia'][i]])
            ep_mudec_obs = np.array([t['epoch_dec_hip'][i], t['epoch_dec_gaia'][i]])
            mualp_obs = np.array([t['pmra_hip'][i], t['pmra_gaia'][i]])
            mualp_obs_err = np.array([t['pmra_hip_error'][i], t['pmra_gaia_error'][i]])
            mudec_obs = np.array([t['pmdec_hip'][i], t['pmdec_gaia'][i]])
            mudec_obs_err = np.array([t['pmdec_hip_error'][i], t['pmdec_gaia_error'][i]])
            for i in range(len(ep_mualp_obs)):
                ep_mualp_obs[i] = self.calendar_to_JD(ep_mualp_obs[i])
                ep_mudec_obs[i] = self.calendar_to_JD(ep_mudec_obs[i])
        except:
            self.have_pmdat = False
        return ep_mualp_obs, ep_mudec_obs, mualp_obs, mudec_obs, mualp_obs_err, mudec_obs_err
    
    def calc_RV_offset(self, walker_idx, step_idx):
        """
            Function to calculate the offset of the observed RV data
        """
        try:
            # calculate the offsets of the RV curves
            assert self.multi_instr
            offset_dic = {}
            for i in np.arange(8, 8 + self.nInst, 1):
                offset = self.extras[walker_idx, step_idx, i]
                offset_dic[i-8] = offset
        except:
            offset = self.extras[walker_idx, step_idx, 8]
        return offset_dic

    def bestfit_orbit(self):
        """
            Function to calculate the most likely orbit
        """
        data = orbit.Data(self.Hip, self.RVfile, self.relAstfile)
        par = orbit.Params(self.chain[self.beststep][0]) # beststep is the best fit orbit
        plx = self.extras[self.beststep[0], self.beststep[1], 0]
        data.custom_epochs(self.epoch)
        model = orbit.Model(data)

        orbit.calc_EA_RPP(data, par, model)
        orbit.calc_offsets(data, par, model, 0) # change 0 into an argument
        orbit.calc_RV(data, par, model)
        
        # most likely orbit for delRA and delDec
        dRAs_G, dDecs_G, dRAs_H1, dDecs_H1, dRAs_H2, dDecs_H2 = model.return_dRA_dDec()
        ratio = -(1. + par.mpri/par.msec)*plx
        dras_ml, ddecs_ml = ratio*dRAs_G, ratio*dDecs_G       #only considering Gaia
        # most likely orbit for RV
        RV_ml = model.return_RVs()
        # most likely orbit for proper motions
        mu_RA_ml, mu_Dec_ml =  model.return_proper_motions(par)
        # most likely orbit for relative separation
        relsep_ml = model.return_relsep()*plx             # relsep in arcsec
        # most likely orbit for position angle
        PA_ml = model.return_PAs()*180/np.pi
        PA_ml = PA_ml % 360.
        # most likely orbit for proper motions
        mualp_ml, mudec_ml = model.return_proper_motions(par)
        mualp_ml, mudec_ml = 1e3*plx*365.25*mualp_ml, 1e3*plx*365.25*mudec_ml   # convert from arcsec/day to mas/yr
        mualp_ml += self.extras[self.beststep[0], self.beststep[1], 1]*1000
        mudec_ml += self.extras[self.beststep[0], self.beststep[1], 2]*1000
        
        # interpolation
        f_RVml = interp1d(self.epoch, RV_ml, fill_value="extrapolate")
        f_relsepml = interp1d(self.epoch, relsep_ml, fill_value="extrapolate")
        f_mualpml = interp1d(self.epoch, mualp_ml, fill_value="extrapolate")
        f_mudecml = interp1d(self.epoch, mudec_ml, fill_value="extrapolate")
        
        # redefine RV_ml
        self.offset_ml = self.calc_RV_offset(self.beststep[0], self.beststep[1])
        
        # find the positions of nodes and periastron
        TA_ml = model.return_TAs(par)
        # When par.arg (omega) is negative, it means the periastron is below the plane
        # of the sky. We can set omega = -par.arg, which is the angle symmetric with
        # respect to the plane of the sky. Although this angle doesn't mean anything, the
        # same algorithm below can be applied to locate the position of nodes.
        omega = abs(par.arg)
        
        idx_node0 = np.where(abs(TA_ml - (np.pi - omega)) == min(abs(TA_ml - (np.pi - omega))))[0]
        idx_node1 = np.where(abs(TA_ml - (-omega)) == min(abs(TA_ml - (-omega))))[0]
        node0 = (dras_ml[idx_node0], ddecs_ml[idx_node0])
        node1 = (dras_ml[idx_node1], ddecs_ml[idx_node1])
        idx_pras = np.where(abs(TA_ml) == min(abs(TA_ml)))[0]
        pras = (dras_ml[idx_pras], ddecs_ml[idx_pras])
        
        return data, dras_ml, ddecs_ml, RV_ml, mu_RA_ml, mu_Dec_ml, relsep_ml, mualp_ml, mudec_ml, f_RVml, f_relsepml, f_mualpml, f_mudecml, PA_ml, TA_ml, node0, node1, idx_pras, pras

    def random_orbits(self):
        """
            Function to calculate more orbits with parameters randomly drawn from the mcmc chian, num_orbits can be specified in the config.ini file
        """
        RV_dic = {}
        dras_dic = {}
        ddecs_dic = {}
        relsep_dic = {}
        PA_dic = {}
        mualp_dic = {}
        mudec_dic = {}
        
        for i in range(self.num_lines):

            # get parameters from one single step of the mcmc chain
            walker_idx = randrange(self.chain.shape[0])
            step_idx = randrange(self.burnin, self.chain.shape[1])
            par = orbit.Params(self.chain[walker_idx, step_idx])
            plx = self.extras[walker_idx, step_idx, 0]

            # calculate and assign variables
            data = self.data
            data.custom_epochs(self.epoch)
            model = orbit.Model(data)

            orbit.calc_EA_RPP(data, par, model)
            orbit.calc_offsets(data, par, model, 0)
            orbit.calc_RV(data, par, model)
            dRAs_G, dDecs_G, dRAs_H1, dDecs_H1, dRAs_H2,  dDecs_H2 = model.return_dRA_dDec()
            ratio = -(1. + par.mpri/par.msec)*plx
            dras, ddecs = ratio*dRAs_G, ratio*dDecs_G # convert from AU to mas
            RV = model.return_RVs()
            relsep = model.return_relsep()*plx
            PA = model.return_PAs()*180/np.pi
            PA = PA % 360.
            mualp, mudec = model.return_proper_motions(par)
            mualp, mudec = 1e3*plx*365.25*mualp, 1e3*plx*365.25*mudec   # convert from arcsec/day to mas/yr
            
            # shift the RV curve wrt data points
            offset = np.sum(self.calc_RV_offset(walker_idx, step_idx)[i] for i in range(self.nInst))/self.nInst
            RV += np.sum(self.offset_ml[i] for i in range(self.nInst))/self.nInst - offset
            # shift curves to the data points
            mualp += self.extras[walker_idx, step_idx, 1]*1000
            mudec += self.extras[walker_idx, step_idx, 2]*1000
        
            cmref = getattr(par, self.cmref)
            RV_dic[cmref] = RV
            dras_dic[cmref] = dras
            ddecs_dic[cmref] = ddecs
            relsep_dic[cmref] = relsep
            PA_dic[cmref] = PA
            mualp_dic[cmref] = mualp
            mudec_dic[cmref] = mudec

        # sort the diconaries in terms of msec/ecc/etc.
        RV_dic = dict(sorted(RV_dic.items(), key=lambda key: key[0])) # if key[1], sort in terms of values
        dras_dic = dict(sorted(dras_dic.items(), key=lambda key: key[0]))
        ddecs_dic = dict(sorted(ddecs_dic.items(), key=lambda key: key[0]))
        relsep_dic = dict(sorted(relsep_dic.items(), key=lambda key: key[0]))
        PA_dic = dict(sorted(PA_dic.items(), key=lambda key: key[0]))
        mualp_dic = dict(sorted(mualp_dic.items(), key=lambda key: key[0]))
        mudec_dic = dict(sorted(mudec_dic.items(), key=lambda key: key[0]))
        dic_keys = list(RV_dic.keys())  # this gives a list of msec from orbits, from small to large
        RV_dic_vals = list(RV_dic.values())
        dras_dic_vals = list(dras_dic.values())
        ddecs_dic_vals = list(ddecs_dic.values())
        relsep_dic_vals = list(relsep_dic.values())
        PA_dic_vals = list(PA_dic.values())
        mualp_dic_vals = list(mualp_dic.values())
        mudec_dic_vals = list(mudec_dic.values())

        return RV_dic, dras_dic, ddecs_dic, relsep_dic, PA_dic, mualp_dic, mudec_dic, dic_keys, RV_dic_vals, dras_dic_vals, ddecs_dic_vals, relsep_dic_vals, PA_dic_vals, mualp_dic_vals, mudec_dic_vals


########################################## Plotting ################################################

# 1. Astrometric plots
# Finalized
# need to fix some issues with scipy method to rotate the labels

    ###############################################################################################
    ############################## plot astrometric orbits ######################
    
    def astrometry(self):
        
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111)

        # plot the num_lines randomly selected curves
        for i in range(self.num_lines):
            ax.plot(self.dras_dic_vals[i], self.ddecs_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.4, linewidth=0.8)

        # plot the most likely orbit
        ax.plot(self.dras_ml, self.ddecs_ml, color='black')

        # plot the relAst data points
        f_drasml = interp1d(self.epoch, self.dras_ml)
        f_ddecsml = interp1d(self.epoch, self.ddecs_ml)

        try:
            assert self.have_reldat == True
            ra_obs = self.relsep_obs * np.sin(self.PA_obs*np.pi /180.)
            dec_obs = self.relsep_obs * np.cos(self.PA_obs*np.pi /180.)
            ax.scatter(ra_obs, dec_obs, s=45, facecolors=self.marker_color, edgecolors='none', zorder=99)
            ax.scatter(ra_obs, dec_obs, s=45, facecolors='none', edgecolors='k', zorder=100)
        except:
            pass

        # plot the predicted positions of the companion star from 1990 to 2030
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        
        epoch_int = []
        for year in self.epoch_calendar:
            epoch_int.append(int(year))
        
        # define some functions
        def calc_linear(x,y):
            x1, x2 = x[0], x[1]
            y1, y2 = y[0], y[1]
            m = (y1- y2)/(x1 - x2)
            b = y1 - m*x1
            return m,b
            
        def array_list(array_num):
            num_list = array_num.tolist() # list
            return num_list
        
        # new method to rotate the labels according to angle of normal to the curve tangent
        # need to install the geomdl package to calculate the tangent line
        for year in self.predicted_ep:
            year = int(year)
            idx = epoch_int.index(year)
            x = self.dras_ml[idx]
            y = self.ddecs_ml[idx]
            r = np.sqrt(x**2 + y**2)
            
            # rotate the labels
            data_list = array_list(np.hstack(([np.vstack(self.dras_ml),np.vstack(self.ddecs_ml)])))

            # Create a BSpline curve instance
            curve = BSpline.Curve()
            curve.degree = 3
            curve.ctrlpts = data_list
            curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))
            curve.delta = 0.01
            curve.evaluate()

            # Evaluate curve tangent at idx/step
            curvetan = []
            ct = operations.tangent(curve, idx/self.steps, normalize=True)  # 1500 = self.steps
            curvetan.append(ct)

            ctrlpts = np.array(curve.ctrlpts)
            ctarr = np.array(curvetan) # Convert tangent list into a NumPy array
            ax.scatter(x, y, s=55, facecolors='none', edgecolors='k', zorder=100)

            x=[ctarr[:, 0, 0] ,ctarr[:, 0, 0]+ctarr[:, 1, 0] ]
            y=[ctarr[:, 0, 1], ctarr[:, 0, 1]+ctarr[:, 1, 1] ]

            m, b = calc_linear(x,y)
            
            #changed
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            xlim = [x0,x1]
            ylim = [y0,y1]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            dd = self.ddecs_ml
            dra = self.dras_ml

            y_intercept = m*xlim[0]+b
            x_intercept = (ylim[0]-b)/m
            angle = np.arctan(np.abs(y_intercept - ylim[0])/np.abs(x_intercept - xlim[0]))*180./np.pi
            angle -= 90

            # calculate the line of semi-major axis
            max_y = dd[(np.where(dra == max(dra)))[0][0]]
            min_y = dd[(np.where(dra == min(dra)))[0][0]]
            max_x = dra[(np.where(dd == max(dd)))[0][0]]
            min_x = dra[(np.where(dd == min(dd)))[0][0]]

            m_semimajor, b_semimajor = calc_linear([min_x,max_x],[min(dd),max(dd)])
            m_semiminor, b_semiminor = calc_linear([min(dra),max(dra)],[min_y,max_y])
            if dra[idx] < (dd[idx]-b_semimajor)/m_semimajor and dd[idx] > m_semiminor*dra[idx] + b_semiminor:
                ax.annotate('  ' + str(year), xy=(ctarr[0][0][0], ctarr[0][0][1]), verticalalignment='bottom', horizontalalignment='left',rotation =-angle[0])
            if dra[idx] < (dd[idx]-b_semimajor)/m_semimajor and dd[idx] < m_semiminor*dra[idx] + b_semiminor:
                ax.annotate('  ' + str(year), xy=(ctarr[0][0][0], ctarr[0][0][1]), verticalalignment='top', horizontalalignment='left',rotation =angle[0])
            elif dra[idx] > (dd[idx]-b_semimajor)/m_semimajor and dd[idx] < m_semiminor*dra[idx] + b_semiminor:
                ax.annotate(str(year)+'  ', xy=(ctarr[0][0][0], ctarr[0][0][1]), verticalalignment='top', horizontalalignment='right',rotation= -angle[0])
            elif dra[idx] > (dd[idx]-b_semimajor)/m_semimajor and dd[idx] > m_semiminor*dra[idx] + b_semiminor:
                ax.annotate(str(year)+'  ', xy=(ctarr[0][0][0], ctarr[0][0][1]), verticalalignment='bottom', horizontalalignment='right', rotation=angle[0])

        # plot line of nodes, periastron and the direction of motion of the companion
        ax.plot([self.node0[0], self.node1[0]], [self.node0[1], self.node1[1]], 'k--',linewidth = 1)
        ax.plot([0, self.pras[0]], [0, self.pras[1]], 'k:')
        # changed, self.idx_pras+1 maybe out of range, changed to -1
        arrow = mpatches.FancyArrowPatch((self.dras_ml[self.idx_pras-1][0], self.ddecs_ml[self.idx_pras-1][0]),(self.dras_ml[self.idx_pras][0], self.ddecs_ml[self.idx_pras][0]),  arrowstyle='->', mutation_scale=25, zorder=100)
        ax.add_patch(arrow)
        ax.plot(0, 0, marker='*',  color='black', markersize=10)
        
        # from advanced plotting settings in config.ini
        if self.set_limit is True:
            ax.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
            ax.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            aspect0 = np.float(self.user_xlim[1])-np.float(self.user_xlim[0])
            aspect1 = np.float(self.user_ylim[1])-np.float(self.user_ylim[0])
            ax.set_aspect(np.abs(aspect0/aspect1))
        if self.show_title is True:
            ax.set_title('Astrometric Orbits')
        if self.add_text is True:
            ax.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
        if self.usecolorbar is True:
            cbar = fig.colorbar(self.sm, ax=ax, fraction=self.colorbar_size, pad=self.colorbar_pad)
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
        plt.tight_layout()
        plt.savefig(os.path.join(self.outputdir,'astrometric_orbit_' + self.title)+'.pdf') # or +'.png'

# 2. RV orbits plot
# FINALIZED

    ################################################################################################
    ########################### plot the RV orbits #######################
    
    def RV(self):
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_axes((0.15, 0.3, 0.8, 0.6))

        # plot the 50 randomly selected curves
        for i in range(self.num_lines):
            ax.plot(self.epoch_calendar, self.RV_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

        # plot the most likely one
        ax.plot(self.epoch_calendar, self.RV_ml, color='black')

        # plot the observed data points (RV & relAst)
        try:
            assert self.multi_instr
            rv_epoch_list = []
            for i in range(self.nInst):
                epoch_obs_Inst = np.zeros(len(self.epoch_obs_dic[i]))
                for j in range(len(self.epoch_obs_dic[i])):
                    epoch_obs_Inst[j] = self.JD_to_calendar(self.epoch_obs_dic[i][j])
                rv_epoch_list.append(epoch_obs_Inst)
                
            for i in range(self.nInst):
                ax.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + self.offset_ml[i], yerr=self.RV_obs_err_dic[i], fmt=self.color_list[i]+'o', ecolor='black', markersize = 1, elinewidth = 0.3, capsize=1, capthick = 0.3, zorder = 200+i, alpha = 0.8)
                ax.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + self.offset_ml[i], s=3, facecolors='none', edgecolors='k', zorder=200+i, alpha = 0.8)
        except:
            ax.plot(self.JD_to_calendar(self.epoch_obs), self.RV_obs, 'ro', markersize=2)
            
        ax.set_xlim(self.start_epoch, self.end_epoch)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        if self.usecolorbar is True:
            cbar = fig.colorbar(self.sm, ax=ax, fraction=self.colorbar_size, pad=self.colorbar_pad)
            cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
            cbar.ax.get_yaxis().labelpad=20
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax.set_xlabel('Epoch (year)', labelpad = 10, fontsize=13)
        ax.set_ylabel('RV (m/s)', fontsize=13)

        plt.tight_layout()
        print("Plotting RV orbits, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'RV_orbit_' + self.title)+'.pdf') # or +'.pdf'



# 3. relative RV plot
# FINALIZED
    ################################################################################################
    ############### plot the relative RV and O-C #####################

    def relRV(self):
        fig = plt.figure(figsize=(5, 6))
        ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
        ax2 = fig.add_axes((0.15, 0.1, 0.8, 0.15))

        # plot the 50 randomly selected curves
        self.f_RVml = interp1d(self.epoch, self.RV_ml)
        RV_OC = self.RV_dic_vals
        for i in range(self.num_lines):
            ax1.plot(self.epoch_calendar, self.RV_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
            for j in range(len(self.epoch)):
                RV_OC[i][j] -= self.f_RVml(self.epoch[j])
            ax2.plot(self.epoch_calendar, RV_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

        # plot the most likely one
        ax1.plot(self.epoch_calendar, self.RV_ml, color='black')
        ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

        try:
            assert self.multi_instr
            rv_epoch_list = []
            for i in range(self.nInst):
                epoch_obs_Inst = np.zeros(len(self.epoch_obs_dic[i]))
                for j in range(len(self.epoch_obs_dic[i])):
                    epoch_obs_Inst[j] = self.JD_to_calendar(self.epoch_obs_dic[i][j])
                rv_epoch_list.append(epoch_obs_Inst)
            
            # plot the observed data points
            datOC_list = []
            # plot which instrument, self.whichInst=1 means 1st instrument, etc.
            whichInst = self.whichInst - 1
            i = whichInst
            if i == whichInst:
                ax1.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + self.offset_ml[i], yerr=self.RV_obs_err_dic[i], fmt=self.color_list[i]+'o', ecolor='black', capsize=3, zorder=299)#, ecolor='black', markersize = 1, elinewidth = 0.3, capsize=1, capthick = 0.3, zorder = 200+i, alpha = 0.8)
                ax1.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + self.offset_ml[i], facecolors='none', edgecolors='k', zorder=300, alpha = 0.8)
                y_err = []
                epoch_list = []
                for j in range(len(self.RV_obs_dic[i])):
                    OC = self.RV_obs_dic[i][j] + self.offset_ml[i] - self.f_RVml(self.epoch_obs_dic[i][j])
                    datOC_list.append(OC)
                    y_err.append(self.RV_obs_err_dic[i][j])
                    epoch_list.append(rv_epoch_list[i][j])
                ax2.errorbar(epoch_list, datOC_list, yerr=y_err, fmt=self.color_list[i]+'o', ecolor='black', capsize=3, zorder=299)
                ax2.scatter(epoch_list, datOC_list, s=45, facecolors='none', edgecolors='k', alpha=0.8, zorder=300)
                
            else:
                for i in range(self.nInst):
                    ax1.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + self.offset_ml[i], yerr=self.RV_obs_err_dic[i], fmt=self.color_list[i]+'o', ecolor='black', markersize = 1, elinewidth = 0.3, capsize=1, capthick = 0.3, zorder = 200+i, alpha = 0.8)
                    ax1.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + self.offset_ml[i], s=3, facecolors='none', edgecolors='k', zorder=200+i, alpha = 0.8)
                    y_err = []
                    epoch_list = []
                    for j in range(len(self.offset_ml)):
                        OC = self.RV_obs_dic[i][j] + self.offset_ml[i] - self.f_RVml(self.epoch_obs_dic[i][j])
                        datOC_list.append(OC)
                        y_err.append(self.RV_obs_err_dic[i][j])
                        epoch_list.append(rv_epoch_list[i][j])
                    ax2.errorbar(epoch_list, datOC_list, yerr=y_err, fmt=self.color_list[i]+'o', ecolor='black', capsize=3)
                    ax2.scatter(epoch_list, datOC_list, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=0.5)
                
            # axes settings
            # ax1
            ax1.get_shared_x_axes().join(ax1, ax2)
            range_ep_obs = max(rv_epoch_list[i]) - min(rv_epoch_list[i])
            range_RV_obs = max(self.RV_obs_dic[i] + self.offset_ml[i]) - min(self.RV_obs_dic[i] + self.offset_ml[i])
            ax1.set_xlim(min(rv_epoch_list[i]) - range_ep_obs/20., max(rv_epoch_list[i]) + range_ep_obs/20.)
            ax1.set_ylim(min(self.RV_obs_dic[i] + self.offset_ml[i]) - range_RV_obs/10., max(self.RV_obs_dic[i] + self.offset_ml[i]) + range_RV_obs/10.)
            ax1.xaxis.set_major_formatter(NullFormatter())
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_ylabel('Relative RV (m/s)', fontsize=13)
            # ax2
            range_datOC = max(datOC_list) - min(datOC_list)
            if np.abs(min(datOC_list)) <= np.abs(max(datOC_list)):
                ax2.set_ylim(-np.abs(max(datOC_list)) - range_datOC/7., max(datOC_list) + range_datOC/7.)
            elif np.abs(min(datOC_list)) > np.abs(max(datOC_list)):
                ax2.set_ylim(min(datOC_list) - range_datOC/7., np.abs(min(datOC_list)) + range_datOC/7.)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (yr)', fontsize=13)
            ax2.set_ylabel('O-C', fontsize=13)
            
            # from advanced plotting settings in config.ini
            if self.set_limit is True:
                ax2.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                ax2.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title is True:
                ax1.set_title('Relative RV vs. Epoch')
            if self.add_text is True:
                ax1.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar is True:
                cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                fig.delaxes(cbar_ax)

        except:
            ax1.errorbar(self.epoch_obs, self.RV_obs, yerr=self.RV_obs_err, fmt='bo', ecolor='black', capsize=3)
            ax1.scatter(self.epoch_obs, self.RV_obs, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            for i in range(len(self.epoch_obs)):
                OC = self.RV_obs[i] + self.offset_ml[i] - self.f_RVml(self.epoch_obs[i])
                datOC_list.append(OC)
                ax2.errorbar(self.epoch_obs[i], OC, yerr=self.RV_obs_err[i], fmt='bo', ecolor='black', capsize=3)
                ax2.scatter(self.epoch_obs[i], OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)

        print("Plotting relative RV, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'relRV_OC_' + self.title)+'.pdf')
################################################################################################



# 4. relative separation plot
# FINALIZED
    ################################################################################################
    ##################### plot Relative Separation vs. Epoch and O-C ############
    
    def relsep(self):
    
        try:
            assert self.have_reldat == True
            fig = plt.figure(figsize=(4.7, 5.5))
            ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
            ax2 = fig.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height
            
            # plot the randomly selected curves
            relsep_OC = self.relsep_dic_vals
            for i in range(self.num_lines):
                ax1.plot(self.epoch_calendar, self.relsep_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                for j in range(len(self.epoch)):
                    relsep_OC[i][j] -= self.f_relsepml(self.epoch[j])
                ax2.plot(self.epoch_calendar, relsep_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                
            # plot the highest likelihood orbit
            ax1.plot(self.epoch_calendar, self.relsep_ml, color='black')
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

            # plot the observed data points
            datOC_list = []
            ep_relAst_obs_calendar = []
            for i in range(len(self.ep_relAst_obs)):
                ep_relAst_obs_calendar.append(self.JD_to_calendar(self.ep_relAst_obs[i]))
            ax1.errorbar(ep_relAst_obs_calendar, self.relsep_obs, yerr=self.relsep_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax1.scatter(ep_relAst_obs_calendar, self.relsep_obs, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
            for i in range(len(self.ep_relAst_obs)):
                dat_OC = self.relsep_obs[i] - self.f_relsepml(self.ep_relAst_obs[i])
                datOC_list.append(dat_OC)
            ax2.errorbar(ep_relAst_obs_calendar, datOC_list, yerr=self.relsep_obs_err[i], color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax2.scatter(ep_relAst_obs_calendar, datOC_list, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
                
            # axes settings
            # ax1
            ax1.get_shared_x_axes().join(ax1, ax2)
            range_eprA_obs = max(ep_relAst_obs_calendar) - min(ep_relAst_obs_calendar)
            range_relsep_obs = max(self.relsep_obs)  - min(self.relsep_obs)
            ax1.set_xlim(min(ep_relAst_obs_calendar) - range_eprA_obs/8., max(ep_relAst_obs_calendar) + range_eprA_obs/8.)
            ax1.set_ylim(min(self.relsep_obs) - range_relsep_obs/2., max(self.relsep_obs) + range_relsep_obs/2.)
            ax1.xaxis.set_major_formatter(NullFormatter())
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_ylabel('Seperation (arcsec)', labelpad=13, fontsize=13)
            # ax2
            range_datOC = max(datOC_list) - min(datOC_list)
            if np.abs(min(datOC_list)) <= np.abs(max(datOC_list)):
                ax2.set_ylim(-np.abs(max(datOC_list)) - range_datOC, max(datOC_list) + range_datOC)
            elif np.abs(min(datOC_list)) > np.abs(max(datOC_list)):
                ax2.set_ylim(min(datOC_list) - range_datOC, np.abs(min(datOC_list)) + range_datOC)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
            ax2.set_ylabel('O-C', labelpad=-2, fontsize=13)
            
            # from advanced plotting settings in config.ini
            if self.set_limit is True:
                ax2.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                ax2.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title is True:
                ax1.set_title('Relative Separation vs. Epoch')
            if self.add_text is True:
                ax1.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar is True:
                cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                fig.delaxes(cbar_ax)
                
        except:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_axes((0.15, 0.1, 0.8, 0.8))
            
            # plot the 50 randomly selected curves
            for i in range(self.num_lines):
                ax.plot(self.epoch, self.relsep_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                
            # plot the most likely one
            ax.plot(self.epoch, self.relsep_ml, color='black')
            
            # manually change the x tick labels from JD to calendar years
            epoch_ticks = np.linspace(self.start_epoch, self.end_epoch, 5)
            epoch_label = np.zeros(len(epoch_ticks))
            for i in range(len(epoch_ticks)):
                epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))
                
            # axes settings
            ax.set_xlim(self.start_epoch, self.end_epoch)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.set_xticks(epoch_ticks)
            ax.set_xticklabels([str(int(i)) for i in epoch_label])
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax.set_ylabel('Seperation (arcsec)')
                
        plt.tight_layout()
        print("Plotting Separation, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'relsep_OC_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
################################################################################################



# 5. position angle plot
# FINALIZED
    ################################################################################################
    ############### plot the position angle and O-C ###############

    def PA(self):
    
        try:
            assert self.have_reldat == True
            fig = plt.figure(figsize=(4.7, 5.5))
            ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
            ax2 = fig.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height

            # plot the randomly selected curves from mcmc
            f_PAml = interp1d(self.epoch, self.PA_ml)
            PA_OC = self.PA_dic_vals
            for i in range(self.num_lines):
                ax1.plot(self.epoch_calendar, self.PA_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                for j in range(len(self.epoch)):
                    PA_OC[i][j] -= f_PAml(self.epoch[j])
                ax2.plot(self.epoch_calendar, PA_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the highest likelihood orbit
            ax1.plot(self.epoch_calendar, self.PA_ml, color='black')
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))
                       
            # plot the observed data points
            datOC_list = []
            ep_relAst_obs_calendar = []
            for i in range(len(self.ep_relAst_obs)):
                ep_relAst_obs_calendar.append(self.JD_to_calendar(self.ep_relAst_obs[i]))
            ax1.errorbar(ep_relAst_obs_calendar, self.PA_obs, yerr=self.PA_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax1.scatter(ep_relAst_obs_calendar, self.PA_obs, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
            for i in range(len(self.ep_relAst_obs)):
                dat_OC = self.PA_obs[i] - f_PAml(self.ep_relAst_obs[i])
                datOC_list.append(dat_OC)
            ax2.errorbar(ep_relAst_obs_calendar, datOC_list, yerr=self.PA_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax2.scatter(ep_relAst_obs_calendar, datOC_list, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
            
            # axes settings
            # ax1
            ax1.get_shared_x_axes().join(ax1, ax2)
            range_eprA_obs = max(ep_relAst_obs_calendar) - min(ep_relAst_obs_calendar)
            range_PA_obs = max(self.PA_obs)  - min(self.PA_obs)
            ax1.set_xlim(min(ep_relAst_obs_calendar) - range_eprA_obs/8., max(ep_relAst_obs_calendar) + range_eprA_obs/8.)
            ax1.set_ylim(min(self.PA_obs) - range_PA_obs/5., max(self.PA_obs) + range_PA_obs/5.)
            ax1.xaxis.set_major_formatter(NullFormatter())
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_ylabel(r'Position Angle ($^{\circ}$)', labelpad=5, fontsize=13)
            # ax2
            range_datOC = max(datOC_list) - min(datOC_list)
            if np.abs(min(datOC_list)) <= np.abs(max(datOC_list)):
                ax2.set_ylim(-np.abs(max(datOC_list)) - range_datOC, max(datOC_list) + range_datOC)
            elif np.abs(min(datOC_list)) > np.abs(max(datOC_list)):
                ax2.set_ylim(min(datOC_list) - range_datOC, np.abs(min(datOC_list)) + range_datOC)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
            ax2.set_ylabel('O-C', labelpad=17, fontsize=13)

            # from advanced plotting settings in config.ini
            if self.set_limit is True:
                ax2.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                ax2.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title is True:
                ax1.set_title('Position angle vs. Epoch')
            if self.add_text is True:
                ax1.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar is True:
                cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                fig.delaxes(cbar_ax)

        except:
            fig = plt.figure(figsize=(5, 5.5))
            ax = fig.add_axes((0.15, 0.1, 0.8, 0.8))

            # plot the 50 randomly selected curves
            for i in range(self.num_lines):
                ax.plot(self.epoch, self.PA_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the most likely one
            ax.plot(self.epoch, self.PA_ml, color='black')

            # manually change the x tick labels from JD to calendar years
            epoch_ticks = np.linspace(self.start_epoch, self.end_epoch, 5)
            epoch_label = np.zeros(len(epoch_ticks))
            for i in range(len(epoch_ticks)):
                epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))

            ax.set_xlim(self.start_epoch, self.end_epoch)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.set_xticks(epoch_ticks)
            ax.set_xticklabels([str(int(i)) for i in epoch_label])
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax.set_title(self.title)
            ax.set_ylabel(r'Position Angle ($^{\circ}$)')
        
        plt.tight_layout()
        print("Plotting Position Angle, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir,'PA_OC_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
################################################################################################




# 6. Proper motion plots
# FINALIZED
    ################################################################################################
    ############### plot the proper motions, together or separately ###############
    
    def proper_motions(self):
        try:
            assert self.have_pmdat == True
            if self.pm_separate is True:
                fig = plt.figure(figsize=(4.7, 5.5))
                ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
                ax2 = fig.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height
                fig1 = plt.figure(figsize=(4.7, 5.5))
                ax3 = fig1.add_axes((0.15, 0.3, 0.8, 0.6))
                ax4 = fig1.add_axes((0.15, 0.12, 0.8, 0.16)) # X, Y, width, height
            else:
                fig = plt.figure(figsize=(11,5.5))
                ax1 = fig.add_axes((0.10, 0.30, 0.35, 0.60))
                ax2 = fig.add_axes((0.10, 0.12, 0.35, 0.15))
                ax3 = fig.add_axes((0.6, 0.30, 0.35, 0.60))
                ax4 = fig.add_axes((0.6, 0.12, 0.35, 0.15))

            # plot the num_lines randomly selected curves
            mualp_OC = self.mualp_dic_vals
            mudec_OC = self.mudec_dic_vals

            for i in range(self.num_lines):
                ax1.plot(self.epoch_calendar, self.mualp_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                ax3.plot(self.epoch_calendar, self.mudec_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                for j in range(len(self.epoch)):
                    mualp_OC[i][j] -= self.f_mualpml(self.epoch[j])
                    mudec_OC[i][j] -= self.f_mudecml(self.epoch[j])
                ax2.plot(self.epoch_calendar, mualp_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                ax4.plot(self.epoch_calendar, mudec_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the highest likelihood curve
            ax1.plot(self.epoch_calendar, self.mualp_ml, color='black')
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch_calendar)), 'k--', dashes=(5, 5))
            ax3.plot(self.epoch_calendar, self.mudec_ml, color='black')
            ax4.plot(self.epoch_calendar, np.zeros(len(self.epoch_calendar)), 'k--', dashes=(5, 5))

            # plot the observed data points
            mualpdatOC_list = []
            mudecdatOC_list = []
            ep_mualp_obs_calendar = []
            ep_mudec_obs_calendar = []
            for i in range(len(self.ep_mualp_obs)):
                ep_mualp_obs_calendar.append(self.JD_to_calendar(self.ep_mualp_obs[i]))
            for i in range(len(self.ep_mudec_obs)):
                ep_mudec_obs_calendar.append(self.JD_to_calendar(self.ep_mudec_obs[i]))
            
            ax1.errorbar(ep_mualp_obs_calendar, self.mualp_obs, yerr=self.mualp_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax1.scatter(ep_mualp_obs_calendar, self.mualp_obs, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)
            ax3.errorbar(ep_mudec_obs_calendar, self.mudec_obs, yerr=self.mualp_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax3.scatter(ep_mudec_obs_calendar, self.mudec_obs, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)
            
            # plot the O-C for observed data points
            for i in range(len(self.ep_mualp_obs)):
                dat_OC = self.mualp_obs[i] - self.f_mualpml(self.ep_mualp_obs[i])
                mualpdatOC_list.append(dat_OC)
            ax2.errorbar(ep_mualp_obs_calendar, mualpdatOC_list, yerr=self.mualp_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax2.scatter(ep_mualp_obs_calendar, mualpdatOC_list, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)
            for i in range(len(self.ep_mudec_obs)):
                dat_OC = self.mudec_obs[i] - self.f_mudecml(self.ep_mudec_obs[i])
                mudecdatOC_list.append(dat_OC)
            ax4.errorbar(ep_mudec_obs_calendar  , mudecdatOC_list, yerr=self.mudec_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax4.scatter(ep_mudec_obs_calendar, mudecdatOC_list, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)

            # axes settings
            # ax1
            ax1.get_shared_x_axes().join(ax1, ax2)
            range_eppm_obs = max(ep_mualp_obs_calendar) - min(ep_mualp_obs_calendar)
            range_mualp_obs = max(self.mualp_obs)  - min(self.mualp_obs)
            ax1.set_xlim(min(ep_mualp_obs_calendar) - range_eppm_obs/8., max(ep_mualp_obs_calendar) + range_eppm_obs/8.)
            ax1.set_ylim(min(self.mualp_obs) - range_mualp_obs/5., max(self.mualp_obs) + range_mualp_obs/5.)
            ax1.set_ylabel(r'$\Delta \mu_{\alpha}$ (mas/yr)', labelpad = 5, fontsize = 13)
            # ax3
            ax3.get_shared_x_axes().join(ax3, ax4)
            range_mudec_obs = max(self.mudec_obs)  - min(self.mudec_obs)
            ax3.set_xlim(min(ep_mudec_obs_calendar) - range_eppm_obs/8., max(ep_mudec_obs_calendar) + range_eppm_obs/8.)
            ax3.set_ylim(min(self.mudec_obs) - range_mudec_obs/5., max(self.mudec_obs) + range_mudec_obs/5.)
            ax3.set_ylabel(r'$\Delta \mu_{\delta}$ (mas/yr)', labelpad = 6, fontsize = 13)
            for ax in [ax1, ax3]:
                ax.xaxis.set_major_formatter(NullFormatter())
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            # ax2
            range_mualpdatOC = max(mualpdatOC_list) - min(mualpdatOC_list)
            if np.abs(min(mualpdatOC_list)) <= np.abs(max(mualpdatOC_list)):
                ax2.set_ylim(-np.abs(max(mualpdatOC_list)) - range_mualpdatOC, max(mualpdatOC_list) + range_mualpdatOC)
            elif np.abs(min(mualpdatOC_list)) > np.abs(max(mualpdatOC_list)):
                ax2.set_ylim(min(mualpdatOC_list) - range_mualpdatOC, np.abs(min(mualpdatOC_list)) + range_mualpdatOC)
            for ax in [ax2,ax4]:
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
                ax.set_xlabel('Epoch (year)', labelpad = 6, fontsize = 13)
                ax.set_ylabel('O-C', labelpad = 8, fontsize = 13)
            # ax4
            range_mudecdatOC = max(mudecdatOC_list) - min(mudecdatOC_list)
            if np.abs(min(mudecdatOC_list)) <= np.abs(max(mudecdatOC_list)):
                ax4.set_ylim(-np.abs(max(mudecdatOC_list)) - range_mudecdatOC, max(mudecdatOC_list) + range_mudecdatOC)
            elif np.abs(min(mudecdatOC_list)) > np.abs(max(mudecdatOC_list)):
                ax4.set_ylim(min(mudecdatOC_list) - range_mudecdatOC, np.abs(min(mudecdatOC_list)) + range_mudecdatOC)
                
            # from advanced plotting settings in config.ini
            if self.set_limit is True:
                for ax in [ax1, ax3]:
                    ax.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                    ax.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title is True:
                ax1.set_title('Right Ascension vs. Epoch')
                ax3.set_title('Declination vs. Epoch')
            if self.add_text is True:
                for ax in [ax1,ax3]:
                    ax.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar is True:
                if self.pm_separate is True:
                    for figure in [fig, fig1]:
                        cbar_ax = figure.add_axes([1.6, 0.16, 0.05, 0.7])
                        cbar = figure.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                        cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                        cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                        figure.delaxes(cbar_ax)
                else:
                    cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                    cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                    cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                    cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                    fig.delaxes(cbar_ax)
            
        except:
            fig = plt.figure(figsize=(11, 5.5))
            ax1 = fig.add_axes((0.10, 0.1, 0.35, 0.77))
            ax2 = fig.add_axes((0.60, 0.1, 0.35, 0.77))

            # plot the num_lines randomly selected curves
            for i in range(self.num_lines):
                ax1.plot(self.epoch, self.mualp_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                ax2.plot(self.epoch, self.mudec_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the most likely one
            ax1.plot(self.epoch, self.mualp_ml, color='black')
            ax2.plot(self.epoch, self.mudec_ml, color='black')

            # manually change the x tick labels from JD to calendar years
            epoch_ticks = np.linspace(self.start_epoch, self.end_epoch, 5)
            epoch_label = np.zeros(len(epoch_ticks))
            for i in range(len(epoch_ticks)):
                epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))

            ax1.set_xlim(self.start_epoch, self.end_epoch)
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.set_xticks(epoch_ticks)
            ax1.set_xticklabels([str(int(i)) for i in epoch_label])
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax1.set_title(self.title)
            ax1.set_xlabel('date (yr)')
            ax1.set_ylabel(r'$\Delta \mu_{\alpha}$ (mas/yr)',labelpad = 1)

            ax2.set_xlim(self.start_epoch, self.end_epoch)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.set_xticks(epoch_ticks)
            ax2.set_xticklabels([str(int(i)) for i in epoch_label])
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax2.set_title(self.title)
            ax2.set_xlabel('date (yr)')
            ax2.set_ylabel(r'$\Delta \mu_{\delta}$ (mas/yr)')
        
        plt.tight_layout()
        print("Plotting Proper Motions, your plot is generated at " + self.outputdir)
        
        if self.pm_separate is True:
            fig.savefig(os.path.join(self.outputdir, 'ProperMotions_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
            fig1.savefig(os.path.join(self.outputdir, 'ProperMotions1_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
        else:
            fig.savefig(os.path.join(self.outputdir, 'Proper_Motions_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
################################################################################################



# 7. Corner plot
# Finalized
    ################################################################################################
    ############### plot a nicer corner plot###############
    
    def plot_corner(self, title_fmt=".3f", **kwargs):
        labels=[r'$\mathrm{M_{pri}\, (M_{\odot})}$', r'$\mathrm{M_{sec}\, (M_{Jup})}$', 'a (AU)', r'$\mathrm{e\, (^{\circ})}$', r'$\mathrm{i\, (^{\circ})}$']
        rcParams["lines.linewidth"] = 1.0
        rcParams["axes.labelpad"] = 80.0
        rcParams["xtick.labelsize"] = 10.0
        rcParams["ytick.labelsize"] = 10.0
        
        burnin = self.burnin
        chain = self.chain
        ndim = chain[:,burnin:,0].flatten().shape[0]
        Mpri = chain[:,burnin:,1].flatten().reshape(ndim,1)                      # in M_{\odot}
        Msec = (chain[:,burnin:,2]*1989/1.898).flatten().reshape(ndim,1)         # in M_{jup}
        Semimajor = chain[:,burnin:,3].flatten().reshape(ndim,1)                       # in AU
        Ecc = (chain[:,burnin:,4]**2 +chain[:,burnin:,5]**2).flatten().reshape(ndim,1)
        #Omega=(np.arctan2(chain[:,burnin:,4],chain[:,burnin:,5])).flatten().reshape(ndim,1)
        Inc = (chain[:,burnin:,6]*180/np.pi).flatten().reshape(ndim,1)
        
        chain =np.hstack([Mpri,Msec,Semimajor,Ecc,Inc])
        
        figure = corner.corner(chain, labels=labels, quantiles=[0.16, 0.5, 0.84], verbose=False, show_titles=True, title_kwargs={"fontsize": 12}, hist_kwargs={"lw":1.}, label_kwargs={"fontsize":15}, xlabcord=(0.5,-0.45), ylabcord=(-0.45,0.5), title_fmt=title_fmt,**kwargs)

        print("Plotting corner plot, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'Corner_' + self.title)+'.pdf')

###################################################################################################
# end of code
