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

"""
    Example:
    OPs = orbit_plots.OrbitPlots(title, Hip, start_ep, end_ep, cmref, num_lines, cm_name, burnin, mcmcfile, RVfile, AstrometryFile, HGCAFile, outputdir):
    OPs.astrometry()
    OPs.RV()
    OPs.relRV()
    OPs.relsep()
    OPs.PA()
    OPs.proper_motions()
    OPs.plot_corner()
"""

class OrbitPlots:

    ############################## Initialize Class ############################################
    def __init__(self, title, Hip, start_ep, end_ep, cmref, num_lines, cm_name, burnin, mcmcfile, RVfile, AstrometryFile, HGCAFile, outputdir):

        self.title = title
        self.Hip = Hip
        self.start_epoch = start_ep
        self.end_epoch = end_ep
        self.cmref = cmref
        self.num_lines = num_lines
        self.cm_name = cm_name
        self.burnin = burnin
        self.MCMCfile = mcmcfile
        self.RVfile = RVfile
        self.relAstfile = AstrometryFile
        self.HGCAFile = HGCAFile
        self.outputdir = outputdir
        
        
        self.cmlabel_dic = {'msec': r'$\mathrm{M_{comp} (M_\odot)}$', 'ecc': 'Eccentricity'}
        self.color_list = ['r', 'g', 'y', 'm', 'c', 'b']
        
        
        ################ load in data #####################
        # define epochs
        self.epoch, self.epoch_JD = self.define_epochs()
        # load mcmc data
        self.tt, self.beststep, self.extras = self.load_mcmc_data()
        # load observed RV data
        self.epoch_obs, self.RV_obs, self.RV_obs_err, self.nInst, self.epoch_obs_dic, self.RV_obs_dic, self.RV_obs_err_dic = self.load_obsRV_data()
        # load relative astrometry data:
        self.epochrA_obs, self.relsep_obs, self.relsep_obs_err, self.PA_obs, self.PA_obs_err = self.load_relAst_data()
        # load HGCA data:
        self.ep_mualp_obs, self.ep_mudec_obs, self.mualp_obs, self.mudec_obs, self.mualp_obs_err, self.mudec_obs_err = self.load_HGCA_data()
        
        ############## calculate orbits ###################
        # might be better to calculate the best fit orbit by passing an argument in the random_orbits(tt[beststep])
        
        # calcualte the best fit orbit
        self.data, self.plx, self.dras_ml, self.ddecs_ml, self.RV_ml, self.mu_RA_ml, self.mu_Dec_ml, self.relsep_ml, self.mualp_ml, self.mudec_ml, self.f_RVml, self.f_relsepml, self.f_mualpml, self.f_mudecml, self.PA_ml, self.TA_ml, self.node0, self.node1, self.pras = self.bestfit_orbit()  # plx in units of arcsecs
        # calculate more random orbits drawn from the mcmc chian
        self.RV_dic, self.dras_dic, self.ddecs_dic, self.relsep_dic, self.PA_dic, self.mualp_dic, self.mudec_dic, self.dic_keys, self.RV_dic_vals, self.dras_dic_vals, self.ddecs_dic_vals, self.relsep_dic_vals, self.PA_dic_vals, self.mualp_dic_vals, self.mudec_dic_vals = self.random_orbits()
        
        ################ set colorbar #####################
        # setup the normalization and the colormap
        self.nValues = np.array(self.dic_keys)
        self.normalize = mcolors.Normalize(vmin=self.nValues.min(), vmax=self.nValues.max())
        self.colormap = getattr(cm, self.cm_name)
        # setup the colorbar
        self.sm = cm.ScalarMappable(norm=self.normalize, cmap=self.colormap)
        self.sm.set_array(self.nValues)

        print("Generating plots for target " + self.title)

    ################################## Define Functions ########################################
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

    def chi_sqr(self, offset, f, epoch_obs, data_obs, data_obs_err): #changed
        """
            A chi-square function for fitting
        """
        chi_sqr = 0
        for i in range(len(epoch_obs)):
            chi_sqr += (f(epoch_obs[i]) - data_obs[i] - offset)**2 / data_obs_err[i]**2
        return chi_sqr

    def define_epochs(self): # need to switch self.epoch and self.epoch_JD
        """
            Function to define a custom range of epochs
        """
        start_epoch = self.calendar_to_JD(self.start_epoch)
        end_epoch = self.calendar_to_JD(self.end_epoch)
        range_epoch = end_epoch - start_epoch
        epoch = np.linspace(start_epoch - 0.1*range_epoch, end_epoch + 0.5*range_epoch, 1500)
        epoch_JD = np.zeros(len(epoch))
        for i in range(len(epoch_JD)):
            epoch_JD[i] = self.JD_to_calendar(epoch[i])
        return epoch, epoch_JD

    def load_mcmc_data(self):
        """
            Function to load in the MCMC chain from fit_orbit
        """
        source = self.MCMCfile.split('_')[0]
        tt, lnp, extras = [fits.open(self.MCMCfile)[i].data for i in range(3)]
        beststep = np.where(lnp == lnp.max())
        return tt, beststep, extras
        
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
            
            epochrA_obs = reldat[:, 0]
            for i in range(len(epochrA_obs)):
                epochrA_obs[i] = self.calendar_to_JD(epochrA_obs[i])
            relsep_obs = reldat[:, 1]
            relsep_obs_err = reldat[:, 2]
            PA_obs = reldat[:, 3]
            PA_obs_err = reldat[:, 4]
        except:
            self.have_reldat = False
        return epochrA_obs, relsep_obs, relsep_obs_err, PA_obs, PA_obs_err
    
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
    
    def calc_RV_offset(self, f_RVml):
        """
            Function to calculate the offset of the observed RV data
        """
        try:
            # calculate the differences of offsets (del_offset), and shift the data according to
            # del_offset, which normalizes the offset for different instruments.
            assert self.multi_instr
            offset_dic = {}
            del_offset_dic = {}
            for i in range(self.nInst):
                result = op.minimize(self.chi_sqr, 70, args=(f_RVml, self.epoch_obs_dic[i], self.RV_obs_dic[i], self.RV_obs_err_dic[i]))
                offset_dic[i] = result['x']
            offset = min(offset_dic.values())[0]
            for i in range(self.nInst):
                del_offset_dic[i] = offset_dic[i] - offset
                self.RV_obs_dic[i] += del_offset_dic[i]
        except:
            result = op.minimize(self.chi_sqr, 70, args=(f_RVml, self.epoch_obs, self.RV_obs, self.RV_obs_err))   # 70 is initial guess
            offset = result['x']
        return offset
    
    def bestfit_orbit(self):
        """
            Function to calculate the most likely orbit
        """
        data = orbit.Data(self.Hip, self.RVfile, self.relAstfile)
        par = orbit.Params(self.tt[self.beststep][0]) # beststep is the best fit orbit
        data.custom_epochs(self.epoch)
        model = orbit.Model(data)

        orbit.calc_EA_RPP(data, par, model)
        orbit.calc_offsets(data, par, model, 0)
        orbit.calc_RV(data, par, model)
        
        # most likely orbit for delRA and delDec
        dRAs_G, dDecs_G, dRAs_H1, dDecs_H1, dRAs_H2, dDecs_H2 = model.return_dRA_dDec()
        ratio = -(1. + par.mpri/par.msec)*data.plx
        dras_ml, ddecs_ml = ratio*dRAs_G, ratio*dDecs_G       #only considering Gaia
        # most likely orbit for RV
        RV_ml = model.return_RVs()
        # most likely orbit for proper motions
        mu_RA_ml, mu_Dec_ml =  model.return_proper_motions(par)
        # most likely orbit for relative separation
        relsep_ml = model.return_relsep()*data.plx             # relsep in arcsec
        # most likely orbit for position angle
        PA_ml = model.return_PAs()*180/np.pi
        for i in range(len(PA_ml)):
            if PA_ml[i] <= 0:
                PA_ml[i] += 360
        # most likely orbit for proper motions
        mualp_ml, mudec_ml = model.return_proper_motions(par)
        mualp_ml, mudec_ml = 1e3*data.plx*365.25*mualp_ml, 1e3*data.plx*365.25*mudec_ml   # convert from arcsec/day to mas/yr
        mualp_ml += self.extras[self.beststep[0], self.beststep[1], 1]*1000
        mudec_ml += self.extras[self.beststep[0], self.beststep[1], 2]*1000
        
        # interpolation
        f_RVml = interp1d(self.epoch, RV_ml, fill_value="extrapolate")
        f_relsepml = interp1d(self.epoch, relsep_ml, fill_value="extrapolate")
        f_mualpml = interp1d(self.epoch, mualp_ml, fill_value="extrapolate")
        f_mudecml = interp1d(self.epoch, mudec_ml, fill_value="extrapolate")
        
        # redefine RV_ml
        offset = self.calc_RV_offset(f_RVml)
        RV_ml -= offset
        
        # find the positions of nodes and periastron
        TA_ml = model.return_TAs(par)
        # When par.arg (omega) is negative, it means the periastron is below the plane
        # of the sky. We can set omega = -par.arg, which is the angle symmetric with
        # respect to the plane of the sky. Although this agnle doesn't mean anything, the
        # same algorithm below can be applied to locate the position of nodes.
        omega = abs(par.arg)
        
        idx_node0 = np.where(abs(TA_ml - (np.pi - omega)) == min(abs(TA_ml - (np.pi - omega))))[0]
        idx_node1 = np.where(abs(TA_ml - (-omega)) == min(abs(TA_ml - (-omega))))[0]
        node0 = (dras_ml[idx_node0], ddecs_ml[idx_node0])
        node1 = (dras_ml[idx_node1], ddecs_ml[idx_node1])
        idx_pras = np.where(abs(TA_ml) == min(abs(TA_ml)))[0]
        pras = (dras_ml[idx_pras], ddecs_ml[idx_pras])
        
        return data, data.plx, dras_ml, ddecs_ml, RV_ml, mu_RA_ml, mu_Dec_ml, relsep_ml, mualp_ml, mudec_ml, f_RVml, f_relsepml, f_mualpml, f_mudecml, PA_ml, TA_ml, node0, node1, pras

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
            walker_idx = randrange(self.tt.shape[0])
            step_idx = randrange(self.burnin, self.tt.shape[1])
            par = orbit.Params(self.tt[walker_idx, step_idx])

            # calculate and assign variables
            data = self.data
            data.custom_epochs(self.epoch)
            model = orbit.Model(data)

            orbit.calc_EA_RPP(data, par, model)
            orbit.calc_offsets(data, par, model, 0)
            orbit.calc_RV(data, par, model)
            dRAs_G, dDecs_G, dRAs_H1, dDecs_H1, dRAs_H2,  dDecs_H2 = model.return_dRA_dDec()
            ratio = -(1. + par.mpri/par.msec)*self.plx
            dras, ddecs = ratio*dRAs_G, ratio*dDecs_G # convert from AU to mas
            RV = model.return_RVs()
            relsep = model.return_relsep()*self.plx
            PA = model.return_PAs()*180/np.pi
            for i in range(len(PA)):
                if PA[i] <= 0:
                    PA[i] += 360
            mualp, mudec = model.return_proper_motions(par)
            mualp, mudec = 1e3*self.plx*365.25*mualp, 1e3*self.plx*365.25*mudec   # convert from arcsec/day to mas/yr
            
            # shift the RV curve wrt data points
            midep_RV_obs = (self.epoch_obs[0] + self.epoch_obs[-1])/2
            RV_ref_val = self.f_RVml(midep_RV_obs)
            f_RV = interp1d(self.epoch, RV, fill_value="extrapolate")
            offset = self.calc_RV_offset(f_RV)
            RV -= offset + (f_RV(midep_RV_obs) - RV_ref_val)
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


    ################################## Plotting ################################################

    ############### plot astrometric orbit ###############

    def astrometry(self):
        #rcParams["axes.labelpad"] = 10.0
        
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        # plot the num_lines randomly selected curves
        for i in range(self.num_lines):
            ax.plot(self.dras_dic_vals[i], self.ddecs_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

        # plot the most likely one
        ax.plot(self.dras_ml, self.ddecs_ml, color='black')

        # plot the relAst data points
        f_drasml = interp1d(self.epoch, self.dras_ml)
        f_ddecsml = interp1d(self.epoch, self.ddecs_ml)

        try:
            assert self.have_reldat == True
            for i in range(len(self.epochrA_obs)):
                ra_exp, dec_exp = f_drasml(self.epochrA_obs[i]), f_ddecsml(self.epochrA_obs[i])
                relsep_exp = self.f_relsepml(self.epochrA_obs[i])
                ra_obs = ra_exp * self.relsep_obs[i] / relsep_exp    # similar triangles
                dec_obs = dec_exp * self.relsep_obs[i] / relsep_exp
                ax.scatter(ra_obs, dec_obs, s=45, facecolors='coral', edgecolors='none', zorder=99)
                ax.scatter(ra_obs, dec_obs, s=45, facecolors='none', edgecolors='k', zorder=100)
        except:
            pass

        # plot the 5 predicted positions of the companion star from 1990 to 2030
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        epoch_int = []
        for year in self.epoch_JD:
            epoch_int.append(int(year))

        text_x_max = 0
        text_x_min = 0
        text_y_max = 0
        text_y_min = 0
        for year in [1990, 2000, 2010, 2020, 2030]:
            idx = epoch_int.index(year)
            x = self.dras_ml[idx]
            y = self.ddecs_ml[idx]
            r = np.sqrt(x**2 + y**2)
            ax.scatter(x, y, s=55, facecolors='none', edgecolors='k', zorder=100)

            # avoid overlapping between text and plot
            if x >= 0 and y >= 0:
                text_x = x*(r + 0.12*(x1 - x0))/r
                text_y = y*(r + 0.03*(y1 - y0))/r
                ax.text(text_x, text_y, str(year), fontsize=10)
            elif x >= 0 and y <= 0:
                text_x = x*(r + 0.12*(x1 - x0))/r
                text_y = y*(r + 0.05*(y1 - y0))/r
                ax.text(text_x, text_y, str(year), fontsize=10)
            elif x <= 0 and y >= 0:
                text_x = x*(r + 0.03*(x1 - x0))/r
                text_y = y*(r + 0.03*(y1 - y0))/r
                ax.text(text_x, text_y, str(year), fontsize=10)
            elif x <= 0 and y <= 0:
                text_x = x*(r + 0.03*(x1 - x0))/r
                text_y = y*(r + 0.05*(y1 - y0))/r
                ax.text(text_x, text_y, str(year), fontsize=10)

            if text_x > text_x_max:
                text_x_max = text_x
            if text_x < text_x_min:
                text_x_min = text_x
            if text_y > text_y_max:
                text_y_max = text_y
            if text_y < text_y_min:
                text_y_min = text_y

        # avoid the text exceeding the box
        if abs(text_x_min - x0) < 0.05*(x1 - x0):
            x0 -= 0.10*(x1 - x0)
        if abs(text_x_max - x1) < 0.05*(x1 - x0):
            x1 += 0.10*(x1 - x0)
        if abs(text_y_min - y0) < 0.05*(y1 - y0):
            y0 -= 0.05*(y1 - y0)
        if abs(text_y_max - y1) < 0.05*(y1 - y0):
            y1 += 0.05*(y1 - y0)

        # plot line of nodes, periastron and the direction of motion of companion star, and label the host star
        ax.plot([self.node0[0], self.node1[0]], [self.node0[1], self.node1[1]], 'k--')
        ax.plot([0, self.pras[0]], [0, self.pras[1]], 'k:')
        arrow = mpatches.FancyArrowPatch((self.dras_ml[self.idx_pras][0], self.ddecs_ml[self.idx_pras][0]), (self.dras_ml[self.idx_pras+1][0], self.ddecs_ml[self.idx_pras+1][0]), arrowstyle='->', mutation_scale=25, zorder=100)
        ax.add_patch(arrow)
        ax.plot(0, 0, marker='*',  color='black')

        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect(abs((x1-x0)/(y1-y0)))
        fig.colorbar(self.sm, ax=ax, label=self.cmlabel_dic[self.cmref])
        ax.invert_xaxis()
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax.set_xlabel(r'$\Delta \alpha$ [arcsec]')
        ax.set_ylabel(r'$\Delta \delta$ [arcsec]')
        ax.set_title(self.title + ' Astrometric Orbits')
        plt.tight_layout()
        print("Plotting Astrometry orbits, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir,'astrometric_orbit_' + self.title))


    ############### plot the RV orbits ###############
    def RV(self):

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        # plot the 50 randomly selected curves
        for i in range(self.num_lines):
            ax.plot(self.epoch_JD, self.RV_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

        # plot the most likely one
        ax.plot(self.epoch_JD, self.RV_ml, color='black')

        # plot the observed data points (RV & relAst)
        try:
            assert self.multi_instr
            for i in range(self.nInst):
                epoch_obs_Inst = np.zeros(len(self.epoch_obs_dic[i]))
                for j in range(len(self.epoch_obs_dic[i])):
                    epoch_obs_Inst[j] = self.JD_to_calendar(self.epoch_obs_dic[i][j])
                ax.plot(epoch_obs_Inst, self.RV_obs_dic[i], self.color_list[i]+'o', markersize=2)
        except:
            ax.plot(self.JD_to_calendar(self.epoch_obs), self.RV_obs, 'ro', markersize=2)

        ax.set_xlim(self.start_epoch, self.end_epoch)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        fig.colorbar(self.sm, ax=ax, label=self.cmlabel_dic[self.cmref])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax.set_xlabel('Date (yr)')
        ax.set_ylabel('RV (m/s)')
        ax.set_title(self.title + ' RV Orbits')
        print("Plotting RV orbits, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'RV_orbit_' + self.title))

    ############### plot the relative RV and O-C ###############

    def relRV(self):

        fig = plt.figure(figsize=(5, 6))
        ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
        ax2 = fig.add_axes((0.15, 0.1, 0.8, 0.15))

        # plot the 50 randomly selected curves
        self.f_RVml = interp1d(self.epoch, self.RV_ml)
        RV_OC = self.RV_dic_vals

        for i in range(self.num_lines):
            ax1.plot(self.epoch, self.RV_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
            for j in range(len(self.epoch)):
                RV_OC[i][j] -= self.f_RVml(self.epoch[j])
            ax2.plot(self.epoch, RV_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

        # plot the most likely one
        ax1.plot(self.epoch, self.RV_ml, color='black')
        ax2.plot(self.epoch, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

        # plot the observed data points
        datOC_list = []
        try:
            assert self.multi_instr
            for i in range(self.nInst):
                ax1.errorbar(self.epoch_obs_dic[i], self.RV_obs_dic[i], yerr=self.RV_obs_err_dic[i], fmt=self.color_list[i]+'o', ecolor='black', capsize=3)
                ax1.scatter(self.epoch_obs_dic[i], self.RV_obs_dic[i], s=45, facecolors='none', edgecolors='k', zorder=100, alpha=0.5)
                for j in range(len(self.epoch_obs_dic[i])):
                    OC = self.RV_obs_dic[i][j] - self.f_RVml(self.epoch_obs_dic[i][j])
                    datOC_list.append(OC)
                    ax2.errorbar(self.epoch_obs_dic[i][j], OC, yerr=self.RV_obs_err_dic[i][j], fmt=self.color_list[i]+'o', ecolor='black', capsize=3)
                    ax2.scatter(self.epoch_obs_dic[i][j], OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=0.5)
        except:
            ax1.errorbar(self.epoch_obs, self.RV_obs, yerr=self.RV_obs_err, fmt='bo', ecolor='black', capsize=3)
            ax1.scatter(self.epoch_obs, self.RV_obs, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            for i in range(len(self.epoch_obs)):
                OC = self.RV_obs[i] - self.f_RVml(self.epoch_obs[i])
                datOC_list.append(OC)
                ax2.errorbar(self.epoch_obs[i], OC, yerr=self.RV_obs_err[i], fmt='bo', ecolor='black', capsize=3)
                ax2.scatter(self.epoch_obs[i], OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)

        # manually change the x tick labels from JD to calendar years
        epoch_ticks = np.linspace(self.epoch_obs[0], self.epoch_obs[-1], 5)
        epoch_label = np.zeros(len(epoch_ticks))
        for i in range(len(epoch_ticks)):
            epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))

        range_ep_obs = max(self.epoch_obs) - min(self.epoch_obs)
        range_RV_obs = max(self.RV_obs) - min(self.RV_obs)
        ax1.set_xlim(min(self.epoch_obs) - range_ep_obs/20., max(self.epoch_obs) + range_ep_obs/20.)
        ax1.set_ylim(min(self.RV_obs) - range_RV_obs/10., max(self.RV_obs) + range_RV_obs/10.)
        ax1.xaxis.set_major_formatter(NullFormatter())
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax1.set_title(self.title)
        ax1.set_ylabel('Relative RV (m/s)')

        range_datOC = max(datOC_list) - min(datOC_list)
        ax2.set_xlim(min(self.epoch_obs) - range_ep_obs/20., max(self.epoch_obs) + range_ep_obs/20.)
        ax2.set_ylim(min(datOC_list) - range_datOC/5., max(datOC_list) + range_datOC/5.)
        ax2.set_xticks(epoch_ticks)
        ax2.set_xticklabels([str(int(i)) for i in epoch_label])
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax2.set_xlabel('Epoch (yr)')
        ax2.set_ylabel('O-C')
        
        print("Plotting relative RV, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'relRV_OC_' + self.title))


    ############### plot the seperation and O-C ###############

    def relsep(self):

        try:
            assert self.have_reldat == True
            fig = plt.figure(figsize=(5, 6))
            ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
            ax2 = fig.add_axes((0.15, 0.1, 0.8, 0.15))

            # plot the 50 randomly selected curves
            relsep_OC = self.relsep_dic_vals

            for i in range(self.num_lines):
                ax1.plot(self.epoch, self.relsep_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                for j in range(len(self.epoch)):
                    relsep_OC[i][j] -= self.f_relsepml(self.epoch[j])
                ax2.plot(self.epoch, relsep_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the most likely one
            ax1.plot(self.epoch, self.relsep_ml, color='black')
            ax2.plot(self.epoch, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

            # plot the observed data points
            datOC_list = []
            ax1.errorbar(self.epochrA_obs, self.relsep_obs, yerr=self.relsep_obs_err, color='coral', fmt='o', ecolor='black', capsize=3)
            ax1.scatter(self.epochrA_obs, self.relsep_obs, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            for i in range(len(self.epochrA_obs)):
                dat_OC = self.relsep_obs[i] - self.f_relsepml(self.epochrA_obs[i])
                datOC_list.append(dat_OC)
                ax2.errorbar(self.epochrA_obs[i], dat_OC, yerr=self.relsep_obs_err[i], color='coral', fmt='o', ecolor='black', capsize=3)
                ax2.scatter(self.epochrA_obs[i], dat_OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)

            # manually change the x tick labels from JD to calendar years
            epoch_ticks = np.linspace(self.epochrA_obs[0], self.epochrA_obs[-1], 5)
            epoch_label = np.zeros(len(epoch_ticks))
            for i in range(len(epoch_ticks)):
                epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))

            self.range_eprA_obs = max(self.epochrA_obs) - min(self.epochrA_obs)
            range_relsep_obs = max(self.relsep_obs)  - min(self.relsep_obs)
            ax1.set_xlim(min(self.epochrA_obs) - self.range_eprA_obs/8., max(self.epochrA_obs) + self.range_eprA_obs/8.)
            ax1.set_ylim(min(self.relsep_obs) - range_relsep_obs/2., max(self.relsep_obs) + range_relsep_obs/2.)
            ax1.xaxis.set_major_formatter(NullFormatter())
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_title(self.title)
            ax1.set_ylabel('Seperation (arcsec)')

            range_datOC = max(datOC_list) - min(datOC_list)
            ax2.set_xlim(min(self.epochrA_obs) - self.range_eprA_obs/8., max(self.epochrA_obs) + self.range_eprA_obs/8.)
            ax2.set_ylim(min(datOC_list) - range_datOC, max(datOC_list) + range_datOC)
            ax2.set_xticks(epoch_ticks)
            ax2.set_xticklabels([str(int(i)) for i in epoch_label])
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (yr)')
            ax2.set_ylabel('O-C')

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

            ax.set_xlim(self.start_epoch, self.end_epoch)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.set_xticks(epoch_ticks)
            ax.set_xticklabels([str(int(i)) for i in epoch_label])
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax.set_title(self.title)
            ax.set_ylabel('Seperation (arcsec)')
        
        print("Plotting Separation, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'relsep_OC_' + self.title))


    ############### plot the position angle and O-C ###############

    def PA(self):

        try:
            assert self.have_reldat == True
            fig = plt.figure(figsize=(5, 6))
            ax1 = fig.add_axes((0.15, 0.3, 0.8, 0.6))
            ax2 = fig.add_axes((0.15, 0.1, 0.8, 0.15))

            # plot the 50 randomly selected curves
            f_PAml = interp1d(self.epoch, self.PA_ml)
            PA_OC = self.PA_dic_vals

            for i in range(self.num_lines):
                ax1.plot(self.epoch, self.PA_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                for j in range(len(self.epoch)):
                    PA_OC[i][j] -= f_PAml(self.epoch[j])
                ax2.plot(self.epoch, PA_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the most likely one
            ax1.plot(self.epoch, self.PA_ml, color='black')
            ax2.plot(self.epoch, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

            # plot the observed data points
            datOC_list = []
            ax1.errorbar(self.epochrA_obs, self.PA_obs, yerr=self.PA_obs_err, color='coral', fmt='o', ecolor='black', capsize=3)
            ax1.scatter(self.epochrA_obs, self.PA_obs, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            for i in range(len(self.epochrA_obs)):
                dat_OC = self.PA_obs[i] - f_PAml(self.epochrA_obs[i])
                datOC_list.append(dat_OC)
                ax2.errorbar(self.epochrA_obs[i], dat_OC, yerr=self.PA_obs_err[i], color='coral', fmt='o', ecolor='black', capsize=3)
                ax2.scatter(self.epochrA_obs[i], dat_OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)

            # manually change the x tick labels from JD to calendar years
            epoch_ticks = np.linspace(self.epochrA_obs[0], self.epochrA_obs[-1], 5)
            epoch_label = np.zeros(len(epoch_ticks))
            for i in range(len(epoch_ticks)):
                epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))

            self.range_eprA_obs = max(self.epochrA_obs) - min(self.epochrA_obs)
            range_PA_obs = max(self.PA_obs)  - min(self.PA_obs)
            ax1.set_xlim(min(self.epochrA_obs) - self.range_eprA_obs/8., max(self.epochrA_obs) + self.range_eprA_obs/8.)
            ax1.set_ylim(min(self.PA_obs) - range_PA_obs/5., max(self.PA_obs) + range_PA_obs/5.)
            ax1.xaxis.set_major_formatter(NullFormatter())
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_title(self.title)
            ax1.set_ylabel(r'Position Angle ($^{\circ}$)')

            range_datOC = max(datOC_list) - min(datOC_list)
            ax2.set_xlim(min(self.epochrA_obs) - self.range_eprA_obs/8., max(self.epochrA_obs) + self.range_eprA_obs/8.)
            ax2.set_ylim(min(datOC_list) - range_datOC, max(datOC_list) + range_datOC)
            ax2.set_xticks(epoch_ticks)
            ax2.set_xticklabels([str(int(i)) for i in epoch_label])
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (yr)')
            ax2.set_ylabel('O-C')

        except:
            fig = plt.figure(figsize=(5, 5))
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
            ax.set_title(self.title)
            ax.set_ylabel(r'Position Angle ($^{\circ}$)')

        print("Plotting Position Angle, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir,'PA_OC_' + self.title))
   
    def plot_corner(self, title_fmt=".4f", **kwargs):
        labels=[r'$\mathrm{M_{pri}}$', r'$\mathrm{M_{sec}}$', 'Sep', 'Ecc', 'Inc']
        rcParams["lines.linewidth"] = 1.0
        rcParams["axes.labelpad"] = 80.0
        rcParams["xtick.labelsize"] = 10.0
        rcParams["ytick.labelsize"] = 10.0
        
        burnin = self.burnin
        tt = self.tt
        ndim = tt[:,burnin:,0].flatten().shape[0]
        Mpri = tt[:,burnin:,1].flatten().reshape(ndim,1)                      # in M_{\odot}
        Msec = (tt[:,burnin:,2]*1989/1.898).flatten().reshape(ndim,1)         # in M_{jup}
        Sep = tt[:,burnin:,3].flatten().reshape(ndim,1)                       # in AU
        Ecc = (tt[:,burnin:,4]**2 +tt[:,burnin:,5]**2).flatten().reshape(ndim,1)
        #Omega=(np.arctan2(tt[:,burnin:,4],tt[:,burnin:,5])).flatten().reshape(ndim,1)
        Inc = (tt[:,burnin:,6]*180/np.pi).flatten().reshape(ndim,1)
        
        chain =np.hstack([Mpri,Msec,Sep,Ecc,Inc])
        
        figure = corner.corner(chain, labels=labels, quantiles=[0.16, 0.5, 0.84], verbose=False, show_titles=True, title_kwargs={"fontsize": 14}, hist_kwargs={"lw":1.}, label_kwargs={"fontsize":14}, xlabcord=(0.5,-0.45), ylabcord=(-0.45,0.5), title_fmt=title_fmt,**kwargs)

        print("Plotting corner plot, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'Corner_' + self.title))
        
        
    def proper_motions(self):

        try:
            assert self.have_pmdat == True
            fig = plt.figure(figsize=(11, 6))
            ax1 = fig.add_axes((0.10, 0.30, 0.35, 0.60))
            ax2 = fig.add_axes((0.10, 0.10, 0.35, 0.15))
            ax3 = fig.add_axes((0.55, 0.30, 0.35, 0.60))
            ax4 = fig.add_axes((0.55, 0.10, 0.35, 0.15))

            # plot the num_lines randomly selected curves
            mualp_OC = self.mualp_dic_vals
            mudec_OC = self.mudec_dic_vals

            for i in range(self.num_lines):
                ax1.plot(self.epoch, self.mualp_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                ax3.plot(self.epoch, self.mudec_dic_vals[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                for j in range(len(self.epoch)):
                    mualp_OC[i][j] -= self.f_mualpml(self.epoch[j])
                    mudec_OC[i][j] -= self.f_mudecml(self.epoch[j])
                ax2.plot(self.epoch, mualp_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)
                ax4.plot(self.epoch, mudec_OC[i], color=self.colormap(self.normalize(self.nValues[i])), alpha=0.3)

            # plot the most likely one
            ax1.plot(self.epoch, self.mualp_ml, color='black')
            ax2.plot(self.epoch, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))
            ax3.plot(self.epoch, self.mudec_ml, color='black')
            ax4.plot(self.epoch, np.zeros(len(self.epoch)), 'k--', dashes=(5, 5))

            # plot the observed data points
            mualpdatOC_list = []
            mudecdatOC_list = []
            ax1.errorbar(self.ep_mualp_obs, self.mualp_obs, yerr=self.mualp_obs_err, color='coral', fmt='o', ecolor='black', capsize=3, zorder=99)
            ax1.scatter(self.ep_mualp_obs, self.mualp_obs, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            ax3.errorbar(self.ep_mudec_obs, self.mudec_obs, yerr=self.mualp_obs_err, color='coral', fmt='o', ecolor='black', capsize=3, zorder=99)
            ax3.scatter(self.ep_mudec_obs, self.mudec_obs, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            for i in range(len(self.ep_mualp_obs)):
                dat_OC = self.mualp_obs[i] - self.f_mualpml(self.ep_mualp_obs[i])
                mualpdatOC_list.append(dat_OC)
            ax2.errorbar(self.ep_mualp_obs, mualpdatOC_list, yerr=self.mualp_obs_err, color='coral', fmt='o', ecolor='black', capsize=3, zorder=99)
            ax2.scatter(self.ep_mualp_obs, mualpdatOC_list, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)
            for i in range(len(self.ep_mudec_obs)):
                dat_OC = self.mudec_obs[i] - self.f_mudecml(self.ep_mudec_obs[i])
                mudecdatOC_list.append(dat_OC)
            ax4.errorbar(self.ep_mudec_obs, mudecdatOC_list, yerr=self.mudec_obs_err, color='coral', fmt='o', ecolor='black', capsize=3, zorder=99)
            ax4.scatter(self.ep_mudec_obs, mudecdatOC_list, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=1)

            # manually change the x tick labels from JD to calendar years
            epoch_ticks = np.linspace(self.ep_mualp_obs[0], self.ep_mualp_obs[-1], 5)
            epoch_label = np.zeros(len(epoch_ticks))
            for i in range(len(epoch_ticks)):
                epoch_label[i] = round(self.JD_to_calendar(epoch_ticks[i]))

            self.range_eppm_obs = max(self.ep_mualp_obs) - min(self.ep_mualp_obs)
            range_mualp_obs = max(self.mualp_obs)  - min(self.mualp_obs)
            ax1.set_xlim(min(self.ep_mualp_obs) - self.range_eppm_obs/8., max(self.ep_mualp_obs) + self.range_eppm_obs/8.)
            ax1.set_ylim(min(self.mualp_obs) - range_mualp_obs/5., max(self.mualp_obs) + range_mualp_obs/5.)
            ax1.xaxis.set_major_formatter(NullFormatter())
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax1.set_title(self.title)
            ax1.set_ylabel(r'$\Delta \mu_{\alpha}$ (mas/yr)')

            range_mudec_obs = max(self.mudec_obs)  - min(self.mudec_obs)
            ax3.set_ylabel(r'$\Delta \mu_{\alpha}$ mas/yr')
            ax3.set_xlim(min(self.ep_mudec_obs) - self.range_eppm_obs/8., max(self.ep_mudec_obs) + self.range_eppm_obs/8.)
            ax3.set_ylim(min(self.mudec_obs) - range_mudec_obs/5., max(self.mudec_obs) + range_mudec_obs/5.)
            ax3.xaxis.set_major_formatter(NullFormatter())
            ax3.xaxis.set_minor_locator(AutoMinorLocator())
            ax3.yaxis.set_minor_locator(AutoMinorLocator())
            ax3.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax3.set_title(self.title)
            ax3.set_ylabel(r'$\Delta \mu_{\delta}$ (mas/yr)')

            range_mualpdatOC = max(mualpdatOC_list) - min(mualpdatOC_list)
            ax2.set_xlim(min(self.ep_mualp_obs) - self.range_eppm_obs/8., max(self.ep_mualp_obs) + self.range_eppm_obs/8.)
            ax2.set_ylim(min(mualpdatOC_list) - range_mualpdatOC, max(mualpdatOC_list) + range_mualpdatOC)
            ax2.set_xticks(epoch_ticks)
            ax2.set_xticklabels([str(int(i)) for i in epoch_label])
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (yr)')
            ax2.set_ylabel('O-C')

            range_mudecdatOC = max(mudecdatOC_list) - min(mudecdatOC_list)
            ax4.set_xlim(min(self.ep_mudec_obs) - self.range_eppm_obs/8., max(self.ep_mudec_obs) + self.range_eppm_obs/8.)
            ax4.set_ylim(min(mudecdatOC_list) - range_mudecdatOC, max(mudecdatOC_list) + range_mudecdatOC)
            ax4.set_xticks(epoch_ticks)
            ax4.set_xticklabels([str(int(i)) for i in epoch_label])
            ax4.xaxis.set_minor_locator(AutoMinorLocator())
            ax4.yaxis.set_minor_locator(AutoMinorLocator())
            ax4.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax4.set_xlabel('Epoch (yr)')
            ax4.set_ylabel('O-C')
            
        except:
            fig = plt.figure(figsize=(11, 5))
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
            ax1.set_title(self.title)
            ax1.set_xlabel('date (yr)')
            ax1.set_ylabel(r'$\Delta \mu_{\alpha}$ (mas/yr)')

            ax2.set_xlim(self.start_epoch, self.end_epoch)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.set_xticks(epoch_ticks)
            ax2.set_xticklabels([str(int(i)) for i in epoch_label])
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_title(self.title)
            ax2.set_xlabel('date (yr)')
            ax2.set_ylabel(r'$\Delta \mu_{\delta}$ (mas/yr)')
        print("Plotting Proper Motions, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'ProperMotions_' + self.title))

        def properMotion_mualp(self):
            pass

        def properMotion_mudec(self):
            pass


    def test(self):
        print("Testing code and debugging, test your code here")
        #plt.plot(self.epoch_JD, self.RV_ml, color='black')
        #plt.show()
        #epoch_obs_each_ins = np.zeros(len(self.epoch_obs_dic[0]))
        #for i in range(len(self.epoch_obs_dic[0])):
        #    epoch_obs_each_ins[i] = self.JD_to_calendar(self.epoch_obs_dic[0][i])
        
        print("exit")
        
