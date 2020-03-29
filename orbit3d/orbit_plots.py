import os
import numpy as np
from random import randrange
from orbit3d import corner_modified
from scipy.interpolate import interp1d
import scipy.optimize as op
from scipy import stats, signal
from orbit3d import orbit
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator

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
# Fix the excepts, or delete them?


class Orbit:

    def __init__(self, OP, step='best', epochs='custom'):

        data = orbit.Data(OP.Hip, OP.RVfile, OP.relAstfile, verbose=False)
        
        if isinstance(epochs, list) or isinstance(epochs, np.ndarray):
            data.custom_epochs(epochs)
        elif epochs == 'custom':
            data.custom_epochs(OP.epoch)
        elif epochs == 'observed':
            pass
        else:
            raise ValueError("epochs should be 'custom', 'observed', or array-like")
        model = orbit.Model(data)

        if step == 'best':
            idx = int(np.amin(np.where(OP.lnp == np.amax(OP.lnp))))
        elif step == int(step):
            idx = step
            
        self.plx = OP.extras[idx, 0]
        for i in range(OP.nplanets):
            par = orbit.Params(OP.chain[idx], iplanet=i, nplanets=OP.nplanets)
            if i == OP.iplanet:
                self.par = par
            orbit.calc_EA_RPP(data, par, model)
            orbit.calc_offsets(data, par, model, OP.iplanet)
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

        self.mu_RA_CM = 1e3*OP.extras[idx, 1]
        self.mu_Dec_CM = 1e3*OP.extras[idx, 2]
        
        self.mu_RA = 1e3*(self.mu_RA*self.plx*365.25 + OP.extras[idx, 1])
        self.mu_Dec = 1e3*(self.mu_Dec*self.plx*365.25 + OP.extras[idx, 2])
        self.offset = OP.calc_RV_offset(idx)

        if OP.cmref == 'msec':
            self.colorpar = self.par.msec*1989/1.898
        elif OP.cmref == 'ecc':
            self.colorpar = self.par.ecc
        
class OrbitPlots:

###################################### Initialize Class ############################################

    def __init__(self):
        pass        

    def start(self):

        self.cmlabel_dic = {'msec': r'$\mathrm{M_{comp} (M_{Jup})}$', 'ecc': 'Eccentricity'}
        self.color_list = ['r', 'g', 'b', 'y', 'c', 'b']
        
        ############################### load in data #######################
        # define epochs
        self.epoch, self.epoch_calendar = self.define_epochs()
        # load mcmc data
        self.chain, self.beststep, self.extras = self.load_mcmc_data()

        self.rand_idx = [randrange(self.chain.shape[0]) for i in range(self.num_orbits)]
        
        # load observed RV data
        self.epoch_obs, self.RV_obs, self.RV_obs_err, self.nInst, self.epoch_obs_dic, self.RV_obs_dic, self.RV_obs_err_dic = self.load_obsRV_data()
        # load relative astrometry data:
        if os.access(self.relAstfile,os.R_OK):
            self.have_reldat = True
            self.ep_relAst_obs, self.relsep_obs, self.relsep_obs_err, self.PA_obs, self.PA_obs_err = self.load_relAst_data()
        else:
            self.have_reldat = False
        # load HGCA data:
        self.ep_mualp_obs, self.ep_mudec_obs, self.mualp_obs, self.mudec_obs, self.mualp_obs_err, self.mudec_obs_err = self.load_HGCA_data()        
        
        ################################ set colorbar ######################
        # setup the normalization and the colormap

        if self.cmref == 'msec':
            vmin = 1989/1.898*stats.scoreatpercentile(self.chain[:, 2], 1)
            vmax = 1989/1.898*stats.scoreatpercentile(self.chain[:, 2], 99)
        elif self.cmref == 'ecc':
            ecc = self.chain[:, 4]**2 + self.chain[:, 5]**2
            vmin = stats.scoreatpercentile(ecc, 1)
            vmax = stats.scoreatpercentile(ecc, 99)
        else:
            raise ValueError("Reference parameter for color should be either 'msec' or 'ecc'")
            
        self.normalize = mcolors.Normalize(vmin=vmin, vmax=vmax)
        self.colormap = getattr(cm, self.color_map)
        self.sm = cm.ScalarMappable(norm=self.normalize, cmap=self.colormap)
        self.sm.set_array(np.linspace(vmin, vmax, 10))
        
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
        epoch = np.linspace(start_epoch - 0.1*range_epoch, end_epoch + 0.5*range_epoch, self.num_steps)
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
        chain = chain[:, self.burnin:, :].reshape(-1, chain.shape[-1])
        self.lnp = lnp[:, self.burnin:].flatten()
        extras = extras[:, self.burnin:, :].reshape(-1, extras.shape[-1])
        beststep = np.where(self.lnp == self.lnp.max())
        self.nplanets = chain.shape[-1]//7
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
            self.RVinst = RVinst
            # Check to see that the column we loaded was an integer
            assert np.all(RVinst == rvdat[:, 3])
            nInst = int(np.amax(rvdat[:, 3]) + 1)
            
            self.multi_instr = True
        except:
            self.multi_instr = False
            nInst = 1
            self.RVinst = (RV_obs*0).astype(int)
    
        idx_dic = {}
        epoch_obs_dic = {}
        RV_obs_dic = {}
        RV_obs_err_dic = {}
        
        for i in range(nInst):
            idx_dic[i] = (np.where(self.RVinst == i)[0][0], np.where(self.RVinst == i)[0][-1])
            epoch_obs_dic[i] = epoch_obs[idx_dic[i][0]: idx_dic[i][-1] + 1]
            RV_obs_dic[i] = RV_obs[idx_dic[i][0]: idx_dic[i][-1] + 1]
            RV_obs_err_dic[i] = RV_obs_err[idx_dic[i][0]: idx_dic[i][-1] + 1]
        return epoch_obs, RV_obs, RV_obs_err, nInst, epoch_obs_dic, RV_obs_dic, RV_obs_err_dic

    def load_relAst_data(self):
        """
            Function to load in the relative astrometry data
        """
        reldat = np.genfromtxt(self.relAstfile, usecols=(1,2,3,4,5), skip_header=1)

        ep_relAst_obs = reldat[:, 0]
        for i in range(len(ep_relAst_obs)):
            ep_relAst_obs[i] = self.calendar_to_JD(ep_relAst_obs[i])
            relsep_obs = reldat[:, 1]
            relsep_obs_err = reldat[:, 2]
            PA_obs = reldat[:, 3]
            PA_obs_err = reldat[:, 4]
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
    
    def calc_RV_offset(self, idx):
        """
            Function to calculate the offset of the observed RV data
        """
        try:
            # calculate the offsets of the RV curves
            assert self.multi_instr
            offset_dic = {}
            for i in np.arange(8, 8 + self.nInst, 1):
                offset = self.extras[idx, i]
                offset_dic[i-8] = offset
        except:
            offset = offset_dic = [self.extras[idx, 8]]
        return offset_dic



########################################## Plotting ################################################

# 1. Astrometric plots
# Finalized
# scipy seems to have issues finding the tangent line of an ecplise of any shape

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
        
        dra = (B*X + G*Y)*(-par.sau)*plx
        ddec = (A*X + F*Y)*(-par.sau)*plx
        
        return dra, ddec

    def astrometry(self):

        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111)

        # plot the num_orbits randomly selected curves
        for i in range(self.num_orbits):

            orb = Orbit(self, step=self.rand_idx[i])
            dra, ddec = self.closed_orbit(orb.par, orb.plx)

            ax.plot(dra, ddec, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.4, linewidth=0.8)

        # plot the most likely orbit
        
        orb_ml = Orbit(self, step='best')
        dra, ddec = self.closed_orbit(orb_ml.par, orb_ml.plx)
        #dd = ddec
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
                ax.annotate('  ' + str(year), xy=(x, y), verticalalignment='bottom', horizontalalignment='left',rotation =-angle)
            # if the point you are labeling has x value larger than the corresponding x value on the semimajor axis (same y), and has y value smaller than the corresponding y value on the semiminor axis (same x), then it's in the fourth quardret and the labels should be aligned to the left rotated by +angle
            elif x < (y-b_semimajor)/m_semimajor and y < m_semiminor*x + b_semiminor:
                ax.annotate('  ' + str(year), xy=(x, y), verticalalignment='top', horizontalalignment='left',rotation =angle)
            # if the point you are labeling has x value smaller than the corresponding x value on the semimajor axis (same y), and has y value larger than the corresponding y value on the semiminor axis (same x), then it's in the second quardret and the labels should be aligned to the right and rotated by -angle
            elif x > (y-b_semimajor)/m_semimajor and y < m_semiminor*x + b_semiminor:
                ax.annotate(str(year)+'  ', xy=(x, y), verticalalignment='top', horizontalalignment='right',rotation= -angle)
            # if the point you are labeling has x value smaller than the corresponding x value on the semimajor axis (same y), and has y value smaller than the corresponding y value on the semiminor axis (same x), then it's in the third quardret and the labels should be aligned to the right and rotated by +angle
            elif x > (y-b_semimajor)/m_semimajor and y > m_semiminor*x + b_semiminor:
                ax.annotate(str(year)+'  ', xy=(x, y), verticalalignment='bottom', horizontalalignment='right', rotation=angle)

        # plot line of nodes, periastron and the direction of motion of the companion
        dra_nd, ddec_nd = self.closed_orbit(orb_ml.par, orb_ml.plx, nodes=True)
        
        ax.plot(dra_nd[2:], ddec_nd[2:], 'k--',linewidth = 1)
        ax.plot([0, dra_nd[1]], [0, ddec_nd[1]], 'k:')
        # add arrow
        arrow = mpatches.FancyArrowPatch((dra_nd[0], ddec_nd[0]), (dra_nd[1], ddec_nd[1]),  arrowstyle='->', mutation_scale=25, zorder=100)

        ax.add_patch(arrow)
        ax.plot(0, 0, marker='*',  color='black', markersize=10)
        
        # from advanced plotting settings in config.ini
        if self.set_limit:
            ax.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
            ax.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            aspect0 = np.float(self.user_xlim[1])-np.float(self.user_xlim[0])
            aspect1 = np.float(self.user_ylim[1])-np.float(self.user_ylim[0])
            ax.set_aspect(np.abs(aspect0/aspect1))
        if self.show_title:
            ax.set_title('Astrometric Orbits')
        if self.add_text:
            ax.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
        if self.usecolorbar:
            if self.colorbar_size < 0 or self.colorbar_pad < 0:
                cbar = fig.colorbar(self.sm, ax=ax, fraction=0.046, pad=0.04)
            else:
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
        rv_epoch_list = []
        for i in range(self.nInst):
            epoch_obs_Inst = np.zeros(len(self.epoch_obs_dic[i]))
            for j in range(len(self.epoch_obs_dic[i])):
                epoch_obs_Inst[j] = self.JD_to_calendar(self.epoch_obs_dic[i][j])
            rv_epoch_list.append(epoch_obs_Inst)

        for i in range(self.nInst):
            ax.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], yerr=self.RV_obs_err_dic[i], fmt=self.color_list[i]+'o', ecolor='black', alpha = 0.8, zorder = 299)
            ax.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], facecolors='none', edgecolors='k', alpha = 0.8, zorder=300)
            
        ax.set_xlim(self.start_epoch, self.end_epoch)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        
        if self.set_limit:
            ax.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
            ax.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
        if self.show_title:
            ax.set_title('Relative RV vs. Epoch')
        if self.add_text:
            ax.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
        if self.usecolorbar:
            cbar_ax = fig.add_axes([1.3, 0.325, 0.03, 0.55])
            if self.colorbar_size < 0 or self.colorbar_pad < 0:
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
            else:
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12) # default value = 12
            cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
            if self.colorbar_size < 0 or self.colorbar_pad < 0:
                cbar.ax.get_yaxis().labelpad=20
            else:
                cbar.ax.get_yaxis().labelpad=self.colorbar_pad
            fig.delaxes(cbar_ax)
            
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax.set_xlabel('Epoch (year)', labelpad = 10, fontsize=13)
        ax.set_ylabel('RV (m/s)', fontsize=13)

        plt.tight_layout()
        print("Plotting RV orbits, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'RV_orbit_' + self.title)+'.pdf',bbox_inches='tight', dpi=200) # or +'.pdf'


# 3. relative RV plot
# FINALIZED
    ################################################################################################
    ############### plot the relative RV and O-C #####################

    def relRV(self):
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
        
        #if self.whichInst == np.str('All'):
        #    print('You have chosen to plot relative RV for all the Instruments')
        all_OC = []
        all_OC_err = []
        
        for i in range(self.nInst):
            plot_this = True
            if not self.whichInst == np.str('All'):
                plot_this = False
                whichInst = np.int(self.whichInst)
                if i + 1 == whichInst and i < self.nInst:
                    plot_this = True
            if not plot_this:
                continue
            
            ax1.errorbar(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], yerr=self.RV_obs_err_dic[i], fmt=self.color_list[i]+'o', ecolor='black', capsize=3, alpha = 0.8, zorder=199+i)#, ecolor='black', markersize = 1, elinewidth = 0.3, capsize=1, capthick = 0.3, zorder = 200+i, alpha = 0.8)
            ax1.scatter(rv_epoch_list[i], self.RV_obs_dic[i] + orb_ml.offset[i], s=45, facecolors='none', edgecolors='k', zorder=200+i, alpha = 0.8)
            
            OC = self.RV_obs_dic[i] + orb_ml_obs.offset[i] - orb_ml_obs.RV[self.RVinst == i]
            y_err = self.RV_obs_err_dic[i]
            all_OC += list(OC)
            all_OC_err += list(y_err)
            
            ax2.errorbar(rv_epoch_list[i], OC, yerr=y_err, fmt=self.color_list[i]+'o', ecolor='black', capsize=3)
            ax2.scatter(rv_epoch_list[i], OC, s=45, facecolors='none', edgecolors='k', zorder=100, alpha=0.5)

            #else:
            #    print('ValueError: Please enter a valid instrument number between 1 and '+np.str(self.nInst)+ ' or enter All to plot the observed data points from all Instruments')
            #    raise SystemExit
        
        # axes settings
        # ax1
        ax1.get_shared_x_axes().join(ax1, ax2)
        range_ep_obs = max(all_RV_eps) - min(all_RV_eps)
        RV_obs_max = max([max(self.RV_obs_dic[i] + orb_ml.offset[i]) for i in range(self.nInst)])
        RV_obs_min = min([min(self.RV_obs_dic[i] + orb_ml.offset[i]) for i in range(self.nInst)])
        range_RV_obs = RV_obs_max - RV_obs_min
        
        ax1.set_xlim(min(all_RV_eps) - range_ep_obs/20., max(all_RV_eps) + range_ep_obs/20.)
        ax1.set_ylim(RV_obs_min - range_RV_obs/10., RV_obs_max + range_RV_obs/10.)
        ax1.xaxis.set_major_formatter(NullFormatter())
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax1.set_ylabel('Relative RV (m/s)', fontsize=13)
        # ax2

        min_OC = np.amin(np.asarray(all_OC) - np.asarray(all_OC_err))
        max_OC = np.amax(np.asarray(all_OC) + np.asarray(all_OC_err))
        min_OC = min(min_OC, -max_OC)
        max_OC = max(max_OC, -min_OC)
        range_OC = max_OC - min_OC
        
        ax2.set_ylim(min_OC - range_OC/7., max_OC + range_OC/7.)

        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())

        ax2.set_xticks(ax2.get_xticks()[::2])
            
        ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
        ax2.set_xlabel('Epoch (yr)', fontsize=13)
        ax2.set_ylabel('O-C', fontsize=13)
        
        # from advanced plotting settings in config.ini
        if self.set_limit:
            ax2.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
            ax2.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
        if self.show_title:
            ax1.set_title('Relative RV vs. Epoch')
        if self.add_text:
            ax1.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
        if self.usecolorbar:
            cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
            if self.colorbar_size < 0 or self.colorbar_pad < 0:
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
            else:
                cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
            cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
            if self.colorbar_size < 0 or self.colorbar_pad < 0:
                cbar.ax.get_yaxis().labelpad=20
            else:
                cbar.ax.get_yaxis().labelpad=self.colorbar_pad
            fig.delaxes(cbar_ax)

        print("Plotting relative RV, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'relRV_OC_' + self.title + '_Inst' + np.str(self.whichInst) +'.pdf'),bbox_inches='tight', dpi=200)
################################################################################################



# 4. relative separation plot
# FINALIZED
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
            ax1.errorbar(ep_relAst_obs_calendar, self.relsep_obs, yerr=self.relsep_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax1.scatter(ep_relAst_obs_calendar, self.relsep_obs, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)

            dat_OC = self.relsep_obs - orb_ml.relsep
            
            ax2.errorbar(ep_relAst_obs_calendar, dat_OC, yerr=self.relsep_obs_err[i], color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax2.scatter(ep_relAst_obs_calendar, dat_OC, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
                
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
            range_datOC = max(dat_OC) - min(dat_OC)
            if np.abs(min(dat_OC)) <= np.abs(max(dat_OC)):
                ax2.set_ylim(-np.abs(max(dat_OC)) - range_datOC, max(dat_OC) + range_datOC)
            elif np.abs(min(dat_OC)) > np.abs(max(dat_OC)):
                ax2.set_ylim(min(dat_OC) - range_datOC, np.abs(min(dat_OC)) + range_datOC)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
            ax2.set_ylabel('O-C', labelpad=-2, fontsize=13)
            
            # from advanced plotting settings in config.ini
            if self.set_limit:
                ax2.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                ax2.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title:
                ax1.set_title('Relative Separation vs. Epoch')
            if self.add_text:
                ax1.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar:
                cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                if self.colorbar_size < 0 or self.colorbar_pad < 0:
                    cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                else:
                    cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                if self.colorbar_size < 0 or self.colorbar_pad < 0:
                    cbar.ax.get_yaxis().labelpad=20
                else:
                    cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
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
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax.set_ylabel('Seperation (arcsec)', fontsize=13)
            ax.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
                
        plt.tight_layout()
        print("Plotting Separation, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'relsep_OC_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
################################################################################################



# 5. position angle plot
# FINALIZED
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
            ax1.errorbar(ep_relAst_obs_calendar, self.PA_obs, yerr=self.PA_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax1.scatter(ep_relAst_obs_calendar, self.PA_obs, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
            dat_OC = self.PA_obs - orb_ml.PA
            
            ax2.errorbar(ep_relAst_obs_calendar, dat_OC, yerr=self.PA_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize=5, zorder=299)
            ax2.scatter(ep_relAst_obs_calendar, dat_OC, s=45, facecolors='none', edgecolors='k', alpha=1, zorder=300)
            
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
            range_datOC = max(dat_OC) - min(dat_OC)
            if np.abs(min(dat_OC)) <= np.abs(max(dat_OC)):
                ax2.set_ylim(-np.abs(max(dat_OC)) - range_datOC, max(dat_OC) + range_datOC)
            elif np.abs(min(dat_OC)) > np.abs(max(dat_OC)):
                ax2.set_ylim(min(dat_OC) - range_datOC, np.abs(min(dat_OC)) + range_datOC)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            ax2.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
            ax2.set_ylabel('O-C', labelpad=17, fontsize=13)

            # from advanced plotting settings in config.ini
            if self.set_limit:
                ax2.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                ax2.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title:
                ax1.set_title('Position angle vs. Epoch')
            if self.add_text:
                ax1.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar:
                cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                if self.colorbar_size < 0 or self.colorbar_pad < 0:
                    cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                else:
                    cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                if self.colorbar_size < 0 or self.colorbar_pad < 0:
                    cbar.ax.get_yaxis().labelpad=20
                else:
                    cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
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
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax.set_title(self.title)
            ax.set_ylabel(r'Position Angle ($^{\circ}$)', fontsize=13)
            ax.set_xlabel('Epoch (year)', labelpad=6, fontsize=13)
        
        plt.tight_layout()
        print("Plotting Position Angle, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir,'PA_OC_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
################################################################################################




# 6. Proper motion plots
# FINALIZED
    ################################################################################################
    ############### plot the proper motions, together or separately ###############
    
    def proper_motions(self):
        if self.have_pmdat:
            if self.pm_separate:
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

                
            # set up the observed data points
            mualpdatOC_list = []
            mudecdatOC_list = []
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


            t1_RA = min(ep_mualp_obs_calendar)
            t2_RA = max(ep_mualp_obs_calendar)
            dt_RA = t2_RA - t1_RA
            
            t1_Dec = min(ep_mudec_obs_calendar)
            t2_Dec = max(ep_mudec_obs_calendar)
            dt_Dec = t2_Dec - t1_Dec
            
            ax1.set_xlim([t1_RA - dt_RA/8., t2_RA + dt_RA/8.])
            ax3.set_xlim([t1_Dec - dt_Dec/8., t2_Dec + dt_Dec/8.])
            ax2.set_xlim([t1_RA - dt_RA/8., t2_RA + dt_RA/8.])
            ax4.set_xlim([t1_Dec - dt_Dec/8., t2_Dec + dt_Dec/8.])
                            
            orb_ml = Orbit(self, 'best')

            minT = min(t1_RA - dt_RA/5., t1_Dec - dt_Dec/5.)
            maxT = max(t2_RA + dt_RA/5., t2_Dec + dt_Dec/5.)            
            T = np.asarray(self.epoch_calendar)
            indx = (T > minT)*(T < maxT) > 0
            
            for i in range(self.num_orbits):

                orb = Orbit(self, step=self.rand_idx[i]) 
                
                ax1.plot(T[indx], orb.mu_RA[indx], color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax3.plot(T[indx], orb.mu_Dec[indx], color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax2.plot(T[indx], orb.mu_RA[indx] - orb_ml.mu_RA[indx], color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax4.plot(T[indx], orb.mu_Dec[indx] - orb_ml.mu_Dec[indx], color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                
            # plot the highest likelihood curve

            ax1.plot(T[indx], orb_ml.mu_RA[indx], color='black')
            ax2.plot(self.epoch_calendar, np.zeros(len(self.epoch_calendar)), 'k--', dashes=(5, 5))
            ax3.plot(T[indx], orb_ml.mu_Dec[indx], color='black')
            ax4.plot(self.epoch_calendar, np.zeros(len(self.epoch_calendar)), 'k--', dashes=(5, 5))
            
            ax1.errorbar(ep_mualp_obs_calendar, self.mualp_obs, yerr=self.mualp_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax1.scatter(ep_mualp_obs_calendar, self.mualp_obs, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)
            ax3.errorbar(ep_mudec_obs_calendar, self.mudec_obs, yerr=self.mualp_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax3.scatter(ep_mudec_obs_calendar, self.mudec_obs, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)
            
            # plot the O-C for observed data points
            dat_OC_RA = self.mualp_obs - orb_ml_RA.mu_RA
            ax2.errorbar(ep_mualp_obs_calendar, dat_OC_RA, yerr=self.mualp_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax2.scatter(ep_mualp_obs_calendar, dat_OC_RA, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)

            dat_OC_Dec = self.mudec_obs - orb_ml_Dec.mu_Dec
            ax4.errorbar(ep_mudec_obs_calendar, dat_OC_Dec, yerr=self.mudec_obs_err, color=self.marker_color, fmt='o', ecolor='black', capsize=3, markersize = 5, zorder = 299)
            ax4.scatter(ep_mudec_obs_calendar, dat_OC_Dec, s=45, facecolors='none', edgecolors='k', zorder = 300, alpha=1)

            # axes settings
            # ax1
            ax1.get_shared_x_axes().join(ax1, ax2)
            ax2.set_xlim([t1_RA - dt_RA/8., t2_RA + dt_RA/8.])

            ax1.set_ylabel(r'$\Delta \mu_{\alpha}$ (mas/yr)', labelpad = 5, fontsize = 13)
            # ax3
            ax3.get_shared_x_axes().join(ax3, ax4)
            ax4.set_xlim([t1_Dec - dt_Dec/8., t2_Dec + dt_Dec/8.])
            range_mudec_obs = max(self.mudec_obs)  - min(self.mudec_obs)
            ax3.set_ylabel(r'$\Delta \mu_{\delta}$ (mas/yr)', labelpad = 6, fontsize = 13)
            for ax in [ax1, ax3]:
                ax.xaxis.set_major_formatter(NullFormatter())
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            # ax2
            
            y1_RA = max(np.abs(dat_OC_RA) + max(self.mualp_obs_err))
            y1_Dec = max(np.abs(dat_OC_Dec) + max(self.mudec_obs_err))
            
            ax2.set_ylim(-1.2*y1_RA, 1.2*y1_RA)
            ax4.set_ylim(-1.2*y1_Dec, 1.2*y1_Dec)
            
            for ax in [ax2,ax4]:
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
                ax.set_xlabel('Epoch (year)', labelpad = 6, fontsize = 13)
                ax.set_ylabel('O-C', labelpad = 8, fontsize = 13)
                
            # from advanced plotting settings in config.ini
            if self.set_limit:
                for ax in [ax1, ax3]:
                    ax.set_xlim(self.calendar_to_JD(np.float(self.user_xlim[0])), self.calendar_to_JD(np.float(self.user_xlim[1])))
                    ax.set_ylim(np.float(self.user_ylim[0]),np.float(self.user_ylim[1]))
            if self.show_title:
                ax1.set_title('Right Ascension vs. Epoch')
                ax3.set_title('Declination vs. Epoch')
            if self.add_text:
                for ax in [ax1,ax3]:
                    ax.text(self.calendar_to_JD(self.x_text),self.y_text, self.text_name, fontsize=15)
            if self.usecolorbar:
                if self.pm_separate:
                    for figure, name in [[fig, 'RA_'], [fig1, 'Dec_']]:
                        cbar_ax = figure.add_axes([1.6, 0.16, 0.05, 0.7])
                        if self.colorbar_size < 0 or self.colorbar_pad < 0:
                            cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                        else:
                            cbar = figure.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                        cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                        if self.colorbar_size < 0 or self.colorbar_pad < 0:
                            cbar.ax.get_yaxis().labelpad=20
                        else:
                            cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                        figure.delaxes(cbar_ax)
                        figure.savefig(os.path.join(self.outputdir, 'ProperMotions_' + name + self.title)+'.pdf',bbox_inches='tight', dpi=200)
                else:
                    cbar_ax = fig.add_axes([1.6, 0.16, 0.05, 0.7])
                    if self.colorbar_size < 0 or self.colorbar_pad < 0:
                        cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=12)
                    else:
                        cbar = fig.colorbar(self.sm, ax=cbar_ax, fraction=self.colorbar_size) # default value = 12
                    cbar.ax.set_ylabel(self.cmlabel_dic[self.cmref], rotation=270, fontsize=13)
                    if self.colorbar_size < 0 or self.colorbar_pad < 0:
                        cbar.ax.get_yaxis().labelpad=20
                    else:
                        cbar.ax.get_yaxis().labelpad=self.colorbar_pad # defult value = 20
                    fig.delaxes(cbar_ax)
            
        else:
            fig = plt.figure(figsize=(11, 5.5))
            ax1 = fig.add_axes((0.10, 0.1, 0.35, 0.77))
            ax2 = fig.add_axes((0.60, 0.1, 0.35, 0.77))

            # plot the num_orbits randomly selected curves
            orb_ml = Orbit(self, 'best')
            
            for i in range(self.num_orbits):
                orb = Orbit(self, step=self.rand_idx[i])

                ax1.plot(self.epoch_calendar, orb.mu_RA, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                ax2.plot(self.epoch_calendar, orb.mu_Dec, color=self.colormap(self.normalize(orb.colorpar)), alpha=0.3)
                
            # plot the most likely one
            ax1.plot(self.epoch_calendar, orb_ml.mu_RA, color='black')
            ax2.plot(self.epoch_calendar, orb_ml.mu_Dec, color='black')

            ax1.set_xlim(self.start_epoch, self.end_epoch)
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax1.set_title(self.title)
            ax1.set_xlabel('date (yr)')
            ax1.set_ylabel(r'$\Delta \mu_{\alpha}$ (mas/yr)',labelpad = 1)

            ax2.set_xlim(self.start_epoch, self.end_epoch)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(direction='in', which='both', left=True, right=True, bottom=True, top=True)
            #ax2.set_title(self.title)
            ax2.set_xlabel('date (yr)')
            ax2.set_ylabel(r'$\Delta \mu_{\delta}$ (mas/yr)')
        
        plt.tight_layout()
        print("Plotting Proper Motions, your plot is generated at " + self.outputdir)
        
        if self.pm_separate and not self.usecolorbar:
            fig.savefig(os.path.join(self.outputdir, 'ProperMotions_RA_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
            fig1.savefig(os.path.join(self.outputdir, 'ProperMotions_Dec_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
        else:
            fig.savefig(os.path.join(self.outputdir, 'Proper_Motions_' + self.title)+'.pdf',bbox_inches='tight', dpi=200)
################################################################################################


    def astrometric_prediction(self, nbins=500):

        # Fetch parameters as ndarrays
        par = self.chain.astype(float)
        mpri = par[:, 1]
        msec = par[:, 2 + 7*self.iplanet]
        sau = par[:, 3 + 7*self.iplanet]
        esino = par[:, 4 + 7*self.iplanet]
        ecoso = par[:, 5 + 7*self.iplanet]
        inc = par[:, 6 + 7*self.iplanet]
        asc = par[:, 7 + 7*self.iplanet]
        lam = par[:, 8 + 7*self.iplanet]
        arg = np.arctan2(esino, ecoso)
        ecc = esino**2 + ecoso**2

        plx = self.extras[:, 0]
        plx = np.reshape(plx, -1)

        # The date we want
        JD_predict = self.calendar_to_JD(self.predicted_ep_ast)

        data = orbit.Data(self.Hip, self.RVfile, self.relAstfile, verbose=False)

        # Solve Kepler's equation in array format given a different
        # eccentricity for each point.  This is the same Newton solver
        # used by radvel.
        
        period = np.sqrt(sau**3/(mpri + msec))*365.25
        MA = (2*np.pi/period*(JD_predict - data.refep) + lam - arg)%(2*np.pi)
        E = MA + np.sign(np.sin(MA))*0.85*ecc
        fi = E - ecc*np.sin(E) - MA

        for i in range(10):
            fip = 1 - ecc*np.cos(E)
            fipp = ecc*np.sin(E)
            fippp = 1 - fip
            d1 = -fi/fip
            d1 = -fi/(fip + d1*fipp/2.)
            d1 = -fi/(fip + d1*fipp/2. + d1**2*fippp/6.)
            E += d1
            fi = E - ecc*np.sin(E) - MA

        # Thiele-Innes constants -> relative separation.
        A = np.cos(arg)*np.cos(asc) - np.sin(arg)*np.sin(asc)*np.cos(inc)
        B = np.cos(arg)*np.sin(asc) + np.sin(arg)*np.cos(asc)*np.cos(inc)
        F = -np.sin(arg)*np.cos(asc) - np.cos(arg)*np.sin(asc)*np.cos(inc)
        G = -np.sin(arg)*np.sin(asc) + np.cos(arg)*np.cos(asc)*np.cos(inc)
        
        X = np.cos(E) - ecc
        Y = np.sin(E)*np.sqrt(1 - ecc**2)
        
        dra = (B*X + G*Y)*(-sau)*plx
        ddec = (A*X + F*Y)*(-sau)*plx

        # Set limits on plot to include basically all of the data.
        
        ramin = stats.scoreatpercentile(dra, 0.1)
        ramax = stats.scoreatpercentile(dra, 99.9)
        decmin = stats.scoreatpercentile(ddec, 0.1)
        decmax = stats.scoreatpercentile(ddec, 99.9)

        diff = max(ramax - ramin, decmax - decmin)*1.7
        xmin = 0.5*(ramin + ramax) - diff/2.
        xmax = xmin + diff
        ymin = 0.5*(decmin + decmax) - diff/2.
        ymax = ymin + diff

        x = np.linspace(xmin, xmax, nbins)
        y = np.linspace(ymin, ymax, nbins)

        # Bin it up, then smooth it.  More points -> less smoothing.
        
        dens = np.histogram2d(dra, ddec, bins=[x, y])[0].T
        
        _x, _y = np.mgrid[-20:21, -20:21]
        window = np.exp(-(_x**2 + _y**2)/20.*len(par)/len(x)**2)
        dens = signal.convolve2d(dens, window, mode='same')
        dens /= np.amax(dens)

        # Make one-, two-, and three-sigma contours.
        
        dens_sorted = np.sort(dens.flatten())
        cdf = np.cumsum(dens_sorted)/np.sum(dens_sorted)
        cdf_func = interp1d(cdf, dens_sorted)
        
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111)
        ax.imshow(dens[::-1], extent=(xmin, xmax, ymin, ymax),
                  interpolation='nearest', aspect=1, cmap=cm.hot_r)

        # Mark the central star if (0, 0) is within the axis limits
        if xmin*xmax < 0 and ymin*ymax < 0:
            ax.plot(0, 0, marker='*', markersize=15, color='c')

        x = 0.5*(x[1:] + x[:-1])
        y = 0.5*(y[1:] + y[:-1])
        levels = [cdf_func(p) for p in [1 - 0.9973, 1 - 0.954, 1 - 0.683]]
        ax.contour(x, y, dens, levels=levels, colors=['k', 'C0', 'b'])
        ax.set_xlabel(r'$\mathrm{\Delta \alpha}$ (arcsec)', fontsize=14)
        ax.set_ylabel(r'$\mathrm{\Delta \delta}$ [arcsec]', fontsize=14)
        ax.text(xmax - 3e-2*(xmax - xmin), ymax - 7e-2*(ymax - ymin),
                'Location of %s, %.1f' % (self.text_name, self.predicted_ep_ast), fontsize=14)
        ax.invert_xaxis()
        plt.tight_layout()
        plt.savefig(os.path.join(self.outputdir,'astrometric_prediction_' + self.title)+'.pdf')

# 7. Corner plot
# find out a way to fix the labels (2 sig figure for errors)...
    ################################################################################################
    ############### plot a nicer corner plot###############
    
    def plot_corner(self, **kwargs):
        labels=[r'$\mathrm{M_{pri}\, (M_{\odot})}$', r'$\mathrm{M_{sec}\, (M_{Jup})}$', 'a (AU)', r'$\mathrm{ecc}$', r'$\mathrm{i\, (^{\circ})}$']
        rcParams["lines.linewidth"] = 1.0
        rcParams["axes.labelpad"] = 80.0
        rcParams["xtick.labelsize"] = 10.0
        rcParams["ytick.labelsize"] = 10.0
        
        #burnin = self.burnin
        chain = self.chain
        ndim = chain[:, 0].flatten().shape[0]
        Mpri = chain[:, 1].flatten().reshape(ndim,1)                      # in M_{\odot}
        Msec = (chain[:,2]*1989/1.898).flatten().reshape(ndim,1)         # in M_{jup}
        Semimajor = chain[:,3].flatten().reshape(ndim,1)                       # in AU
        Ecc = (chain[:,4]**2 +chain[:,5]**2).flatten().reshape(ndim,1)
        #Omega=(np.arctan2(chain[:,burnin:,4],chain[:,burnin:,5])).flatten().reshape(ndim,1)
        Inc = (chain[:,6]*180/np.pi).flatten().reshape(ndim,1)
        
        chain =np.hstack([Mpri,Msec,Semimajor,Ecc,Inc])
        
        # in corner_modified, the error is modified to keep 2 significant figures
        figure = corner_modified.corner(chain, labels=labels, quantiles=[0.16, 0.5, 0.84], verbose=False, show_titles=True, title_kwargs={"fontsize": 12}, hist_kwargs={"lw":1.}, label_kwargs={"fontsize":15}, xlabcord=(0.5,-0.45), ylabcord=(-0.45,0.5),  **kwargs)

        print("Plotting Corner plot, your plot is generated at " + self.outputdir)
        plt.savefig(os.path.join(self.outputdir, 'Corner_' + self.title)+'.pdf')

###################################################################################################
# end of code
