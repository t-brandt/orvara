#!/usr/bin/env python
"""
Plotting code. The run function is the console entry point,
accessed by calling plot_orbit from the command line.
"""

from __future__ import print_function
import numpy as np
import pandas as pd
import os
import time
import emcee, corner
import scipy.optimize as op
from random import randrange
from scipy.interpolate import interp1d
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator
from configparser import ConfigParser
from htof.main import Astrometry
from orbit3d import orbit
#from orbit3d.config import parse_args
import argparse
from configparser import ConfigParser

def initialize_data(config):
    """
        initialize the ConfigParser object
    """
    # load in items from the ConfigParser object
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    RVFile = config.get('data_paths', 'RVFile')
    AstrometryFile = config.get('data_paths', 'AstrometryFile')
    GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback=None)
    Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback=None)
    Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback=None)
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)

    data = orbit.Data(HipID, RVFile, AstrometryFile)
    
    

    if use_epoch_astrometry:
        Gaia_fitter = Astrometry('GaiaDR2', '%06d' % (HipID), GaiaDataDir,
                                 central_epoch_ra=data.epRA_G,
                                 central_epoch_dec=data.epDec_G,
                                 central_epoch_fmt='frac_year')
        Hip2_fitter = Astrometry('Hip2', '%06d' % (HipID), Hip2DataDir,
                                 central_epoch_ra=data.epRA_H,
                                 central_epoch_dec=data.epDec_H,
                                 central_epoch_fmt='frac_year')
        Hip1_fitter = Astrometry('Hip1', '%06d' % (HipID), Hip1DataDir,
                                 central_epoch_ra=data.epRA_H,
                                 central_epoch_dec=data.epDec_H,
                                 central_epoch_fmt='frac_year')
        # instantiate C versions of the astrometric fitter which are must faster than HTOF's Astrometry
        hip1_fast_fitter = orbit.AstrometricFitter(Hip1_fitter)
        hip2_fast_fitter = orbit.AstrometricFitter(Hip2_fitter)
        gaia_fast_fitter = orbit.AstrometricFitter(Gaia_fitter)

        data = orbit.Data(HipID, RVFile, AstrometryFile, use_epoch_astrometry,
                          epochs_Hip1=Hip1_fitter.data.julian_day_epoch(),
                          epochs_Hip2=Hip2_fitter.data.julian_day_epoch(),
                          epochs_Gaia=Gaia_fitter.data.julian_day_epoch())
    else:
        hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter = None, None, None

    return data, hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter



#def  initialize_plot_options():
    
    
########## define some useful functions
def JD_to_calendar(JD):
    """
        Julian to calendar date conversion
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
    M = f - 1 - 12*int(f/14)
    Y = d - 4715 - int((7 + M)/10)
    year = Y + M/12 + D/365.25
    return year

def calendar_to_JD(year, M=1, D=1):
    """
        Calendar to Julian date conversion
    """
    if M <= 2:
        y = year - 1
        m = M + 12
    else:
        y = year
        m = M
    if year <= 1583:
        B = -2
    else:
        B = int(y/400) - int(y/100)
    UT = 0
    JD = int(365.25*y) + int(30.6001*(m+1)) + B + 1720996.5  + D + UT/24
    return JD


def chi_sqr(offset, epoch, obs, obs_err):
    """
        define a chi-square function for fitting
    """
    chi_sqr = 0
    for i in range(len(epoch)):
        chi_sqr += (func(epoch[i]) - obs[i] - offset)**2 / obs_err[i]**2
    return chi_sqr

# plot the proper motions

def plot_proper_motions(number_orbits , mu_RA_or_mu_Dec):

    rcParams['xtick.major.width']=1
    rcParams['xtick.major.size']=4
    rcParams['xtick.minor.width']=0.5
    rcParams['xtick.minor.size']=2
    rcParams['xtick.direction'] = "in"
    rcParams['ytick.direction'] = "in"

    fig = plt.figure(figsize=(5, 6))
    ax1 = fig.add_axes((0.15, 0.3, 0.7, 0.5))
    ax2 = fig.add_axes((0.15, 0.1, 0.7, 0.15))

    ratio = ((1. + params.mpri/params.msec))*1000.
    
    # plotting mu_RA or plotting mu_Dec
    if mu_RA_or_mu_Dec == 'mu_RA':
        # proper motions data points from paper
        mu = mu_RA
        x = [1991.0 , 2015.5]
        y = [0.6, -0.38]
        y_err = [0.49, 0.05]
        func = interp1d(ep_jd, mu_RA)
        offset_init_guess = -1.5
    else:
        mu = mu_Dec
        x = [1991.25 , 2015.75]
        y = [1.41, -0.90]
        y_err = [0.40, 0.]
        func = interp1d(ep_jd, mu_Dec)
        offset_init_guess = 1.5
    
    # shift the best fit curve to fit the data and calculate the offset with scipy
    optimize = op.minimize(chi_sqr, offset_init_guess, args=(x, y, y_err))
    offset_best = optimize['x']
    
    # plot the most likely proper motion curve
    ax1.plot(ep_jd, mu*ratio + offset_best,color= 'k', linewidth = 2, zorder = 100)
    
    # plot more orbits drawn from mcmc posteriors
    for i in range(number_orbits):

        i_walker = randrange(tt.shape[0])
        i_step = randrange(1000, tt.shape[1])
        params = orbit.Params(tt[i_walker, i_step])

        data.custom_epochs(ep)
        model = orbit.Model(data)

        orbit.calc_EA_RPP(data, params, model)
        orbit.calc_offsets(data, params, model, 0)
        orbit.calc_RV(data, params, model)
        mu_ra, mu_dec =  model.return_proper_motions(params)
        
        if mu_RA_or_mu_Dec == 'mu_RA':
            i_mu = mu_ra
        else:
            i_mu = mu_dec
           
        for i in range(len(ep)):
            ep_jd[i] = JD_to_calendar(ep[i])

        # calculating the offsets of each individual curve and shift the midpoints of the curves to align with the best fit curve
        d = {'epoch': ep_jd, 'mu': mu}
        df = pd.DataFrame(data=d)
        mid_pt = df[int(len(ep)/2):int(len(ep)/2)+1]['mu']
        
        i_d = {'epoch': ep_jd, 'i_mu': i_mu*ratio}
        i_df = pd.DataFrame(data=i_d)
        i_mid_pt = i_df[int(len(ep)/2):int(len(ep)/2)+1]['i_mu']
        
        i_offset = i_mid_pt[int(len(ep)/2)] - mid_pt[int(len(ep)/2)]
        #define the colormap
        cmap = plt.cm.cubehelix
        mu_y = i_mu*ratio - i_offset
        
        # plotting
        ax1.plot(ep_jd, mu_y, c= cmap((params.msec*1989/1.898 - 34.)/(42-34.)))
        ax2.plot(ep_jd, np.zeros(len(ep)), 'k--', dashes=(5, 5))
        for j in range(len(ep)):
                mu_y[j] -= (func(ep_jd[j])*ratio + offset_best)
        ax2.plot(ep_jd, mu_y, c =cmap((params.msec*1989/1.898 - 34.)/(42-34.)) , alpha=0.3)
        ax2.scatter(x, y - (func(x)*ratio + offset_best),zorder = 10000)

    ax1.errorbar(x, y, yerr= y_err ,fmt='o', ecolor='k', capthick=3,capsize=4,zorder=1000)
    ax1.set_ylabel(r'$\mathrm{\Delta}$' + mu_RA_or_mu_Dec+ r'$\mathrm{ \, (mas \, yr^{-1})}$')
    ax1.set_xlim(np.min(ep_jd),np.max(ep_jd))
    ax1.minorticks_on()
    #ax1.text(2010,1.5,'Gl 758B', fontsize =16)
    ax2.set_xlim(np.min(ep_jd),np.max(ep_jd))
    ax2.set_ylim(-2,2)
    ax2.set_xlabel("Epoch")
    ax2.set_ylabel("O-C")
    ax2.minorticks_on()
    plt.show()
    
def run():
    """
    Initialize and make plots
    """
    
    start_time = time.time()
    
    # this following function should be in the config.py file
    def parse_args():
        parser = argparse.ArgumentParser(description='Plot astrometry or RV orbits and other relavant plots. Required arguments are shown with [].')
        parser.add_argument("--output-dir", required=True,
                            help="Directory within which to save the plots.")
        parser.add_argument("--config-file", required=True,
                            help="Path to the planet-specific configuration file.")
        args = parser.parse_args()
        return args
    
    args = parse_args()
    config = ConfigParser()
    config.read(args.config_file)
    
    
    
    
    #nplanets = config.getint('mcmc_settings', 'nplanets')
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)
   
    #define plot settings
    plot_settings = {}
    plot_corner = config.getboolean('plotting', 'Corner_plot', fallback=False)
    plot_rv = config.getboolean('plotting', 'RV_orbits_plot', fallback=False)
    plot_astr = config.getboolean('plotting', 'Astrometry_orbits_plot', fallback=False)
    plot_rel_sep = config.getboolean('plotting', 'Relative_separation', fallback=False)
    plot_proper_motions = config.getboolean('plotting', 'Proper_motion_plot', fallback=False)
    plot_colorbar = config.getboolean('plotting', 'colorbar', fallback=False)
    #color_map = config.getboolean('plotting', 'colormap', fallback='plt.cm.cubehelix')
    
    #set arguments for plotting
    
  
  
    print(plot_corner)
  
    
    
    # initialize the data
    data, H1f, H2f, Gf = initialize_data(config)
    
    #read in the mcmc chains
    filename = os.path.join(args.output_dir, 'HIP%d_chain%03d.fits' % (HipID, 1))
    source = filename.split('_')[0]
    tt, lnp, extras = [fits.open(filename)[i].data for i in range(3)]
    nsteps = 50*tt.shape[1]
    beststep = np.where(lnp==lnp.max())
    
    # calculate the best fit or the most likely orbit
    params = orbit.Params(tt[beststep][0])
    ep = np.linspace(calendar_to_JD(1983), calendar_to_JD(2023) , 1000)
    data.custom_epochs(ep)
    model = orbit.Model(data)

    orbit.calc_EA_RPP(data, params, model)
    orbit.calc_offsets(data, params, model, 0)
    orbit.calc_RV(data, params, model)
    
    mu_RA, mu_Dec =  model.return_proper_motions(params)

    #plot_proper_motions(50 , 'mu_Dec')
    
    
if __name__ == "__main__":
    run()
