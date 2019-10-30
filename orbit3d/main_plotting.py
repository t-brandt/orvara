#!/usr/bin/env python
"""
Plotting code. The run function is the console entry point,
accessed by calling plot_orbit from the command line.

Example:

plot_orbit --output-dir ./orbit3d --config-file config.ini

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
from orbit3d import orbit_plots          #import yunlin's plotting package


def initialize_plot_options(config):
    """
        initialize the user defined plotting options from config.ini
    """

    target = config.get('plotting', 'target')
    nplanets = config.getint('mcmc_settings', 'nplanets')
    RVFile = config.get('data_paths', 'RVFile')
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    AstrometryFile = config.get('data_paths', 'AstrometryFile')
    GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback=None)
    Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback=None)
    Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback=None)
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)
     
    # colorbar settings
    plot_colorbar = config.getboolean('plotting', 'colorbar', fallback=False)
    color_map = config.get('plotting', 'colormap', fallback= 'viridis')
     
    # custom range of epochs
    start_ep = config.getfloat('plotting', 'start_epoch', fallback=0)
    end_ep = config.getfloat('plotting', 'end_epoch', fallback=0)
     
    # how many orbits
    num_orbits = config.getint('plotting', 'num_orbits', fallback = 50)
    cm_ref = config.get('plotting', 'reference')
    
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

    #read in the mcmc chains
    MCMCFile = os.path.join(args.output_dir, 'HIP%d_chain%03d.fits' % (HipID, 1))
    source = MCMCFile.split('_')[0]
    tt, lnp, extras = [fits.open(MCMCFile)[i].data for i in range(3)]
    nsteps = 50*tt.shape[1]
    beststep = np.where(lnp==lnp.max())
   
    # initialize the OP object
    OP = orbit_plots.OrbitPlots(target, HipID, (start_ep, end_ep), cm_ref, num_orbits, color_map, MCMCFile, RVFile, AstrometryFile, use_epoch_astrometry, args.output_dir)

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

    return OP , hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter



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
    
    # initialize the OP class object
    OPs, H1f, H2f, Gf = initialize_plot_options(config)
    
    # which plot
    plot_settings = {}
    burnin = config.getint('plotting', 'burnin', fallback=0)
    plot_corner = config.getboolean('plotting', 'Corner_plot', fallback=False)
    plot_astr = config.getboolean('plotting', 'Astrometry_orbits_plot', fallback=False)
    plot_rv = config.getboolean('plotting', 'RV_orbits_plot', fallback=False)
    plot_rel_rv = config.getboolean('plotting', 'Relative_RV_plot', fallback=False)
    plot_rel_sep = config.getboolean('plotting', 'Relative_separation_plot', fallback=False)
    plot_position_angle = config.getboolean('plotting', 'Position_angle_plot', fallback=False)
    plot_proper_motions = config.getboolean('plotting', 'Proper_motion_plot', fallback=False)

   
    if plot_corner is True:
        OPs.plot_corner(burnin)
        
    if plot_astr is True:
        OPs.astrometry()
    
    if plot_rv is True:
        Ops.RV()
    
    if plot_rel_rv is True:
        OPs.relRV()
    
    if plot_rel_sep is True:
        OPs.relsep()
        
    if plot_position_angle is True:
        OPs.PA()
    
    if plot_proper_motions is True:
        OPs.proper_motions()
    
if __name__ == "__main__":
    run()
