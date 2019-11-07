#!/usr/bin/env python
"""
Plotting code. The run function is the console entry point,
accessed by calling plot_orbit from the command line.

Example:

plot_orbit --output-dir ./Plots --config-file config_HD4747.ini

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
from orbit3d.config import parse_args_plotting
import argparse
from configparser import ConfigParser
from orbit3d import orbit_plots         #import orbit_plots plotting package


def initialize_plot_options(config):
    """
        initialize the user defined plotting options from config.ini
    """

    target = config.get('plotting', 'target')
    nplanets = config.getint('mcmc_settings', 'nplanets')
    RVFile = config.get('data_paths', 'RVFile')
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    AstrometryFile = config.get('data_paths', 'AstrometryFile', fallback=None)
    GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback=None)
    Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback=None)
    Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback=None)
    HGCAFile = config.get('data_paths', 'HGCAFile', fallback=None)
    
    # colorbar settings
    plot_colorbar = config.getboolean('plotting', 'colorbar', fallback=False)
    color_map = config.get('plotting', 'colormap', fallback= 'viridis')
     
    # custom range of epochs
    start_ep = config.getfloat('plotting', 'start_epoch', fallback=0)
    end_ep = config.getfloat('plotting', 'end_epoch', fallback=0)
     
    # how many orbits
    num_orbits = config.getint('plotting', 'num_orbits', fallback = 50)
    cm_ref = config.get('plotting', 'reference')
    
    args = parse_args_plotting()

    #read in the mcmc chains
    burnin = config.getint('plotting', 'burnin', fallback=0)
    MCMCFile = config.get('plotting', 'McmcDataDir', fallback=None)
    source = MCMCFile.split('_')[0]
    tt, lnp, extras = [fits.open(MCMCFile)[i].data for i in range(3)]
    nsteps = 50*tt.shape[1]
    beststep = np.where(lnp==lnp.max())
   
    # initialize the OP object
    OP = orbit_plots.OrbitPlots(target, HipID, start_ep, end_ep, cm_ref, num_orbits, color_map, burnin, MCMCFile, RVFile, AstrometryFile, HGCAFile, args.output_dir)
    return OP


def run():
    """
    Initialize and make plots
    """
    
    start_time = time.time()
   
    args = parse_args_plotting()
    config = ConfigParser()
    config.read(args.config_file)
    
    # initialize the OP class object
    OPs = initialize_plot_options(config)
    
    # which plot
    plot_settings = {}
    burnin = config.getint('plotting', 'burnin', fallback=0)
    plot_astr = config.getboolean('plotting', 'Astrometry_orbits_plot', fallback=False)
    plot_rv = config.getboolean('plotting', 'RV_orbits_plot', fallback=False)
    plot_rel_rv = config.getboolean('plotting', 'Relative_RV_plot', fallback=False)
    plot_rel_sep = config.getboolean('plotting', 'Relative_separation_plot', fallback=False)
    plot_position_angle = config.getboolean('plotting', 'Position_angle_plot', fallback=False)
    plot_proper_motions = config.getboolean('plotting', 'Proper_motion_plot', fallback=False)
    plot_corner = config.getboolean('plotting', 'Corner_plot', fallback=False)
    test_code = config.getboolean('plotting', 'test_code', fallback=False)

        
    if plot_astr is True:
        OPs.astrometry()
    
    if plot_rv is True:
        OPs.RV()
    
    if plot_rel_rv is True:
        OPs.relRV()
    
    if plot_rel_sep is True:
        OPs.relsep()
        
    if plot_position_angle is True:
        OPs.PA()
    
    if plot_proper_motions is True:
        OPs.proper_motions()
    
    if plot_corner is True:
        OPs.plot_corner()
        
    if test_code is True:
        OPs.test()
    
if __name__ == "__main__":
    run()
