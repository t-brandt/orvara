#!/usr/bin/env python

### Finalized version
"""
Plotting code. The run function is the console entry point,
accessed by calling plot_orbit from the command line.

Example:

plot_orbit --output-dir ./Plots --config-file ./orbit3d/tests/config_Gl758.ini

"""

from __future__ import print_function
import numpy as np
import pandas as pd
import os
import time
import emcee, corner
import scipy.optimize as op
from random import randrange
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
from orbit3d import orbit_plots         # import orbit_plots plotting package


def initialize_plot_options(config):
    """
        initialize the user defined plotting options from config.ini
    """

    # target information
    target = config.get('plotting', 'target')
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    nplanets = config.getint('mcmc_settings', 'nplanets')
    
    # read data
    RVFile = config.get('data_paths', 'RVFile')
    AstrometryFile = config.get('data_paths', 'AstrometryFile', fallback=None)
    GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback=None)
    Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback=None)
    Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback=None)
    HGCAFile = config.get('data_paths', 'HGCAFile', fallback=None)
    
    #read in the mcmc chains
    burnin = config.getint('plotting', 'burnin', fallback=0)
    MCMCFile = config.get('plotting', 'McmcDataDir', fallback=None)
    
    # colorbar settings
    use_colorbar = config.getboolean('plotting', 'use_colorbar', fallback=False)
    color_map = config.get('plotting', 'colormap', fallback= 'viridis')
    cm_ref = config.get('plotting', 'reference')
    colorbar_size = config.getfloat('plotting', 'fraction', fallback=0.046)
    colorbar_pad = config.getfloat('plotting', 'pad', fallback=0.04)
    
    # customized range of epochs
    start_ep = config.getfloat('plotting', 'start_epoch', fallback=0)
    end_ep = config.getfloat('plotting', 'end_epoch', fallback=0)
    
    # predicted epoch positions
    predict_ep = config.get('plotting', 'predicted_years', fallback=0).split(",")
    
    # how many random orbits
    num_orbits = config.getint('plotting', 'num_orbits', fallback = 50)
    
    # step size
    num_steps = config.getint('plotting', 'num_steps', fallback = 1000)
    
    # plot axes settings
    set_limit = config.getboolean('plotting', 'set_limit', fallback=True)
    xlim = config.get('plotting', 'xlim', fallback=None).split(",")
    ylim = config.get('plotting', 'ylim', fallback=None).split(",")
    
    # show or not show title, add a text on the plot
    show_title = config.getboolean('plotting', 'show_title', fallback=True)
    add_text = config.getboolean('plotting', 'add_text', fallback=False)
    text_name = config.get('plotting', 'text_name', fallback=None)
    x_text = config.getfloat('plotting', 'x_text', fallback=None)
    y_text = config.getfloat('plotting', 'y_text', fallback=None)
    
    # marker settings
    marker_color = config.get('plotting', 'marker_color', fallback= 'coral')
    
    # plot the two proper motion plots separately or together
    separate_pm_plots = config.getboolean('plotting', 'Proper_motion_separate_plots', fallback=False)
    
    args = parse_args_plotting()

    # initialize the OP object
    OP = orbit_plots.OrbitPlots(target, HipID, start_ep, end_ep, predict_ep, cm_ref, num_orbits, color_map, use_colorbar, colorbar_size, colorbar_pad, marker_color, burnin, set_limit, xlim, ylim, show_title, add_text, text_name, x_text, y_text, num_steps, MCMCFile, RVFile, AstrometryFile, HGCAFile, args.output_dir, separate_pm_plots)
    return OP


def run():
    """
    Initialize OP and make plots
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
    
if __name__ == "__main__":
    run()
