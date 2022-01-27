#!/usr/bin/env python

### Finalized version
"""
Plotting code. The run function is the console entry point,
accessed by calling plot_orbit from the command line.

Example:

plot_orbit --output-dir ./Plots --config-file ./orvara/tests/config_Gl758.ini

"""

from __future__ import print_function
import numpy as np
import pandas as pd
import os
import time
#import emcee, corner
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
from orvara import orbit
from orvara.config import parse_args_plotting
import argparse
from configparser import ConfigParser
from orvara import orbit_plots         # import orbit_plots plotting package


def initialize_plot_options(config):
    """
        initialize the user defined plotting options from config.ini
    """

    OP = orbit_plots.OrbitPlots()

    # target information
    OP.target = OP.title = config.get('plotting', 'target', fallback='')
    OP.Hip = config.getint('data_paths', 'HipID', fallback=0)
    OP.nplanets = config.getint('mcmc_settings', 'nplanets')

    # read data
    OP.RVfile = config.get('data_paths', 'RVFile', fallback=None)
    OP.relAstfile = config.get('data_paths', 'AstrometryFile', fallback=None)
    OP.GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback=None)
    OP.Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback=None)
    OP.Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback=None)
    OP.HGCAFile = config.get('data_paths', 'HGCAFile', fallback=None)
    
    #read in the mcmc chains
    OP.burnin = config.getint('plotting', 'burnin', fallback=0)
    OP.MCMCfile = config.get('plotting', 'McmcDataFile', fallback=None)
    
    # colorbar settings
    OP.usecolorbar = config.getboolean('plotting', 'use_colorbar', fallback=False)
    OP.color_map = config.get('plotting', 'colormap', fallback= 'viridis')
    OP.cmref = config.get('plotting', 'reference', fallback='msec_jup')

    # which planet to plot?
    OP.iplanet = config.getint('plotting', 'iplanet', fallback=0)
    
    # customized range of epochs
    OP.start_epoch = config.getfloat('plotting', 'start_epoch', fallback=1950)
    OP.end_epoch = config.getfloat('plotting', 'end_epoch', fallback=2030)
    
    # predicted epoch positions
    OP.predicted_ep = config.get('plotting', 'predicted_years', fallback=('2010,2020')).split(",")
    OP.predicted_ep_ast = config.getfloat('plotting', 'position_predict', fallback=2000)
    # how many random orbits
    OP.num_orbits = config.getint('plotting', 'num_orbits', fallback = 50)
    
    # step size
    OP.num_steps = config.getint('plotting', 'num_steps', fallback = 1000)
    
    # plot axes settings
    OP.set_limit = config.getboolean('plotting', 'set_limit', fallback=False)
    OP.user_xlim = config.get('plotting', 'xlim', fallback=None).split(",")
    OP.user_ylim = config.get('plotting', 'ylim', fallback=None).split(",")
    
    # show or not show title, add a text on the plot
    OP.show_title = config.getboolean('plotting', 'show_title', fallback=False)
    OP.add_text = config.getboolean('plotting', 'add_text', fallback=False)
    OP.text_name = config.get('plotting', 'text_name', fallback=None)
    OP.x_text = config.getfloat('plotting', 'x_text', fallback=None)
    OP.y_text = config.getfloat('plotting', 'y_text', fallback=None)
    
    # marker settings
    OP.marker_color = config.get('plotting', 'marker_color', fallback= 'coral')
    
    # plot which instrument for the RV plot, starting from 1,2 ... n
    OP.whichInst = config.get('plotting', 'RV_Instrument', fallback=False)
    
    # plot the two proper motion plots separately or together
    OP.pm_separate = config.getboolean('plotting', 'Proper_motion_separate_plots', fallback=False)
    
    #save data
    OP.save_params = config.getboolean('save_results', 'save_params', fallback=True)
    OP.err_margin = config.get('save_results', 'err_margin', fallback= ('0.16, 0.5, 0.84')).split(",")
    
    args = parse_args_plotting()
    OP.outputdir = args.output_dir

    # initialize the OP object
    OP.start()
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
    plot_astr_pred = config.getboolean('plotting', 'Astrometric_prediction_plot', fallback=False)
    plot_astr_topdown = config.getboolean('plotting', 'Astrometry_topdown_plot', fallback=False)
    plot_rv_full = config.getboolean('plotting', 'RV_orbits_plot', fallback=False)
    plot_rv = config.getboolean('plotting', 'RV_plot', fallback=False)
    plot_rel_sep = config.getboolean('plotting', 'Relative_separation_plot', fallback=False)
    plot_position_angle = config.getboolean('plotting', 'Position_angle_plot', fallback=False)
    plot_proper_motions = config.getboolean('plotting', 'Proper_motion_plot', fallback=False)
    plot_corner = config.getboolean('plotting', 'Corner_plot', fallback=False)
    save_params = config.getboolean('save_results', 'save_params', fallback=True)
    checkconv = config.getboolean('plotting', 'check_convergence', fallback=False)
    
    if checkconv:
        OPs.plot_chains()
    if plot_astr:
        OPs.astrometry()
    if plot_astr_pred:
        OPs.astrometric_prediction()
    if plot_astr_topdown:
        OPs.astrometry_topdown()
    if plot_rv_full:
        OPs.RV_fullorbit()
    if plot_rv:
        OPs.RV()
    if plot_rel_sep:
        OPs.relsep()
    if plot_position_angle:
        OPs.PA()
    if plot_proper_motions:
        OPs.proper_motions()
    if plot_corner:
        OPs.plot_corner()
    if save_params:
        OPs.save_data()

    
if __name__ == "__main__":
    run()
