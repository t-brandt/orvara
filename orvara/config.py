import argparse
from configparser import ConfigParser
import os


def parse_args():
    parser = argparse.ArgumentParser(description='Fit a Keplerian orbit or orbits using orvara.')
    parser.add_argument("config_file", metavar='configfile', 
                        help="Path to the planet-specific configuration file.")
    parser.add_argument("--output-dir",
                        help="Directory within which to save the MCMC results, default current directory.")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.getcwd()

    return args

def parse_args_plotting():
    parser = argparse.ArgumentParser(description='Plot astrometry or RV orbits and other relavant plots.')
    parser.add_argument("config_file", metavar='configfile',
                        help="Path to the planet-specific configuration file.")
    parser.add_argument("--output-dir", 
                        help="Directory within which to save the plots, default current directory.")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.getcwd()

    return args
