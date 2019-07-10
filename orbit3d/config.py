import argparse
from configparser import ConfigParser


def parse_args():
    parser = argparse.ArgumentParser(description='Fit an orbit. Required arguments are shown with [].')
    parser.add_argument("--output-dir", required=True,
                        help="Directory within which to save the MCMC results.")
    parser.add_argument("--config-file", required=True,
                        help="Path to the planet-specific configuration file.")
    args = parser.parse_args()
    return args
