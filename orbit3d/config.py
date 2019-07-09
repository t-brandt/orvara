import argparse


def parse_config(config_path, delimiter='='):
    with open(config_path) as f:
        keys_vals = [line.split(delimiter) for line in f.read().splitlines()]
    return dict(keys_vals)


def parse_args():
    parser = argparse.ArgumentParser(description='Fit an orbit. Required arguments are shown without [].')
    parser.add_argument("--nstep", required=False, default=100, metavar='', type=int,
                        help="Number of MCMC steps per walker")
    parser.add_argument("--nthreads", required=False, default=1, metavar='', type=int,
                        help="Number of threads to use for parallelization")
    parser.add_argument("-a", "--use-epoch-astrometry", action="store_true",
                        required=False, default=False,
                        help="Whether or not to use intermediate astrometry data")
    parser.add_argument("--output-dir", required=True,
                        help="Directory within which to save the MCMC results.")
    parser.add_argument("--config-file", required=True,
                        help="Path to the planet-specific configuration file.")
    parser.add_argument("--ntemps", required=False, default=5, metavar='', type=int,
                        help="number of MCMC temperatures.")
    parser.add_argument("--nwalkers", required=False, default=100, metavar='', type=int,
                        help="number of MCMC walkers.")
    parser.add_argument("--nplanets", required=False, default=1, metavar='', type=int,
                        help="Assumed number of planets in the system.")
    parser.add_argument("--start-file", required=False, default=None, metavar='',
                        help="Filepath for the orbit initial conditions.")
    args = parser.parse_args()
    return args