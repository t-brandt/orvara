from configparser import ConfigParser


def test_example_config_file():
    parser = ConfigParser()
    parser.read('orbit3d/tests/config.ini')
    assert len(parser.items('data_paths')) == 7
    assert len(parser.items('mcmc_settings')) == 6
    assert parser.getint('mcmc_settings', 'nthreads') == 2
    assert parser.getboolean('mcmc_settings', 'use_epoch_astrometry')
