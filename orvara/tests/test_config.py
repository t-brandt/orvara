from configparser import ConfigParser


def test_example_config_file():
    parser = ConfigParser()
    parser.read('orvara/tests/config.ini')
    assert len(parser.items('data_paths')) == 8
    assert len(parser.items('mcmc_settings')) == 6
    assert parser.getint('mcmc_settings', 'nthreads') == 1
    assert parser.getboolean('mcmc_settings', 'use_epoch_astrometry')
