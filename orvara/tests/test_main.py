import numpy as np
import pytest
import tempfile
import mock
import os
import time

from orvara.main import set_initial_parameters, run, get_priors
from orvara.tests.utils import FakeArguments
from astropy.io import fits
from configparser import ConfigParser


def test_set_initial_parameters():
    ntemps, nplanets, nwalkers = np.random.randint(1, 30), np.random.randint(1, 30), np.random.randint(1, 30)
    params = set_initial_parameters('none', ntemps, nplanets, nwalkers)
    assert np.isclose(params[:, 0, 0].size, ntemps)
    assert np.isclose(params[0, :, 0].size, nwalkers)


def test_get_priors():
    config = ConfigParser()
    config.read('orvara/tests/config_with_secondary_priors.ini')
    priors = get_priors(config)
    assert priors['m_secondary0'] == 1
    assert priors['m_secondary0_sig'] == 2
    assert priors['m_secondary7'] == 3
    assert priors['m_secondary7_sig'] == 4


@pytest.mark.integration
@mock.patch('orvara.main.parse_args')
def test_run(fake_args):
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_args.return_value = FakeArguments('orvara/tests/config.ini', tmp_dir)
        run()
        assert True

@pytest.mark.integration
@mock.patch('orvara.main.parse_args')
def test_run_with_secondary_priors(fake_args):
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_args.return_value = FakeArguments('orvara/tests/config_with_secondary_priors.ini', tmp_dir)
        run()
        assert True

@pytest.mark.e2e
@mock.patch('orvara.main.parse_args')
def test_converges_to_accurate_values(fake_args):
    """
    A test of a fit to HIP3850, verifying that the fit yields the same values as an earlier fit to HIP3850.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_args.return_value = FakeArguments('orvara/tests/diagnostic_config.ini', tmp_dir)
        tt = run()[1].data
        i = -1  # walker index.
        burn = 250  # number of burn in steps to discard
        rv_jitter = np.mean(tt['jitter'][i, burn:])
        rv_jitter_err = np.std(tt['jitter'][i, burn:])
        companion_jup_mass = np.mean(tt['msec0'][i, burn:]*1989/1.898)
        companion_mass_stderr = np.std(tt['msec0'][i, burn:]*1989/1.898)
        separation_AU = np.mean(tt['sau0'][i, burn:])
        separation_stderr = np.std(tt['sau0'][i, burn:])
        eccentricity = np.mean(tt['esino0'][i, burn:]**2 + tt['ecoso0'][i, burn:]**2)
        eccentricity_stderr = np.std(tt['esino0'][i, burn:]**2 + tt['ecoso0'][i, burn:]**2)
        inclination_deg = np.mean(tt['inc0'][i, burn:]*180/np.pi)
        inclination_err = np.std(tt['inc0'][i, burn:]*180/np.pi)

        expected_1_sigma_errors = [0.6282, 2.9215, 0.44668, 0.0030392, 2.3431]
        expected_values = [4.9378, 67.04218, 10.189, 0.73568, 49.89184]
        values = [rv_jitter, companion_jup_mass, separation_AU, eccentricity, inclination_deg]
        errors = [rv_jitter_err, companion_mass_stderr, separation_stderr,
                  eccentricity_stderr, inclination_err]
        for value, expected, sigma in zip(values, expected_values, expected_1_sigma_errors):
            assert np.isclose(value, expected, atol=3 * sigma)
        assert np.allclose(errors, expected_1_sigma_errors, rtol=.5)


@pytest.mark.e2e
@mock.patch('orvara.main.parse_args')
def test_constraints_improve_with_fake_7parameter_dr3_data(fake_args):
    """
    A test of a fit to HIP3850, similar to test_converges_to_accurate_values, however:
    We include here fake DR3 acceleration data (6th and 7th parameters), that are highly precise.
    The orbital constraints should improve markedly. This test verifies that they do indeed improve.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        # TODO give this config file a very good start file, then run the chain for like 20k steps.
        #  will be faster
        fake_args.return_value = FakeArguments('orvara/tests/diagnostic_config_fake_fulldr3.ini', tmp_dir)
        tt = run()[1].data
