import numpy as np
import pytest
import tempfile
import mock
import os

from orbit3d.main import set_initial_parameters, run
from orbit3d.tests.utils import FakeArguments
from astropy.io import fits


def test_set_initial_parameters():
    ntemps, nplanets, nwalkers = np.random.randint(1, 30), np.random.randint(1, 30), np.random.randint(1, 30)
    params = set_initial_parameters('none', ntemps, nplanets, nwalkers)
    assert np.isclose(params[:, 0, 0].size, ntemps)
    assert np.isclose(params[0, :, 0].size, nwalkers)


@pytest.mark.integration
@mock.patch('orbit3d.main.parse_args')
def test_run(fake_args):
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_args.return_value = FakeArguments('orbit3d/tests/config.ini', tmp_dir)
        run()
        assert True


@pytest.mark.e2e
@mock.patch('orbit3d.main.parse_args')
def test_converges_to_accurate_values(fake_args):
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_args.return_value = FakeArguments('orbit3d/tests/diagnostic_config.ini', tmp_dir)
        run()
        # load file and check params
        file = 'HIP3850_chain000.fits'
        tt = fits.open(os.path.join(tmp_dir, file))[0].data
        i = -1  # walker index.
        nsteps = 50 * tt.shape[1]
        burn = 250  # number of burn in steps to discard
        rv_jitter = np.mean(tt[i, burn:, 0])
        rv_jitter_err = np.std(tt[i, burn:, 0])
        companion_jup_mass = np.mean(tt[i, burn:, 2]*1989/1.898)
        companion_mass_stderr = np.std(tt[i, burn:, 2]*1989/1.898)
        separation_AU = np.mean(tt[i, burn:, 3])
        separation_stderr = np.std(tt[i, burn:, 3])
        eccentricity = np.mean(tt[i, burn:, 4]**2 + tt[i, burn:, 5]**2)
        eccentricity_stderr = np.std(tt[i, burn:, 4]**2 + tt[i, burn:, 5]**2)
        inclination_deg = np.mean(tt[i, burn:, 6]*180/np.pi)
        inclination_err = np.std(tt[i, burn:, 6]*180/np.pi)

        expected_1_sigma_errors = [0.11231, 2.9215, 0.44668, 0.0030392, 2.3431]
        expected_values = [1.37637, 67.04218, 10.189, 0.73568, 49.89184]
        values = [rv_jitter, companion_jup_mass, separation_AU, eccentricity, inclination_deg]
        errors = [rv_jitter_err, companion_mass_stderr, separation_stderr,
                  eccentricity_stderr, inclination_err]
        for value, expected, sigma in zip(values, expected_values, expected_1_sigma_errors):
            assert np.isclose(value, expected, atol=3 * sigma)
        assert np.allclose(errors, expected_1_sigma_errors, rtol=.5)
