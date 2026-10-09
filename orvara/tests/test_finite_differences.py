"""
Unit tests for the finite-difference extraction of proper motion, acceleration
and jerk in ``orvara.orbit.calc_PMs_no_epoch_astrometry``.

The model samples the simulated sky path at a few epochs around the Gaia
central epoch and recovers the polynomial coefficients Gaia would report.  Feed
it an exact polynomial and the recovered coefficients must be exact.
"""
import os
import tempfile

import numpy as np
import pytest

from orvara import orbit
from orvara.config import load_config, nss_settings
from orvara.gaia_accel import gen_fit

HERE = os.path.dirname(os.path.abspath(__file__))
GAIA_ACCEL_DIR = os.path.join(HERE, 'data', 'gaia_accel')

DT_G = 2.76          # Gaia DR3 baseline in years


def _data_for(config_name, tmp_dir):
    path = os.path.join(GAIA_ACCEL_DIR, config_name)
    settings = nss_settings(load_config(path), path=path)
    settings['output_fits'] = os.path.join(tmp_dir, 'target.fits')
    gen_fit.write_nss_fits(settings)
    return orbit.Data(settings['hip_id'], settings['output_fits'], '', '',
                      gaia_mission_length_yrs=DT_G, verbose=False)


def _cubic(t, p, v, a, j):
    """f(t) = p + v t + (a/2) t^2 + (j/6) t^3, so f'=v, f''=a, f'''=j at t=0."""
    return p + v * t + 0.5 * a * t ** 2 + j * t ** 3 / 6.0


@pytest.fixture(scope='module')
def data9():
    with tempfile.TemporaryDirectory() as tmp:
        yield _data_for('gaia_accel_prep_Acc9.ini', tmp)


@pytest.fixture(scope='module')
def data7():
    with tempfile.TemporaryDirectory() as tmp:
        yield _data_for('gaia_accel_prep_J2043.ini', tmp)


def test_npar9_recovers_mean_pm_acceleration_and_jerk_exactly(data9):
    assert data9.gaia_npar == 9
    assert data9.nGaia == 8
    h = DT_G / 2.0
    times = np.array([-h, -h / 3.0, h / 3.0, h])

    p_ra, v_ra, a_ra, j_ra = 0.7, -3.1, 0.42, -0.017
    p_dec, v_dec, a_dec, j_dec = -1.3, 2.4, -0.31, 0.026

    model = orbit.Model(data9)
    model.set_dRA_dDec(
        dRA_G=list(_cubic(times, p_ra, v_ra, a_ra, j_ra)) + [0.0] * 4,
        dDec_G=[0.0] * 4 + list(_cubic(times, p_dec, v_dec, a_dec, j_dec)),
        dRA_H1=[0.0] * data9.nHip1, dDec_H1=[0.0] * data9.nHip1)
    orbit.calc_PMs_no_epoch_astrometry(data9, model)

    # The DR3 NSS "pmra" is the mean proper motion over the window, which for
    # a cubic is v + j h^2 / 6 -- not the instantaneous v at the centre.
    mean_v_ra = v_ra + j_ra * h ** 2 / 6.0
    mean_v_dec = v_dec + j_dec * h ** 2 / 6.0
    assert np.isclose(model.pmra_G, mean_v_ra)
    assert not np.isclose(model.pmra_G, v_ra)
    assert np.isclose(model.accra_G, a_ra)
    assert np.isclose(model.jerkra_G, j_ra)
    assert np.isclose(model.pmdec_G, mean_v_dec)
    assert np.isclose(model.accdec_G, a_dec)
    assert np.isclose(model.jerkdec_G, j_dec)


def test_npar9_centroid_is_the_time_average_over_the_baseline(data9):
    """
    pmra_HG is a positional difference over ~25 yr, so the Gaia-side position
    must be the time-average of the path over the mission -- p + (a/2) h^2/3 --
    not the plain mean of the four samples, which would give p + (a/2) 5h^2/9.
    """
    h = DT_G / 2.0
    times = np.array([-h, -h / 3.0, h / 3.0, h])
    p, v, a, j = 0.0, 1.7, 0.9, 0.05

    model = orbit.Model(data9)
    model.set_dRA_dDec(dRA_G=list(_cubic(times, p, v, a, j)) + [0.0] * 4,
                       dDec_G=[0.0] * 8,
                       dRA_H1=[0.0] * data9.nHip1, dDec_H1=[0.0] * data9.nHip1)
    orbit.calc_PMs_no_epoch_astrometry(data9, model)

    baseline = data9.epRA_G - data9.epRA_H
    recovered_RA_G = model.pmra_HG * baseline
    exact_time_average = p + 0.5 * a * h ** 2 / 3.0
    plain_mean = np.mean(_cubic(times, p, v, a, j))

    assert np.isclose(recovered_RA_G, exact_time_average)
    assert not np.isclose(recovered_RA_G, plain_mean)


def test_npar7_recovers_pm_and_acceleration_and_leaves_jerk_zero(data7):
    assert data7.gaia_npar == 7
    assert data7.nGaia == 6
    h = DT_G / 2.0
    times = np.array([-h, 0.0, h])

    v_ra, a_ra = -2.2, 0.35
    v_dec, a_dec = 1.1, -0.44

    model = orbit.Model(data7)
    model.set_dRA_dDec(
        dRA_G=list(_cubic(times, 0.0, v_ra, a_ra, 0.0)) + [0.0] * 3,
        dDec_G=[0.0] * 3 + list(_cubic(times, 0.0, v_dec, a_dec, 0.0)),
        dRA_H1=[0.0] * data7.nHip1, dDec_H1=[0.0] * data7.nHip1)
    orbit.calc_PMs_no_epoch_astrometry(data7, model)

    assert np.isclose(model.pmra_G, v_ra)
    assert np.isclose(model.accra_G, a_ra)
    assert np.isclose(model.pmdec_G, v_dec)
    assert np.isclose(model.accdec_G, a_dec)
    # three points cannot constrain a cubic term
    assert model.jerkra_G == 0
    assert model.jerkdec_G == 0


def test_set_dRA_dDec_rejects_wrong_lengths(data9):
    model = orbit.Model(data9)
    with pytest.raises(ValueError, match='nGaia'):
        model.set_dRA_dDec(dRA_G=[0.0] * 3)
