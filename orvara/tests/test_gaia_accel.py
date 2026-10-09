"""
Tests for orvara.gaia_accel: building a one-star HGCA-format file from a Gaia
DR3 NSS acceleration solution, and reading it back through orvara.orbit.Data.
"""
import os
import tempfile

import numpy as np
import pytest
from astropy.io import fits

from orvara import orbit
from orvara.config import ConfigError, load_config, nss_settings
from orvara.gaia_accel import gen_fit

HERE = os.path.dirname(os.path.abspath(__file__))
GAIA_ACCEL_DIR = os.path.join(HERE, 'data', 'gaia_accel')
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
HGCA_EDR3 = os.path.join(REPO_ROOT, 'HGCA_vEDR3.fits')


def _settings(path):
    return nss_settings(load_config(path), path=path)


# ---------------------------------------------------------------------------
# corr_vec
# ---------------------------------------------------------------------------

def test_corr_vec_to_matrix_is_symmetric_with_unit_diagonal():
    npar = 7
    vec = np.arange(npar * (npar - 1) // 2) / 100.
    matrix = gen_fit.corr_vec_to_matrix(vec, npar)
    assert matrix.shape == (npar, npar)
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 1.0)


def test_corr_vec_fill_order_matches_gaia():
    """Gaia stores (0,1), (0,2), (1,2), (0,3), ... -- upper triangle by column."""
    matrix = gen_fit.corr_vec_to_matrix([0.1, 0.2, 0.3], 3)
    assert np.isclose(matrix[0, 1], 0.1)
    assert np.isclose(matrix[0, 2], 0.2)
    assert np.isclose(matrix[1, 2], 0.3)


def test_corr_vec_wrong_length_raises():
    with pytest.raises(ValueError, match='corr_vec'):
        gen_fit.corr_vec_to_matrix(np.zeros(20), 7)


# ---------------------------------------------------------------------------
# build_row
# ---------------------------------------------------------------------------

def test_gaia_only_row_has_the_columns_orvara_reads():
    settings = _settings(os.path.join(GAIA_ACCEL_DIR, 'gaia_accel_prep_J2043.ini'))
    row = gen_fit.build_row(settings)
    for name in ('gaia_npar', 'accra_gaia', 'accdec_gaia', 'accra_gaia_error',
                 'accdec_gaia_error', 'accra_accdec_gaia', 'parallax_gaia',
                 'pmra_gaia', 'epoch_ra_gaia', 'epoch_ra_hip'):
        assert name in row, name
    assert row['gaia_npar'] == 7
    # No Hipparcos data: those rows must not constrain anything.
    assert row['pmra_hip_error'] == gen_fit.BIG_PM_ERROR
    assert row['pmra_hg_error'] == gen_fit.BIG_PM_ERROR


def test_error_inflation_is_applied_to_acceleration_only():
    settings = _settings(os.path.join(GAIA_ACCEL_DIR, 'gaia_accel_prep_J2043.ini'))
    inflate = settings['error_inflation_factor']
    row = gen_fit.build_row(settings)
    assert np.isclose(row['accra_gaia_error'], settings['accra_err'] * inflate)
    assert np.isclose(row['accdec_gaia_error'], settings['accdec_err'] * inflate)
    # the Gaia proper-motion errors are the catalogue values, untouched
    assert np.isclose(row['pmra_gaia_error'], settings['pmra_err'])


def test_npar9_row_carries_the_jerk_terms():
    settings = _settings(os.path.join(GAIA_ACCEL_DIR, 'gaia_accel_prep_Acc9.ini'))
    row = gen_fit.build_row(settings)
    assert row['gaia_npar'] == 9
    for name in ('jerkra_gaia', 'jerkdec_gaia', 'jerkra_gaia_error',
                 'jerkdec_gaia_error', 'jerkra_jerkdec_gaia'):
        assert name in row, name
    assert np.isclose(row['jerkra_gaia'], settings['deriv_accra'])
    assert np.isclose(row['jerkra_gaia_error'],
                      settings['deriv_accra_err'] * settings['error_inflation_factor'])


@pytest.mark.parametrize('config_name', ['gaia_accel_prep_J2043.ini',
                                         'gaia_accel_prep_Acc9.ini'])
def test_gaia_window_is_centred_on_the_nss_reference_epoch(config_name):
    """
    The NSS acceleration model is written in t = epoch - 2016.0 and its proper
    motion is the mean over that window (Gaia DR3 documentation Sect. 7.2.4),
    so orvara's sample window must be centred on 2016.0.  The HGCA-style
    "uncorrelated position epoch", 2016 - corr(ra, pmra) sig_ra / sig_pmra, is
    only right for a linear (5-parameter) source; used here it would bias the
    model proper motion by acc * dt.
    """
    settings = _settings(os.path.join(GAIA_ACCEL_DIR, config_name))
    row = gen_fit.build_row(settings)
    assert row['epoch_ra_gaia'] == gen_fit.GAIA_REF_EPOCH == 2016.0
    assert row['epoch_dec_gaia'] == gen_fit.GAIA_REF_EPOCH
    # and the window the fit will use is the NSS model's 2 x 517.5 d
    assert np.isclose(gen_fit.GAIA_BASELINE_YR, 2 * 517.5 / 365.25)


# ---------------------------------------------------------------------------
# Writing the FITS file, and reading it back through orvara
# ---------------------------------------------------------------------------

def _build(config_name, tmp_dir, **overrides):
    settings = _settings(os.path.join(GAIA_ACCEL_DIR, config_name))
    settings['output_fits'] = os.path.join(tmp_dir, 'target.fits')
    settings.update(overrides)
    gen_fit.write_nss_fits(settings)
    return settings


def test_written_file_reproduces_the_reference_output():
    """The npar=7 example must still produce what it produced before."""
    reference = os.path.join(GAIA_ACCEL_DIR, 'J2043_gaia_accel.fits')
    with tempfile.TemporaryDirectory() as tmp:
        settings = _build('gaia_accel_prep_J2043.ini', tmp)
        with fits.open(settings['output_fits']) as new, fits.open(reference) as old:
            new_row, old_row = new[1].data[0], old[1].data[0]
            for name in ('hip_id', 'gaia_npar', 'accra_gaia', 'accdec_gaia',
                         'accra_gaia_error', 'accdec_gaia_error',
                         'accra_accdec_gaia', 'pmra_gaia', 'pmdec_gaia',
                         'parallax_gaia', 'pmra_hip_error', 'pmra_hg_error'):
                assert np.isclose(float(new_row[name]), float(old_row[name]),
                                  rtol=1e-6), name
            # Deliberate difference from the reference file: the window is now
            # centred on the NSS reference epoch rather than the 5-parameter
            # "uncorrelated position" epoch the old generator computed.
            assert float(new_row['epoch_ra_gaia']) == 2016.0
            assert not np.isclose(float(old_row['epoch_ra_gaia']), 2016.0)


@pytest.mark.parametrize('config_name,npar', [
    ('gaia_accel_prep_J2043.ini', 7),
    ('gaia_accel_prep_Acc9.ini', 9),
])
def test_orvara_data_reads_back_the_nss_solution(config_name, npar):
    with tempfile.TemporaryDirectory() as tmp:
        settings = _build(config_name, tmp)
        data = orbit.Data(settings['hip_id'], settings['output_fits'], '', '',
                          verbose=False)
        assert data.use_abs_ast == 1
        assert data.gaia_npar == npar
        # orvara stores astrometry in arcsec, the file is in mas.
        assert np.isclose(data.accra_G, 1e-3 * settings['accra'])
        assert np.isclose(data.accdec_G, 1e-3 * settings['accdec'])
        if npar == 9:
            assert np.isclose(data.jerkra_G, 1e-3 * settings['deriv_accra'])
            assert np.isclose(data.jerkdec_G, 1e-3 * settings['deriv_accdec'])


def test_npar9_uses_four_gaia_sample_points_without_epoch_astrometry():
    """
    A cubic (acceleration + jerk) needs four sample points per coordinate, so
    the no-epoch-astrometry path must switch from nGaia = 6 to nGaia = 8.
    """
    with tempfile.TemporaryDirectory() as tmp:
        seven = _build('gaia_accel_prep_J2043.ini', tmp)
        data7 = orbit.Data(seven['hip_id'], seven['output_fits'], '', '', verbose=False)
        assert data7.nGaia == 6
    with tempfile.TemporaryDirectory() as tmp:
        nine = _build('gaia_accel_prep_Acc9.ini', tmp)
        data9 = orbit.Data(nine['hip_id'], nine['output_fits'], '', '', verbose=False)
        assert data9.nGaia == 8


def test_jerk_inverse_covariance_is_zero_for_npar7():
    """A 7-parameter source must incur no jerk penalty."""
    with tempfile.TemporaryDirectory() as tmp:
        settings = _build('gaia_accel_prep_J2043.ini', tmp)
        data = orbit.Data(settings['hip_id'], settings['output_fits'], '', '',
                          verbose=False)
        assert np.allclose(data.Cinv_G_jerk_terms, 0)
        assert not np.allclose(data.Cinv_G_acc_terms, 0)


def test_hip_and_gaia_mode_keeps_the_real_hipparcos_row():
    if not os.path.exists(HGCA_EDR3):
        pytest.skip('HGCA_vEDR3.fits is not present in the repository root.')
    hip_id = 3850
    with fits.open(HGCA_EDR3, memmap=False) as hdul:
        catalogue = hdul[1].data
        reference = catalogue[catalogue['hip_id'] == hip_id][0]
        expected_pmra_hip = float(reference['pmra_hip'])
        expected_pmra_hip_error = float(reference['pmra_hip_error'])

    with tempfile.TemporaryDirectory() as tmp:
        settings = _build('gaia_accel_prep_Acc9.ini', tmp,
                          mode='hip_and_gaia', hip_id=hip_id,
                          hgca_template=HGCA_EDR3)
        with fits.open(settings['output_fits']) as hdul:
            row = hdul[1].data[0]
        assert np.isclose(float(row['pmra_hip']), expected_pmra_hip)
        assert np.isclose(float(row['pmra_hip_error']), expected_pmra_hip_error)
        # the real Hipparcos errors are small, unlike the gaia_only placeholders
        assert float(row['pmra_hip_error']) < gen_fit.BIG_PM_ERROR
        # ... and the NSS terms are attached to that same row
        assert int(row['gaia_npar']) == 9
        assert np.isclose(float(row['accra_gaia']), settings['accra'])


def test_hip_and_gaia_mode_rejects_a_star_that_is_not_in_the_hgca():
    if not os.path.exists(HGCA_EDR3):
        pytest.skip('HGCA_vEDR3.fits is not present in the repository root.')
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.raises(ConfigError, match='not in'):
            _build('gaia_accel_prep_Acc9.ini', tmp, mode='hip_and_gaia',
                   hip_id=999999, hgca_template=HGCA_EDR3)


# ---------------------------------------------------------------------------
# The plotting layer reads the same file
# ---------------------------------------------------------------------------

def _orbit_plots_for(fits_path, hip_id):
    """
    OrbitPlots sets its attributes from the config rather than in __init__, so
    the HGCA loader can be exercised on its own.
    """
    import matplotlib
    matplotlib.use('Agg')
    from orvara.orbit_plots import OrbitPlots

    plots = OrbitPlots()
    plots.Hip = hip_id
    plots.HGCAFile = fits_path
    return plots


@pytest.mark.parametrize('config_name,npar,has_jerk', [
    ('gaia_accel_prep_J2043.ini', 7, False),
    ('gaia_accel_prep_Acc9.ini', 9, True),
])
def test_plotting_loads_the_nss_terms(config_name, npar, has_jerk):
    with tempfile.TemporaryDirectory() as tmp:
        settings = _build(config_name, tmp)
        plots = _orbit_plots_for(settings['output_fits'], settings['hip_id'])
        result = plots.load_HGCA_data()

        # load_HGCA_data swallows exceptions and falls back to have_pmdat =
        # False, so a regression here would silently drop the proper-motion
        # panels rather than raise.
        assert plots.have_pmdat is True
        assert all(part is not None for part in result)
        # gaia_only: the placeholder Hipparcos point (1e7 mas/yr error) must
        # not reach the plot, or it sets the y-range to +/- 1e7.
        ep_ra, ep_dec, pm_ra, pm_dec, err_ra, err_dec = result
        assert len(pm_ra) == 1 and len(pm_dec) == 1
        assert err_ra[0] < 1.0 and err_dec[0] < 1.0
        assert plots.gaia_npar == npar
        assert np.isclose(plots.ep_gaia_acc, 2016.0)
        assert np.isclose(plots.dt_gaia, 2 * 517.5 / 365.25)
        assert plots.have_nss_acc is True
        assert np.isclose(plots.nss_accra, settings['accra'])
        assert np.isclose(plots.nss_accdec, settings['accdec'])
        assert plots.have_nss_jerk is has_jerk
        if has_jerk:
            assert np.isclose(plots.nss_jerkra, settings['deriv_accra'])
            assert np.isclose(plots.nss_jerkdec, settings['deriv_accdec'])


def test_plotting_reports_a_missing_star_instead_of_crashing():
    with tempfile.TemporaryDirectory() as tmp:
        settings = _build('gaia_accel_prep_J2043.ini', tmp)
        plots = _orbit_plots_for(settings['output_fits'], 999999)
        assert plots.load_HGCA_data() == [None] * 6
        assert plots.have_pmdat is False


def test_output_fits_is_relative_to_the_working_directory(tmp_path, monkeypatch):
    """
    fit_orbit opens [data_paths] HGCAFile relative to the cwd, so
    orvara-gaia-accel must write [gaia_accel] output_fits relative to the cwd
    too -- otherwise "keep the two identical" is wrong whenever the config
    lives in a subdirectory.
    """
    import shutil
    sub = tmp_path / 'sub'
    sub.mkdir()
    shutil.copy(os.path.join(GAIA_ACCEL_DIR, 'gaia_accel_prep_J2043.ini'),
                sub / 'prep.ini')
    monkeypatch.chdir(tmp_path)
    written = gen_fit.build_from_config(str(sub / 'prep.ini'))
    # the legacy example names its output J2043_gaia_accel.fits
    assert os.path.abspath(written['fits']) == str(tmp_path / 'J2043_gaia_accel.fits')
    assert (tmp_path / 'J2043_gaia_accel.fits').exists()
    assert not (sub / 'J2043_gaia_accel.fits').exists()
