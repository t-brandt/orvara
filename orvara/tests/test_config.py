import os
import tempfile

import numpy as np
import pytest

from orvara.config import (ConfigError, load_config, nss_settings, validate_config,
                           write_template, template_cli, TEMPLATES)

HERE = os.path.dirname(os.path.abspath(__file__))
GAIA_ACCEL_DIR = os.path.join(HERE, 'data', 'gaia_accel')


def _write(text):
    handle = tempfile.NamedTemporaryFile('w', suffix='.ini', delete=False)
    handle.write(text)
    handle.close()
    return handle.name


def test_example_config_file():
    config = load_config(os.path.join(HERE, 'config.ini'))
    assert len(config.items('data_paths')) == 8
    assert len(config.items('mcmc_settings')) == 6
    assert config.getint('mcmc_settings', 'nthreads') == 1
    assert config.getboolean('mcmc_settings', 'use_epoch_astrometry')


@pytest.mark.parametrize('name', TEMPLATES)
def test_shipped_templates_parse(name):
    with tempfile.TemporaryDirectory() as tmp:
        path = write_template(name, os.path.join(tmp, f'{name}.ini'))
        # The templates are deliberately blank, so validation must not choke on
        # empty values -- only on wrong ones.
        config = load_config(path)
        assert config.has_section('data_paths')
        if name != 'fit':
            assert config.has_section('gaia_accel')


def test_template_cli_writes_and_checks(capsys):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, 'acc9.ini')
        assert template_cli(['--template', 'acc9', '-o', out]) == 0
        assert os.path.exists(out)
        assert template_cli(['--check', os.path.join(HERE, 'config.ini')]) == 0
    assert 'OK' in capsys.readouterr().out


def test_unknown_template_rejected():
    with pytest.raises(ConfigError):
        write_template('not_a_template')


def test_missing_file():
    with pytest.raises(ConfigError):
        load_config('/no/such/config.ini')


def test_bad_type_is_reported_with_section_and_key():
    path = _write('[data_paths]\nHGCAFile = x.fits\n'
                  '[mcmc_settings]\nnwalkers = many\n')
    with pytest.raises(ConfigError) as info:
        load_config(path)
    message = str(info.value)
    assert 'mcmc_settings' in message and 'nwalkers' in message


def test_missing_required_key_is_reported_when_the_fit_needs_it():
    path = _write('[mcmc_settings]\nnwalkers = 10\n')
    # required-ness is scoped to the command: only fit_orbit needs HGCAFile,
    # so a prep-only file is not an error by itself.
    load_config(path)
    with pytest.raises(ConfigError) as info:
        load_config(path, require_fit=True)
    assert 'HGCAFile' in str(info.value)


def test_all_problems_reported_at_once():
    path = _write('[data_paths]\nHGCAFile = x.fits\n'
                  '[mcmc_settings]\nnwalkers = many\nntemps = lots\n')
    with pytest.raises(ConfigError) as info:
        load_config(path)
    assert '2 problem(s)' in str(info.value)


def test_typo_in_key_warns_and_suggests():
    path = _write('[data_paths]\nHGCAFile = x.fits\n'
                  '[mcmc_settings]\nnwalker = 10\n')
    with pytest.warns(RuntimeWarning, match='nwalkers'):
        load_config(path)


def test_strict_turns_unknown_keys_into_errors():
    path = _write('[data_paths]\nHGCAFile = x.fits\n'
                  '[mcmc_settings]\nnwalker = 10\n')
    with pytest.raises(ConfigError):
        load_config(path, strict=True)


def test_epoch_astrometry_needs_data_directories():
    path = _write('[data_paths]\nHGCAFile = x.fits\n'
                  '[mcmc_settings]\nuse_epoch_astrometry = True\n')
    with pytest.raises(ConfigError) as info:
        load_config(path)
    assert 'GaiaDataDir' in str(info.value)


# ---------------------------------------------------------------------------
# The Gaia NSS block
# ---------------------------------------------------------------------------

def test_legacy_prep_ini_is_still_understood():
    """The old gaia_accel_prep_AccN.ini layout must keep working unchanged."""
    path = os.path.join(GAIA_ACCEL_DIR, 'gaia_accel_prep_J2043.ini')
    settings = nss_settings(load_config(path), path=path)
    assert settings['npar'] == 7
    assert settings['hip_id'] == 200000
    assert settings['mode'] == 'gaia_only'
    assert settings['corr_vec'].size == 21
    assert np.isclose(settings['accra'], -2.45177558385822)
    assert np.isclose(settings['error_inflation_factor'], 2)


def test_legacy_acc9_prep_ini():
    path = os.path.join(GAIA_ACCEL_DIR, 'gaia_accel_prep_Acc9.ini')
    settings = nss_settings(load_config(path), path=path)
    assert settings['npar'] == 9
    assert settings['corr_vec'].size == 36
    assert np.isclose(settings['deriv_accra'], -0.01384330)
    assert np.isclose(settings['deriv_accdec'], 0.02039483)


def test_corr_vec_length_is_checked_against_npar():
    path = _write('[gaia_accel]\nnpar = 9\nhip_id = 200000\n'
                  'output_fits = out.fits\n'
                  'ra = 1\nra_err = 1\ndec = 1\ndec_err = 1\n'
                  'parallax = 1\nparallax_err = 1\n'
                  'pmra = 1\npmra_err = 1\npmdec = 1\npmdec_err = 1\n'
                  'accra = 1\naccra_err = 1\naccdec = 1\naccdec_err = 1\n'
                  'deriv_accra = 1\nderiv_accra_err = 1\n'
                  'deriv_accdec = 1\nderiv_accdec_err = 1\n'
                  'corr_vec = ' + ' '.join(['0.1'] * 21) + '\n')
    with pytest.raises(ConfigError) as info:
        nss_settings(load_config(path), path=path)
    assert 'corr_vec has 21 entries' in str(info.value)


def test_npar9_requires_the_jerk_terms():
    path = _write('[gaia_accel]\nnpar = 9\nhip_id = 200000\n'
                  'output_fits = out.fits\n'
                  'ra = 1\nra_err = 1\ndec = 1\ndec_err = 1\n'
                  'parallax = 1\nparallax_err = 1\n'
                  'pmra = 1\npmra_err = 1\npmdec = 1\npmdec_err = 1\n'
                  'accra = 1\naccra_err = 1\naccdec = 1\naccdec_err = 1\n'
                  'corr_vec = ' + ' '.join(['0.1'] * 36) + '\n')
    with pytest.raises(ConfigError) as info:
        nss_settings(load_config(path), path=path)
    assert 'deriv_accra' in str(info.value)


def test_unsupported_npar_is_rejected():
    path = _write('[gaia_accel]\nnpar = 5\nhip_id = 1\noutput_fits = out.fits\n'
                  'corr_vec = 0.1\naccra = 1\naccra_err = 1\n'
                  'accdec = 1\naccdec_err = 1\n')
    with pytest.raises(ConfigError) as info:
        nss_settings(load_config(path), path=path)
    assert 'npar' in str(info.value)


def test_mismatched_hgcafile_and_output_fits_warns():
    path = _write('[data_paths]\nHGCAFile = one.fits\n'
                  '[gaia_accel]\noutput_fits = other.fits\n')
    with pytest.warns(RuntimeWarning, match='output_fits'):
        validate_config(load_config(path, validate=False), path=path)


def test_no_nss_block_is_a_clear_message():
    path = os.path.join(HERE, 'config.ini')
    with pytest.raises(ConfigError) as info:
        nss_settings(load_config(path), path=path)
    assert 'no Gaia NSS block' in str(info.value)


def test_keys_are_case_insensitive_as_they_always_were():
    """
    orvara configs have always been read with ConfigParser's default
    lower-casing, so `hipid = 3850` and `HipID = 3850` are the same key.  A
    revision that preserved case silently turned the former into HipID = 0.
    """
    path = _write('[data_paths]\nhipid = 3850\nhgcafile = x.fits\n'
                  '[mcmc_settings]\nNPLANETS = 1\nnstep = 10\n')
    config = load_config(path, require_fit=True)
    assert config.getint('data_paths', 'HipID') == 3850
    assert config.get('data_paths', 'HGCAFile') == 'x.fits'
    assert config.getint('mcmc_settings', 'nplanets') == 1


def test_fit_requires_nplanets_and_nstep():
    """fit_orbit reads these without a fallback, so they must be validated."""
    path = _write('[data_paths]\nHGCAFile = x.fits\n[mcmc_settings]\nnwalkers = 10\n')
    with pytest.raises(ConfigError) as info:
        load_config(path, require_fit=True)
    assert 'nplanets' in str(info.value) and 'nstep' in str(info.value)
