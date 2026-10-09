"""[data_paths] relast_error_scale: a single factor on the relative-astrometry
uncertainties, applied when the file is read and recorded in the fit log."""
import os
import numpy as np
import pytest

from orvara import orbit
from orvara.config import load_config
from orvara.main import _log_error_conventions

HERE = os.path.dirname(os.path.abspath(__file__))
HGCA = os.path.join(HERE, '..', '..', 'HGCA_vEDR3.fits')
RELAST = os.path.join(HERE, 'data', 'Gl758_relAST.txt')


def _data(scale):
    return orbit.Data(95319, HGCA, '', RELAST, relast_error_scale=scale, verbose=False)


def test_scale_multiplies_sep_and_pa_errors():
    d1, d2 = _data(1.0), _data(2.5)
    assert d1.nAst == d2.nAst > 0
    assert np.allclose(np.asarray(d2.relsep_err), 2.5 * np.asarray(d1.relsep_err))
    assert np.allclose(np.asarray(d2.PA_err), 2.5 * np.asarray(d1.PA_err))
    assert d2.relast_error_scale == 2.5


def test_default_is_unity_and_nonpositive_is_rejected():
    d = orbit.Data(95319, HGCA, '', RELAST, verbose=False)
    assert d.relast_error_scale == 1.0
    with pytest.raises(ValueError):
        _data(0.0)


def test_key_is_accepted_by_the_validator(tmp_path):
    cfg = open(os.path.join(HERE, 'config.ini')).read()
    cfg = cfg.replace('[data_paths]', '[data_paths]\nrelast_error_scale = 1.3', 1)
    f = tmp_path / 'c.ini'
    f.write_text(cfg)
    with _nowarn():
        c = load_config(str(f), require_fit=True, strict=True)
    assert c.getfloat('data_paths', 'relast_error_scale') == 1.3


class _nowarn:
    """Fail the test if load_config emits any warning (an unknown key would)."""
    def __enter__(self):
        import warnings
        self._w = warnings.catch_warnings(record=True)
        self._rec = self._w.__enter__()
        warnings.simplefilter('always')
        return self
    def __exit__(self, *a):
        self._w.__exit__(*a)
        assert not self._rec, [str(w.message) for w in self._rec]


def test_fit_log_records_both_conventions(capsys):
    _log_error_conventions(HGCA, RELAST, 1.37)
    out = capsys.readouterr().out
    assert 'scaled by 1.3700' in out          # the relative-astrometry scale
    assert 'Gl758_relAST.txt' in out
    # HGCA_vEDR3.fits carries no ERRINFL keyword, so no inflation line for it
    assert 'inflated' not in out
