"""
Smoke tests for the plotting entry point.

``plot_orbit`` was untested, which is how a config without ``xlim`` came to
crash it and how a NumPy 2 incompatibility in ``load_HGCA_data`` could have
silently dropped the proper-motion panels.  These tests run the real entry
point on the tutorial chain and check that files come out.
"""
import os
import shutil
import tempfile
from configparser import ConfigParser

import matplotlib
import mock
import pytest

matplotlib.use('Agg')

from orvara.main_plotting import run                      # noqa: E402
from orvara.tests.utils import FakeArguments              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CHAIN = os.path.join('orvara', 'tests', 'chains',
                     'HD159062B_chain_TUTORIAL_BACKUP.fits')


def _config_pointing_at_the_tutorial_chain(tmp_dir, **plotting):
    """Copy config_HD159062.ini, aim it at the chain that is actually shipped."""
    config = ConfigParser()
    config.optionxform = str
    config.read(os.path.join(HERE, 'config_HD159062.ini'))
    config.set('plotting', 'McmcDataFile', CHAIN)
    config.set('plotting', 'num_orbits', '3')
    config.set('plotting', 'num_steps', '200')
    config.set('plotting', 'check_convergence', 'True')   # exercises plot_chains too
    for name in ('Astrometry_orbits_plot', 'Astrometric_prediction_plot',
                 'RV_orbits_plot', 'RV_plot', 'Relative_separation_plot',
                 'Position_angle_plot', 'Proper_motion_plot', 'Corner_plot'):
        config.set('plotting', name, 'False')
    for name, value in plotting.items():
        config.set('plotting', name, value)
    path = os.path.join(tmp_dir, 'plot.ini')
    with open(path, 'w') as handle:
        config.write(handle)
    return path


@pytest.mark.integration
@mock.patch('orvara.main_plotting.parse_args_plotting')
def test_every_plot_type_runs(fake_args):
    """
    Every figure plot_orbit can make, on the tutorial chain.  Each of these has
    broken on a current Matplotlib or NumPy at some point (ContourSet.collections
    removed in 3.10, np.float, length-1 arrays in pull_chain_params, ...).
    """
    if not os.path.exists(CHAIN):
        pytest.skip('run this from the repository root; the tutorial chain is '
                    'addressed by a relative path.')
    with tempfile.TemporaryDirectory() as tmp:
        path = _config_pointing_at_the_tutorial_chain(
            tmp, Astrometry_orbits_plot='True', Astrometric_prediction_plot='True',
            RV_orbits_plot='True', RV_plot='True', Relative_separation_plot='True',
            Position_angle_plot='True', Proper_motion_plot='True', Corner_plot='True')
        fake_args.return_value = FakeArguments(path, tmp)
        run()
        produced = sorted(name for name in os.listdir(tmp) if name.endswith('.pdf'))
        expected = ('astrometric_orbit', 'astrometric_prediction', 'RV_orbit', 'RV_OC',
                    'relsep_OC', 'PA_OC', 'Proper_Motions', 'Corner')
        for stem in expected:
            assert any(name.startswith(stem) for name in produced), \
                f'no {stem}* PDF produced; got {produced}'


@pytest.mark.integration
@mock.patch('orvara.main_plotting.parse_args_plotting')
def test_relative_separation_plot_runs_without_xlim(fake_args):
    """A config with no xlim/ylim used to crash in initialize_plot_options."""
    if not os.path.exists(CHAIN):
        pytest.skip('run this from the repository root.')
    with tempfile.TemporaryDirectory() as tmp:
        path = _config_pointing_at_the_tutorial_chain(
            tmp, Relative_separation_plot='True')
        config = ConfigParser()
        config.optionxform = str
        config.read(path)
        config.remove_option('plotting', 'xlim')
        config.remove_option('plotting', 'ylim')
        with open(path, 'w') as handle:
            config.write(handle)

        fake_args.return_value = FakeArguments(path, tmp)
        run()
        produced = [name for name in os.listdir(tmp) if name.endswith('.pdf')]
        assert produced, 'plot_orbit produced no PDF'
