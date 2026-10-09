"""
Tests for the model acceleration and jerk used in the proper-motion figures.

Two things are checked:

1. ``calc_offsets`` is what populates the Gaia sky-path offsets.  Calling
   ``calc_EA_RPP`` alone leaves them at zero, which silently made the plotted
   model acceleration identically zero and the O-C equal to the full measured
   acceleration.

2. A closure test of the forward model: fit a synthetic single-star Gaia NSS
   source, then re-derive the acceleration (and jerk) from the posterior sky
   path through the plotting code and confirm it reproduces the NSS values the
   fit was constrained by.
"""
import os
import tempfile

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from orvara import orbit                                       # noqa: E402
from orvara.config import load_config, nss_settings            # noqa: E402
from orvara.gaia_accel import gen_fit                          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
GAIA_ACCEL_DIR = os.path.join(HERE, 'data', 'gaia_accel')
DT_G = 2.76

# mpri, msec, sau, esino, ecoso, inc, asc, lam, jitter
THETA = np.array([1.0, 0.1, 10.0, 0.4, 0.3, 1.0, 1.0, 1.0, 0.5])


def _nss_file(config_name, tmp_dir):
    path = os.path.join(GAIA_ACCEL_DIR, config_name)
    settings = nss_settings(load_config(path), path=path)
    settings['output_fits'] = os.path.join(tmp_dir, 'target.fits')
    gen_fit.write_nss_fits(settings)
    return settings


def test_calc_offsets_is_what_populates_the_gaia_offsets():
    with tempfile.TemporaryDirectory() as tmp:
        settings = _nss_file('gaia_accel_prep_Acc9.ini', tmp)
        data = orbit.Data(settings['hip_id'], settings['output_fits'], '', '',
                          gaia_mission_length_yrs=DT_G, verbose=False)
        h = DT_G / 2.0
        epochs = 2455197.5 + (np.array([-h, -h / 3, h / 3, h])
                              + settings['ra'] * 0.0 + 6.0) * 365.25
        data.custom_epochs(epochs)
        model = orbit.Model(data)
        par = orbit.Params(THETA, 0, 1, data.nInst, 1)

        orbit.calc_EA_RPP(data, par, model)
        dRA_G, dDec_G = model.return_dRA_dDec()[:2]
        assert np.allclose(dRA_G, 0), 'calc_EA_RPP should not fill dRA_G'
        assert np.allclose(dDec_G, 0)

        orbit.calc_offsets(data, par, model, 0)
        dRA_G, dDec_G = model.return_dRA_dDec()[:2]
        assert not np.allclose(dRA_G, 0), 'calc_offsets must fill dRA_G'
        assert not np.allclose(dDec_G, 0)


FIT_CONFIG = """
[data_paths]
HipID = {hip_id}
HGCAFile = {fits}
RVFile =
AstrometryFile =
start_file = None

[mcmc_settings]
ntemps = 4
nwalkers = 20
nplanets = 1
nstep = 400
nthreads = 1
thin = 10
use_epoch_astrometry = False

[priors_settings]
mpri = 1.0
mpri_sig = 0.1
minjitter = 1e-3
maxjitter = 100

[plotting]
McmcDataFile = {chain}
target = TEST
burnin = 5
iplanet = 0
start_epoch = 1990.
end_epoch = 2030.
num_orbits = 5
num_steps = 100
check_convergence = False
Proper_motion_plot = True
Astrometry_orbits_plot = False
Astrometric_prediction_plot = False
RV_orbits_plot = False
RV_plot = False
Relative_separation_plot = False
Position_angle_plot = False
Corner_plot = False
use_colorbar = True
colormap = viridis
reference = msec_solar
marker_color = blue
set_limit = False
show_title = False
add_text = False
"""


@pytest.mark.integration
def test_model_acceleration_reproduces_the_nss_solution_it_was_fit_to():
    """
    Closure test of the acceleration/jerk forward model.

    The fit is constrained by the NSS acceleration (and, here, jerk).  Running
    the sky path of the resulting posterior back through the plotting code must
    return those same coefficients -- if it does not, either the likelihood or
    the plotted model is wrong.  Before ``calc_offsets`` was called, this
    returned exactly zero.
    """
    import mock
    from orvara.main import run
    from orvara.main_plotting import initialize_plot_options
    from orvara.tests.utils import FakeArguments

    with tempfile.TemporaryDirectory() as tmp:
        settings = _nss_file('gaia_accel_prep_Acc9.ini', tmp)
        config_path = os.path.join(tmp, 'fit.ini')
        with open(config_path, 'w') as handle:
            handle.write(FIT_CONFIG.format(hip_id=settings['hip_id'],
                                           fits=settings['output_fits'],
                                           chain=''))
        with mock.patch('orvara.main.parse_args') as fake_args:
            fake_args.return_value = FakeArguments(config_path, tmp)
            run()
        chain = [name for name in os.listdir(tmp) if name.endswith('.fits')
                 and 'chain' in name]
        assert chain, 'the fit produced no chain'

        with open(config_path, 'w') as handle:
            handle.write(FIT_CONFIG.format(hip_id=settings['hip_id'],
                                           fits=settings['output_fits'],
                                           chain=os.path.join(tmp, chain[0])))
        # initialize_plot_options reads the parsed command line for output_dir
        with mock.patch('orvara.main_plotting.parse_args_plotting') as fake_plot_args:
            fake_plot_args.return_value = FakeArguments(config_path, tmp)
            plots = initialize_plot_options(
                load_config(config_path, require_fit=True))
        best = int(np.where(plots.chain['lnp'] == np.amax(plots.chain['lnp']))[0][0])
        accra, accdec, jerkra, jerkdec = plots.model_accel_jerk(best)

        assert plots.gaia_npar == 9
        assert accra != 0 and accdec != 0, 'model acceleration is identically zero'
        # the fit is driven by these numbers, so it must reproduce them
        assert abs(accra - plots.nss_accra) < 3 * plots.nss_accra_err
        assert abs(accdec - plots.nss_accdec) < 3 * plots.nss_accdec_err
        assert abs(jerkra - plots.nss_jerkra) < 3 * plots.nss_jerkra_err
        assert abs(jerkdec - plots.nss_jerkdec) < 3 * plots.nss_jerkdec_err
