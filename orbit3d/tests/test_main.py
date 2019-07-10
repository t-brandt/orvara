import numpy as np

from orbit3d.main import set_initial_parameters


def test_set_initial_parameters():
    ntemps, nplanets, nwalkers = np.random.randint(1, 30), np.random.randint(1, 30), np.random.randint(1, 30)
    params = set_initial_parameters('none', ntemps, nplanets, nwalkers)
    assert np.isclose(params[:, 0, 0].size, ntemps)
    assert np.isclose(params[0, :, 0].size, nwalkers)
