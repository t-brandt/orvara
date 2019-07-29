import numpy as np
import pytest
import tempfile
import mock

from orbit3d.main import set_initial_parameters, run
from orbit3d.tests.utils import FakeArguments


def test_set_initial_parameters():
    ntemps, nplanets, nwalkers = np.random.randint(1, 30), np.random.randint(1, 30), np.random.randint(1, 30)
    params = set_initial_parameters('none', ntemps, nplanets, nwalkers)
    assert np.isclose(params[:, 0, 0].size, ntemps)
    assert np.isclose(params[0, :, 0].size, nwalkers)


@pytest.mark.e2e
@mock.patch('orbit3d.main.parse_args')
def test_run(fake_args):
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_args.return_value = FakeArguments('orbit3d/tests/config.ini', tmp_dir)
        run()
        assert True