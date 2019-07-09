from orbit3d.config import parse_config


def test_parse_config():
    params = parse_config('orbit3d/tests/config.txt')
    assert len(params) == 6
    assert params['HipID'] == '95319'
