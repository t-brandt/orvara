def parse_config(config_path, delimiter='='):
    with open(config_path) as f:
        keys_vals = [line.split(delimiter) for line in f.read().splitlines()]
    return dict(keys_vals)
