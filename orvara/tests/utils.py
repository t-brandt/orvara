class FakeArguments(object):
    """Stand-in for the parsed command line, used by mock.patch in the tests."""

    def __init__(self, config_file='', output_dir='', no_validate=False, strict=False):
        self.config_file = config_file
        self.output_dir = output_dir
        # Mirror the real argparse namespace so main.run()/main_plotting.run()
        # can read these without getattr fallbacks.
        self.no_validate = no_validate
        self.strict = strict
