orbit3d
===============

This repo contains orbit3d, the package for fitting orbits of exoplanets.


Installation
------------
orbit3d is built by running :code:`python setup.py build_ext --inplace` or :code:`pip install -e .`
while in the the root directory of this repo. The :code:`-e` flag builds the Cython modules. HTOF is a requirement
for this package. Install it using :code:`pip install git+https://www.github.com/gmbrandt/HTOF` or by following
the installation directions for that repo.

Configuration
-------------
First, assign the appropriate file directories and settings inside of a config.ini file. See the example config.ini file in
:code:`orbit3d/tests/data/config.ini`. If you are using relative astrometry, you must
give paths for :code:`GaiaDataDir`, :code:`Hip1DataDir`, and :code:`Hip2DataDir`. Those are the paths
to the intermediate data for GaiaDR2, the original Hipparcos data reduction, and the second Hipparcos data reduction.

Usage
-----
After setting paths and MCMC (markov-chain monte-carlo)  settings in a config.ini file,
you fit an orbit by running the following from the command line (while in the root directory of the repo).

.. code-block:: bash

    fit_orbit --output-dir /path/to/output --config-file path/to/config.ini

One can set the number of threads in the config.ini file via :code:`nthreads`. Note that the built-in parallelization
is poor. It is better to set nthreads to 1 then simply run multiple instances of orbit3d
on separate cores. One can set the initial conditions of the orbit via the config.ini file.
You can access the help menu with the --help flag as follows.

.. code-block:: bash

    fit_orbit --help

Examples
--------
To run a quick test using the test data and test config.txt in orbit3d/tests, I would cd
to the root directory of orbit3d, then run the following

.. code-block:: bash

    fit_orbit --output-dir ~/Downloads --config-file orbit3d/tests/config.ini

This will create a .fits file in the downloads folder. The MCMC should terminate in less than
one second because of the short number of steps indicated in the example config file.

License
-------

...