orbit3d
===============

This repo contains orbit3d, the package for fitting orbits of exoplanets.


Installation
------------
orbit3d is built by running :code:`pip install -e .` while in the the root directory
of this repo. The :code:`-e` flag builds the Cython modules. HTOF is a requirement
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

The output of the MCMC is a .fits file and is contained within your given output directory. The output file
contains three .fits extensions with all the MCMC parameters sampled every 50 steps.
That output file is formatted as follows:

HDU0: Parameters. 3d-array of shape (nwalkers,  nsteps/50, nparameters) with nparameters=2+7*nplanets. E.g.
HDU0[10, 40, :] will be the parameters of walker 10 at step 2000 (50 * 40).
Parameters are in order of 0, 1, 2,...: (RV_jitter, ...)

HDU1: Log likelyhood. 2d-array of shape (nwalkers,  nsteps/50) which is the log likelyhood for each set
of parameters. E.g. HDU1[10, 40] will be the log likelyhood for the paremeters given
by HDU0[10, 40, :]. Note that this likelyhood includes matrix determinants; it isn't just chisq.

HDU2: ?? containing the following 8 best-fit values (in the following order) along the 3rd axis:

::

    1. best-fit parallax
    2. best-fit center-of-mass RA proper motion
    3. best-fit center-of-mass Dec proper motion
    4. chi squared of relative separations
    5. chi squared of position angles
    6. chi squared of Hipparcos proper motions
    7. chi squared of Hipparcos-> Gaia mean proper motions
    8. chi squared of Gaia proper motions

If you want an overall astrometric chi squared, you would add the values from items (6), (7), and (8) above.
There are effectively four measurements since the mean proper motion of the system was fit (values (2) and (3)).

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