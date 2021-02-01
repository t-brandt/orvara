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

HDU0: Parameters:
~~~~~~~~~~~~~~~~~
This extension is a 3d-array of shape (nwalkers,  nsteps/50, nparameters) with nparameters=2+7*nplanets. E.g.
HDU0[10, 40, :] will be the parameters of walker 10 at step 2000 (50 * 40).
Parameters are in order of 0, 1, 2,...:

0. (e.g. fits.open(chain)[0].data[0]) RV jitter. Note that by jitter here, we do not mean irreducible RV scatter on top of the RV error bars that is due to things like stellar convection. By jitter, we mean extra scatter that could be due to *any* source: error underestimation, convective scatter, etc...
1. (e.g. fits.open(chain)[0].data[1]) Primary mass (Msun)
2. Secondary mass (Msun)
3. Semi major axis (A.U.)
4. sqrt(e) * sin(omega), where e is the eccentricity and omega is the argument of periastron in radians
5. sqrt(e) * cos(omega), where e is the eccentricity and omega is the argument of periastron in radians
6. inclination in radians
7. Position angle of the ascending node
8. Mean longitude at the reference epoch. The reference epoch is currently (always) BJD 2455197.50

Then parameters 2-8 repeat for any additional companions, e.g.

9. companion 2 (Tertiary) mass (Msun)
10. Semi major axis of companion 2 (A.U.)
11. sqrt(e) * sin(omega) of companion 2, where e is the eccentricity and omega is the argument of periastron in radians
12. sqrt(e) * cos(omega) of companion 2, where e is the eccentricity and omega is the argument of periastron in radians
13. inclination in radians of companion 2
14. Position angle of the ascending node of companion 2
15. Mean longitude at the reference epoch of companion 2

and so forth for any additional companions.

HDU1: Log likelyhood:
~~~~~~~~~~~~~~~~~~~~~
2d-array of shape (nwalkers,  nsteps/50) which is the log likelyhood for each set
of parameters. E.g. HDU1[10, 40] will be the log likelyhood for the paremeters given
by HDU0[10, 40, :]. Note that this likelyhood includes matrix determinants; it isn't just the chisquared.


For example, one can use this extension to compare the likelyhoods of the best orbits if a certain posterior is multimodal.
Assume that the marginalized posterior in PA is multimodal, with a mode at a value > 180 degrees, and
a mode at a value less than 180 degrees. The following code would print the likelyhoods of the maximum likelyhood orbits
at each of those two modes.


.. code-block:: python

    tt = fits.open('my/path/chain.fits')[0].data
    logl = fits.open('my/path/chain.fits')[1].data
    pa_data = (tt[:,:,7]*180/np.pi) % 360
    ls180 = pa_data < 180
    g180 = pa_data >= 180

    print(np.max(logl[ls180]))  # the max log likelyhood of all orbits with PA of ascending nodes < 180 degrees
    print(np.max(logl[g180]))  # the max log likelyhood of all orbits with PA of ascending nodes > 180 degrees

HDU2:
~~~~~

This extension contains the following the chains of the fit (and derived) parameters.
It is a 3d array of shape (nwalkers, nsteps, 8 + nRV_inst) where nRV_inst is the number of
rv instruments in the fit. That nRV_inst scaling is because the last rows in this 3d array are the radial
velocity offsets for each instrument.

The arrays in these extensions should be treated just like the chains in HDU0.

1. Parallax
2. center-of-mass RA* (right ascension times cos delta) proper motion
3. center-of-mass Dec (declination or delta) proper motion
4. formal chi squared of the fit to the relative separations
5. formal chi squared of the fit to the position angles
6. formal chi squared of the fit to the Hipparcos proper motions
7. formal chi squared of the fit to the Hipparcos-Gaia mean proper motions (from the HGCA)
8. formal chi squared of the fit to the Gaia proper motions (from the HGCA)
9. RV offset for instrument labelled 0 in the input data files
10. RV offset for instrument labelled 1 in the input data files
11. RV offset etc..

Note that if you have no RV instruments, HDU2 will only have length 8 along the last column.

If you want an overall absolute astrometric chi squared, you would add the values from items (6), (7), and (8) above.
There are effectively four measurements since the mean proper motion of the system was fit (values (2) and (3)).

Examples
--------
To run a quick test using the test data and test config.txt in orbit3d/tests, I would cd
to the root directory of orbit3d, then run the following

.. code-block:: bash

    fit_orbit --output-dir ~/Downloads --config-file orbit3d/tests/config.ini

This will create a .fits file in the downloads folder. The MCMC should terminate in less than
one second because of the short number of steps indicated in the example config file.

The end-to-end tests in test_e2e check that the code is converging to previously accepted
values for HIP3850. If you wanted to run the code yourself on this test case and
check the results yourself against those in misc/Diagnostic_plots.ipynb, you can run:

.. code-block:: bash

    fit_orbit --output-dir ~/Downloads --config-file orbit3d/tests/diagnostic_config.ini

The diagnostic_config.ini has the same parameters as those used to create the plots in
Diagnostic_plots.ipynb

Plotting Examples
-----------------

You can access the help menu with the --help flag as follows.

.. code-block:: bash

    plot_orbit --help

To plot orbits, run the plot_orbit command from the root directory, for example

.. code-block:: bash

    plot_orbit --output-dir ./plots --config-file orbit3d/tests/config_HD4747.ini

License
-------

...
