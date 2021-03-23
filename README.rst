orbit3d
===============

This repo contains orbit3d, the package for fitting orbits of exoplanets.


Installation
------------
To install, we will git clone the repo. Run
:code:`git clone https://github.com/t-brandt/orbit3d`
Then:
:code:`cd orbit3d`
Then finally:
:code:`pip install -e .`
orbit3d is built by running :code:`pip install -e .` while in the the root directory
of this repo. The :code:`-e` flag builds the Cython modules.

verifying
~~~~~~~~~

Cd to the root directory of the repo (if you are not already there). Run:
:code:`pytest -sv`
This will run a small suite of tests. This should take about 5 minutes. If any of the tests fail, something
is wrong with the install. Try running :code:`pip install -e .` . If the issue persists, please submit an issue ticket!

Configuration
-------------
First, assign the appropriate file directories and settings inside of a config.ini file. See the example config.ini file in
:code:`orbit3d/tests/data/config.ini`. If you are using relative astrometry, you must
give paths for :code:`GaiaDataDir`, :code:`Hip1DataDir`, and :code:`Hip2DataDir`. Those are the paths
to the intermediate data for GaiaDR2, the original Hipparcos data reduction, and the second Hipparcos data reduction.
Note: if your Hip2 intermediate data come from the DVD, you will want to point to the 'resrec' folder. This should be e.g.:
Hip2_DVD_Book/IntermediateData/resrec

Also in the config.ini file, we recommend always using absolute paths to the various files
(e.g. /home/username/Documents/Hip2_DVD_Book/IntermediateData/resrec). Relative paths do work, but first fits to any source
should use absolute paths until you confirm that all the input data work.


The config file
~~~~~~~~~~~~~~~
Example configuration files can be found in orbit3d/tests/ e.g. orbit3d/tests/config.ini.

This one in particular looks like:

[data_paths]
# Hipparcos ID of the star in question. This is used for fetching it's intermediate astrometry.

HipID = 95319

# The file containing the radial velocity time series for the star.

RVFile = orbit3d/tests/data/Gl758_RV.dat

# The Hipparcos Gaia Catalog

HGCAFile = HGCA_vDR2_corrected.fits

# The file containing the relative astrometry for the star.

AstrometryFile = orbit3d/tests/data/Gl758_relAST.txt

# The path to all the Gaia DR2 intermediate data

GaiaDataDir = orbit3d/tests/data/gaia

# The path to all the Hipparcos (original reduction) intermediate data

Hip1DataDir = orbit3d/tests/data/hip1

# The path to all the Hipparcos (second reduction) intermediate data

Hip2DataDir = orbit3d/tests/data/hip2

# the file path to the initial conditions to the orbit. Set to None for default guess.

start_file = None

[mcmc_settings]

# number of temperatures to use in the parallel tempering chain

ntemps = 5

# number of walkers. Each walker will have ntemps number of chains.

nwalkers = 100

# number of planets to fit.

nplanets = 1

# number of steps contained in each chain

nstep = 100

# number of threads to use with emcee. Note this built-in parellelization is poor.

nthreads = 2

# True if you want to use the epoch astrometry in GaiaDataDir, Hip1DataDir etc... False if not.

use_epoch_astrometry = True

Notes: We recommend starting with use_epoch_astrometry = False. If this fails, then there is something
wrong with the RVFile, HGCAFile, or (relative) AstrometryFile. If that chain finishes fine, then set use_epoch_astrometry = True

Setting priors
~~~~~~~~~~~~~~
Adding Gaussian mass priors are supported. As well, uniform priors on the RV jitter are supported.  To add a guassian mass prior on the primary, you will want to add the following
section to your configuration file:

[priors_settings]

mpri = 1

mpri_sig = inf

minjitter = 1e-5

maxjitter = 1e3


minjitter and maxjitter are the lower and upper bounds respectively for the Uniform prior. E.g.
If you wanted to set a Uniform prior on jitter between 1 and 300 meters/second (and zero outside that prior), you would
set minjitter = 1 and maxjitter = 300.

mpri is the mean of the gaussian prior on the primary mass. mpri_sig is the standard deviation
of that distribution. E.g. for a gaussian prior on the primary mass of 1 solar mass and 0.1 solar mass deviation, you would
set mpri = 1 and mpri_sig = 0.1.

Leaving mpri_sig = inf will turn off the prior.

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

0. Parallax
1. center-of-mass RA* (right ascension times cos delta) proper motion
2. center-of-mass Dec (declination or delta) proper motion
3. formal chi squared of the fit to the relative separations
4. formal chi squared of the fit to the position angles
5. formal chi squared of the fit to the Hipparcos proper motions
6. formal chi squared of the fit to the Hipparcos-Gaia mean proper motions (from the HGCA)
7. formal chi squared of the fit to the Gaia proper motions (from the HGCA)
8. RV offset for instrument labelled 0 in the input data files
9. RV offset for instrument labelled 1 in the input data files
10. RV offset etc..

Note that if you have no RV instruments, HDU2 will only have length 8 along the last column.

If you want an overall absolute astrometric chi squared, you would add the values from items (6), (7), and (8) above.
There are effectively four measurements since the mean proper motion of the system was fit (values (2) and (3)).

For instance, displaying hdu2[:, :, 0] will show all the walkers for the parallax chain (however this parameter
is marginalized over in orbit3d, it is not fit). numpy.mean(hdu2[:, burn:, 0]), numpy.std(hdu2[:, burn:, 0])
would give the mean and standard deviation of the parallax (with burn = some integer that is the number of steps/thinning factor
that you are discarding as burn in)

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

Usage
-----
Once a .fits file from the output of the MCMC is generated, you can produce several plots of 
an orbit by running the following in the command line in the root directory of the repo. To do
this, specify the path to the directory containing the .fits MCMC output file. 

.. code-block:: bash

    plot_orbit --output-dir /path/to/output --config-file path/to/config.ini
    
You can access the help menu with the --help flag as follows.

.. code-block:: bash

    plot_orbit --help

Main plots orvara can produce:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Astrometry orbit of the companion
2. Radial Velocity (RV) orbit
3. Relative RV orbit
4. Relative separation of the two companions
5. Position angle between the two companions
6. Astrometric acceleration or proper motion fit to Hipparcos-Gaia Astrometry

To generate any of these plots, simply set the corresponding parameters under the
[plotting section] in the config.ini file to a boolean variable True. If False, 
a plot would not be produced. Here, for 1. Astrometry orbit plots, you can modify the
predicted_years parameter to plot random predicted epoch positions on the Astrometry plot.
For 2. RV orbit of the companion, you can choose to plot a specific instrument (by name) or
all of the RV instruments by changing the Relative_RV_Instrument parameter to either the
name of the instrument or All. For 6. Proper motion plots, you can plot the proper motions
in RA and DEC in one plot (Proper_motion_separate_plots = False) or 
two (Proper_motion_separate_plots = True). In general, you can also set a customized range of
epochs you want to plot, as well as number of orbits sampled from the poserior distributions
and the resolution (step size). 

Other outputs:
~~~~~~~~~~~~~~
In addition to the six plots, you can check convergence of fitted parameters in
the HDU0 extention by setting the parameter check_convergence to True. You can define
the length of the burn-in phase, note that the parameters are sampled every 50 steps. And you can 
save the results from the fitted and infered parameters from the HDU1 extention
with save_params = True in the [save_results] section, with an option of setting 
the sigma percentages for the errors. 

Color bar settings:
~~~~~~~~~~~~~~~~~~~
User has the options of showing an error bar via use_colorbar = False or True, setting a colormap from 
matplotlib list of colormaps, and a reference scheme for the colorbar. Three reference schemes
are avaliable: the eccentricity as ecc, the secondary companion in jupiter mass as msec_jup and
the secondary companion in solar mass as msec_solar.

Multiple Keplerian orbit fits:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the case of a 3-body or multiple-body fit, you can plot the results for each companion 
by setting iplanet to the corresponding companion ID used in the fitting. 
iplanet starts from 0 for the innermost companion.


Examples
--------

To plot orbits, run a quick test with the plot_orbit command from the root directory, for example

.. code-block:: bash

    plot_orbit --output-dir ./plots --config-file orbit3d/tests/config_HD4747.ini

Then, plot your MCMC chains by pointing to the paths for the configuration file following -config-file and 
the output directory for the plots following -output-dir.

    
Contribution Guidelines
-----------------------
We encourage contributions to orbit3d. The workflow for contributing is the following.

First time contributers:
 * Fork the repository
 * Checkout a new branch for your feature or bug fix.
 * Make your changes to that branch.
 * When you are ready to submit a pull request into the main orbit3d branch (currently called master), run :code:`pytest -sv` to make sure that the required tests pass.
 * If the tests pass, submit your pull request.
 * One approving administrator review is required to approve a pull request.

Users who are invited to be collaborators on the repo:
The same as above, except there is no need to fork the repository once you accept your invite!


License
-------

...
