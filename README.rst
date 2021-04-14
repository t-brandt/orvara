orvara
===============

This repo contains orvara, the package for fitting orbits of exoplanets.


Installation
------------
To install, we will git clone the repo. Run
:code:`git clone https://github.com/t-brandt/orvara`
Then:
:code:`cd orvara`
Then finally:
:code:`pip install -e .`
orvara is built by running :code:`pip install -e .` while in the the root directory
of this repo. The :code:`-e` flag builds the Cython modules.

verifying
~~~~~~~~~

cd to the root directory of the repo (if you are not already there). Run:
:code:`pytest -sv`
This will run a small suite of tests. This should take about 5 minutes. If any of the tests fail, something
is wrong with the install. Try running :code:`pip install -e .` . If the issue persists, please submit an issue ticket!

Configuration
-------------
First, assign the appropriate file directories and settings inside of a config.ini file. See the example config.ini file in
:code:`orvara/tests/data/config.ini`. If you are using relative astrometry, you must
give paths for :code:`GaiaDataDir`, :code:`Hip1DataDir`, and :code:`Hip2DataDir`. Those are the paths
to the intermediate data for GaiaDR2, the original Hipparcos data reduction, and the second Hipparcos data reduction.
Note: if your Hip2 intermediate data come from the DVD, you will want to point to the 'resrec' folder. This should be e.g.:
Hip2_DVD_Book/IntermediateData/resrec

Also in the config.ini file, we recommend always using absolute paths to the various files
(e.g. /home/username/Documents/Hip2_DVD_Book/IntermediateData/resrec). Relative paths do work, but first fits to any source
should use absolute paths until you confirm that all the input data work.


The config file
~~~~~~~~~~~~~~~
Example configuration files can be found in orvara/tests/ e.g. orvara/tests/config.ini.

This one in particular looks like:

[data_paths]
# Hipparcos ID of the star in question. This is used for fetching it's intermediate astrometry.

HipID = 95319

# The file containing the radial velocity time series for the star.

RVFile = orvara/tests/data/Gl758_RV.dat

# The Hipparcos Gaia Catalog

HGCAFile = HGCA_vDR2_corrected.fits

# The file containing the relative astrometry for the star.

AstrometryFile = orvara/tests/data/Gl758_relAST.txt

# The path to all the Gaia DR2 intermediate data

GaiaDataDir = orvara/tests/data/gaia

# The path to all the Hipparcos (original reduction) intermediate data

Hip1DataDir = orvara/tests/data/hip1

# The path to all the Hipparcos (second reduction) intermediate data

Hip2DataDir = orvara/tests/data/hip2

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
is poor. It is better to set nthreads to 1 then simply run multiple instances of orvara
on separate cores. One can set the initial conditions of the orbit via the config.ini file.
You can access the help menu with the --help flag as follows.

.. code-block:: bash

    fit_orbit --help

The output of the MCMC is a .fits file and is contained within your given output directory. The output file
contains two .fits extensions: an empty one, and a fits table with all the MCMC parameters sample.

HDU0: empty
~~~~~~~~~~~~~~~~~
The first extension is empty for table data.

HDU1: table
~~~~~~~~~~~~~~~~~~~~~
This is a fits table object.  Each table column is of shape (nwalkers, nsteps/thin) where thin is the thinning used in the configuration file (default 50, to save every 50th step).  You may access a column by, e.g.,

lnlike = hdulist[1].data['lnp']

The column names and descriptions are:

'mpri' : Primary mass (Solar masses)
'msec0' : Secondary mass of the first (index 0) companion, Solar masses
'sau0' : Semimajor axis of the first companion, Solar masses
'esino0' : sqrt(ecc)*sin(omega) for the first companion
'ecoso0' : sqrt(ecc)*cos(omega) for the first companion
'inc0' : inclination (radians) for the first companion
'asc0' : PA of the ascending node (radians) for the first companion
'lam0' : Mean longitude at reference epoch for the first companion

If there is more than one companion, then there are additional fields with, e.g., 'msec1', 'msec2', etc.

'jitter' : log RV jitter (jitter in m/s is 10**(0.5*hdulist[1].data['jitter']))
'jitter0' : log RV jitter for instrument 0 
Note that 'jitter0', 'jitter1', etc. are present and 'jitter' is not if using one jitter per instrument.  The default is to use the same jitter for all instruments.  In this case 'jitter' is present but 'jitter0', 'jitter1', etc. are not.

'lnp' : natural log of the (unnormalized) probability.  Note that this includes matrix determinants and is not simply chi squared.

'plx_ML' : maximum likelihood (ML) parallax at this chain step
'pmra_ML' : ML proper motion in RA at this chain step
'pmdec_ML' : ML proper motion in Dec at this chain step
'chisq_sep' : The chi squared in separation at the ML parallax at this chain step
'chisq_PA' : The chi squared in position angle at this chain step
'chisq_H' : The chi squared for the two Hipparcos proper motions
'chisq_HG' : The chi squared for the two long-term Hipparcos-Gaia proper motions
'chisq_G' : The chi squared for the two Gaia proper motions
'RV_ZP_0_ML' : The ML zero point (barycenter RV) for instrument 0

There will be an 'RV_ZP_1_ML' for instrument 1, etc., up to the number of RV instruments.  

If you want an overall absolute astrometric chi squared, you would add the values from items 'chisq_H', 'chisq_HG', and 'chisq_G' above.
There are effectively four measurements since the mean proper motion of the system was fit ('pmra_ML' and 'pmdec_ML').

For instance, displaying hdulist[1].data['plx_ML'] will show all the walkers for the parallax chain (however this parameter
is marginalized over in orvara, it is not fit). numpy.mean(hdulist[1].data['plx_ML'][:, burn:]), numpy.std(hdulist[1].data['plx_ML'][:, burn:])
would give the mean and standard deviation of the parallax (with burn = some integer that is the number of steps/thinning factor
that you are discarding as burn in)

One can use the 'lnp' column to compare the likelihoods of the best orbits if a certain posterior is multimodal.
Assume that the marginalized posterior in PA is multimodal, with a mode at a value > 180 degrees, and
a mode at a value less than 180 degrees. The following code would print the likelihoods of the maximum likelihood orbits
at each of those two modes.


.. code-block:: python

    tt = fits.open('my/path/chain.fits')[1].data
    logl = tt['lnp']
    pa_data = (tt['asc0']*180/np.pi) % 360
    ls180 = pa_data < 180
    g180 = pa_data >= 180

    print(np.max(logl[ls180]))  # the max log likelyhood of all orbits with PA of ascending nodes < 180 degrees
    print(np.max(logl[g180]))  # the max log likelyhood of all orbits with PA of ascending nodes > 180 degrees


Examples
--------
To run a quick test using the test data and test config.txt in orvara/tests, I would cd
to the root directory of orvara, then run the following

.. code-block:: bash

    fit_orbit --output-dir ~/Downloads --config-file orvara/tests/config.ini

This will create a .fits file in the downloads folder. The MCMC should terminate in less than
one second because of the short number of steps indicated in the example config file.

The end-to-end tests in test_e2e check that the code is converging to previously accepted
values for HIP3850. If you wanted to run the code yourself on this test case and
check the results yourself against those in misc/Diagnostic_plots.ipynb, you can run:

.. code-block:: bash

    fit_orbit --output-dir ~/Downloads --config-file orvara/tests/diagnostic_config.ini

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

Main plots orvara is configured to produce from the orbital fit:
~~~~~~~~~~~~~~~~~
1. Astrometry orbit of the companion
2. Radial Velocity (RV) orbit
3. Relative RV orbit
4. Relative separation of the two companions
5. Position angle between the two companions
6. Astrometric acceleration or proper motion fit to Hipparocs-Gaia Astrometry

To generate any of these plots, simply set the correspondig parameters under the 
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
~~~~~~~~~~~~~~~~~
In addition to the six plots, you can check convergence of fitted parameters in
the HDU0 extention by setting the parameter check_convergence to True. You can define
the length of the burn-in phase, note that the parameters are sampled every 50 steps. And you can 
save the results from the fitted and infered parameters from the HDU1 extention
with save_params = True in the [save_results] section, with an option of setting 
the sigma percentages for the errors. 

Color bar settings:
~~~~~~~~~~~~~~~~~
User has the options of showing an error bar via use_colorbar = False or True, setting a colormap from 
matplotlib list of colormaps, and a reference scheme for the colorbar. Three reference schemes
are avaliable: the eccentricity as ecc, the secondary companion in jupiter mass as msec_jup and
the secondary companion in solar mass as msec_solar.

Multiple Keplerian orbit fits:
~~~~~~~~~~~~~~~~~
In the case of a 3-body or multiple-body fit, you can plot the results for each companion 
by setting iplanet to the corresponding companion ID used in the fitting. 
iplanet starts from 0 for the innermost companion.


Examples
--------

To plot orbits, run a quick test with the plot_orbit command from the root directory, for example

.. code-block:: bash

    plot_orbit --output-dir ./plots --config-file orvara/tests/config_HD4747.ini

Then, plot your MCMC chains by pointing to the paths for the configuration file following -config-file and 
the output directory for the plots following -output-dir.

    
Contribution Guidelines
-----------------------
We encourage contributions to orvara. The workflow for contributing is the following.

First time contributers:
 * Fork the repository
 * Checkout a new branch for your feature or bug fix.
 * Make your changes to that branch.
 * When you are ready to submit a pull request into the main orvara branch (currently called master), run :code:`pytest -sv` to make sure that the required tests pass.
 * If the tests pass, submit your pull request.
 * One approving administrator review is required to approve a pull request.

Users who are invited to be collaborators on the repo:
The same as above, except there is no need to fork the repository once you accept your invite!


License
-------

...
