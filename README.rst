======
orvara
======
This repo contains orvara, the package for fitting orbits of bright stars and
their faint companions (exoplanets, brown dwarfs, white dwarfs, and low-mass
stars).


This build additionally supports Gaia DR3 non-single-star (NSS) astrometric
acceleration solutions -- both 7-parameter (acceleration) and 9-parameter
(acceleration + jerk) -- with or without epoch astrometry. See
`Gaia NSS acceleration sources`_ below.

Requirements
============
Python 3.10-3.13, NumPy >= 1.24 (NumPy 2.x supported), SciPy >= 1.9,
Astropy >= 5.1, pandas >= 1.5, Matplotlib >= 3.6, emcee >= 3.1 and
`htof <https://github.com/gmbrandt/HTOF>`_ >= 1.2.0.

ptemcee is **not** a dependency: a patched copy is bundled in
``orvara/_vendor/ptemcee`` because the PyPI release calls ``np.float``, which
NumPy removed in 1.24. See that package's docstring.

This repository is a fork of https://github.com/t-brandt/orvara at commit 66349265, maintained at https://github.com/Realaqe/orvara.

Installation
============
Cython is declared as a build requirement, so a modern pip or uv installs it
automatically:

.. code-block:: bash

    git clone <this repo>
    cd orvara
    pip install -e ".[test]"

If you use an older pip that does not honour ``pyproject.toml`` build
requirements, install Cython first (``pip install "cython>=3"``).

htof must be installed too; if you are using the matching modernised htof, do
``pip install -e ../htof`` before installing orvara.

Verifying
---------
cd to the root directory of the repo (if you are not already there). Run:

.. code-block:: bash

    pytest -sv

This runs the full suite: unit tests plus end-to-end MCMC fits. Expect roughly
10-20 minutes, most of it in the ``e2e`` fits. To skip those and get a fast
check of the install (about 30 seconds):

.. code-block:: bash

    pytest -sv -m "not e2e"

To run just the Gaia-acceleration and configuration tests:

.. code-block:: bash

    pytest -sv orvara/tests/test_gaia_accel.py orvara/tests/test_config.py \
              orvara/tests/test_finite_differences.py

If any test fails, something is wrong with the install. If you recently switched
branches, delete the generated ``orvara/orbit.c`` and the compiled
``orvara/orbit.*.so`` and run ``pip install -e .`` again. If the issue persists,
please submit an issue ticket on github!

You can also check a configuration file without running anything:

.. code-block:: bash

    orvara-config --check my_target.ini

Gaia NSS acceleration sources
=============================
Gaia DR3 publishes ``gaiadr3.nss_acceleration_astro`` solutions for stars whose
astrometry is not consistent with uniform motion. They come in two flavours:

- **Acceleration7** (``npar = 7``): position, parallax, proper motion,
  acceleration.
- **Acceleration9** (``npar = 9``): the above plus the derivative of the
  acceleration, i.e. jerk.

orvara reads absolute astrometry from an HGCA-format table, so the workflow is
to turn one NSS solution into a one-star HGCA-format FITS file.

Everything lives in a single config file. Start from a template:

.. code-block:: bash

    orvara-config --template acc9 -o my_target.ini    # or acc7, hip_and_gaia

Fill in the ``[gaia_accel]`` section with the values from the archive, set
``[data_paths] HGCAFile`` to the same path as ``[gaia_accel] output_fits``, then:

.. code-block:: bash

    orvara-gaia-accel my_target.ini    # writes the one-star FITS
    fit_orbit my_target.ini            # fits the orbit
    plot_orbit my_target.ini           # makes the plots

Two modes
---------
``mode = gaia_only``
    The star has no Hipparcos data. The Hipparcos and Hipparcos-Gaia rows are
    filled with the Gaia proper motion and given errors of 1e7 mas/yr, so they
    contribute nothing to chi-squared and only the NSS acceleration (and jerk)
    constrains the fit. ``hip_id`` is a synthetic integer of your choosing;
    200000 and above is safe, since real Hipparcos IDs stop at 120404.

``mode = hip_and_gaia``
    The star *is* in the HGCA. Its real Hipparcos and Hipparcos-Gaia proper
    motions are copied from ``hgca_template`` and the NSS acceleration and jerk
    terms are attached to that same row, so the long Hipparcos-Gaia baseline and
    the within-Gaia curvature constrain the orbit together. ``hip_id`` is the
    star's real Hipparcos number.

What the generator writes
-------------------------
The output table carries the 35 standard HGCA columns plus:

===========================  =========  ==================================
Column                       npar       Meaning
===========================  =========  ==================================
``gaia_npar``                7, 9       Which NSS solution this is
``accra_gaia``               7, 9       RA acceleration, mas/yr^2
``accdec_gaia``              7, 9       Dec acceleration, mas/yr^2
``accra_gaia_error``         7, 9       ...and its inflated 1-sigma error
``accdec_gaia_error``        7, 9
``accra_accdec_gaia``        7, 9       Their correlation, from ``corr_vec``
``jerkra_gaia``              9          RA jerk, mas/yr^3
``jerkdec_gaia``             9          Dec jerk, mas/yr^3
``jerkra_gaia_error``        9
``jerkdec_gaia_error``       9
``jerkra_jerkdec_gaia``      9          Their correlation
===========================  =========  ==================================

``epoch_ra_gaia`` and ``epoch_dec_gaia`` are written as **2016.0**, the NSS
reference epoch. The DR3 acceleration model is written in ``t = epoch - 2016.0``
and its position and proper-motion parameters are *means* over the observation
window ``2016.0 +/- 517.5 d`` (Gaia DR3 documentation, Sect. 7.2.4, Eq. 7.5), so
that is the window orvara has to average over. The HGCA's convention for
5-parameter sources -- placing the Gaia epoch where position and proper motion
decorrelate, ``2016 - corr(ra, pmra) sigma_ra / sigma_pmra`` -- must **not** be
used for an NSS source: shifting the window by ``dt`` biases the model proper
motion by ``acc * dt``, tens of sigma for a typical accelerator. A Gaia-only fit
absorbs that into the free barycentric proper motion; a ``hip_and_gaia`` fit,
where Hipparcos pins the barycentric motion, cannot.

``error_inflation_factor`` (default 2 in the templates) multiplies the NSS
acceleration and jerk errors only; the Gaia proper-motion errors are left as
published.

Epoch astrometry, or not
------------------------
With ``use_epoch_astrometry = True`` the Gaia intermediate astrometry is fitted
with htof, at polynomial degree 1, 2 or 3 for ``npar`` 5, 7 or 9.

With ``use_epoch_astrometry = False`` orvara simulates the Gaia solution by
sampling the model sky path over a window of length ``T`` centred on
``epoch_ra_gaia`` / ``epoch_dec_gaia``. For a one-star NSS file ``T`` is the NSS
model's own ``2 x 517.5 d = 2.834 yr`` (34 months); for a full HGCA catalogue it
is the 2.76 yr the HGCA assumes.

- ``npar`` 5 or 7: three points per coordinate at ``-T/2, 0, +T/2``. The
  proper motion is ``(f(+T/2) - f(-T/2)) / T`` -- the mean proper motion over
  the window, exactly what the NSS model reports -- and the acceleration is the
  quadratic difference. Jerk stays zero and carries no penalty.
- ``npar`` 9: four points per coordinate at ``-T/2, -T/6, +T/6, +T/2``. The
  proper motion is again the mean over the window, ``(f(+T/2) - f(-T/2)) / T``
  (the NSS ``pmra`` for an Acceleration9 solution is defined as the mean, not the
  value at 2016.0; the two differ by ``jerk * (T/2)^2 / 6``); the acceleration at
  2016.0 and the jerk come from four-point stencils exact for a cubic. The
  Gaia-side position that enters the Hipparcos-Gaia proper motion is Simpson's
  3/8 rule over the four points, the time-average of the path, matching the
  Simpson-rule average the three-point scheme uses.

Earlier versions raised a ``RuntimeError`` for 9-parameter sources with
``use_epoch_astrometry = False``. That restriction is gone.

corr_vec
--------
``corr_vec`` must have 21 entries for ``npar = 7`` and 36 for ``npar = 9``,
comma- or whitespace-separated. Gaia stores the off-diagonal correlations once,
in the order ``(0,1), (0,2), (1,2), (0,3), ...`` -- the upper triangle read
column by column. The archive column is zero-padded to a fixed length; a longer
vector is truncated to the first ``npar*(npar-1)/2`` entries. A shorter one is
an error, and ``orvara-config --check`` will tell you so before you run anything.


Orbit fitting
=============

Configuration
-------------
First, assign the appropriate file directories and settings inside of a
config.ini file. Example configuration files can be found in
``orvara/tests/`` e.g. ``orvara/tests/config.ini``, and annotated blank
templates come from ``orvara-config --template fit``.

One file holds everything for a target: the fit settings and, for a Gaia NSS
source, the ``[gaia_accel]`` block described above. ``fit_orbit``,
``plot_orbit`` and ``orvara-gaia-accel`` all read the same file.

The file is validated before anything runs. A wrong type, an unusable
combination (``use_epoch_astrometry = True`` with no IAD directories), a missing
required key, or an ``output_fits`` that does not match ``HGCAFile`` are all
reported together, up front, naming the section and key -- rather than surfacing
as a ``KeyError`` or a silently-applied default some way into a run. A key that
is not recognised gets a warning with a "did you mean" suggestion, since it is
almost always a typo. To check a file without running a fit:

.. code-block:: bash

    orvara-config --check my_target.ini

``fit_orbit`` and ``plot_orbit`` accept ``--strict`` (unknown sections and keys
become errors) and ``--no-validate`` (skip the checks entirely).

Epoch Astrometry
~~~~~~~~~~~~~~~~
If you are using epoch astrometry (observational epochs and scan angles), you
must give paths for ``GaiaDataDir``, ``Hip1DataDir``, and
``Hip2DataDir``. Those are the paths of the intermediate epoch and scan
angles for Gaia/GOST in csv format, the intermediate astrometry from the
original Hipparcos data reduction, and from the second Hipparcos data
reduction.

We recommend starting with ``use_epoch_astrometry = False``. If this
fails, then there is something wrong with the RVFile, HGCAFile, or (relative)
AstrometryFile. If that chain finishes fine, then set ``use_epoch_astrometry = True``.


It is preferred if your data come from the Java Tool (published in 2021). You will want to point orvara to the directory
that contains the H00, H01, etc. folders. It is also OK to have all the IAD in one folder, in which case, point orvara
to that large folder.
For more information about the Java Tool IAD, see:
https://www.cosmos.esa.int/web/hipparcos/interactive-data-access and https://www.cosmos.esa.int/web/hipparcos/catalogues
you can download the data there as a zip file. See https://arxiv.org/abs/2109.06761 Brandt et al. 2021 for
the explanation as to why the Java Tool data are preferred. We recommend only using the Java Tool IAD. If you only use
the Java Tool IAD, then ignore the following three warnings.

DVD IAD Warning 1: If your Hip2 intermediate data come from the DVD, you will
want to point to the 'resrec' folder. This should be e.g.:
``Hip2_DVD_Book/IntermediateData/resrec``. You can keep the original
directory structure as it is.

DVD IAD Warning 2: If you have fits of both the Java Tool and DVD IAD, then
it is imperative that you keep those different IAD in different directories. The DVD IAD and Java tool IAD have
slightly different data formats. Orvara will abort if you have more than one IAD
file with the same hipparcos id. (E.g., a folder with both H003850.d and HIP003850.d)

DVD IAD Warning 3: As described in https://arxiv.org/abs/2109.06761 Brandt et al. 2021,
there are roughly 6600 sources whose DVD IAD is suspected to be corrupted. This corruption is automatically
fixed by htof if you use the Java Tool IAD.

Setting priors
~~~~~~~~~~~~~~
Adding Gaussian mass priors are supported. As well, bounded log-uniform priors
on the RV jitter are supported. To add a Gaussian mass prior of 1 ± 0.1 M_sun
on the primary, you will want to add the following section to your
configuration file:

.. code-block:: bash

    [priors_settings]
    mpri = 1
    mpri_sig = 0.1
    minjitter = 1e-5
    maxjitter = 1e3

``minjitter`` and ``maxjitter`` are the lower and upper bounds respectively for
the log-uniform prior. E.g. If you wanted to set a log-uniform prior on jitter
between 1 and 300 meters/second (and zero outside that prior), you would set
``minjitter = 1`` and ``maxjitter = 300``.

``mpri`` is the mean of the gaussian prior on the primary mass. ``mpri_sig`` is
the standard deviation of that distribution. E.g. for our example of a gaussian
prior on the primary mass of 1 solar mass and 0.1 solar mass deviation, you set
``mpri = 1`` and ``mpri_sig = 0.1``.

Leaving ``mpri_sig = inf`` will turn off the prior and default to a 1/m prior.

Fits without using HGCA (hipparcos and gaia) absolute astrometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can fit an object that does not exist in the HGCA. Simply comment out the ``HipID`` line
in the configuration file. E.g.,:

.. code-block:: text

    . more stuff above...
    .
    # In that case you must supply a parallax and parallax_error in [priors_settings].
    #HipID = 3850
    #^ # note that we have commented the hip id out!!!
    .
    . more stuff below...


And set a parallax and parallax error in the priors section. So:

.. code-block:: text

    . more stuff above...
    .
    [priors_settings]

    parallax=53.05263
    parallax_error=2.818e-02
    .
    . more stuff below...

And that is it! The fit will proceed as normal, but gaia and hipparcos
epoch astrometry *will not* constrain the orbit.


Fits without subsets of data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can do fits with orvara to any subset of the three data sources:
1. Radial velocities
2. Absolute astometry (via the HGCA)
3. Relative astrometry.

E.g., you can do fits without the absolute astrometry (so just 1. and 3.), or a fit to purely radial velocities. Etc.
Below are descriptions with how to do fits to a source without a specific set of data.

1. To do a fit without radial velocity data: In the config file, set ``RVFile = None``
2. To do a fit without Absolute astrometry: In the config file, set ``HipID = None``
3. To do a fit without Relative Astrometry: ``AstrometryFile = None``

So if I wanted to do a fit to just radial velocities, I would set ``HipID = None`` and ``AstrometryFile = None``, and
``RVFile = path/to/the_rv_data.dat``.

Starting conditions
~~~~~~~~~~~~~~~~~~~
You can set the initial conditions of the orbit (starting parameters of the
walkers of the MCMC chain) via a starting file specified in the config.ini
file. There is an example starting file provided in ``orvara/tests/start_Gl758.dat``.

Usage
-----
After setting paths and MCMC (markov-chain monte-carlo)  settings in a
config.ini file, you fit an orbit by running the following from the command
line.

.. code-block:: bash

    fit_orbit --output-dir /path/to/output path/to/config.ini

If you do not specify an output directory using ``--output-dir``, then
orvara will write its output files to the current working directory. One can
set the number of threads in the config.ini file via ``nthreads``. Note
that the built-in parallelization is poor, and that parallelization might not
work at all on some systems. If you get an error when running the code check to
see if it goes away when setting ``nthreads=1``. It is often better to set
``nthreads`` to ``1`` then simply run multiple instances of orvara on separate
cores.


You can access the help menu with the --help flag as follows.

.. code-block:: bash

    fit_orbit --help


Input Data
~~~~~~~~~~

Orvara can fit orbits to a combination of input data. Accepted are (any combination of): relative radial velocities
(the radial velocity of the planet minus the radial velocity of the star), radial velocities, relative astrometry.

See orvara/tests/data for examples of these input data files. Be sure that the units are consistent
with those presented in the examples.

Orvara will automatically include absolute astrometric constraints when available. If the star is not inside of the
Hipparcos Gaia Catalog of Accelerations, these data will be excluded in the fit.

Output
~~~~~~
The output of the MCMC is a .fits file contained within your given output
directory (default current working directory). The name of the file will be
given by the ``target`` field within the ``[plotting]`` field of the
config file, with ``_chain000.fits`` appended (the integer will be
incremented to the next free filename). If this field is empty, the backup
target name will follow the form ``HIP%d_chain%03d.fits``.

The output file contains two .fits extensions: an empty one, and a fits table
with all the MCMC parameters at each step of the chain, together with the log
of the (unnormalized) posterior probability and some additional quantities.

HDU0: empty
"""""""""""
The first extension is empty for table data.

HDU1: table
"""""""""""
This is a fits table object. Each table column is of shape ``(nwalkers,
nsteps/thin)`` where ``thin`` is the thinning used in the configuration file (default
50, to save every 50th step). You may access a column by, e.g.,

``lnlike = hdulist[1].data['lnp']``

The column names and descriptions are:

    'mpri' : Primary mass (Solar masses)

    'msec0' : Secondary mass of the first (index 0) companion, Solar masses

    'sau0' : Semimajor axis of the first companion, AU

    'esino0' : sqrt(ecc)*sin(omega) for the first companion

    'ecoso0' : sqrt(ecc)*cos(omega) for the first companion

    'inc0' : inclination (radians) for the first companion

    'asc0' : PA of the ascending node (radians) for the first companion

    'lam0' : Mean longitude at reference epoch for the first companion

If there is more than one companion, then there are additional fields with,
e.g., 'msec1', 'msec2', etc.

    'jitter' : log RV jitter in m/s

    'jitter0' : log RV jitter for instrument 0

Note that 'jitter0', 'jitter1', etc. are present and 'jitter' is not if using
one jitter per instrument. The default is to use the same jitter for all
instruments. In this case 'jitter' is present but 'jitter0', 'jitter1', etc.
are not.

    'lnp' : natural log of the (unnormalized) probability. Note that this includes matrix determinants and is not simply chi squared.

    'plx_ML' : maximum likelihood (ML) parallax at this chain step

    'pmra_ML' : ML proper motion in RA at this chain step

    'pmdec_ML' : ML proper motion in Dec at this chain step

    'chisq_sep' : The chi squared in separation at the ML parallax at this chain step

    'chisq_PA' : The chi squared in position angle at this chain step

    'chisq_H' : The chi squared for the two Hipparcos proper motions

    'chisq_HG' : The chi squared for the two long-term Hipparcos-Gaia proper motions

    'chisq_G' : The chi squared for the two Gaia proper motions

    'chisq_relRV' : The chi squared for the relative RV data (RV planet - RV star). This will be zero if no relative RV data is used. It will be "None" if the chain you are plotting is from an older version of orvara which did not have relRV data support.

    'RV_ZP_0_ML' : The ML zero point (barycenter RV) for instrument 0

There will be an 'RV_ZP_1_ML' for instrument 1, etc., up to the number of RV instruments.

NOTE: chisq_RV, the chisquared of the fit to the RV data *is not saved*. This is for the following reason. Orvara
uses a jitter term such that the reduced chisquared ( equal to (chisquared of RV fit)/(number of RV data points) ) is 1. So
if you want chisq_RV, it is simply a value very close to the number of RV data points.

If you want an overall absolute astrometric chi squared, you would add the
values from items 'chisq_H', 'chisq_HG', and 'chisq_G' above. There are
effectively four measurements since the mean proper motion of the system was
fit ('pmra_ML' and 'pmdec_ML').

For instance, displaying ``hdulist[1].data['plx_ML']`` will show all the
walkers for the parallax chain (however this parameter is marginalized over in
orvara, it is not fit). ``numpy.mean(hdulist[1].data['plx_ML'][:, burn:])``
and ``numpy.std(hdulist[1].data['plx_ML'][:, burn:])`` would give the mean
and standard deviation of the maximum likelihood parallax (with burn = some
integer that is the number of steps/thinning factor that you are discarding as
burn in). You should add the measurement error of the parallax in quadrature
with this particular uncertainty.

One can use the 'lnp' column to compare the likelihoods of the best orbits if a
certain posterior is multimodal. Assume that the marginalized posterior in PA
is multimodal, with a mode at a value > 180 degrees, and a mode at a value less
than 180 degrees. The following code would print the likelihoods of the maximum
likelihood orbits at each of those two modes.


.. code-block:: python

    tt = fits.open('my/path/chain.fits')[1].data
    logl = tt['lnp']
    pa_data = (tt['asc0']*180/np.pi) % 360
    ls180 = pa_data < 180
    g180 = pa_data >= 180

    print(np.max(logl[ls180]))  # the max log likelihood of all orbits with PA of ascending nodes < 180 degrees
    print(np.max(logl[g180]))  # the max log likelihood of all orbits with PA of ascending nodes > 180 degrees


Examples
--------
To run a quick test using the test data and test config.ini in orvara/tests,
you could cd to the root directory of orvara, then run the following

.. code-block:: bash

    fit_orbit orvara/tests/config.ini --output-dir ~/Downloads

This will create a .fits file in the Downloads folder. The MCMC should
terminate in less than one second because of the short number of steps
indicated in the example config file.

The end-to-end tests in test_e2e check that the code is converging to
previously accepted values for HIP3850. If you wanted to run the code yourself
on this test case and check the results yourself against those in
misc/Diagnostic_plots.ipynb, you can run:

.. code-block:: bash

    fit_orbit orvara/tests/diagnostic_config.ini --output-dir ~/Downloads

The diagnostic_config.ini has the same parameters as those used to create the
plots in Diagnostic_plots.ipynb

Plotting
========

Usage
-----
Once a .fits file from the output of the MCMC is generated, you can produce
several plots of an orbit by running the following in the command line in the
root directory of the repo. To do this, specify the path to the .fits MCMC
output file within the configuration file.

.. code-block:: bash

    plot_orbit path/to/config.ini --output-dir /path/to/output

You can access the help menu with the --help flag as follows.

.. code-block:: bash

    plot_orbit --help

Configuration
-------------
Main plots orvara is configured to produce from the orbital fit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Astrometry orbit of the companion
2. Radial Velocity (RV) orbit over an extended time baseline
3. RV orbit over the observed baseline with O-C
4. Relative separation of the two companions
5. Position angle between the two companions
6. Astrometric acceleration or proper motion fit to Hipparcos-Gaia Astrometry
7. A density plot showing the predicted position at a chosen epoch
8. An on-sky predicted position of the planet, at any given epoch. A plot is produced, and the predicted ra, dec (in mas, with uncertainties) is printed to console.


To generate any of these plots, simply set the corresponding parameters under
the ``[plotting section]`` in the config.ini file to a boolean variable ``True``. If
``False``, a plot would not be produced. 

Here, for 1. Astrometry orbit plots, you can modify the ``predicted_years``
parameter to plot random predicted epoch positions on the Astrometry plot. 

For 2. RV orbit of the companion, you can choose to plot a specific instrument
(by name) or all of the RV instruments by changing the ``Relative_RV_Instrument``
parameter to either the name of the instrument or ``All``. 

For 6. Proper motion plots, you can plot the proper motions
in RA and DEC in one plot (``Proper_motion_separate_plots = False``) or two
(``Proper_motion_separate_plots = True``). In general, you can also set a
customized range of epochs you want to plot, as well as number of orbits
sampled from the posterior distributions and the resolution (step size).

Other outputs
~~~~~~~~~~~~~
In addition to the six plots, you can check convergence of fitted parameters in
the HDU1 extension by setting the parameter ``check_convergence`` to ``True``. You can
define the length of the burn-in phase, note that the parameters are sampled
every ``thin`` steps (as set in the configuration file; default 50). And
you can save the results from the fitted and inferred parameters from the HDU1
extension with ``save_params = True`` in the ``[save_results]`` section, with
an option of setting the quantiles for the uncertainties.

Color bar settings
~~~~~~~~~~~~~~~~~~
To color-code orbits with a key, choose a colormap from the matplotlib list of
colormaps (default viridis) and a reference scheme for the colorbar. Three
reference schemes are available: the eccentricity as ``ecc``, the secondary
companion in jupiter mass as ``msec_jup`` and the secondary companion in solar mass
as ``msec_solar``. Use ``use_colorbar`` to toggle the colorbar key on and off
with ``True`` or ``False``.

Multiple Keplerian orbit fits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the case of a 3-body or multiple-body fit, you can plot the results for each
companion by setting ``iplanet`` to the corresponding companion ID used in the
fitting. ``iplanet`` starts from ``0``.


Examples
--------
To plot orbits, run the plot_orbit command. We can use HD4747 as an example.
First, run the orvara fit, then the plotting, e.g.:

.. code-block:: bash

    fit_orbit --output-dir ./orvara/tests/chains orvara/tests/config_HD4747.ini

    plot_orbit --output-dir ./orvara/tests/plots orvara/tests/config_HD4747.ini

plot_orbit will plot the chain given in the ``McmcDataFile`` in the config
file. If you want the results from a different chain (e.g. a second, different
fit of the same star) you must update this parameter to the new file.


Contribution Guidelines
=======================
We encourage contributions to orvara. The workflow for contributing is the following.

First time contributers:
 * Fork the repository
 * Checkout a new branch for your feature or bug fix.
 * Make your changes to that branch.
 * When you are ready to submit a pull request into the main orvara branch
   (currently called master), run ``pytest -sv`` to make sure that the
   required tests pass.
 * If the tests pass, submit your pull request.
 * One approving administrator review is required to approve a pull request.

Users who are invited to be collaborators on the repo:
 * The same as above, except there is no need to fork the repository once you
   accept your invite!

Citation
========
Placeholder FIXME

License
=======
BSD 3-clause license
