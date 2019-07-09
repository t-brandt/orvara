htof
===============

This repo contains orbit3d, the package for fitting orbits of exoplanets.


Installation
------------
orbit3d is built by running :code:`python setup.py build_ext --inplace` while in the the root directory of this repo.

Configuration
-------------
First, assign the appropriate file directories inside of config.py. If you are using relative astrometry, you must
give paths for :code:`GaiaDataDir`, :code:`Hip1DataDir`, and :code:`Hip2DataDir`. Those are the paths
to the intermediate data for GaiaDR2, the original Hipparcos data reduction, and the second Hipparcos data reduction.

Usage
-----
After setting paths in the config.py file, you fit an orbit by running the following from the command line

.. code-block:: python

    python orbitfit_3d.py -output-dir /path/to/output

The number of MCMC (markov-chain monte-carlo) walkers, temperatures, and steps can be set with the appropriate arguments.
For example:

.. code-block:: python

    python orbitfit_3d.py --ntemps 10 --nwalkers 10 --nplanets 1 --nthreads 2

We have also specified the number of planets in the star-system and the number of threads to
parallelize to via the nthreads and nplanets keywords. Not that the built-in parallelization is poor. It is better
to set nthreads to 1 then simply run multiple instances of orbit3d on separate cores. You can access the help menu
with the --help flag.

.. code-block:: python

    python orbitfit_3d.py --help

License
-------

To be decided.