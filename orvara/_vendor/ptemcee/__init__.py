"""
Vendored copy of ptemcee 1.0.0 (parallel-tempered emcee).

ptemcee -- Copyright (c) 2016 Will Vousden, MIT licence.
Upstream: https://github.com/willvousden/ptemcee

Why this is bundled rather than installed from PyPI: ptemcee has had no release
since 2019 and ``ptemcee.Sampler`` calls ``np.float``, which NumPy removed in
1.24.  Installing ptemcee therefore gives a package that imports cleanly and
then raises ``AttributeError: module 'numpy' has no attribute 'float'`` in the
middle of a long MCMC run.  Only those six ``np.float`` -> ``float`` changes
(plus raw-string docstrings, to silence Python 3.12 SyntaxWarnings) have been
made; the sampling algorithm is untouched.

Only ``sampler`` and ``util`` are bundled.  ptemcee's ``InterruptiblePool`` and
``MPIPool`` are not, because orvara does not use them: pass any object with a
``map`` method as ``pool=`` if you want multiprocessing.
"""

from orvara._vendor.ptemcee.sampler import Sampler, default_beta_ladder   # noqa: F401
from orvara._vendor.ptemcee import util                                   # noqa: F401

__version__ = '1.0.0+orvara'
__all__ = ['Sampler', 'default_beta_ladder', 'util']
