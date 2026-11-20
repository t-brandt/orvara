"""
orvara -- Orbits from Radial Velocity, Absolute, and/or Relative Astrometry.

This build adds support for Gaia DR3 non-single-star (NSS) astrometric
acceleration solutions: 7-parameter (acceleration) and 9-parameter
(acceleration + jerk) fits, with or without epoch astrometry.
"""

try:
    from importlib.metadata import version as _version, PackageNotFoundError
    try:
        __version__ = _version("orvara")
    except PackageNotFoundError:          # running from a source checkout
        __version__ = "1.2.1"
except ImportError:                       # pragma: no cover
    __version__ = "1.2.1"

__all__ = ["__version__"]
