"""
Gaia DR3 non-single-star (NSS) astrometric acceleration support for orvara.

The Gaia archive publishes ``gaiadr3.nss_acceleration_astro`` solutions with
either 7 parameters (position, parallax, proper motion, acceleration) or 9
parameters (the same plus the derivative of the acceleration, i.e. jerk).
orvara consumes absolute astrometry through an HGCA-format table, so this
subpackage turns one NSS solution into a one-star HGCA-format FITS file that
``fit_orbit`` can read.

Two modes:

``gaia_only``
    No Hipparcos data exist for the target.  The Hipparcos and Hip-Gaia rows
    are filled with the Gaia proper motion and given effectively infinite
    errors, so only the Gaia NSS acceleration (and jerk) constrain the fit.

``hip_and_gaia``
    The target *is* in the HGCA.  Its real Hipparcos and Hip-Gaia proper
    motions are copied across from the full catalogue and the NSS acceleration
    and jerk terms are attached to them, so both constrain the fit together.

See :func:`orvara.gaia_accel.gen_fit.write_nss_fits`.
"""

from orvara.gaia_accel.gen_fit import (              # noqa: F401
    GAIA_REF_EPOCH,
    GAIA_BASELINE_YR,
    HIP_REF_EPOCH,
    corr_vec_to_matrix,
    build_row,
    write_nss_fits,
    build_from_config,
)

__all__ = ['GAIA_REF_EPOCH', 'GAIA_BASELINE_YR', 'HIP_REF_EPOCH',
           'corr_vec_to_matrix', 'build_row', 'write_nss_fits',
           'build_from_config']
