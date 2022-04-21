1.1.2 (2022-04-25)
------------------
- Added the ability to constrain orbits with a relative RV measurement (e.g. that between beta pic A and beta
pic b), by providing a relative RV file (`relRVFile=relative_rv_filename.dat`) in the config.ini. This is implemented
programmatically so that it will work for multiple companions (3, 4 etc body fits), i.e., it won't crash.
However the logic to compute the relative RV is technically only valid in detail for 2-body systems, or multi-planet
systems where the primary is substantially more massive than any companion. See the note within `def calc_relRV` inside
of orbit.pyx.
- the chisquared of the relative RV fit is now saved (in the same way as the hipparcos/gaia chisquareds).

1.1.1 (2022-04-21)
------------------
- Added support for the Java Tool IAD (colloquially called hip21).

1.1.0 (2022-04-21)
------------------
- Added Seven and nine parameter fits in preparation for Gaia DR3, with tests in test_main.py.
- Added informative comments to the log likelihood function.
- Changed the astrometric fitting routine to work natively in units of years, instead of days. This
improved the condition number of the matrix that is solved via SVD.
- Bumped required version of HTOF to 1.1.0. This is needed for the fits with parallax.

1.0.5 (2021-12-06)
------------------
- Fixed a bug in the astrometric orbit plotting.

1.0.4 (2021-05-14)
------------------
- Fixed units bug in periastron time calculation.

1.0.3 (2021-05-14)
------------------
- Another small tweak to secondary priors: ensure 1/M if a given prior is not set.

1.0.2 (2021-05-13)
------------------
- Fixed a memory leak when fitting multiple companions, fixed prior on secondary masses to fully undo a 1/M prior.

1.0.1 (2021-04-28)
------------------
- Fixed the diagnostic test in test_main.py so that it will work on windows. File i/o failed on the
test due to https://github.com/astropy/astropy/issues/7404 

1.0.0 (2021-04-04)
------------------
- Changed format of chain files, incompatible with previous versions.  Support for stars not in HGCA, ability to set one jitter per instrument.  Various small changes and bug fixes.  Renamed repo and orbit3d directory to orvara.

0.1.5 (2021-02-26)
------------------
- Added the ability to place mass priors on the secondary, tertiary, etc.. companions.

0.1.4 (2021-02-26)
------------------
- No changes, incorrect tag.

0.1.3 (2021-02-21)
------------------
- Added documentation.

0.1.2 (2021-01-28)
------------------
- Fixed bug where priors on the primary mass were applied to the primary mass + masses of interior companions.

0.1.1 (2019-08-15)
------------------
- Refactored code and added documentation. Bumped requirement
for HTOF to 0.1.1.

0.1.0 (prior to 2019-08-15)
---------------------------
- Initial release.
