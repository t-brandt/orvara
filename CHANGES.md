1.1.0 (2022-XX-XX)
------------------
- Added Seven and nine parameter fits in preparation for Gaia DR3, with tests in test_main.py.
- Added informative comments to the log likelihood function.

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
