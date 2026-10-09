1.2.1 (2026-09-18)
------------------

- **`[data_paths] relast_error_scale`**: one factor that multiplies every
  separation and position-angle uncertainty in `AstrometryFile` when it is
  read (default 1). It replaces the practice of rewriting the relative-astrometry
  file by hand between passes; the value is printed at the start of the fit and
  kept, with the rest of the configuration, in the chain-file header. The
  recommended use is iterative: set it to `sqrt(chi2_rel / 2 N_rel)` of the
  previous fit's maximum-likelihood step, cumulatively, until the ratio is
  about unity (see USAGE.md, Section 7).
- `fit_orbit` now logs the two error conventions that enter the likelihood at
  start-up: the NSS acceleration/jerk inflation recorded in the one-star HGCA
  file (`ERRINFL`, written by `orvara-gaia-accel`) and the relative-astrometry
  scale above.
- `Data.relsep_err` and `Data.PA_err` are readable from Python (they were
  private Cython attributes), so the scaling is testable; four tests added
  (`test_relast_error_scale.py`). 63 tests pass.
- Verified on 2026-10-05: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1, Astropy 8.0.1.
  59 unit tests and 4 end-to-end tests pass. Sdist and wheel pass `twine check` and
  install into a fresh Python 3.13 environment. The 05231-3142 fit of An et al.
  reproduces from the release install.

1.2.0 (2026-08-28)
------------------

### 2026-09-11 review pass

An adversarial review of the 1.2.0 changes against the Gaia DR3 documentation
(Sect. 7.2.4, Eq. 7.5, the NSS acceleration model) turned up one significant
problem and several smaller ones. All are fixed here; 59 tests pass.

- **The Gaia window was centred at the wrong epoch for NSS sources.** The
  generator (inherited from the original `gaia_accel_gen_fit.py`) wrote
  `epoch_ra_gaia = 2016 - corr(ra,pmra) sigma_ra / sigma_pmra`, the HGCA's
  "uncorrelated position epoch" for a 5-parameter source. The DR3 acceleration
  model is written in `t = epoch - 2016.0`, its proper motion is the *mean* over
  `2016.0 +/- 517.5 d`, and its acceleration and jerk are the derivatives at
  2016.0. Centring the model window elsewhere biases the model proper motion by
  `acc * dt` -- for the shipped J2043 example, (-0.36, +1.69) mas/yr against a
  0.027 mas/yr error. A Gaia-only fit absorbs this into the free barycentric
  proper motion, so `gaia_only` results were not affected beyond a negligible
  `jerk * dt` shift in the acceleration; a `hip_and_gaia` fit, where Hipparcos
  pins the barycentric motion, could not absorb it and would have shown it as
  an inflated chi^2_G. The window is now centred on 2016.0.
- **The window length for NSS sources is the model's own.** `identify_hgca_version`
  returns `2 x 517.5 d = 2.834 yr` (34 months) for a one-star NSS file, instead
  of the 2.76 yr the HGCA_vEDR3 catalogue assumes. The difference is under 3 %
  and changes nothing measurable; it makes the model average over the window
  the catalogue parameters are defined on.
- **The 9-parameter proper motion is now the mean over the window**,
  `(f(+T/2) - f(-T/2)) / T`, which is how the NSS `pmra` is defined for an
  Acceleration9 solution, rather than the instantaneous derivative at 2016.0.
  The two differ by `jerk (T/2)^2 / 6`, ~0.2-0.3 sigma for a typical jerk.
- Configuration keys are case-insensitive again (the 1.2.0 validator preserved
  case, which silently turned a lower-case `hipid = 3850` into `HipID = 0`).
  `nplanets` and `nstep` are now required for a fit, since `fit_orbit` reads
  them without a fallback.
- `orvara-gaia-accel` resolves `output_fits` relative to the working directory,
  as `fit_orbit` resolves `HGCAFile`, so the two paths can be kept identical.
- The proper-motion figure drops the `gaia_only` placeholder Hipparcos point
  (its 1e7 mas/yr error bar set the y-range to +/- 1e7) and evaluates the Dec
  window about `epoch_dec_gaia`.
- The `gaia_npar not in (5, 7, 9)` check was inside a `try/except` that
  swallowed it; moved out.
- Packaging: the sdist now ships `orbit.pyx` (and not the generated `orbit.c`),
  so a wheel can be built from it; licence metadata uses the SPDX form. Verified
  by building sdist + wheel and installing from the sdist on Python 3.10, 3.12
  and 3.13.
- A `conftest.py` explains, instead of a bare ImportError, that the test suite
  runs against an editable install.


This release merges the two Gaia-acceleration branches (`orvara_gaia_accel_Acc9`
and `orvara_hip_and_gaia`) into one package, and brings the whole thing up to
Python 3.10-3.13 with NumPy 2.x.

### Gaia DR3 NSS acceleration support

- Gaia 9-parameter (acceleration + jerk) solutions are now fitted without epoch
  astrometry. `calc_PMs_no_epoch_astrometry` samples the Gaia baseline at four
  points (-T/2, -T/6, +T/6, +T/2) instead of three when `gaia_npar == 9`, and
  recovers proper motion, acceleration and jerk with finite-difference formulas
  that are exact for a cubic. 5- and 7-parameter sources keep the original
  3-point quadratic scheme, and their jerk stays zero.
- The old `RuntimeError` refusing 9-parameter sources with
  `use_epoch_astrometry = False` is gone, because that case now works.
- **Bug fix.** In the 9-parameter branch the Gaia-side position used for the
  Hipparcos-Gaia proper motion was the plain mean of the four samples, which
  weights the acceleration term by `5h^2/9`. It is now Simpson's 3/8 rule,
  `(f0 + 3 f1 + 3 f2 + f3)/8`, giving the correct time-average `p + (a/2)h^2/3`
  and matching what the 3-point branch computes for 5- and 7-parameter sources.
- `orvara.gaia_accel` (new subpackage, was the standalone `gaia_acceleration/`
  directory) builds the one-star HGCA-format FITS file from a Gaia NSS solution,
  via the `orvara-gaia-accel` console script. Two modes:
  - `gaia_only` -- no Hipparcos data; the Hip and Hip-Gaia rows get the Gaia
    proper motion with effectively infinite errors, so only the NSS solution
    constrains the fit.
  - `hip_and_gaia` -- the star is in the HGCA; its real Hipparcos and Hip-Gaia
    proper motions are copied from the full catalogue and the NSS acceleration
    and jerk terms are attached to them.
- The generator no longer needs `HGCA_vEDR3.fits` as a column template in
  `gaia_only` mode: the HGCA schema is defined in code. The template is still
  used, and required, for `hip_and_gaia`.
- Optional RV and relative-astrometry CSVs are written both as extra FITS
  extensions and as the whitespace-delimited text files `fit_orbit` reads.
- `identify_hgca_version` replaces a bare `try/except` around the HGCA version
  sniffing. A one-star NSS file (no HIP 1 row, but a `gaia_npar` column) is now
  recognised as the normal case and reported quietly; a file that is neither
  gets a warning that says what is wrong.

### One validated configuration file

- Everything for a target lives in one `.ini`: the fit settings and, in a new
  `[gaia_accel]` section, the NSS solution. `fit_orbit`, `plot_orbit` and
  `orvara-gaia-accel` all read the same file.
- `orvara.config` gained a schema. Wrong types, bad values, missing required
  keys, `use_epoch_astrometry = True` without the IAD directories, and an
  `output_fits` that does not match `HGCAFile` are reported together, up front,
  with section and key named -- instead of surfacing as a `KeyError` or a
  silently-applied default halfway through a run. Unknown keys are reported as
  warnings with a "did you mean" suggestion; `--strict` makes them errors.
- The legacy `[paths]/[star]/[advanced]/[gaia_solution]` prep files still work
  unchanged.
- New console script `orvara-config`: `--template {fit,acc7,acc9,hip_and_gaia}`
  writes an annotated starter config, `--check FILE` validates an existing one.
- `fit_orbit` and `plot_orbit` accept `--strict` and `--no-validate`.

### Python 3.10-3.13 and NumPy 2.x

- **Critical.** `orvara.orbit.Data` assigned length-1 arrays from the HGCA table
  to C `double`/`int` attributes. NumPy 2.0 removed that conversion, so every
  fit raised `TypeError: only 0-dimensional arrays can be converted to Python
  scalars`. The matching row is now taken as a record, giving NumPy scalars.
- Removed the `np.float` / `np.int` / `np.str` aliases that had been
  reintroduced in `orbit_plots.py` and `corner_modified.py` (removed in NumPy
  1.24).
- `ax.get_shared_x_axes().join(...)` replaced with `Axes.sharex` (removed in
  Matplotlib 3.8).
- `pkg_resources` replaced with `importlib.metadata` (`pkg_resources` is not
  installed by default on Python 3.12+).
- ptemcee is now bundled in `orvara/_vendor/ptemcee` (MIT, Will Vousden). The
  PyPI release has had no update since 2019 and calls `np.float`, so it imports
  cleanly and then dies in the middle of a long MCMC. Only the six
  `np.float` -> `float` changes were made; the algorithm is untouched. ptemcee
  is no longer a dependency.
- The dead `emcee.PTSampler` fallback path (emcee 2) is now selected by feature
  detection rather than by catching an exception.
- Packaging moved from `distutils`/`setup.py` to `pyproject.toml` with a
  setuptools + Cython 3 build. `python_requires >= 3.10`.
- Every bare `except:` in the package is now `except Exception:`.

### Other bug fixes found while modernising

- **The acceleration O-C panel plotted a model acceleration of exactly zero.**
  The panel built its own `Data`/`Model` at three epochs around the Gaia central
  epoch and called `calc_EA_RPP`, but not `calc_offsets` -- and `calc_offsets` is
  what writes `dRA_G`/`dDec_G`. The offsets therefore stayed at the zeros set in
  `Model.__init__`, so the plotted model acceleration was identically zero and
  the O-C equalled the full measured acceleration. The likelihood was never
  affected; this was a plotting fault only. The panel now goes through the new
  `OrbitPlots.model_accel_jerk`, which samples the sky path at the same epochs
  the likelihood uses and reduces them with the same finite-difference
  coefficients -- three points for a 5- or 7-parameter source, four for a
  9-parameter one -- so the panel and chi^2_G describe the same quantity, and a
  9-parameter source now gets a model jerk as well as a measured one.

- `orbit_plots.load_HGCA_data` did `int(np.where(...)[0])`, which NumPy 2
  rejects. Because the whole block sits inside a `try/except`, the failure was
  silent: every proper-motion panel would simply have disappeared from the
  plots. It now indexes the match properly and raises a named `KeyError` when
  the star is genuinely absent.
- `plot_orbit` crashed with `'NoneType' object has no attribute 'split'` on any
  config that omitted `xlim`/`ylim` — they were parsed unconditionally with
  `fallback=None` even though they are only used when `set_limit = True`. They
  now default to empty, and `set_limit = True` without them is a clear error.
- `ErfaWarning` is imported from `erfa`, not `astropy.utils.exceptions`, which
  deprecated the re-export in 6.1.
- `make_header` truncates long values to fit an 80-character FITS card and no
  longer leaves the config file handle open. This matters now that `corr_vec`
  (several hundred characters) can appear in the config.
- `Data.Cinv_G_acc_terms` and `Data.Cinv_G_jerk_terms` are public, so the
  acceleration and jerk weighting can be inspected from Python.
- Added `orvara/tests/__init__.py`, which was missing, so `orvara.tests` is a
  real package rather than an implicit namespace one.

### Tests

- `orvara/tests/test_finite_differences.py`: feeds an exact cubic to
  `calc_PMs_no_epoch_astrometry` and checks the recovered proper motion,
  acceleration and jerk, plus the centroid fix above. `Model.set_dRA_dDec` was
  added to make that path testable.
- `orvara/tests/test_gaia_accel.py`: `corr_vec` handling, both build modes, a
  regression check against the previously generated `J2043_gaia_accel.fits`,
  and a round trip through `orvara.orbit.Data`.
- `orvara/tests/test_config.py`: schema validation, the shipped templates, and
  backwards compatibility with the legacy prep files.
- `orvara/tests/test_plotting.py`: `plot_orbit` had no tests at all. It now has
  smoke tests that run the real entry point on the tutorial chain, including
  the no-`xlim` case above.
- `orvara/tests/test_accel_forward_model.py`: asserts that `calc_offsets` is
  what populates the Gaia offsets (the bug above), and closes the loop on the
  forward model -- fit a synthetic Gaia-only Acceleration9 source, then
  re-derive the acceleration and jerk from the posterior sky path and confirm
  they reproduce the NSS coefficients the fit was constrained by. On the shipped
  Acceleration9 example a 5 x 50 x 6000-step chain returns posterior medians of
  (-0.614, +0.939) mas/yr^2 against a measured (-0.613 +/- 0.064,
  +0.941 +/- 0.056) and (-0.011, +0.021) mas/yr^3 against (-0.014 +/- 0.008,
  +0.020 +/- 0.007), with chi^2_G = 1.0 at the maximum-likelihood step; the
  Acceleration7 example returns (-2.456, +7.714) against (-2.452 +/- 0.139,
  +7.710 +/- 0.101).

### Results are unchanged

The fit was checked against the pre-modernisation code (`orvara_hip_and_gaia`
on Python 3.10 / NumPy 1.23 / Astropy 5.2 / ptemcee from PyPI), running
`diagnostic_config_fake_fulldr3_7pfit.ini` in both:

| quantity | before | after |
|---|---|---|
| RV jitter (m/s)   | 5.00 ± 0.62      | 4.99 ± 0.62      |
| M_comp (M_Jup)    | 66.06 ± 1.56     | 66.19 ± 1.55     |
| a (AU)            | 9.984 ± 0.099    | 9.992 ± 0.098    |
| eccentricity      | 0.73433 ± 0.00134| 0.73429 ± 0.00133|
| inclination (deg) | 48.72 ± 0.37     | 48.75 ± 0.37     |
| chi2_H / chi2_G   | 0.95 / 4.17      | 0.97 / 4.15      |
| ln P              | -147.86 ± 2.20   | -147.87 ± 2.17   |

Agreement is within the Monte Carlo scatter of a single chain.

One caution for anyone repackaging this: do **not** set `boundscheck`,
`wraparound` or `cdivision` globally in `cythonize()`. `orbit.pyx` applies
`@cython.boundscheck(False)` / `@cython.wraparound(False)` per function, where
the author established they are safe; applying them module-wide makes the RV
jitter run away (5 m/s to ~270 m/s) and biases every orbital element. `setup.py`
sets only `language_level`, and carries a comment saying so.

1.1.4 (2022-08-23)
------------------
- Changed the internals of the astrometric prediction code to use the orvara internals. This is 
more robust and works for more than just 2-body systems.
- A position (ra, dec) and uncertainty is now printed to standard out with the astrometric plot. 

1.1.3 (2022-05-24)
------------------
- Added the ability to constrain orbits with a relative RV measurement (e.g. that between beta pic A and beta
pic b), by providing a relative RV file (`relRVFile=relative_rv_filename.dat`) in the config.ini. This is implemented
programmatically so that it will work for multiple companions (3, 4 etc body fits), i.e., it won't crash.
However the logic to compute the relative RV is technically only valid in detail for 2-body systems, or multi-planet
systems where the primary is substantially more massive than any companion. See the note within `def calc_relRV` inside
of orbit.pyx.
- the chisquared of the relative RV fit is now saved (in the same way as the hipparcos/gaia chisquareds).

1.1.2 (2022-04-28)
------------------
- Bugfix for plotting proper motions. When the internals changed to mas/yr with 1.1.1, the plotting 
was still done in units of mas/day.

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
