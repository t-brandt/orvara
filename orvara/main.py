#!/usr/bin/env python
"""
Orbit fitting code. The run function is the console entry point,
accessed by calling fit_orbit from the command line.
"""

from __future__ import print_function
import numpy as np
import os
import time
import emcee
# ptemcee is bundled (orvara/_vendor/ptemcee); see that package's __init__
# for why the PyPI release cannot be used with NumPy >= 1.24.
from orvara._vendor.ptemcee import Sampler as PTSampler
from astropy.io import fits
from htof.main import Astrometry
from astropy.time import Time
from astropy.coordinates import Angle
import sys
import re
import random
import warnings
from orvara import orbit
from orvara.config import parse_args, load_config
from orvara.format_fits import make_header, pack_cols
from orvara.gaia_accel.gen_fit import GAIA_BASELINE_YR as GAIA_NSS_BASELINE_YR

_loglkwargs = {}

def set_initial_parameters(start_file, ntemps, nplanets, nwalkers, njit=1,
                           minjit=-20, maxjit=20):
    
    par0 = np.ones((ntemps, nwalkers, 2 + 7 * nplanets))

    if start_file.lower() == 'none':
        mpri = 1
        jit = 0.5
        sau = 10
        esino = 0.5
        ecoso = 0.5
        inc = 1
        asc = 1
        lam = 1
        msec = 0.1

        sig = np.ones((ntemps, nwalkers, 2 + 7 * nplanets))*0.5    
        init = [jit, mpri]
        for i in range(nplanets):
             init += [msec, sau, esino, ecoso, inc, asc, lam]
        par0 *= np.asarray(init)

    else:
        init, sig = np.loadtxt(start_file).T
        init[0] = max(2*np.log10(init[0]), minjit) # Convert from m/s to units used in code
        try:
            par0 *= init
        except Exception:
            raise ValueError('Starting file %s has the wrong format/length.' % (start_file))

    #######################################################################
    # Introduce scatter.  Normal in most parameters, lognormal in
    # mass and semimajor axis.
    #######################################################################
    
    scatter = sig*np.random.randn(np.prod(par0.shape)).reshape(par0.shape)
    par0 += scatter
    par0[..., 2::7] = (par0[..., 2::7] - scatter[..., 2::7])*np.exp(scatter[..., 2::7])
    par0[..., 3::7] = (par0[..., 3::7] - scatter[..., 3::7])*np.exp(scatter[..., 3::7])

    # Ensure that values are within allowable ranges.
    
    bounds = [[0, minjit, maxjit],   # jitter
              [1, 1e-4, 1e3],        # mpri (Solar masses)
              [2, 1e-5, 1e3],        # msec (Solar masses)
              [3, 1e-5, 2e5],        # semimajor axis (AU)
              [6, 1e-5, np.pi],      # inclination (radians)
              [7, -np.pi, 3*np.pi],  # longitude of ascending node (rad)
              [8, -np.pi, 3*np.pi]]  # long at ref epoch (rad)
        
    for i in range(len(bounds)):
        j, minval, maxval = bounds[i]
        if j <= 1:
            par0[..., j][par0[..., j] < minval] = minval
            par0[..., j][par0[..., j] > maxval] = maxval
        else:
            par0[..., j::7][par0[..., j::7] < minval] = minval
            par0[..., j::7][par0[..., j::7] > maxval] = maxval
            
    # Eccentricity is a special case.  Cap at 0.99.
    ecc = par0[..., 4::7]**2 + par0[..., 5::7]**2
    fac = np.ones(ecc.shape)
    fac[ecc > 0.99] = np.sqrt(0.99)/np.sqrt(ecc[ecc > 0.99])
    par0[..., 4::7] *= fac
    par0[..., 5::7] *= fac

    # Move jitter to the end, add (shuffled) realizations if needed.
    par0_jitlast = np.zeros((ntemps, nwalkers, par0.shape[-1] + njit - 1))
    par0_jitlast[..., :-njit] = par0[..., 1:]
    for i in range(njit):
        random.shuffle(par0[..., 0].T)
        par0_jitlast[..., -1 - i] = par0[..., 0]

    return par0_jitlast


# Reference epochs of HIP 1 in the released HGCA catalogues.  These are the
# fingerprints used to tell which Gaia data release an HGCA file was built from.
_HGCA_EPOCH_SIGNATURES = (
    # (epoch_ra_gaia of HIP 1, HGCA version label, effective Gaia baseline [yr])
    (2015.60211565, 'GaiaDR2', 1.75),
    (2015.92749023, 'GaiaeDR3', 2.76),
)
# A one-star file built by orvara.gaia_accel from a DR3 NSS acceleration
# solution contains no HIP 1.  Its parameters are means over the observation
# window the NSS model assumes, +/- 517.5 d about 2016.0 (Gaia DR3 documentation
# Sect. 7.2.4), so the model has to average over that same window.
_DEFAULT_HGCA_VERSION = 'GaiaeDR3'
_NSS_GAIA_BASELINE_YRS = GAIA_NSS_BASELINE_YR          # 2 x 517.5 d = 2.834 yr
_DEFAULT_GAIA_BASELINE_YRS = 2.76                      # plain HGCA_vEDR3 files


def identify_hgca_version(HGCAFile):
    """
    Work out which Gaia data release an HGCA-style file corresponds to.

    Returns
    -------
    (version, gaia_mission_length_yrs) : (str, float)
        ``version`` is the label htof expects ('GaiaDR2' or 'GaiaeDR3');
        ``gaia_mission_length_yrs`` is the effective astrometric baseline.

    Notes
    -----
    The released HGCA catalogues are identified from the Gaia reference epoch of
    HIP 1.  A one-star file produced by ``orvara.gaia_accel`` (a Gaia-only NSS
    acceleration solution) contains no HIP 1 row; that is expected, not an
    error, and such files are always on the DR3/eDR3 reference frame.
    """
    with fits.open(HGCAFile, memmap=False) as hdul:
        data = hdul[1].data
        names = list(data.names)
        hip1 = data[data['hip_id'] == 1]

    if len(hip1) == 0:
        if 'gaia_npar' in names:
            # A single-star Gaia NSS acceleration file. Expected; say so quietly.
            print(f"No HIP 1 row in '{HGCAFile}'; treating it as a one-star Gaia "
                  f"NSS file on the {_DEFAULT_HGCA_VERSION} reference frame "
                  f"(NSS window {_NSS_GAIA_BASELINE_YRS:.3f} yr).")
            return _DEFAULT_HGCA_VERSION, _NSS_GAIA_BASELINE_YRS
        else:
            warnings.warn(
                f"HIP 1 is missing from the HGCA file '{HGCAFile}' and the file "
                f"carries no 'gaia_npar' column, so the Gaia data release cannot "
                f"be determined. Assuming {_DEFAULT_HGCA_VERSION} with a "
                f"{_DEFAULT_GAIA_BASELINE_YRS} yr baseline. If you meant to use a "
                f"full HGCA catalogue (HGCA_vEDR3.fits / HGCA_vDR2_corrected.fits) "
                f"this indicates a setup problem.",
                RuntimeWarning)
        return _DEFAULT_HGCA_VERSION, _DEFAULT_GAIA_BASELINE_YRS

    epoch_ra_gaia = float(hip1['epoch_ra_gaia'][0])
    for reference_epoch, version, baseline in _HGCA_EPOCH_SIGNATURES:
        if np.isclose(epoch_ra_gaia, reference_epoch):
            return version, baseline

    raise ValueError(
        f"Cannot match '{HGCAFile}' to DR2, eDR3 or DR3: the Gaia RA epoch of "
        f"HIP 1 is {epoch_ra_gaia}, which matches none of "
        f"{[e for e, _, _ in _HGCA_EPOCH_SIGNATURES]}.")


def _log_error_conventions(HGCAFile, AstrometryFile, relast_error_scale):
    """Print the two error conventions that enter the likelihood, so that a
    fit's log records them: the NSS acceleration/jerk error inflation written
    into a one-star HGCA file by orvara-gaia-accel (header ERRINFL), and the
    scaling applied to the relative-astrometry uncertainties."""
    inflation = None
    try:
        with fits.open(HGCAFile) as hdul:
            inflation = hdul[0].header.get('ERRINFL')
    except Exception:
        pass
    if inflation is not None:
        print('Gaia NSS acceleration/jerk uncertainties inflated by %s (from %s)'
              % (inflation, os.path.basename(HGCAFile)))
    if AstrometryFile:
        print('Relative-astrometry uncertainties in %s scaled by %.4f'
              % (os.path.basename(AstrometryFile), relast_error_scale))


def initialize_data(config, companion_gaia):
    # load in items from the ConfigParser object
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    HGCAFile = config.get('data_paths', 'HGCAFile')
    if not os.path.exists(HGCAFile):
        raise FileNotFoundError(f'No HGCA file found at {HGCAFile}')

    HGCAVersion, gaia_mission_length_yrs = identify_hgca_version(HGCAFile)

    RVFile = config.get('data_paths', 'RVFile', fallback='')
    relRVFile = config.get('data_paths', 'relRVFile', fallback='')
    AstrometryFile = config.get('data_paths', 'AstrometryFile', fallback='')
    GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback='')
    Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback='')
    Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback='')
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)
    relast_error_scale = config.getfloat('data_paths', 'relast_error_scale', fallback=1.0)
    if relast_error_scale <= 0:
        raise ValueError('[data_paths] relast_error_scale must be positive.')
    _log_error_conventions(HGCAFile, AstrometryFile, relast_error_scale)
    data = orbit.Data(HipID, HGCAFile, RVFile, AstrometryFile, relRVFile=relRVFile, companion_gaia=companion_gaia,
                      gaia_mission_length_yrs=gaia_mission_length_yrs,
                      relast_error_scale=relast_error_scale)
    if use_epoch_astrometry and data.use_abs_ast == 1:
        gaia_fit_degree = {5: 1, 7: 2, 9: 3}[data.gaia_npar]
        Gaia_fitter = Astrometry(HGCAVersion, '%06d' % (HipID), GaiaDataDir,
                                 central_epoch_ra=data.epRA_G,
                                 central_epoch_dec=data.epDec_G,
                                 format='jyear', fit_degree=gaia_fit_degree,
                                 use_parallax=True, use_catalog_parallax_factors=True)
        Hip2_fitter = Astrometry('Hip2or21', '%06d' % (HipID), Hip2DataDir,
                                 central_epoch_ra=data.epRA_H,
                                 central_epoch_dec=data.epDec_H,
                                 format='jyear', use_parallax=True, use_catalog_parallax_factors=True)
        Hip1_fitter = Astrometry('Hip1', '%06d' % (HipID), Hip1DataDir,
                                 central_epoch_ra=data.epRA_H,
                                 central_epoch_dec=data.epDec_H,
                                 format='jyear', use_parallax=True, use_catalog_parallax_factors=True)
        hip1_fast_fitter = orbit.AstrometricFitter(Hip1_fitter)
        hip2_fast_fitter = orbit.AstrometricFitter(Hip2_fitter)
        gaia_fast_fitter = orbit.AstrometricFitter(Gaia_fitter)
        data = orbit.Data(HipID, HGCAFile, RVFile, AstrometryFile, 
                          use_epoch_astrometry=use_epoch_astrometry,
                          relRVFile=relRVFile,
                          epochs_Hip1=Hip1_fitter.data.julian_day_epoch(),
                          epochs_Hip2=Hip2_fitter.data.julian_day_epoch(),
                          epochs_Gaia=Gaia_fitter.data.julian_day_epoch(),
                          companion_gaia=companion_gaia, verbose=False,
                          relast_error_scale=relast_error_scale)
    elif data.use_abs_ast == 1:
        hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter = None, None, None
    else:
        hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter = None, None, None
        try:
            data.plx = 1e-3*config.getfloat('priors_settings', 'parallax')
            data.plx_err = 1e-3*config.getfloat('priors_settings', 'parallax_error')
        except Exception:
            raise RuntimeError("Cannot load absolute astrometry. Please supply a prior "
                               "for parallax and parallax_error in the priors_settings area"
                               " of the configuration file.")
    return data, hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter


def lnprob(theta, returninfo=False, RVoffsets=False, use_epoch_astrometry=False,
           data=None, nplanets=1, H1f=None, H2f=None, Gf=None, priors=None, 
           njitters=1):
    """
    Log likelihood function for the joint parameters
    :param theta:
    :param returninfo:
    :param use_epoch_astrometry:
    :param data:
    :param nplanets:
    :param H1f:
    :param H2f:
    :param Gf:
    :return:
    """
    model = orbit.Model(data)
    lnp = 0
    for i in range(nplanets):
        # Note that params.mpri is really the mass contained in the primary + companions inner to the current planet.
        # params.mpri_true is the real mass of the primary. So params.mpri should really be renamed params.interior_mass
        params = orbit.Params(theta, i, nplanets, data.nInst, njitters)
        lnp = lnp + orbit.lnprior(params, minjit=priors['minjit'],
                                  maxjit=priors['maxjit'])

        if not np.isfinite(lnp):
            model.free()
            params.free()
            return -np.inf

        orbit.calc_EA_RPP(data, params, model)
        orbit.calc_RV(data, params, model)
        if data.n_rel_RV > 0:
            # slightly more performant if we ignore the relative RV calculation call entirely if
            # we do not have any relative RVs to calculate.
            orbit.calc_relRV(data, params, model, i)
        orbit.calc_offsets(data, params, model, i)

        chisq_sec = (params.msec - priors[f'm_secondary{i}'])**2/priors[f'm_secondary{i}_sig']**2
        if chisq_sec > 0:  # If chisq_sec > 0 then a prior was set on the mass of the secondary
            # undo the default 1/m prior that is set in orbit.pyx (see lnprior())
            lnp = lnp - 0.5*chisq_sec + np.log(params.msec)
        # else, don't change lnp at all.

        # Free params if we need to cycle through the next companion
        if i < nplanets - 1:
            params.free()

    if use_epoch_astrometry and data.use_abs_ast:
        orbit.calc_PMs_epoch_astrometry(data, model, H1f, H2f, Gf)
    elif data.use_abs_ast:
        orbit.calc_PMs_no_epoch_astrometry(data, model)

    if returninfo:
        return orbit.calcL(data, params, model, chisq_resids=True, RVoffsets=RVoffsets)

    if np.isfinite(priors['mpri_sig']):
        return lnp - 0.5*(params.mpri_true - priors['mpri'])**2/priors['mpri_sig']**2 + orbit.calcL(data, params, model)
    else:
        return lnp - np.log(params.mpri_true) + orbit.calcL(data, params, model)

    
def return_one(theta):
    return 1.


def avoid_pickle_lnprob(theta):
    global _loglkwargs
    return lnprob(theta, **_loglkwargs)


def get_priors(config):
    priors = {}
    priors['mpri'] = config.getfloat('priors_settings', 'mpri', fallback=1.)
    priors['mpri_sig'] = config.getfloat('priors_settings', 'mpri_sig', fallback=np.inf)
    # priors on the masses of the companions (labelled 0 though 9). Limit to 10 planet systems.
    for i in range(10):
        priors[f'm_secondary{i}'] = config.getfloat('priors_settings', f'm_secondary{i}', fallback=1.)
        priors[f'm_secondary{i}_sig'] = config.getfloat('priors_settings', f'm_secondary{i}_sig', fallback=np.inf)
    # priors on the RV jitter. Converting from m/s to orvara internal units.
    priors['minjit'] = config.getfloat('priors_settings', 'minjitter', fallback = 1e-5)
    priors['minjit'] = max(priors['minjit'], 1e-20) # effectively zero, but we need the log
    priors['minjit'] = 2*np.log10(priors['minjit'])
    priors['maxjit'] = config.getfloat('priors_settings', 'maxjitter', fallback = 1e3)
    priors['maxjit'] = 2*np.log10(priors['maxjit'])
    if priors['maxjit'] < priors['minjit']:
        raise ValueError("Requested maximum jitter < minimum jitter")
    return priors


def get_gaia_catalog_companion(config):
    companion_gaia = {}
    companion_gaia['ID'] = config.getint('secondary_gaia', 'companion_ID', fallback=-1)
    companion_gaia['pmra'] = config.getfloat('secondary_gaia', 'pmra', fallback=0)
    companion_gaia['pmdec'] = config.getfloat('secondary_gaia', 'pmdec', fallback=0)
    companion_gaia['e_pmra'] = config.getfloat('secondary_gaia', 'epmra', fallback=1)
    companion_gaia['e_pmdec'] = config.getfloat('secondary_gaia', 'epmdec', fallback=1)
    companion_gaia['corr_pmra_pmdec'] = config.getfloat('secondary_gaia', 'corr_pmra_pmdec', fallback=0)
    return companion_gaia


def run():
    """
    Initialize and run the MCMC sampler.
    """
    global _loglkwargs
    start_time = time.time()

    args = parse_args()
    if not os.path.exists(args.output_dir):
        raise FileNotFoundError(f'No output_dir found at {args.output_dir}')
    # load_config raises a single ConfigError listing every problem it found,
    # rather than letting a typo surface later as a silent default.
    config = load_config(args.config_file, require_fit=True, strict=getattr(args, 'strict', False),
                         validate=not getattr(args, 'no_validate', False))

    # set the mcmc parameters
    nwalkers = config.getint('mcmc_settings', 'nwalkers', fallback=100)
    ntemps = config.getint('mcmc_settings', 'ntemps', fallback=10)
    nplanets = config.getint('mcmc_settings', 'nplanets')
    jit_per_inst = config.getboolean('mcmc_settings', 'jit_per_inst', fallback=False)
    nstep = config.getint('mcmc_settings', 'nstep')
    thin = config.getint('mcmc_settings', 'thin', fallback=50)
    nthreads = config.getint('mcmc_settings', 'nthreads', fallback=1)
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    start_file = config.get('data_paths', 'start_file', fallback='none')
    starname = config.get('plotting', 'target', fallback='HIP%d' % (HipID))
    #
    priors = get_priors(config)
    companion_gaia = get_gaia_catalog_companion(config)
    
    # Save configuration file in FITS header format
    header = make_header(args.config_file)

    data, H1f, H2f, Gf = initialize_data(config, companion_gaia)

    # set initial conditions
    if jit_per_inst:
        njit = data.nInst
    else:
        njit = 1
    par0 = set_initial_parameters(start_file, ntemps, nplanets, nwalkers, 
                                  njit=njit, minjit=priors['minjit'], 
                                  maxjit=priors['maxjit'])
    ndim = par0[0, 0, :].size

    # set arguments for emcee PTSampler and the log-likelyhood (lnprob)
    samplekwargs = {'thin': thin}
    loglkwargs = {'returninfo': False, 'use_epoch_astrometry': use_epoch_astrometry,
                  'data': data, 'nplanets': nplanets, 'H1f': H1f, 'H2f': H2f, 'Gf': Gf, 'priors': priors, 'njitters': njit}
    _loglkwargs = loglkwargs
    # run sampler without feeding it loglkwargs directly, since loglkwargs contains non-picklable C objects.

    # emcee dropped PTSampler in 3.0; use it only on a legacy emcee 2 install,
    # otherwise fall back to the bundled ptemcee.
    legacy_pt = getattr(emcee, 'PTSampler', None)
    if legacy_pt is not None:
        use_ptemcee = False
        sample0 = legacy_pt(ntemps, nwalkers, ndim, avoid_pickle_lnprob,
                            return_one, threads=nthreads)
        print('Using emcee.PTSampler (emcee %s).' % emcee.__version__)
    else:
        use_ptemcee = True
        sample0 = PTSampler(ntemps=ntemps, nwalkers=nwalkers, dim=ndim,
                            logl=avoid_pickle_lnprob, logp=return_one,
                            threads=nthreads)
        print('Using the bundled ptemcee sampler.')

    print("Running MCMC ... ")
    #sample0.run_mcmc(par0, nstep, **samplekwargs)
    #add a progress bar
    width = 30
    N = min(100, nstep//thin)
    n_taken = 0
    sys.stdout.write("[{0}]  {1}%".format(' ' * width, 0))
    for ipct in range(N):
        dn = (((nstep*(ipct + 1))//N - n_taken)//thin)*thin
        n_taken += dn
        if ipct == 0:
            sample0.run_mcmc(par0, dn, **samplekwargs)
        else:
            # Continue from last step
            sample0.run_mcmc(sample0.chain[..., -1, :], dn, **samplekwargs)
        n = int((width+1) * float(ipct + 1) / N)
        sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
        sys.stdout.write("%3d%%" % (int(100*(ipct + 1)/N)))
    sys.stdout.write("\n")
        
    print('Total Time: %.0f seconds' % (time.time() - start_time))
    print("Mean acceptance fraction (cold chain): {0:.6f}".format(np.mean(sample0.acceptance_fraction[0, :])))
    # save data
    if not use_ptemcee:
        shape = sample0.lnprobability[0].shape
    else:
        shape = sample0.logprobability[0].shape
    parfit = np.zeros((shape[0], shape[1], 9 + data.nInst))

    loglkwargs['returninfo'] = True
    loglkwargs['RVoffsets'] = True
    
    for i in range(shape[0]):
        for j in range(shape[1]):
            res, RVoffsets = lnprob(sample0.chain[0][i, j], **loglkwargs)
            parfit[i, j, :9] = [res.plx_best, res.pmra_best, res.pmdec_best,
                                res.chisq_sep, res.chisq_PA,
                                res.chisq_H, res.chisq_HG, res.chisq_G, res.chisq_relRV]
            
            if data.nInst > 0:
                parfit[i, j, 9:] = RVoffsets

    colnames = ['mpri']
    units = ['msun']
    for i in range(nplanets):
        colnames += [s + '%d' % (i) for s in ['msec', 'sau', 'esino', 'ecoso',
                                              'inc', 'asc', 'lam']]
        units += ['msun', 'au', '', '', 'radians', 'radians', 'radians']

    if njit == 1:
        colnames += ['jitter']
        units += ['m/s']
        sample0.chain[0][..., -1] = 10**(0.5*sample0.chain[0][..., -1])
    else:
        for i in range(njit):
            colnames += ['jitter%d' % (i)]
            units += ['m/s']
        sample0.chain[0][..., -njit:] = 10**(0.5*sample0.chain[0][..., -njit:])

    colnames += ['lnp']
    colnames += ['plx_ML', 'pmra_ML', 'pmdec_ML', 'chisq_sep', 
                 'chisq_PA', 'chisq_H', 'chisq_HG', 'chisq_G', 'chisq_relRV']
    units += ['', 'arcsec', 'arcsec/yr', 'arcsec/yr', '', '', '', '', '', '']
    colnames += ['RV_ZP_%d_ML' % (i) for i in range(data.nInst)]
    units += ['m/s' for i in range(data.nInst)]

    out = fits.HDUList(fits.PrimaryHDU(None, header))

    if not use_ptemcee:
        lnp = sample0.lnprobability[0]
    else:
        lnp = sample0.logprobability[0]

    out.append(pack_cols(sample0.chain[0], lnp, parfit, colnames, units))

    for i in range(1000):
        filename = os.path.join(args.output_dir, '%s_chain%03d.fits' % (starname, i))
        if not os.path.isfile(filename):
            print('Writing output to {0}'.format(filename))
            out.writeto(filename, overwrite=False)
            break

    return out


def main():
    """Console entry point ``fit_orbit``: run the fit and exit 0 on success."""
    run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
