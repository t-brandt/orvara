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
from ptemcee import Sampler as PTSampler
from configparser import ConfigParser
from astropy.io import fits
from htof.main import Astrometry
from astropy.time import Time
import sys
import re
import random
from orvara import orbit
from orvara.config import parse_args
from orvara.format_fits import make_header, pack_cols
import pkg_resources

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
        except:
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
              [2, 1e-4, 1e3],        # msec (Solar masses)
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


def initialize_data(config, companion_gaia):
    # load in items from the ConfigParser object
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    HGCAFile = config.get('data_paths', 'HGCAFile')
    if not os.path.exists(HGCAFile):
        raise FileNotFoundError(f'No HGCA file found at {HGCAFile}')
    HGCAVersion = config.get('data_paths', 'HGCAVersion', fallback='GaiaDR2')
    if HGCAVersion.lower() != 'GaiaDR2'.lower() and HGCAVersion.lower() != 'GaiaeDR3'.lower():
        raise ValueError('HGCAVersion in the config file is not GaiaDR2 and is not GaiaeDR3')
    RVFile = config.get('data_paths', 'RVFile')
    AstrometryFile = config.get('data_paths', 'AstrometryFile')
    GaiaDataDir = config.get('data_paths', 'GaiaDataDir', fallback=None)
    Hip2DataDir = config.get('data_paths', 'Hip2DataDir', fallback=None)
    Hip1DataDir = config.get('data_paths', 'Hip1DataDir', fallback=None)
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)

    data = orbit.Data(HipID, HGCAFile, RVFile, AstrometryFile, companion_gaia=companion_gaia)
    if use_epoch_astrometry and data.use_abs_ast == 1:
        to_jd = lambda x: Time(x, format='decimalyear').jd + 0.5
        Gaia_fitter = Astrometry(HGCAVersion, '%06d' % (HipID), GaiaDataDir,
                                 central_epoch_ra=to_jd(data.epRA_G),
                                 central_epoch_dec=to_jd(data.epDec_G),
                                 format='jd', normed=False)
        Hip2_fitter = Astrometry('Hip2', '%06d' % (HipID), Hip2DataDir,
                                 central_epoch_ra=to_jd(data.epRA_H),
                                 central_epoch_dec=to_jd(data.epDec_H),
                                 format='jd', normed=False)
        Hip1_fitter = Astrometry('Hip1', '%06d' % (HipID), Hip1DataDir,
                                 central_epoch_ra=to_jd(data.epRA_H),
                                 central_epoch_dec=to_jd(data.epDec_H),
                                 format='jd', normed=False)
        # instantiate C versions of the astrometric fitter which are much faster than HTOF's Astrometry
        #print(Gaia_fitter.fitter.astrometric_solution_vector_components['ra'])
        hip1_fast_fitter = orbit.AstrometricFitter(Hip1_fitter)
        hip2_fast_fitter = orbit.AstrometricFitter(Hip2_fitter)
        gaia_fast_fitter = orbit.AstrometricFitter(Gaia_fitter)

        data = orbit.Data(HipID, HGCAFile, RVFile, AstrometryFile, 
                          use_epoch_astrometry,
                          epochs_Hip1=Hip1_fitter.data.julian_day_epoch(),
                          epochs_Hip2=Hip2_fitter.data.julian_day_epoch(),
                          epochs_Gaia=Gaia_fitter.data.julian_day_epoch(),
                          companion_gaia=companion_gaia, verbose=False)
    elif data.use_abs_ast == 1:
        hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter = None, None, None
    else:
        hip1_fast_fitter, hip2_fast_fitter, gaia_fast_fitter = None, None, None
        try:
            data.plx = 1e-3*config.getfloat('priors_settings', 'parallax')
            data.plx_err = 1e-3*config.getfloat('priors_settings', 'parallax_error')
        except:
            print("Cannot load absolute astrometry.")
            print("Please supply a prior for parallax and parallax_error")
            print("in the priors_settings area of the configuration file.")

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
        orbit.calc_offsets(data, params, model, i)

        if priors is not None:
            lnp = lnp - 0.5*(params.msec - priors.get(f'm_secondary{i}', 1))**2/priors.get(f'm_secondary{i}_sig', np.inf)**2

    if use_epoch_astrometry and data.use_abs_ast:
        orbit.calc_PMs_epoch_astrometry(data, model, H1f, H2f, Gf)
    elif data.use_abs_ast:
        orbit.calc_PMs_no_epoch_astrometry(data, model)

    if returninfo:
        return orbit.calcL(data, params, model, chisq_resids=True, RVoffsets=RVoffsets)

    if priors is not None:
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
    config = ConfigParser()
    if not os.path.exists(args.config_file):
        raise FileNotFoundError(f'No config file found at {args.config_file}')
    if not os.path.exists(args.output_dir):
        raise FileNotFoundError(f'No output_dir found at {args.output_dir}')
    config.read(args.config_file)

    # set the mcmc parameters
    nwalkers = config.getint('mcmc_settings', 'nwalkers')
    ntemps = config.getint('mcmc_settings', 'ntemps')
    nplanets = config.getint('mcmc_settings', 'nplanets')
    jit_per_inst = config.getboolean('mcmc_settings', 'jit_per_inst', fallback=False)
    nstep = config.getint('mcmc_settings', 'nstep')
    thin = config.getint('mcmc_settings', 'thin', fallback=50)
    nthreads = config.getint('mcmc_settings', 'nthreads')
    use_epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry', fallback=False)
    HipID = config.getint('data_paths', 'HipID', fallback=0)
    start_file = config.get('data_paths', 'start_file', fallback='none')
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

    try:
        use_ptemcee = False
        sample0 = emcee.PTSampler(ntemps, nwalkers, ndim, avoid_pickle_lnprob, return_one, threads=nthreads)
        print('Using emcee.PTSampler.')
    except:
        use_ptemcee = True
        sample0 = PTSampler(ntemps=ntemps, nwalkers=nwalkers, dim=ndim,
                            logl=avoid_pickle_lnprob, logp=return_one,
                            threads=nthreads)
        print('Using ptemcee.')

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
    parfit = np.zeros((shape[0], shape[1], 8 + data.nInst))

    loglkwargs['returninfo'] = True
    loglkwargs['RVoffsets'] = True
    
    for i in range(shape[0]):
        for j in range(shape[1]):
            res, RVoffsets = lnprob(sample0.chain[0][i, j], **loglkwargs)
            parfit[i, j, :8] = [res.plx_best, res.pmra_best, res.pmdec_best,
                                res.chisq_sep, res.chisq_PA,
                                res.chisq_H, res.chisq_HG, res.chisq_G]
            
            if data.nInst > 0:
                parfit[i, j, 8:] = RVoffsets

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
                 'chisq_PA', 'chisq_H', 'chisq_HG', 'chisq_G']
    units += ['', 'arcsec', 'arcsec/yr', 'arcsec/yr', '', '', '', '', '']
    colnames += ['RV_ZP_%d_ML' % (i) for i in range(data.nInst)]
    units += ['m/s' for i in range(data.nInst)]

    out = fits.HDUList(fits.PrimaryHDU(None, header))

    if not use_ptemcee:
        lnp = sample0.lnprobability[0]
    else:
        lnp = sample0.logprobability[0]

    out.append(pack_cols(sample0.chain[0], lnp, parfit, colnames, units))

    for i in range(1000):
        filename = os.path.join(args.output_dir, 'HIP%d_chain%03d.fits' % (HipID, i))
        if not os.path.isfile(filename):
            print('Writing output to {0}'.format(filename))
            out.writeto(filename, overwrite=False)
            break


if __name__ == "__main__":
    run()
