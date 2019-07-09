#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import os
import time
import emcee
from astropy.io import fits
import orbit
from htof.main import Astrometry

# Number of MCMC steps per walker
nstep = 2000

# Number of threads to use for parallelization.  
nthreads = 2

######################################################################
HIP = 95319
RVfile = 'Gl758_RV.dat'
relAstfile = 'Gl758_relAST.txt'
use_epoch_astrometry = False
Gaia_intermediate_data = '/home/tbrandt/data/GaiaDR2IntermediateData/'
Hip1_intermediate_data = '/home/tbrandt/data/hipparcosOriginalIntermediateData/'
Hip2_intermediate_data = '/home/tbrandt/data/Hip2/IntermediateData/resrec/'
output_dir = '/home/tbrandt/chains/'
#startfile = 'HD4747_PTinititer1.fits2'
startfile = None
nwalkers = 100
ntemps = 5
nplanets = 1
######################################################################

######################################################################
# Garbage initial guess
######################################################################

if startfile is None:
    mpri = 1
    jit = 0.5
    sau = 10
    esino = 0.5
    ecoso = 0.5
    inc = 1
    asc = 1
    lam = 1
    msec = 0.1
    
    par0 = np.ones((ntemps, 100, 2 + 7*nplanets))
    init = [jit, mpri]
    for i in range(nplanets):
        init += [msec, sau, esino, ecoso, inc, asc, lam]
    par0 *= np.asarray(init)
    par0 *= 2**(np.random.rand(np.prod(par0.shape)).reshape(par0.shape) - 0.5)
    
else:

    #################################################################
    # read in the starting positions for the walkers. The next four
    # lines remove parallax and RV zero point from the optimization,
    # change semimajor axis from arcseconds to AU, and bring the
    # number of temperatures used for parallel tempering down to
    # ntemps.
    #################################################################

    par0 = fits.open(startfile)[0].data
    par0[:, :, 8] = par0[:, :, 9]
    par0[:, :, 9] = par0[:, :, 10]
    par0[:, :, 0] /= par0[:, :, 9]
    par0 = par0[:ntemps, :, :-2]

ntemps = par0[:, 0, 0].size
nwalkers = par0[0, :, 0].size
ndim = par0[0, 0, :].size

######################################################################
# Load in data
######################################################################

data = orbit.Data(HIP, RVfile, relAstfile)

if use_epoch_astrometry:
    Gaia_fitter = Astrometry('GaiaDR2', '%06d' % (HIP), Gaia_intermediate_data,
                             central_epoch_ra=data.epRA_G,
                             central_epoch_dec=data.epDec_G,
                             central_epoch_fmt='frac_year')
    Hip2_fitter = Astrometry('Hip2', '%06d' % (HIP), Hip2_intermediate_data,
                             central_epoch_ra=data.epRA_H,
                             central_epoch_dec=data.epDec_H,
                             central_epoch_fmt='frac_year')
    Hip1_fitter = Astrometry('Hip1', '%06d' % (HIP), Hip1_intermediate_data,
                             central_epoch_ra=data.epRA_H,
                             central_epoch_dec=data.epDec_H,
                             central_epoch_fmt='frac_year')
    
    H1f = orbit.AstrometricFitter(Hip1_fitter)
    H2f = orbit.AstrometricFitter(Hip2_fitter)
    Gf = orbit.AstrometricFitter(Gaia_fitter)

    data = orbit.Data(HIP, RVfile, relAstfile, use_epoch_astrometry,
                      epochs_Hip1=Hip1_fitter.data.julian_day_epoch(),
                      epochs_Hip2=Hip2_fitter.data.julian_day_epoch(),
                      epochs_Gaia=Gaia_fitter.data.julian_day_epoch())

######################################################################
# define likelihood function for joint parameters
######################################################################
 
def lnprob(theta, returninfo=False):
    
    model = orbit.Model(data)

    for i in range(nplanets):
        params = orbit.Params(theta, i, nplanets)
        
        if not np.isfinite(orbit.lnprior(params)):
            model.free()
            return -np.inf

        orbit.calc_EA_RPP(data, params, model)
        orbit.calc_RV(data, params, model)
        orbit.calc_offsets(data, params, model, i)
    
    if use_epoch_astrometry:
        orbit.calc_PMs_epoch_astrometry(data, model, H1f, H2f, Gf)
    else:
        orbit.calc_PMs_no_epoch_astrometry(data, model)

    if returninfo:
        return orbit.calcL(data, params, model, chisq_resids=True)
        
    return orbit.lnprior(params) + orbit.calcL(data, params, model)

def return_one(theta):
    return 1.

######################################################################
# Initialize and run sampler
######################################################################

from emcee import PTSampler
kwargs = {'thin': 50 }
start_time = time.time()

sample0 = emcee.PTSampler(ntemps, nwalkers, ndim, lnprob, return_one, threads=nthreads)
sample0.run_mcmc(par0, nstep, **kwargs)

print('Total Time: %.2f' % (time.time() - start_time))
print("Mean acceptance fraction (cold chain): {0:.6f}".format(np.mean(sample0.acceptance_fraction[0,:])))

shape = sample0.lnprobability[0].shape
parfit = np.zeros((shape[0], shape[1], 8))
for i in range(shape[0]):
    for j in range(shape[1]):
        res = lnprob(sample0.chain[0][i, j], returninfo=True)
        parfit[i, j] = [res.plx_best, res.pmra_best, res.pmdec_best,
                        res.chisq_sep, res.chisq_PA,
                        res.chisq_H, res.chisq_HG, res.chisq_G]

out = fits.HDUList(fits.PrimaryHDU(sample0.chain[0].astype(np.float32)))
out.append(fits.PrimaryHDU(sample0.lnprobability[0].astype(np.float32)))
out.append(fits.PrimaryHDU(parfit.astype(np.float32)))
for i in range(1000):
    filename = os.path.join(output_dir, 'HIP%d_chain%03d.fits' % (HIP, i))
    if not os.path.isfile(filename):
        out.writeto(filename, overwrite=False)
        exit()
