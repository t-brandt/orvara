#!/usr/bin/env python
"""
Plotting the proper motions
"""

from __future__ import print_function
import orbit
import numpy as np
from random import randrange
from scipy.interpolate import interp1d
import scipy.optimize as op
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rcParams
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator
from PyAstronomy import pyasl
import emcee, corner

######### load in data ##############

HIP = 95319
RVfile = './tests/data/Gl758_RV.dat'
relAstfile = './tests/data/Gl758_relAST.txt'
use_epoch_astrometry = False
Gaia_intermediate_data = '/Users/Sophie/Downloads/Gaia/'
Hip1_intermediate_data = '/Users/Sophie/hipparcosOriginalIntermediateData/'
Hip2_intermediate_data = '/Users/Sophie/Downloads/Hip2/IntermediateData/resrec/'
output_dir = './output_plotting/'

# read in mcmc chains
path = '/Users/Sophie/git/orbit3d/orbit3d/'
file = 'HIP95319_chain001.fits'
source = file.split('_')[0]
tt, lnp, extras = [fits.open(path+file)[i].data for i in range(3)]
nsteps = 50*tt.shape[1]
beststep = np.where(lnp==lnp.max())


######### define some useful functions ##############
def JD_to_calendar(JD):
    a = int(JD + 0.5)
    if a < 2299161:
        c = a + 1524
    else:
        b = int((a - 1867216.25)/36524.25)
        c = a + b - int(b/4) + 1525
    d = int((c - 122.1)/365.25)
    e = int(365.25*d)
    f = int((c - e)/30.6001)

    D = c - e - int(30.6001*f) + ((JD + 0.5) - int(JD + 0.5))
    M = f - 1 - 12*int(f/14)
    Y = d - 4715 - int((7 + M)/10)
    year = Y + M/12 + D/365.25
    return year

def calendar_to_JD(year, M=1, D=1):
    if M <= 2:
        y = year - 1
        m = M + 12
    else:
        y = year
        m = M
    if year <= 1583: # more precisely should be less than or equal to 10/4/1582
        B = -2
    else:
        B = int(y/400) - int(y/100)
    UT = 0
    JD = int(365.25*y) + int(30.6001*(m+1)) + B + 1720996.5  + D + UT/24
    return JD
    

def chi_sqr(offset, epoch, ra_dec_obs, ra_dec_obs_err):
    chi_sqr = 0
    for i in range(len(epoch)):
        chi_sqr += (f_mu(epoch[i]) - ra_dec_obs[i] - offset)**2 / ra_dec_obs_err[i]**2
    return chi_sqr






#initiate the data
data = orbit.Data(HIP, RVfile, relAstfile)

# calculate the most likely RV & astrometric orbits
params = orbit.Params(tt[beststep][0])
period = params.per
ep = np.linspace(calendar_to_JD(1983), calendar_to_JD(2023) , 1000)
data.custom_epochs(ep)
model = orbit.Model(data)

orbit.calc_EA_RPP(data, params, model)
orbit.calc_offsets(data, params, model, 0)
orbit.calc_RV(data, params, model)

mu_RA, mu_Dec =  model.return_proper_motions(params)


def plot_proper_motions(number_orbits , flag):
    rcParams['xtick.major.width']=1
    rcParams['xtick.major.size']=4
    rcParams['xtick.minor.width']=0.5
    rcParams['xtick.minor.size']=2
    rcParams['xtick.direction'] = "in"
    rcParams['ytick.direction'] = "in"

    fig = plt.figure(figsize=(5, 6))
    ax1 = fig.add_axes((0.15, 0.3, 0.7, 0.5))
    ax2 = fig.add_axes((0.15, 0.1, 0.7, 0.15))

    ratio = ((1. + params.mpri/params.msec))*1000.

    if flag == 'mu_RA':
        x_Ra = [1991.0 , 2015.5]
        y_Ra = [0.6, -0.38]
        y_Ra_err = [0.49, 0.05]

        def chi_sqr(offset, epoch, ra_dec_obs, ra_dec_obs_err):
                chi_sqr = 0
                for i in range(len(epoch)):
                    chi_sqr += (f_mu(epoch[i]) - ra_dec_obs[i] - offset)**2 / ra_dec_obs_err[i]**2
                return chi_sqr
            
        f_ra = interp1d(ep_jd, mu_RA)

        result = op.minimize(chi_sqr, 0, args=(x_Ra, y_Ra, y_Ra_err))
        #offset0 = result['x']
        offset = - 1.56
        ax1.plot(ep_jd, mu_RA*ratio +offset ,color= 'k', linewidth = 2,zorder = 1000000)

        for i in range(40):

            # get parameters from one single step of the mcmc chain
            walker = randrange(tt.shape[0])
            step = randrange(1000, tt.shape[1])
            params = orbit.Params(tt[walker, step])

            # calculate and assign variables
            data.custom_epochs(ep)
            model = orbit.Model(data)

            orbit.calc_EA_RPP(data, params, model)
            orbit.calc_offsets(data, params, model, 0)
            orbit.calc_RV(data, params, model)
            mu_ra, mu_dec =  model.return_proper_motions(params)
            
            for i in range(len(ep)):
                ep_jd[i] = JD_to_calendar(ep[i])


            import pandas as pd
            d = {'ep_calendar': ep_jd, 'mu_Ra': mu_ra}
            df = pd.DataFrame(data=d)
            mid_pt = df[int(len(ep)/2):int(len(ep)/2)+1]['mu_Ra']
            
            di = {'ep_calendar': ep_jd, 'mu_ra': mu_ra*ratio}
            dfi = pd.DataFrame(data=di)
            mid_pti = dfi[int(len(ep)/2):int(len(ep)/2)+1]['mu_ra']
            
            mu_offset = mid_pti[int(len(ep)/2)] - mid_pt[int(len(ep)/2)]
            #print(mid_pt[int(len(ep)/2)],mid_pti[int(len(ep)/2)],mu_offset)
            cmap = plt.cm.cubehelix
            mu_y = mu_ra*ratio - mu_offset
            ax1.plot(ep_jd, mu_y, c= cmap((params.msec*1989/1.898 - 34.)/(42-34.)))
            
            ax2.plot(ep_jd, np.zeros(len(ep)), 'k--', dashes=(5, 5))
            for j in range(len(ep)):
                    mu_y[j] -= (f_ra(ep_jd[j])*ratio + offset)
            ax2.plot(ep_jd, mu_y, c =cmap((params.msec*1989/1.898 - 34.)/(42-34.)) , alpha=0.3)
            ax2.scatter(x_Ra, y_Ra - (f_ra(x_Ra)*ratio + offset),zorder = 10000)

        ax1.errorbar(x_Ra, y_Ra,yerr= y_Ra_err ,fmt='o', ecolor='k', capthick=3,capsize=4,zorder=1000)
        ax1.set_xlim(np.min(ep_jd),np.max(ep_jd))
        ax2.set_xlim(np.min(ep_jd),np.max(ep_jd))
        ax2.set_ylim(-2,2)
        ax1.set_ylabel(r'$\mathrm{\Delta \mu_{RA} \, (mas \, yr^{-1})}$')
        ax2.set_xlabel("Epoch")
        ax2.set_ylabel("O-C")
        ax1.minorticks_on()
        ax2.minorticks_on()
        ax1.text(2010,1.0,'Gl 758B', fontsize =16)
    
    
    if flag == 'mu_Dec':
        x_Dec = [1991.25 , 2015.75]
        y_Dec = [1.41, -0.90]
        y_Dec_err = [0.40, 0.]
            
        f_mu = interp1d(ep_jd, mu_Dec)

        result = op.minimize(chi_sqr, 1.53, args=(x_Dec, y_Dec, y_Dec_err))
        offset = result['x']

        ax1.plot(ep_jd, mu_Dec*ratio + offset,color= 'k', linewidth = 2,zorder = 100)

        for i in range(number_orbits:

            walker = randrange(tt.shape[0])
            step = randrange(1000, tt.shape[1])
            params = orbit.Params(tt[walker, step])

            data.custom_epochs(ep)
            model = orbit.Model(data)

            orbit.calc_EA_RPP(data, params, model)
            orbit.calc_offsets(data, params, model, 0)
            orbit.calc_RV(data, params, model)
            mu_ra, mu_dec =  model.return_proper_motions(params)
            
            for i in range(len(ep)):
                ep_jd[i] = JD_to_calendar(ep[i])

            import pandas as pd
            d = {'ep_calendar': ep_jd, 'mu_Dec': mu_Dec}
            df = pd.DataFrame(data=d)
            mid_pt = df[int(len(ep)/2):int(len(ep)/2)+1]['mu_Dec']
            
            di = {'ep_calendar': ep_jd, 'mu_dec': mu_dec*ratio}
            dfi = pd.DataFrame(data=di)
            mid_pti = dfi[int(len(ep)/2):int(len(ep)/2)+1]['mu_dec']
            
            mu_offset = mid_pti[int(len(ep)/2)] - mid_pt[int(len(ep)/2)]
            #print(mid_pt[int(len(ep)/2)],mid_pti[int(len(ep)/2)],mu_offset)
            cmap = plt.cm.cubehelix
            mu_y = mu_dec*ratio - mu_offset
            ax1.plot(ep_jd, mu_y, c= cmap((params.msec*1989/1.898 - 34.)/(42-34.)))
            
            ax2.plot(ep_jd, np.zeros(len(ep)), 'k--', dashes=(5, 5))
            for j in range(len(ep)):
                    mu_y[j] -= (f_mu(ep_jd[j])*ratio + offset)
            ax2.plot(ep_jd, mu_y, c =cmap((params.msec*1989/1.898 - 34.)/(42-34.)) , alpha=0.3)
            ax2.scatter(x_Dec, y_Dec - (f_mu(x_Dec)*ratio + offset),zorder = 10000)

        ax1.errorbar(x_Dec, y_Dec,yerr= y_Dec_err ,fmt='o', ecolor='k', capthick=3,capsize=4,zorder=1000)
        ax1.set_ylabel(r'$\mathrm{\Delta \mu_{Dec} \, (mas \, yr^{-1})}$')
        ax1.set_xlim(np.min(ep_jd),np.max(ep_jd))
        ax1.minorticks_on()
        ax2.minorticks_on()
        ax1.text(2010,1.5,'Gl 758B', fontsize =16)
        
        ax2.set_xlim(np.min(ep_jd),np.max(ep_jd))
        ax2.set_ylim(-2,2)
        ax2.set_xlabel("Epoch")
        ax2.set_ylabel("O-C")

