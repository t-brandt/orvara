import numpy as np
from astropy.io import fits
import re
import pkg_resources
import warnings


def make_header(config_file):

    warnings.filterwarnings(action='ignore', 
                            category=fits.verify.VerifyWarning,
                            module=r'astropy.io.fits')

    header = fits.PrimaryHDU().header
    
    version = pkg_resources.get_distribution("orvara").version
    header.append(('version', version, 'Code version'))
    for line in open(config_file):
        line = line[:-1]
        try:
            if '=' in line and not line.startswith('#'):
                keys = re.split('[ ]*=[ ]*', line)
                # hierarch and continue cards are incompatible.  Hacky fix.     
                if len(keys[0]) > 8 and len(keys[0]) + len(keys[1]) + 13 > 80:
                    n_over = len(keys[0]) + len(keys[1]) + 13 - 80
                    header.append((keys[0], keys[1][:-n_over]), end=True)
                else:
                    header.append((keys[0], keys[1]), end=True)
            elif not line.startswith('#'):
                header.append(('comment', line[:80]), end=True)
        except:
            continue
    return header


def pack_cols(chains, lnp, parfit, names, units):

    if chains.shape[-1] + 1 + parfit.shape[-1] != len(names):
        raise ValueError('Number of names for fits table does not match number of data fields')
    
    n = 0
    fmt = '%dD' % (chains.shape[1])
    cols = []
    
    for i in range(chains.shape[-1]):
        cols += [fits.Column(name=names[n], format=fmt, array=chains[..., i], unit=units[n])]
        n += 1

    cols += [fits.Column(name=names[n], format=fmt, array=lnp, unit=units[n])]
    n += 1

    for i in range(parfit.shape[-1]):
        cols += [fits.Column(name=names[n], format=fmt, array=parfit[..., i], unit=units[n])]
        n += 1

    return fits.BinTableHDU.from_columns(cols)


def pull_chain_params(columns, step=0, lnp_name='lnp'):

    params = []
        
    for i in range(len(columns)):
        if columns[i].name == lnp_name:
            break
        params += [float(columns[i].array[step])]

    return params


def burnin_chain(columns, burnin=0, reshape=True):

    if columns[0].array.shape[1] <= burnin:
        raise ValueError("Cannot use a burnin length longer than the chain.")
        
    newcols = []
    for i in range(len(columns)):
        if reshape:
            newcols += [fits.Column(name=columns[i].name,
                                    format=columns[i].format[-1],
                                    array=columns[i].array[:, burnin:].flatten())]
        else:
            arr = columns[i].array[:, burnin:]
            fmt = '%d' % (arr.shape[1]) + columns[i].format[-1]
            newcols += [fits.Column(name=columns[i].name,
                                    format=fmt,
                                    array=arr)]
            
    return fits.BinTableHDU.from_columns(newcols).data












