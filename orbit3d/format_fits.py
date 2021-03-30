import numpy as np
from astropy.io import fits
import re
import pkg_resources


def make_header(config_file):

    header = fits.PrimaryHDU().header
    
    version = pkg_resources.get_distribution("orbit3d").version
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


def pack_cols(chains, lnp, parfit, names):

    if chains.shape[-1] + 1 + parfit.shape[-1] != len(names):
        raise ValueError('Number of names for fits table does not match number of data fields')
    
    n = 0
    fmt = '%dD' % (chains.shape[1])
    cols = []
    
    for i in range(chains.shape[-1]):
        cols += [fits.Column(name=names[n], format=fmt, array=chains[..., i])]
        n += 1

    cols += [fits.Column(name=names[n], format=fmt, array=lnp)]
    n += 1

    for i in range(parfit.shape[-1]):
        cols += [fits.Column(name=names[n], format=fmt, array=parfit[..., i])]
        n += 1

    return fits.BinTableHDU.from_columns(cols)


def pull_chain_params(columns, walker=0, step=0, lnp_name='lnp'):

    params = []
    if step == 'best':
        walker, step = np.where(columns[lnp_name] == np.amax(columns[lnp_name]))
        
    for i in range(len(columns)):
        if columns[i].name == lnp_name:
            break
        params += [float(columns[i][walker, step])]

    return params














