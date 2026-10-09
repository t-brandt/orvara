import numpy as np
from astropy.io import fits
import warnings

from orvara import __version__


#: A FITS card is 80 characters; HIERARCH and CONTINUE cards cannot be combined,
#: so a long key plus a long value has to be truncated rather than continued.
_CARD_LENGTH = 80
_CARD_OVERHEAD = 13


def make_header(config_file):
    """
    Record the configuration file in the primary header of the chain file.

    Long values are truncated to fit an 80-character card.  With the unified
    config this matters more than it used to: ``corr_vec`` alone is several
    hundred characters.
    """
    warnings.filterwarnings(action='ignore',
                            category=fits.verify.VerifyWarning,
                            module=r'astropy.io.fits')

    header = fits.PrimaryHDU().header
    header.append(('version', __version__, 'Code version'))

    with open(config_file) as handle:
        for line in handle:
            line = line.rstrip('\n')
            stripped = line.strip()
            if not stripped or stripped.startswith(('#', ';')):
                continue
            try:
                if '=' in line:
                    key, _, value = line.partition('=')
                    key, value = key.strip(), value.strip()
                    budget = _CARD_LENGTH - len(key) - _CARD_OVERHEAD
                    if budget < 1:
                        # the key alone will not fit; keep it as a comment
                        header.append(('comment', line[:_CARD_LENGTH]), end=True)
                    else:
                        header.append((key, value[:budget]), end=True)
                else:
                    header.append(('comment', line[:_CARD_LENGTH]), end=True)
            except Exception:
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
        value = np.asarray(columns[i].array[step])
        # `step` may be a plain index into a flattened chain (0-d result) or,
        # as in plot_chains, a fancy index into the 2-D chain that returns a
        # length-1 array; NumPy >= 2 refuses float() on the latter.
        params += [float(value.reshape(-1)[0])]

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












