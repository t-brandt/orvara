#!/usr/bin/env python
"""
Build a one-star HGCA-format FITS file from a Gaia DR3 NSS acceleration solution.

Command line::

    orvara-gaia-accel my_target.ini          # unified or legacy prep .ini

then point ``[data_paths] HGCAFile`` at the file it writes (the unified config
does that for you, since both come from ``[gaia_accel] output_fits``).

The output table carries the 35 standard HGCA columns plus the NSS extras that
``orvara.orbit.Data`` looks for:

    gaia_npar, accra_gaia, accdec_gaia, accra_gaia_error, accdec_gaia_error,
    accra_accdec_gaia
    (npar = 9 only) jerkra_gaia, jerkdec_gaia, jerkra_gaia_error,
    jerkdec_gaia_error, jerkra_jerkdec_gaia

Units follow the HGCA convention: proper motions and their errors in mas/yr,
accelerations in mas/yr^2, jerks in mas/yr^3, epochs in decimal years.
"""

import argparse
import os

import numpy as np
from astropy.io import fits

from orvara.config import load_config, nss_settings, ConfigError

#: Gaia DR3 reference epoch (decimal year).  The NSS acceleration model is
#: written in t = epoch - 2016.0, and its position and proper-motion parameters
#: are *means* over the observation span, so this is the epoch at which
#: orvara must centre its sample window -- not the "uncorrelated position
#: epoch" used by the HGCA for 5-parameter sources.
GAIA_REF_EPOCH = 2016.0
#: Hipparcos reference epoch (decimal year).
HIP_REF_EPOCH = 1991.25
#: Half-span of the observation window assumed by the DR3 NSS acceleration
#: model (Gaia DR3 documentation, Sect. 7.2.4, Eq. 7.5: "DeltaT = 517.5 days").
#: The model's position and proper motion are averages over [-DeltaT, +DeltaT].
GAIA_NSS_HALF_SPAN_DAYS = 517.5
#: The corresponding full baseline in years, used by the likelihood and the
#: diagnostic plots to reproduce those averages.
GAIA_BASELINE_YR = 2.0 * GAIA_NSS_HALF_SPAN_DAYS / 365.25
#: Proper-motion error assigned to the synthetic Hipparcos entry in gaia_only
#: mode.  Large enough that those rows contribute nothing to the likelihood.
BIG_PM_ERROR = 1e7

#: The 35 columns of the released HGCA catalogues, in order, with FITS formats.
HGCA_BASE_COLUMNS = [
    ('hip_id', 'J'), ('gaia_source_id', 'K'),
    ('gaia_ra', 'D'), ('gaia_dec', 'D'),
    ('radial_velocity', 'E'), ('radial_velocity_error', 'E'),
    ('radial_velocity_source', '8A'),
    ('parallax_gaia', 'E'), ('parallax_gaia_error', 'E'),
    ('pmra_gaia', 'E'), ('pmdec_gaia', 'E'),
    ('pmra_gaia_error', 'E'), ('pmdec_gaia_error', 'E'), ('pmra_pmdec_gaia', 'E'),
    ('pmra_hg', 'E'), ('pmdec_hg', 'E'),
    ('pmra_hg_error', 'E'), ('pmdec_hg_error', 'E'), ('pmra_pmdec_hg', 'E'),
    ('pmra_hip', 'E'), ('pmdec_hip', 'E'),
    ('pmra_hip_error', 'E'), ('pmdec_hip_error', 'E'), ('pmra_pmdec_hip', 'E'),
    ('epoch_ra_gaia', 'D'), ('epoch_dec_gaia', 'D'),
    ('epoch_ra_hip', 'D'), ('epoch_dec_hip', 'D'),
    ('crosscal_pmra_hip', 'E'), ('crosscal_pmdec_hip', 'E'),
    ('crosscal_pmra_hg', 'E'), ('crosscal_pmdec_hg', 'E'),
    ('nonlinear_dpmra', 'E'), ('nonlinear_dpmdec', 'E'),
    ('chisq', 'E'),
]

#: Extra columns that carry the NSS solution into orvara.
NSS_ACC_COLUMNS = [
    ('gaia_npar', 'J'),
    ('accra_gaia', 'E'), ('accdec_gaia', 'E'),
    ('accra_gaia_error', 'E'), ('accdec_gaia_error', 'E'),
    ('accra_accdec_gaia', 'E'),
]
NSS_JERK_COLUMNS = [
    ('jerkra_gaia', 'E'), ('jerkdec_gaia', 'E'),
    ('jerkra_gaia_error', 'E'), ('jerkdec_gaia_error', 'E'),
    ('jerkra_jerkdec_gaia', 'E'),
]


def corr_vec_to_matrix(vec, npar):
    """
    Rebuild the full ``npar`` x ``npar`` correlation matrix from Gaia's
    ``corr_vec``.

    Gaia stores the off-diagonal correlations once, in the order

        (0,1), (0,2), (1,2), (0,3), (1,3), (2,3), ...

    i.e. the upper triangle read column by column -- which is the same
    sequence as the lower triangle read row by row.

    :param vec: 1-D array of ``npar*(npar-1)/2`` correlations.
    :param npar: 7 or 9.
    :returns: symmetric (npar, npar) correlation matrix with unit diagonal.
    """
    vec = np.asarray(vec, dtype=float).ravel()
    needed = npar * (npar - 1) // 2
    if vec.size != needed:
        raise ValueError(f'corr_vec has {vec.size} entries; npar={npar} needs '
                         f'exactly {needed}.')
    matrix = np.eye(npar, dtype=float)
    k = 0
    for j in range(1, npar):
        for i in range(j):
            matrix[i, j] = matrix[j, i] = vec[k]
            k += 1
    return matrix


def _nss_param_index(npar):
    names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec', 'accra', 'accdec']
    if npar == 9:
        names += ['jerkra', 'jerkdec']
    return {name: i for i, name in enumerate(names)}


def _hgca_row_for(template_path, hip_id):
    """
    Return the HGCA row for ``hip_id`` as a plain dict, or None if absent.

    :raises FileNotFoundError: if the template catalogue does not exist.
    """
    if not os.path.exists(template_path):
        raise FileNotFoundError(
            f'HGCA template {template_path!r} not found. Point '
            f'[gaia_accel] hgca_template at HGCA_vEDR3.fits.')
    with fits.open(template_path, memmap=False) as hdul:
        data = hdul[1].data
        names = list(data.names)
        rows = data[data['hip_id'] == int(hip_id)]
        if len(rows) == 0:
            return None
        row = rows[0]
        out = {}
        for name in names:
            value = row[name]
            if isinstance(value, (bytes, np.bytes_)):
                value = value.decode('utf-8', 'replace').strip()
            elif isinstance(value, str):
                value = value.strip()
            else:
                value = value.item() if hasattr(value, 'item') else value
            out[name] = value
        return out


def build_row(settings):
    """
    Turn a validated ``[gaia_accel]`` block into one HGCA-format table row.

    :param settings: dict from :func:`orvara.config.nss_settings`.
    :returns: dict mapping column name -> Python scalar.
    """
    npar = settings['npar']
    corr = corr_vec_to_matrix(settings['corr_vec'], npar)
    index = _nss_param_index(npar)
    inflate = settings['error_inflation_factor']

    pmra, pmdec = settings['pmra'], settings['pmdec']

    # For a 5-parameter (linear) source the HGCA places the Gaia epoch where
    # position and proper motion decorrelate, 2016.0 - corr(ra,pmra) sig_ra /
    # sig_pmra.  That is *wrong* for an NSS acceleration solution: its proper
    # motion is the mean over the span centred on 2016.0, and its acceleration
    # and jerk are the derivatives at 2016.0 (Gaia DR3 documentation Sect.
    # 7.2.4).  Shifting the window by dt biases the model proper motion by
    # acc * dt -- tens of sigma for a typical accelerator -- which a Gaia-only
    # fit absorbs into the free barycentric proper motion but a combined
    # Hipparcos + Gaia fit cannot.  So the window is centred on 2016.0.
    epoch_ra_gaia = GAIA_REF_EPOCH
    epoch_dec_gaia = GAIA_REF_EPOCH

    row = {
        'hip_id': int(settings['hip_id']),
        'gaia_source_id': int(settings['source_id']),
        'gaia_ra': settings['ra'],
        'gaia_dec': settings['dec'],
        'radial_velocity': settings['rv'],
        'radial_velocity_error': settings['rv_err'],
        'radial_velocity_source': settings['rv_source'][:8],
        'parallax_gaia': settings['parallax'],
        'parallax_gaia_error': settings['parallax_err'],
        'pmra_gaia': pmra,
        'pmdec_gaia': pmdec,
        'pmra_gaia_error': settings['pmra_err'],
        'pmdec_gaia_error': settings['pmdec_err'],
        'pmra_pmdec_gaia': corr[index['pmra'], index['pmdec']],
        'epoch_ra_gaia': float(epoch_ra_gaia),
        'epoch_dec_gaia': float(epoch_dec_gaia),
        'gaia_npar': npar,
        'accra_gaia': settings['accra'],
        'accdec_gaia': settings['accdec'],
        'accra_gaia_error': settings['accra_err'] * inflate,
        'accdec_gaia_error': settings['accdec_err'] * inflate,
        'accra_accdec_gaia': corr[index['accra'], index['accdec']],
    }

    if npar == 9:
        row.update({
            'jerkra_gaia': settings['deriv_accra'],
            'jerkdec_gaia': settings['deriv_accdec'],
            'jerkra_gaia_error': settings['deriv_accra_err'] * inflate,
            'jerkdec_gaia_error': settings['deriv_accdec_err'] * inflate,
            'jerkra_jerkdec_gaia': corr[index['jerkra'], index['jerkdec']],
        })

    if settings['mode'] == 'hip_and_gaia':
        hip = _hgca_row_for(settings['hgca_template'], settings['hip_id'])
        if hip is None:
            raise ConfigError(
                f"[gaia_accel] mode = hip_and_gaia, but HIP {settings['hip_id']} is "
                f"not in {settings['hgca_template']!r}. Use mode = gaia_only for a "
                f"star with no Hipparcos data, or check hip_id.")
        # Keep the real Hipparcos and Hip-Gaia proper motions and epochs; the
        # Gaia columns and the NSS terms above take precedence over them.
        for name in ('pmra_hip', 'pmdec_hip', 'pmra_hip_error', 'pmdec_hip_error',
                     'pmra_pmdec_hip', 'pmra_hg', 'pmdec_hg', 'pmra_hg_error',
                     'pmdec_hg_error', 'pmra_pmdec_hg', 'epoch_ra_hip',
                     'epoch_dec_hip', 'crosscal_pmra_hip', 'crosscal_pmdec_hip',
                     'crosscal_pmra_hg', 'crosscal_pmdec_hg', 'nonlinear_dpmra',
                     'nonlinear_dpmdec', 'chisq'):
            if name in hip:
                row[name] = hip[name]
    else:
        # gaia_only: no Hipparcos measurement exists.  Fill the Hip and Hip-Gaia
        # rows with the Gaia proper motion and give them enormous errors, so
        # they add nothing to chi^2 and only the NSS solution drives the fit.
        row.update({
            'pmra_hip': pmra, 'pmdec_hip': pmdec,
            'pmra_hip_error': BIG_PM_ERROR, 'pmdec_hip_error': BIG_PM_ERROR,
            'pmra_pmdec_hip': 0.0,
            'pmra_hg': pmra, 'pmdec_hg': pmdec,
            'pmra_hg_error': BIG_PM_ERROR, 'pmdec_hg_error': BIG_PM_ERROR,
            'pmra_pmdec_hg': 0.0,
            'epoch_ra_hip': HIP_REF_EPOCH, 'epoch_dec_hip': HIP_REF_EPOCH,
            'crosscal_pmra_hip': pmra, 'crosscal_pmdec_hip': pmdec,
            'crosscal_pmra_hg': pmra, 'crosscal_pmdec_hg': pmdec,
            'nonlinear_dpmra': pmra * 1.0e-4, 'nonlinear_dpmdec': pmdec * 1.0e-4,
            'chisq': 0.0,
        })
    return row


def _column_spec(npar):
    spec = list(HGCA_BASE_COLUMNS) + list(NSS_ACC_COLUMNS)
    if npar == 9:
        spec += list(NSS_JERK_COLUMNS)
    return spec


def _blank(fmt):
    if fmt.endswith('A'):
        return ''
    if fmt in ('J', 'K', 'I', 'B', 'L'):
        return 0
    return np.nan


def _hgca_hdu(row, npar):
    columns = []
    for name, fmt in _column_spec(npar):
        value = row.get(name, _blank(fmt))
        if fmt.endswith('A'):
            array = np.array([str(value)])
        elif fmt in ('J', 'I', 'B', 'L'):
            array = np.array([int(value)], dtype=np.int32)
        elif fmt == 'K':
            array = np.array([int(value)], dtype=np.int64)
        else:
            array = np.array([float(value)], dtype=float)
        columns.append(fits.Column(name=name, format=fmt, array=array))
    hdu = fits.BinTableHDU.from_columns(columns)
    hdu.name = 'HGCA'
    return hdu


def _read_table(path, required, kind):
    """Read a whitespace- or comma-separated table and check its columns."""
    from astropy.table import Table
    fmt = 'ascii.csv' if path.lower().endswith(('.csv', '.tsv')) else 'ascii'
    table = Table.read(path, format=fmt)
    missing = [name for name in required if name not in table.colnames]
    if missing:
        raise ConfigError(f'{kind} file {path!r} is missing column(s) {missing}. '
                          f'Expected columns: {list(required)}.')
    return table


def _rv_hdu(table):
    epoch = np.asarray(table['epoch'], dtype=float)
    instrument = (np.asarray(table['instrument'], dtype=np.int32)
                  if 'instrument' in table.colnames
                  else np.zeros(len(epoch), dtype=np.int32))
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='epoch', format='D', array=epoch, unit='JD'),
        fits.Column(name='rv', format='D',
                    array=np.asarray(table['rv'], dtype=float), unit='m/s'),
        fits.Column(name='rv_err', format='D',
                    array=np.asarray(table['rv_err'], dtype=float), unit='m/s'),
        fits.Column(name='instrument', format='J', array=instrument),
    ])
    hdu.name = 'RV'
    return hdu


def _relast_hdu(table):
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='epoch', format='D',
                    array=np.asarray(table['epoch'], dtype=float), unit='JD'),
        fits.Column(name='sep', format='D',
                    array=np.asarray(table['sep'], dtype=float), unit='mas'),
        fits.Column(name='sep_err', format='D',
                    array=np.asarray(table['sep_err'], dtype=float), unit='mas'),
        fits.Column(name='pa', format='D',
                    array=np.asarray(table['pa'], dtype=float), unit='deg'),
        fits.Column(name='pa_err', format='D',
                    array=np.asarray(table['pa_err'], dtype=float), unit='deg'),
    ])
    hdu.name = 'RELAST'
    return hdu


def write_nss_fits(settings, overwrite=True, write_data_files=True):
    """
    Write the one-star HGCA-format FITS file described by ``settings``.

    :param settings: dict from :func:`orvara.config.nss_settings`.
    :param overwrite: replace an existing output file.
    :param write_data_files: also convert any ``rv_file`` / ``relast_file`` CSVs
        into the whitespace-delimited text files ``fit_orbit`` reads, written
        next to the output FITS.
    :returns: dict with ``fits`` and, when written, ``rv`` and ``relast`` paths.
    """
    npar = settings['npar']
    row = build_row(settings)

    header = fits.Header()
    header['TARGET'] = (str(settings['hip_id']), 'hip_id used by orvara')
    header['SOURCEID'] = (int(settings['source_id']), 'Gaia DR3 source_id')
    header['NSSNPAR'] = (npar, 'NSS solution: 7=Accel7, 9=Accel9')
    header['NSSMODE'] = (settings['mode'], 'gaia_only or hip_and_gaia')
    header['ERRINFL'] = (settings['error_inflation_factor'],
                         'Factor applied to NSS accel/jerk errors')
    header['GAIA_T'] = (GAIA_REF_EPOCH, 'Gaia DR3 reference epoch [yr]')
    header['GAIA_DT'] = (GAIA_BASELINE_YR, 'Gaia astrometric baseline [yr]')

    hdulist = fits.HDUList([fits.PrimaryHDU(header=header), _hgca_hdu(row, npar)])

    written = {}
    outpath = settings['output_fits']
    outdir = os.path.dirname(os.path.abspath(outpath))
    os.makedirs(outdir, exist_ok=True)

    if settings.get('rv_file'):
        table = _read_table(settings['rv_file'], ('epoch', 'rv', 'rv_err'), 'RV')
        hdulist.append(_rv_hdu(table))
        if write_data_files:
            path = os.path.join(outdir, f"{os.path.splitext(os.path.basename(outpath))[0]}_RV.dat")
            columns = ['epoch', 'rv', 'rv_err']
            if 'instrument' in table.colnames:
                columns.append('instrument')
            np.savetxt(path, np.column_stack([np.asarray(table[c], dtype=float)
                                              for c in columns]),
                       fmt='%.8f', header=' '.join(columns))
            written['rv'] = path

    if settings.get('relast_file'):
        table = _read_table(settings['relast_file'],
                            ('epoch', 'sep', 'sep_err', 'pa', 'pa_err'),
                            'relative astrometry')
        hdulist.append(_relast_hdu(table))
        if write_data_files:
            path = os.path.join(outdir, f"{os.path.splitext(os.path.basename(outpath))[0]}_relAST.txt")
            columns = ['epoch', 'sep', 'sep_err', 'pa', 'pa_err']
            np.savetxt(path, np.column_stack([np.asarray(table[c], dtype=float)
                                              for c in columns]),
                       fmt='%.8f', header=' '.join(columns))
            written['relast'] = path

    hdulist.writeto(outpath, overwrite=overwrite)
    written['fits'] = outpath

    print(f'Wrote {outpath}')
    print(f'  hip_id            : {row["hip_id"]}  (use this as [data_paths] HipID)')
    print(f'  Gaia source_id    : {row["gaia_source_id"]}')
    print(f'  NSS solution      : {npar}-parameter, mode = {settings["mode"]}')
    print(f'  Gaia epoch        : {row["epoch_ra_gaia"]:.4f}  (NSS reference epoch; '
          f'window +/- {GAIA_NSS_HALF_SPAN_DAYS} d)')
    print(f'  acceleration      : accra = {row["accra_gaia"]:.6f} +/- '
          f'{row["accra_gaia_error"]:.6f}, accdec = {row["accdec_gaia"]:.6f} +/- '
          f'{row["accdec_gaia_error"]:.6f}  [mas/yr^2]')
    if npar == 9:
        print(f'  jerk              : jerkra = {row["jerkra_gaia"]:.6f} +/- '
              f'{row["jerkra_gaia_error"]:.6f}, jerkdec = {row["jerkdec_gaia"]:.6f} '
              f'+/- {row["jerkdec_gaia_error"]:.6f}  [mas/yr^3]')
    for kind in ('rv', 'relast'):
        if kind in written:
            print(f'  {kind:<18}: {written[kind]}')
    return written


def build_from_config(config_path, overwrite=True, write_data_files=True):
    """
    Read ``config_path``, validate its Gaia NSS block and write the FITS file.

    Accepts the unified config as well as the legacy
    ``gaia_accel_prep_AccN.ini`` layout.

    :returns: dict of written paths (see :func:`write_nss_fits`).
    """
    config = load_config(config_path, require_nss=True)
    settings = nss_settings(config, path=config_path)
    if not settings['enabled']:
        raise ConfigError(f'[gaia_accel] enabled = False in {config_path!r}; '
                          f'nothing to build.')
    # Relative paths are taken from the working directory, exactly as fit_orbit
    # and plot_orbit take [data_paths] HGCAFile -- so "keep output_fits and
    # HGCAFile identical" means what it says.  As a convenience only, an
    # *input* that is missing from the cwd but present next to the config file
    # is picked up from there.
    base = os.path.dirname(os.path.abspath(config_path))
    for key in ('hgca_template', 'rv_file', 'relast_file'):
        value = settings.get(key)
        if value and not os.path.isabs(value) and not os.path.exists(value):
            candidate = os.path.join(base, value)
            if os.path.exists(candidate):
                settings[key] = candidate
    return write_nss_fits(settings, overwrite=overwrite,
                          write_data_files=write_data_files)


def main(argv=None):
    """Console entry point ``orvara-gaia-accel``."""
    parser = argparse.ArgumentParser(
        prog='orvara-gaia-accel',
        description='Build a one-star HGCA-format FITS file from a Gaia DR3 NSS '
                    'astrometric acceleration solution, for use with fit_orbit.')
    parser.add_argument('config_file', metavar='configfile', nargs='?',
                        help='Unified orvara config (with a [gaia_accel] section) '
                             'or a legacy gaia_accel_prep_AccN.ini file.')
    parser.add_argument('--config', dest='config_flag', metavar='CONFIG',
                        help=argparse.SUPPRESS)          # accepted for compatibility
    parser.add_argument('--no-overwrite', action='store_true',
                        help='Fail instead of replacing an existing output file.')
    parser.add_argument('--no-data-files', action='store_true',
                        help='Do not convert the RV / relative-astrometry CSVs into '
                             'the text files fit_orbit reads.')
    args = parser.parse_args(argv)
    config_file = args.config_file or args.config_flag
    if config_file is None:
        parser.error('a configuration file is required')

    build_from_config(config_file,
                      overwrite=not args.no_overwrite,
                      write_data_files=not args.no_data_files)
    return 0


if __name__ == '__main__':          # pragma: no cover
    raise SystemExit(main())
