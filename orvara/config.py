"""
Configuration handling for orvara.

Everything a run needs lives in **one** ``.ini`` file: the orbit-fit settings
(``[data_paths]``, ``[mcmc_settings]``, ``[priors_settings]``, ``[plotting]``,
...) *and*, when the target is a Gaia DR3 non-single-star (NSS) acceleration
source, the NSS solution itself (``[gaia_accel]``).  ``fit_orbit``,
``plot_orbit`` and ``orvara-gaia-accel`` all read the same file.

The module provides

* :func:`load_config`  -- read + validate a config file, returning a
  ``ConfigParser``.  Problems are reported all at once with the section, key
  and reason, instead of surfacing later as a ``KeyError`` or a silent default.
* :func:`nss_settings` -- normalised view of the Gaia NSS block.  Accepts both
  the unified ``[gaia_accel]`` section and the older
  ``[paths]/[star]/[advanced]/[gaia_solution]`` prep files, so existing configs
  keep working.
* :func:`write_template` / ``orvara-config`` -- emit an annotated starter config.

Backwards compatibility: :func:`parse_args` and :func:`parse_args_plotting`
keep the same command line as before.
"""

import argparse
import difflib
import os
import re
import warnings
from configparser import ConfigParser

import numpy as np

__all__ = [
    'ConfigError', 'load_config', 'validate_config', 'nss_settings',
    'write_template', 'template_cli', 'parse_args', 'parse_args_plotting',
    'TEMPLATES',
]


class ConfigError(ValueError):
    """Raised when a configuration file is missing or malformed."""


class Key:
    """One recognised configuration key."""

    __slots__ = ('name', 'kind', 'required', 'choices', 'doc')

    def __init__(self, name, kind='str', required=False, choices=None, doc=''):
        self.name = name
        self.kind = kind
        self.required = required
        self.choices = choices
        self.doc = doc


_BOOLS = {'1': True, 'yes': True, 'true': True, 'on': True,
          '0': False, 'no': False, 'false': False, 'off': False}


def _coerce(kind, raw):
    """Convert a raw string to ``kind``.  Raises ValueError with a clear message."""
    text = raw.strip()
    if kind in ('str', 'path'):
        return text
    if text == '':
        raise ValueError('value is empty')
    if kind == 'int':
        return int(text)
    if kind == 'float':
        return float(text)
    if kind == 'bool':
        try:
            return _BOOLS[text.lower()]
        except KeyError:
            raise ValueError(f"{text!r} is not a boolean (use True/False, yes/no, 1/0)")
    if kind == 'floats':
        tokens = [tok for tok in re.split(r'[,\s]+', text) if tok]
        return np.array([float(tok) for tok in tokens], dtype=float)
    raise AssertionError(f'unknown kind {kind!r}')          # pragma: no cover


#: The orbit-fit half of the config.  Section name -> list of Key.
FIT_SCHEMA = {
    'data_paths': [
        Key('HipID', 'int', doc='Hipparcos ID of the target (or the synthetic ID '
                                'assigned to a Gaia-only NSS source).'),
        Key('HGCAFile', 'path', required=True,
            doc='HGCA catalogue, or the one-star file made by orvara-gaia-accel.'),
        Key('RVFile', 'path'),
        Key('relRVFile', 'path'),
        Key('AstrometryFile', 'path'),
        Key('relast_error_scale', 'float',
            doc='Multiplies every separation and position-angle uncertainty in '
                'AstrometryFile (default 1). Set to sqrt(chi2_rel / 2 N_rel) of a '
                'previous fit, cumulatively, to bring the relative-astrometry '
                'residuals to unit variance.'),
        Key('GaiaDataDir', 'path', doc='Gaia epoch astrometry (GOST) directory.'),
        Key('Hip1DataDir', 'path'),
        Key('Hip2DataDir', 'path'),
        Key('start_file', 'path'),
    ],
    'mcmc_settings': [
        Key('ntemps', 'int'),
        Key('nwalkers', 'int'),
        Key('nplanets', 'int', required=True),
        Key('nstep', 'int', required=True),
        Key('nthreads', 'int'),
        Key('thin', 'int'),
        Key('jit_per_inst', 'bool'),
        Key('use_epoch_astrometry', 'bool',
            doc='Fit the Hipparcos/Gaia intermediate astrometry with htof.'),
    ],
    'priors_settings': [
        Key('mpri', 'float'),
        Key('mpri_sig', 'float'),
        Key('minjitter', 'float'),
        Key('maxjitter', 'float'),
        Key('parallax', 'float',
            doc='mas; required when no absolute astrometry is used.'),
        Key('parallax_error', 'float'),
    ],
    'secondary_gaia': [
        Key('companion_ID', 'int'),
        Key('pmra', 'float'), Key('pmdec', 'float'),
        Key('epmra', 'float'), Key('epmdec', 'float'),
        Key('corr_pmra_pmdec', 'float'),
    ],
    'plotting': [
        Key('McmcDataFile', 'path'), Key('target', 'str'),
        Key('start_epoch', 'float'), Key('end_epoch', 'float'),
        Key('num_orbits', 'int'), Key('num_steps', 'int'),
        Key('burnin', 'int'), Key('iplanet', 'int'),
        Key('check_convergence', 'bool'),
        Key('Relative_separation_plot', 'bool'),
        Key('Position_angle_plot', 'bool'),
        Key('Proper_motion_plot', 'bool'),
        Key('Proper_motion_separate_plots', 'bool'),
        Key('Astrometric_prediction_plot', 'bool'),
        Key('Astrometry_orbits_plot', 'bool'),
        Key('RV_orbits_plot', 'bool'), Key('RV_plot', 'bool'),
        Key('Corner_plot', 'bool'),
        Key('RV_Instrument', 'str'),
        Key('position_predict', 'float'),
        Key('predicted_years', 'str'),
        Key('reference', 'str'),
        Key('show_title', 'bool'), Key('add_text', 'bool'),
        Key('text_name', 'str'), Key('x_text', 'float'), Key('y_text', 'float'),
        Key('marker_color', 'str'), Key('colormap', 'str'),
        Key('use_colorbar', 'bool'),
        Key('set_limit', 'bool'), Key('xlim', 'str'), Key('ylim', 'str'),
    ],
    'save_results': [
        Key('save_params', 'bool'),
        Key('err_margin', 'str'),
    ],
}

#: Keys accepted in ``priors_settings`` that carry an index, e.g. m_secondary0.
_INDEXED_PRIOR = re.compile(r'^m_secondary\d(_sig)?$')   # get_priors reads 0-9 only


#: The Gaia NSS half of the config.
NSS_SCHEMA = {
    'gaia_accel': [
        Key('enabled', 'bool'),
        Key('mode', 'str', choices=('gaia_only', 'hip_and_gaia'),
            doc="'gaia_only': synthetic Hipparcos entry with huge errors, so only "
                "the Gaia NSS solution constrains the fit.  'hip_and_gaia': copy "
                "the real HGCA row for hip_id and attach the NSS terms to it."),
        Key('hgca_template', 'path',
            doc='Full HGCA catalogue used as the column template / Hip source.'),
        Key('output_fits', 'path', required=True,
            doc='One-star file to write; use the same path as data_paths.HGCAFile.'),
        Key('source_id', 'int', doc='Gaia DR3 source_id.'),
        Key('hip_id', 'int', required=True,
            doc='Hipparcos ID.  For gaia_only pick an unused integer (e.g. 200000+).'),
        Key('npar', 'int', required=True, choices=(7, 9),
            doc='7 = Acceleration7, 9 = Acceleration9 (acceleration + jerk).'),
        Key('error_inflation_factor', 'float',
            doc='Multiplies the NSS acceleration/jerk errors. 2 is the suggested default.'),
        Key('ra', 'float'), Key('ra_err', 'float'),
        Key('dec', 'float'), Key('dec_err', 'float'),
        Key('parallax', 'float'), Key('parallax_err', 'float'),
        Key('pmra', 'float'), Key('pmra_err', 'float'),
        Key('pmdec', 'float'), Key('pmdec_err', 'float'),
        Key('rv', 'float'), Key('rv_err', 'float'), Key('rv_source', 'str'),
        Key('accra', 'float', required=True), Key('accra_err', 'float', required=True),
        Key('accdec', 'float', required=True), Key('accdec_err', 'float', required=True),
        Key('deriv_accra', 'float'), Key('deriv_accra_err', 'float'),
        Key('deriv_accdec', 'float'), Key('deriv_accdec_err', 'float'),
        Key('corr_vec', 'floats', required=True,
            doc='Gaia NSS corr_vec: 21 numbers for npar=7, 36 for npar=9.'),
    ],
    'gaia_accel_rv': [
        Key('file', 'path',
            doc='CSV with epoch(JD), rv(m/s), rv_err(m/s)[, instrument].'),
    ],
    'gaia_accel_relast': [
        Key('file', 'path',
            doc='CSV with epoch(JD), sep(mas), pa(deg), sep_err, pa_err.'),
    ],
}

SCHEMA = {**FIT_SCHEMA, **NSS_SCHEMA}

# Legacy prep-file layout -> unified layout.  (old_section, old_key) -> new_key
_LEGACY_MAP = {
    ('paths', 'hgca_template'): 'hgca_template',
    ('paths', 'output_fits'): 'output_fits',
    ('star', 'source_id'): 'source_id',
    ('star', 'hip_id'): 'hip_id',
    ('advanced', 'error_inflation_factor'): 'error_inflation_factor',
    ('advanced', 'mode'): 'mode',
}
_LEGACY_SECTIONS = {'paths', 'star', 'advanced', 'gaia_solution', 'rv', 'relast'}


def _read(path):
    if not os.path.exists(path):
        raise ConfigError(f'No configuration file found at {path!r}.')
    # Keys are case-insensitive, as they always were in orvara: ConfigParser
    # lower-cases them on read and on lookup, so `hipid` and `HipID` are the
    # same key.  (An earlier revision preserved case, which silently turned
    # `hipid = 3850` into HipID = 0.)
    config = ConfigParser(inline_comment_prefixes=('#', ';'))
    if not config.read(path):
        raise ConfigError(f'Could not parse the configuration file {path!r}.')
    return config


def _known_keys(section):
    return {key.name.lower(): key for key in SCHEMA.get(section, [])}


def _suggest(name, candidates):
    match = difflib.get_close_matches(name.lower(), list(candidates), n=1, cutoff=0.75)
    return f" (did you mean {match[0]!r}?)" if match else ''


def validate_config(config, path='<config>', require_fit=False,
                    require_nss=False, strict=False):
    """
    Check ``config`` against the schema.

    Type errors, bad choices, missing required keys and unusable combinations
    are collected and raised together as a single :class:`ConfigError`.
    Unrecognised sections and keys are reported as warnings -- they are almost
    always typos -- or as errors when ``strict`` is True.

    Whether a key is *required* depends on what the file is about to be used
    for, so a prep-only file or a half-filled template is not an error by
    itself.  ``require_fit`` demands the keys ``fit_orbit`` needs;
    ``require_nss`` demands the ``[gaia_accel]`` block ``orvara-gaia-accel``
    needs.  Type and consistency checks always run.

    :returns: the same ``config`` object, for chaining.
    """
    errors, unknown = [], []

    for section in config.sections():
        low = section.lower()
        if low in _LEGACY_SECTIONS:
            continue                                  # handled by nss_settings
        if low not in SCHEMA:
            unknown.append(f'unknown section [{section}]{_suggest(section, SCHEMA)}')
            continue
        known = _known_keys(low)
        for name, raw in config.items(section):
            key = known.get(name.lower())
            if key is None:
                if low == 'priors_settings' and _INDEXED_PRIOR.match(name):
                    key = Key(name, 'float')
                else:
                    unknown.append(f'[{section}] unknown key {name!r}'
                                   f'{_suggest(name, known)}')
                    continue
            if raw.strip() == '':
                continue                              # blank == "not supplied"
            try:
                value = _coerce(key.kind, raw)
            except ValueError as exc:
                errors.append(f'[{section}] {name} = {raw.strip()!r}: {exc} '
                              f'(expected {key.kind})')
                continue
            if key.choices is not None and value not in key.choices:
                errors.append(f'[{section}] {name} = {value!r} is not one of '
                              f'{list(key.choices)}')

    if require_fit:
        for section, keys in FIT_SCHEMA.items():
            for key in keys:
                if key.required and \
                        config.get(section, key.name, fallback='').strip() == '':
                    errors.append(f'[{section}] {key.name} is required but '
                                  f'missing or empty.')

    try:
        epoch_astrometry = config.getboolean('mcmc_settings', 'use_epoch_astrometry',
                                             fallback=False)
    except ValueError:
        epoch_astrometry = False          # already reported as a type error above
    if epoch_astrometry:
        for name in ('GaiaDataDir', 'Hip1DataDir', 'Hip2DataDir'):
            value = config.get('data_paths', name, fallback='').strip()
            if not value:
                errors.append(f'[mcmc_settings] use_epoch_astrometry = True requires '
                              f'[data_paths] {name} to be set.')
            elif not os.path.isdir(value):
                errors.append(f'[data_paths] {name} = {value!r} is not a directory.')

    hgca = config.get('data_paths', 'HGCAFile', fallback='').strip()
    output_fits = config.get('gaia_accel', 'output_fits', fallback='').strip()
    if hgca and output_fits and os.path.abspath(hgca) != os.path.abspath(output_fits):
        warnings.warn(
            f'[data_paths] HGCAFile ({hgca}) and [gaia_accel] output_fits '
            f'({output_fits}) point at different files. orvara-gaia-accel writes '
            f'output_fits and fit_orbit reads HGCAFile, so unless that is '
            f'deliberate the fit will not see the NSS solution you just built.',
            RuntimeWarning)

    if require_nss and not (config.has_section('gaia_accel')
                            or config.has_section('gaia_solution')):
        errors.append('No [gaia_accel] section: this command needs a Gaia NSS block. '
                      'Run `orvara-config --template acc9` for a starting point.')

    if unknown:
        if strict:
            errors.extend(unknown)
        else:
            warnings.warn(f'In {path}:\n  ' + '\n  '.join(unknown), RuntimeWarning)

    if errors:
        raise ConfigError(f'{len(errors)} problem(s) in {path}:\n  ' + '\n  '.join(errors))
    return config


def load_config(path, require_fit=False, require_nss=False, strict=False,
                validate=True):
    """Read ``path`` and validate it.  Returns a ``ConfigParser``."""
    config = _read(path)
    if validate:
        validate_config(config, path=path, require_fit=require_fit,
                        require_nss=require_nss, strict=strict)
    return config


#: Parameter order of the Gaia DR3 NSS acceleration solutions.
NSS_PARAM_ORDER = {
    7: ['ra', 'dec', 'parallax', 'pmra', 'pmdec', 'accra', 'accdec'],
    9: ['ra', 'dec', 'parallax', 'pmra', 'pmdec', 'accra', 'accdec',
        'deriv_accra', 'deriv_accdec'],
}


def _gather_nss_raw(config):
    """Collect the NSS keys from either the unified or the legacy layout."""
    raw = {}
    for (section, key), new_key in _LEGACY_MAP.items():          # legacy first
        for candidate in config.sections():
            if candidate.lower() == section:
                for name, value in config.items(candidate):
                    if name.lower() == key:
                        raw[new_key] = value
    for candidate in config.sections():
        if candidate.lower() == 'gaia_solution':
            for name, value in config.items(candidate):
                raw[name.lower()] = value
    if config.has_section('gaia_accel'):                          # unified wins
        for name, value in config.items('gaia_accel'):
            raw[name.lower()] = value

    files = {}
    for section, alias in (('rv', 'rv'), ('gaia_accel_rv', 'rv'),
                           ('relast', 'relast'), ('gaia_accel_relast', 'relast')):
        if config.has_section(section):
            value = config.get(section, 'file', fallback='').strip()
            if value:
                files[alias] = value
    return raw, files


def nss_settings(config, path='<config>'):
    """
    Normalise and validate the Gaia NSS acceleration block.

    Accepts the unified ``[gaia_accel]`` section as well as the legacy
    ``[paths]/[star]/[advanced]/[gaia_solution]`` prep files.

    :returns: dict of fully-typed, validated values.
    :raises ConfigError: listing every problem found, not just the first.
    """
    raw, files = _gather_nss_raw(config)
    if not raw:
        raise ConfigError(
            f'{path} contains no Gaia NSS block. Add a [gaia_accel] section '
            f'(see `orvara-config --template acc9`).')

    errors = []

    def get(name, kind, required=False, default=None):
        text = raw.get(name, '')
        if text is None or str(text).strip() == '':
            if required:
                errors.append(f'[gaia_accel] {name} is required but missing or empty.')
            return default
        try:
            return _coerce(kind, str(text))
        except ValueError as exc:
            errors.append(f'[gaia_accel] {name} = {str(text).strip()!r}: {exc}')
            return default

    npar = get('npar', 'int', required=True)
    if npar not in (7, 9):
        if npar is not None:
            errors.append(f'[gaia_accel] npar = {npar!r}; only 7 (Acceleration7) and '
                          f'9 (Acceleration9) NSS solutions can be prepared.')
        npar = 7

    mode = (raw.get('mode') or 'gaia_only').strip().lower()
    if mode not in ('gaia_only', 'hip_and_gaia'):
        errors.append(f"[gaia_accel] mode = {mode!r} is not one of "
                      f"['gaia_only', 'hip_and_gaia']")
        mode = 'gaia_only'

    out = {
        'mode': mode,
        'enabled': get('enabled', 'bool', default=True),
        'npar': npar,
        'hgca_template': get('hgca_template', 'path', default='HGCA_vEDR3.fits'),
        'output_fits': get('output_fits', 'path', required=True, default=''),
        'source_id': get('source_id', 'int', default=0),
        'hip_id': get('hip_id', 'int', required=True, default=0),
        'error_inflation_factor': get('error_inflation_factor', 'float', default=1.0),
        'rv': get('rv', 'float', default=np.nan),
        'rv_err': get('rv_err', 'float', default=np.nan),
        'rv_source': (raw.get('rv_source') or '').strip(),
        'rv_file': files.get('rv'),
        'relast_file': files.get('relast'),
    }

    for name in ('ra', 'dec', 'parallax', 'pmra', 'pmdec', 'accra', 'accdec'):
        out[name] = get(name, 'float', required=True)
        out[name + '_err'] = get(name + '_err', 'float', required=True)

    for name in ('deriv_accra', 'deriv_accdec'):
        if npar == 9:
            out[name] = get(name, 'float', required=True)
            out[name + '_err'] = get(name + '_err', 'float', required=True)
        else:
            out[name] = 0.0
            out[name + '_err'] = 0.0

    corr = get('corr_vec', 'floats', required=True)
    expected = npar * (npar - 1) // 2
    if corr is not None:
        if corr.size < expected:
            errors.append(
                f'[gaia_accel] corr_vec has {corr.size} entries but npar={npar} '
                f'needs at least {expected}. Copy the full corr_vec from the '
                f'gaiadr3.nss_acceleration_astro table.')
        else:
            corr = corr[:expected]   # the archive column is zero-padded
    out['corr_vec'] = corr

    if out['error_inflation_factor'] is not None and out['error_inflation_factor'] <= 0:
        errors.append('[gaia_accel] error_inflation_factor must be positive.')
    if mode == 'hip_and_gaia' and not out['hgca_template']:
        errors.append('[gaia_accel] mode = hip_and_gaia needs hgca_template to point '
                      'at the full HGCA catalogue that holds the Hipparcos data.')
    for name in ('rv_file', 'relast_file'):
        value = out[name]
        if value and not os.path.exists(value):
            errors.append(f'[gaia_accel] {name} = {value!r} does not exist.')

    if errors:
        raise ConfigError(f'{len(errors)} problem(s) in the Gaia NSS block of '
                          f'{path}:\n  ' + '\n  '.join(errors))
    return out


def _template_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'gaia_accel', 'templates')


TEMPLATES = ('fit', 'acc7', 'acc9', 'hip_and_gaia')


def write_template(name, dest=None):
    """
    Write the annotated starter config ``name`` (one of :data:`TEMPLATES`).

    :param dest: output path; defaults to ``<name>_config.ini`` in the cwd.
    :returns: the path written.
    """
    if name not in TEMPLATES:
        raise ConfigError(f'Unknown template {name!r}; choose from {list(TEMPLATES)}.')
    source = os.path.join(_template_dir(), f'{name}.ini')
    if not os.path.exists(source):                     # pragma: no cover
        raise ConfigError(f'Template file {source!r} is missing from the install.')
    if dest is None:
        dest = os.path.join(os.getcwd(), f'{name}_config.ini')
    with open(source) as handle:
        text = handle.read()
    with open(dest, 'w') as handle:
        handle.write(text)
    return dest


def template_cli(argv=None):
    """Console entry point ``orvara-config``."""
    parser = argparse.ArgumentParser(
        prog='orvara-config',
        description='Write an annotated orvara configuration file, or check an '
                    'existing one.')
    parser.add_argument('--template', choices=TEMPLATES,
                        help='Which starter config to write. "fit" is an orbit fit '
                             'against the full HGCA; "acc7"/"acc9" add a Gaia NSS '
                             'acceleration block; "hip_and_gaia" combines the real '
                             'Hipparcos data with the NSS solution.')
    parser.add_argument('-o', '--output', help='Where to write it.')
    parser.add_argument('--check', metavar='CONFIG',
                        help='Validate an existing configuration file and exit.')
    parser.add_argument('--strict', action='store_true',
                        help='With --check, treat unknown sections/keys as errors.')
    args = parser.parse_args(argv)

    if args.check:
        # Structural and type checks first: these are always errors.
        config = load_config(args.check, strict=args.strict)
        print(f'{args.check}: structure and types OK.')

        # Then report what the file is ready for.  A prep-only file or a
        # half-filled template is incomplete, not wrong, so say so rather than
        # failing.
        try:
            validate_config(config, path=args.check, require_fit=True,
                            strict=args.strict)
            print('  ready for fit_orbit / plot_orbit: yes')
        except ConfigError as exc:
            print('  ready for fit_orbit / plot_orbit: not yet --'
                  f'{str(exc).split(":", 1)[1]}')
        try:
            nss_settings(config, path=args.check)
            print('  ready for orvara-gaia-accel: yes')
        except ConfigError as exc:
            if 'contains no Gaia NSS block' in str(exc):
                print('  ready for orvara-gaia-accel: no [gaia_accel] section '
                      '(only needed for a Gaia NSS source)')
            else:
                print('  ready for orvara-gaia-accel: not yet --'
                      f'{str(exc).split(":", 1)[1]}')
        return 0

    if not args.template:
        parser.error('give either --template or --check')
    print(f'Wrote {args.template} template to '
          f'{write_template(args.template, args.output)}')
    return 0


def parse_args():
    parser = argparse.ArgumentParser(
        description='Fit a Keplerian orbit or orbits using orvara.')
    parser.add_argument('config_file', metavar='configfile',
                        help='Path to the planet-specific configuration file.')
    parser.add_argument('--output-dir',
                        help='Directory within which to save the MCMC results, '
                             'default current directory.')
    parser.add_argument('--no-validate', action='store_true',
                        help='Skip configuration validation (not recommended).')
    parser.add_argument('--strict', action='store_true',
                        help='Treat unknown configuration sections/keys as errors.')
    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = os.getcwd()
    return args


def parse_args_plotting():
    parser = argparse.ArgumentParser(
        description='Plot astrometry or RV orbits and other relevant plots.')
    parser.add_argument('config_file', metavar='configfile',
                        help='Path to the planet-specific configuration file.')
    parser.add_argument('--output-dir',
                        help='Directory within which to save the plots, '
                             'default current directory.')
    parser.add_argument('--no-validate', action='store_true',
                        help='Skip configuration validation (not recommended).')
    parser.add_argument('--strict', action='store_true',
                        help='Treat unknown configuration sections/keys as errors.')
    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = os.getcwd()
    return args
