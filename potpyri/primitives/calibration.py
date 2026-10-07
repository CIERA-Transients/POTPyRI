"""Master bias, dark, and flat creation from pipeline file tables.

Orchestrates grouping by CalType and calling instrument-specific
create_bias/create_dark/create_flat. Authors: Charlie Kilpatrick.
"""
from potpyri._version import __version__

import os
import time
import numpy as np


def _emit(log, msg, level='error'):
    """Write msg to log at level, or print when no logger is provided."""
    if log:
        getattr(log, level)(msg)
    else:
        print(msg)


def _master_exists(tel, paths, cal_name):
    """True if a night or packaged-caldb master exists at cal_name."""
    found = os.path.exists(cal_name) or os.path.exists(cal_name + '.fz')
    if not found and hasattr(tel, '_find_cal_file'):
        path, _hdu = tel._find_cal_file(paths, os.path.basename(cal_name))
        found = path is not None
    return found


def _warn_missing_science_masters(file_table, tel, paths, nmin_images, log,
        cal_kind, get_name, name_cols):
    """Error for each science setup that has no master calibration file.

    name_cols are file_table columns passed to get_name after paths
    (e.g. ('Amp', 'Binning') or ('Filter', 'Amp', 'Binning')).
    """
    kwds = tel.filetype_keywords
    if 'SCIENCE' not in kwds:
        return
    sci_match = tel.match_type_keywords(kwds['SCIENCE'], file_table)
    if not (len(file_table) and np.any(sci_match)):
        return
    sci = file_table[sci_match]
    seen = set()
    for i in range(len(sci)):
        vals = tuple(str(sci[col][i]) for col in name_cols)
        if vals in seen:
            continue
        seen.add(vals)
        cal_name = get_name(paths, *vals)
        if _master_exists(tel, paths, cal_name):
            continue
        setup = ', '.join(
            f'{col.lower()}={val}' for col, val in zip(name_cols, vals)
        ).replace('binning=', 'bin=')
        _emit(log, (
            f'No master {cal_kind} for science setup {setup}. '
            f'Need at least {nmin_images} {cal_kind} frames for this detector '
            f'configuration (expected {cal_name}). Science in this '
            f'setup will be skipped; other setups will still be reduced.'
        ))


def _not_enough_frames(cal_kind, cal_type, n_found, nmin_images, log):
    _emit(log, (
        f'Not enough {cal_kind} frames for {cal_type} '
        f'({n_found} < {nmin_images}); science with this '
        f'setup will be skipped.'
    ))


def do_bias(file_table, tel, paths, nmin_images=3, log=None):
    """Build master bias frames from file_table; skip if instrument has no bias.

    Parameters
    ----------
    file_table : astropy.table.Table
        File list from sort_files (Type, CalType, File, Amp, Binning).
    tel : Instrument
        Instrument instance (bias, match_type_keywords, create_bias, get_mbias_name).
    paths : dict
        Paths dict from options.add_paths.
    nmin_images : int, optional
        Minimum images per CalType to build master. Default is 3.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
        Master bias FITS written to paths. Missing setups are logged and
        skipped; the pipeline continues with remaining targets.
    """
    # Exit if telescope does not require bias
    if not tel.bias:
        return(None)

    kwds = tel.filetype_keywords
    bias_match = tel.match_type_keywords(kwds['BIAS'], file_table)
    bias_table = file_table[bias_match]

    bias_num = 0
    for cal_type in np.unique(bias_table['CalType']):
        mask = bias_table['CalType']==cal_type
        cal_table = bias_table[mask]

        # Skip if cal_table does not have enough images
        if len(cal_table)<nmin_images:
            _not_enough_frames('bias', cal_type, len(cal_table), nmin_images, log)
            continue
        else:
            if log: 
                log.info(f'Generating bias image with {len(cal_table)} images.')
            else:
                print(f'Generating bias image with {len(cal_table)} images.')
        
        bias_num += 1
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        bias_name = tel.get_mbias_name(paths, amp, binn)
        
        if os.path.exists(bias_name):
            if log: log.info(f'Master bias {bias_name} exists.')
        else:
            if log: log.info(f'Master bias is being created...')
            t1 = time.time()
            if log: log.info('Processing bias files.')
            tel.create_bias(cal_table['File'], amp, binn, paths, 
                log=log)
            t2 = time.time()

            if log: log.info(f'Master bias creation completed in {t2-t1} sec')

    if bias_num==0:
        _emit(log, (
            'No usable master bias could be built (no bias frames, or fewer '
            'than the minimum per detector setup). Science that requires a '
            'bias will be skipped; other setups will still be reduced.'
        ))

    _warn_missing_science_masters(
        file_table, tel, paths, nmin_images, log, 'bias',
        tel.get_mbias_name, ('Amp', 'Binning'),
    )

def do_dark(file_table, tel, paths, nmin_images=3, log=None):
    """Build master dark frames from file_table; skip if instrument has no dark.

    Parameters
    ----------
    file_table : astropy.table.Table
        File list from sort_files.
    tel : Instrument
        Instrument instance (dark, bias, load_bias, create_dark, get_mdark_name).
    paths : dict
        Paths dict from options.add_paths.
    nmin_images : int, optional
        Minimum images per CalType to build master. Default is 3.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
        Master dark FITS written to paths. Missing setups are logged and
        skipped when the instrument requires darks.
    """
    # Exit if telescope does not require dark
    if not tel.dark:
        return(None)

    kwds = tel.filetype_keywords
    dark_match = tel.match_type_keywords(kwds['DARK'], file_table)
    dark_table = file_table[dark_match]

    dark_num = 0
    for cal_type in np.unique(dark_table['CalType']):
        mask = dark_table['CalType']==cal_type
        cal_table = dark_table[mask]

        # Skip if cal_table does not have enough images
        if len(cal_table)<nmin_images:
            _not_enough_frames('dark', cal_type, len(cal_table), nmin_images, log)
            continue
        else:
            if log: 
                log.info(f'Generating dark image with {len(cal_table)} images.')
            else:
                print(f'Generating dark image with {len(cal_table)} images.')
            
        exp = cal_table['Exp'][0]
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        dark_name = tel.get_mdark_name(paths, amp, binn)

        if os.path.exists(dark_name):
            if log: log.info(f'Master dark {dark_name} exists.')
            dark_num += 1
        else:
            if log: log.info(f'Master dark is being created...')
            mbias = None
            if tel.bias:
                if log: log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(paths, amp, binn)
                except Exception as e:
                    _emit(log, (
                        f'Cannot build master dark for exposure={exp}, '
                        f'amp={amp}, bin={binn}: no master bias ({e}). '
                        f'Science in this setup will be skipped; other '
                        f'setups will still be reduced.'
                    ))
                    continue

            t1 = time.time()
            tel.create_dark(cal_table['File'], amp, binn,
                paths, mbias=mbias, log=log)
            t2 = time.time()
            if log: log.info(f'Master dark creation completed in {t2-t1} sec.')
            dark_num += 1

    if dark_num==0:
        _emit(log, (
            'No usable master dark could be built (no dark frames, or fewer '
            'than the minimum per detector setup). Science that requires a '
            'dark will be skipped; other setups will still be reduced.'
        ))

    _warn_missing_science_masters(
        file_table, tel, paths, nmin_images, log, 'dark',
        tel.get_mdark_name, ('Amp', 'Binning'),
    )

def do_flat(file_table, tel, paths, nmin_images=3, log=None):
    """Build master flat frames from file_table; skip if instrument has no flat.

    Parameters
    ----------
    file_table : astropy.table.Table
        File list from sort_files.
    tel : Instrument
        Instrument instance (flat, match_type_keywords, create_flat, get_mflat_name).
    paths : dict
        Paths dict from options.add_paths.
    nmin_images : int, optional
        Minimum images per (filter, amp, binning) to build master. Default is 3.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
        Master flat FITS written to paths. Missing setups are logged and
        skipped when the instrument requires flats. Night cals stay in
        paths['cal']; packaged caldb flats still satisfy a science setup.
    """
    # Exit if telescope does not require flats
    if not tel.flat:
        return(None)

    kwds = tel.filetype_keywords
    flat_match = tel.match_type_keywords(kwds['FLAT'], file_table)
    flat_table = file_table[flat_match]

    # Do not rewrite paths['cal']: night biases/darks live in red/cals.
    # Missing flats are loaded later from red/cals, then packaged caldb.
    flat_num = 0
    for cal_type in np.unique(flat_table['CalType']):
        mask = flat_table['CalType']==cal_type
        cal_table = flat_table[mask]

        # Skip if cal_table does not have enough images
        if len(cal_table)<nmin_images:
            _not_enough_frames('flat', cal_type, len(cal_table), nmin_images, log)
            continue
        else:
            if log: 
                log.info(f'Generating flat image with {len(cal_table)} images.')
            else:
                print(f'Generating flat image with {len(cal_table)} images.')

        fil = cal_table['Filter'][0]
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        is_science = np.any([f=='SCIENCE' for f in cal_table['Type']])

        flat_name = tel.get_mflat_name(paths, fil, amp, binn)

        if os.path.exists(flat_name):
            if log: log.info(f'Master flat {flat_name} exists.')
            flat_num += 1
        else:
            if log: log.info(f'Master flat is being created...')
            mbias = None
            mdark = None
            if tel.bias:
                if log: log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(paths, amp, binn)
                except Exception as e:
                    _emit(log, (
                        f'Cannot build master flat for filter={fil}, '
                        f'amp={amp}, bin={binn}: no master bias ({e}). '
                        f'Science in this setup will be skipped; other '
                        f'setups will still be reduced.'
                    ))
                    continue

            if tel.dark:
                if log: log.info('Loading master dark.')
                try:
                    mdark = tel.load_dark(paths, amp, binn)
                except Exception as e:
                    _emit(log, (
                        f'Cannot build master flat for filter={fil}, '
                        f'amp={amp}, bin={binn}: no master dark ({e}). '
                        f'Science in this setup will be skipped; other '
                        f'setups will still be reduced.'
                    ))
                    continue

            t1 = time.time()
            tel.create_flat(cal_table['File'], fil, amp, binn,
                paths, mbias=mbias, mdark=mdark, is_science=is_science, 
                log=log)
            t2 = time.time()            
            if log: log.info(f'Master flat creation completed in {t2-t1} sec')
            flat_num += 1

    if flat_num==0:
        # Night flats are optional when a packaged caldb master exists.
        _emit(log, (
            'No night master flat was built (no flat frames, or fewer than '
            'the minimum per detector setup). Science will use a packaged '
            'caldb master if one exists; setups without a night or caldb '
            'flat will be skipped.'
        ), level='error')

    _warn_missing_science_masters(
        file_table, tel, paths, nmin_images, log, 'flat',
        tel.get_mflat_name, ('Filter', 'Amp', 'Binning'),
    )

