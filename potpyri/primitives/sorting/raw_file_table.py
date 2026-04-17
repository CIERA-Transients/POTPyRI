"""Build the fixed-width file list table from classified raw FITS paths."""
from __future__ import annotations

import gzip
import logging
import os
import shutil
import sys
import time
import zlib

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table

from potpyri._version import __version__

from .header_classification import (
    is_bad,
    is_bias,
    is_dark,
    is_flat,
    is_science,
    is_spec,
)


def build_raw_file_table(files, file_list, tel, paths, incl_bad=False, log=None):
    """Classify files by type (science, flat, bias, dark) and write file list table.

    Reads headers, applies instrument keyword rules, and writes a fixed-width
    file list.

    Parameters
    ----------
    files : list of str
        Paths to raw FITS files.
    file_list : str
        Path to write fixed-width file list.
    tel : Instrument
        Instrument instance.
    paths : dict
        Paths dict (raw, data, bad, work).
    incl_bad : bool, optional
        If True, include bad files in table. Default is False.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    astropy.table.Table
        Table with Target, Filter, Type, CalType, File, Exp, Time, etc.
    """

    t_start = time.time()

    if log:
        log.info(f'Running build_raw_file_table version: {__version__}')
    else:
        print(f'Running build_raw_file_table version: {__version__}')

    bad_num = 0
    spec_num = 0
    sci_num = 0
    bias_num = 0
    dark_num = 0
    flat_num = 0

    target = ""
    fil = ""
    amp = ""
    binn = ""
    exp = ""
    file_time = 0.0

    params = ('File','Target','TargType','Filter','Amp','Binning','Exp','Type',
        'CalType','Time')
    dtypes = ('S','S','S','S','S','S','S','S','S','float64')
    file_table = Table(names=params, dtype=dtypes)

    for i, f in enumerate(sorted(files)):
        try:
            with fits.open(f, mode='readonly') as file_open:
                ext = tel.raw_header_ext
                hdr = file_open[ext].header

                # Extend header to first extension?
                if tel.extend_header:
                    if len(file_open)>ext+1:
                        extra_hdr = file_open[ext+1].header
                        for key in extra_hdr.keys():
                            if key not in hdr.keys():
                                value = extra_hdr[key]
                                if isinstance(value, (str, int, float, complex,
                                    bool, np.floating, np.integer, np.bool_)):
                                    hdr[key] = value

                check_data = file_open[ext].data
                file_open._verify()
        except IndexError:
            if log:
                log.error(f'Moving file {f} to bad due to error opening file.')
            else:
                print(f'Moving file {f} to bad due to error opening file.')
            moved_path = paths['bad']
            if os.path.dirname(f)!=moved_path: shutil.move(f, paths['bad'])
            continue
        except (TypeError, gzip.BadGzipFile, zlib.error, OSError):
            if log:
                log.error(f'Moving file {f} to bad due to corrupted data.')
            else:
                print(f'Moving file {f} to bad due to corrupted data.')
            moved_path = paths['bad']
            if os.path.dirname(f)!=moved_path: shutil.move(f, paths['bad'])
            continue

        try:
            target = tel.get_target(hdr)
            fil = tel.get_filter(hdr)
            amp = tel.get_ampl(hdr)
            binn = tel.get_binning(hdr)
            exp = str(tel.get_exptime(hdr))
            file_time = tel.get_time(hdr)

            if (is_bad(hdr, tel) and not is_bias(hdr, tel) and
                not is_dark(hdr, tel) and not is_flat(hdr, tel)):
                file_type = 'BAD'
                moved_path = paths['bad']
                bad_num += 1
            elif (is_spec(hdr, tel) and not is_bias(hdr, tel)):
                file_type = 'SPEC'
                moved_path = paths['bad']
                spec_num += 1
            elif is_flat(hdr, tel):
                file_type = 'FLAT'
                moved_path = paths['raw']
                flat_num += 1
            elif is_bias(hdr, tel):
                file_type = 'BIAS'
                moved_path = paths['raw']
                bias_num += 1
            elif is_dark(hdr, tel):
                file_type = 'DARK'
                moved_path = paths['raw']
                dark_num += 1
            elif is_science(hdr, tel):
                file_type = 'SCIENCE'
                moved_path = paths['raw']
                sci_num += 1
            else:
                file_type = 'BAD'
                moved_path = paths['bad']
                bad_num += 1
        except Exception as e:
            if log:
                log.error(f'Moving file {f} to bad due to error: {e}')
            else:
                print(f'Moving file {f} to bad due to error: {e}')
            file_type = 'BAD'
            moved_path = paths['bad']
            bad_num += 1

        if os.path.dirname(f)!=moved_path:
            newfile = os.path.join(moved_path, os.path.basename(f))
            if not os.path.exists(newfile):
                shutil.move(f, moved_path)
            else:
                if log:
                    log.info(f'Removing existing file: {f}')
                else:
                    print(f'Removing existing file: {f}')
                os.remove(f)

        if file_type=='BIAS':
            # Bias only depends on amplifier and bin mode
            cal_type = f'{amp}_{binn}'
            target = 'BIAS'
            targ_type = cal_type
        elif file_type=='DARK':
            # Dark depends on exposure, amplifier, and bin
            cal_type = f'{exp}_{amp}_{binn}'
            target = 'DARK'
            targ_type = cal_type
        elif file_type=='FLAT':
            # Flat depends on filter, amplifier, and bin
            cal_type = f'{fil}_{amp}_{binn}'
            target = 'FLAT'
            targ_type = cal_type
        elif file_type=='SCIENCE':
            # Science is grouped by target, filter, amplifier, bin
            cal_type = f'{fil}_{amp}_{binn}'
            targ_type = f'{target}_{fil}_{amp}_{binn}'
        else:
            cal_type = ''
            targ_type = ''

        currfile = os.path.join(moved_path, os.path.basename(f))
        if not os.path.exists(currfile):
            if log:
                log.critical(f'Lost track of file {f}->{currfile}')
            else:
                print(f'Lost track of file {f}->{currfile}')
            logging.shutdown()
            sys.exit(-1)

        if log:
            log.info(f'File {i+1}/{len(files)}: {currfile} is {file_type},{target},{fil}')
        else:
            print(f'File {i+1}/{len(files)}: {currfile} is {file_type},{target},{fil}')

        if (file_type!='BAD' and file_type!='SPEC') or incl_bad:
            file_table.add_row((currfile,target,targ_type,fil,amp,binn,exp,
                file_type,cal_type,file_time))

    file_table.sort(['Type','Target','CalType','File'])

    if len(file_table)>0:
        ascii.write(file_table, file_list, format='fixed_width',
            formats={'Time':'%5.6f'}, overwrite=True)
    else:
        if log: log.critical('No good files were ingested')

    if sci_num>0 and log: log.info(f'{sci_num} imaging science files found.')
    if bias_num>0 and log: log.info(f'{bias_num} bias files found.')
    if flat_num>0 and log: log.info(f'{flat_num} flat files found.')
    if dark_num>0 and log: log.info(f'{dark_num} dark files found.')
    if bad_num>0 and log: log.info(f'{bad_num} bad files found and removed from reduction.')
    if spec_num>0 and log: log.info(f'{spec_num} spectroscopy files found and removed from reduction.')

    t_end = time.time()
    if log: log.info(f'build_raw_file_table ran in {t_end-t_start} sec')

    return(file_table)
