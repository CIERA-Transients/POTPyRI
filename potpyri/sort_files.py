#!/usr/bin/env python

"Function to sort files for main_pipeline."
"Authors: Owen Eskandari, Kerry Paterson"

__version__ = "2.7" #last updated 18/08/2021

from astropy.io import fits, ascii
from astropy.table import Table
import os
import time
import shutil
import glob
import numpy as np
import gzip
import zlib
import logging
import sys
import re

def is_bad(hdr, tel):
    keywords = tel.bad_keywords
    values = tel.bad_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    bad = np.any([bool(re.search(v, str(hdr[k]).lower())) 
        for k,v in zip(keywords,values)])

    return(bad)

def is_spec(hdr, tel):
    keywords = tel.spec_keywords
    values = tel.spec_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    spec = np.all([bool(re.search(v, str(hdr[k]).lower())) 
        for k,v in zip(keywords,values)])

    return(spec)

def is_flat(hdr, tel):
    keywords = tel.flat_keywords
    values = tel.flat_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    flat = np.all([bool(re.search(v, str(hdr[k]).lower())) 
        for k,v in zip(keywords,values)])

    return(flat)

def is_dark(hdr, tel):
    keywords = tel.dark_keywords
    values = tel.dark_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    dark = np.all([bool(re.search(v, str(hdr[k]).lower())) 
        for k,v in zip(keywords,values)])

    return(dark)

def is_bias(hdr, tel):
    keywords = tel.bias_keywords
    values = tel.bias_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    bias = np.all([bool(re.search(v, str(hdr[k]).lower())) 
        for k,v in zip(keywords,values)])
    
    return(bias)

def is_science(hdr, tel):
    keywords = tel.science_keywords
    values = tel.science_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    science = np.all([bool(re.search(v, str(hdr[k]).lower())) 
        for k,v in zip(keywords,values)])

    # Check minimum exposure time
    if tel.min_exptime:
        exptime = tel.get_exptime(hdr)
        if exptime < tel.min_exptime:
            return(False)

    return(science)

# Overall method to handle files:
def handle_files(file_list, tel, paths, incl_bad=False, proc=None, 
    no_redo=False, log=None):

    file_table = None

    # Always regenerate file list from existing data in data, raw, and bad
    if os.path.exists(file_list): 
        if no_redo:
            file_table = Table.read(file_list, format='ascii.ecsv')
            return(file_table)
        else:
            os.remove(file_list)

    files = glob.glob(os.path.join(paths['raw'], tel.raw_format(proc)))+\
            glob.glob(os.path.join(paths['data'], tel.raw_format(proc)))+\
            glob.glob(os.path.join(paths['bad'], tel.raw_format(proc)))


    if log: log.info('Sorting files and creating file lists.')
    if len(files)!=0:
        log.info(f'{len(files)} files found.')
        file_table = sort_files(files, file_list, tel, paths, incl_bad=incl_bad, 
            log=log)
    else:
        log.critical('No files found, please check data path and rerun.')
        logging.shutdown()
        sys.exit(-1)

    return(file_table)

# Sort the calibration files:
def sort_files(files, file_list, tel, paths, incl_bad=False, log=None):

    '''

    Function used to sort a list of files into a dictionary of files sorted by filter.

    Parameters
    ----------

    :param files: list (string)
        List of strings of files (path should be included).

    :param tel: Telescope
        Telescope with parameters and methods to sort data.

    :param incl_bad: bool, optional

    :param log: log, optional
        Overview log used to write the object and date observed (if ``date`` parameter is not ``None``).
        If no log is inputted, information is printed out instead of being written to ``log2``.
        Default is ``None``.


    Returns
    -------

    :return: python dictionary
        Dictionary of files. Key is the filter of the file, values are the file names themselves.

    '''

    t_start = time.time()
    
    if log: log.info(f'Running sort_files version: {__version__}')

    bad_num = 0
    spec_num = 0
    sci_num = 0
    bias_num = 0
    dark_num = 0
    flat_num = 0

    params = ('File','Target','TargType','Filter','Amp','Binning','Exp','Type',
        'CalType','Time')
    dtypes = ('S','S','S','S','S','S','S','S','S','float64')
    file_table = Table(names=params, dtype=dtypes)

    for i, f in enumerate(sorted(files)):
        try:
            file_open = fits.open(f, mode='readonly')
            ext = tel.raw_header_ext
            hdr = file_open[ext].header
            check_data = file_open[ext].data
            file_open._verify()
        except IndexError:
            if log: log.error(f'Moving file {f} to bad due to error opening file.')
            moved_path = paths['bad']
            if os.path.dirname(f)!=moved_path: shutil.move(f, paths['bad'])
            continue
        except (TypeError, gzip.BadGzipFile, zlib.error):
            if log: log.error(f'Moving file {f} to bad due to corrupted data.')
            moved_path = paths['bad']
            if os.path.dirname(f)!=moved_path: shutil.move(f, paths['bad'])
            continue

        # Edit raw data headers if necessary
        # TODO - adjust method to work only on individual files
        tel.edit_raw_headers([f], log=log)

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
            elif is_spec(hdr, tel):
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
            elif is_science(hdr, tel):
                file_type = 'SCIENCE'
                moved_path = paths['raw']
                sci_num += 1
            elif is_dark(hdr, tel):
                file_type = 'DARK'
                moved_path = paths['raw']
                dark_num += 1
            else:
                file_type = 'BAD'
                moved_path = paths['bad']
                bad_num += 1
        except Exception as e:
            if log: log.error(f'Moving file {f} to bad due to error: {e}')
            file_type = 'BAD'
            moved_path = paths['bad']
            bad_num += 1

        if os.path.dirname(f)!=moved_path:
            newfile = os.path.join(moved_path, os.path.basename(f))
            if not os.path.exists(newfile):
                shutil.move(f, moved_path)
            else:
                if log: log.info(f'Removing existing file: {f}')
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
            if log: log.critical(f'Lost track of file {f}->{currfile}')
            logging.shutdown()
            sys.exit(-1)

        if log: log.info(f'File {currfile} is {file_type}')

        if (file_type!='BAD' and file_type!='SPEC') or incl_bad:
            file_table.add_row((currfile,target,targ_type,fil,amp,binn,exp,
                file_type,cal_type,file_time))

    file_table.sort(['Type','Target','CalType','File'])

    ascii.write(file_table, file_list, format='fixed_width',
        formats={'Time':'%5.6f'})

    if sci_num>0 and log: log.info(f'{sci_num} imaging science files found.')
    if bias_num>0 and log: log.info(f'{bias_num} bias files found.')
    if flat_num>0 and log: log.info(f'{flat_num} flat files found.')
    if dark_num>0 and log: log.info(f'{dark_num} dark files found.')
    if bad_num>0 and log: log.info(f'{bad_num} bad files found and removed from reduction.')
    if spec_num>0 and log: log.info(f'{spec_num} spectroscopy files found and removed from reduction.')

    t_end = time.time()
    log.info(f'Sort_files ran in {t_end-t_start} sec')

    return file_table
