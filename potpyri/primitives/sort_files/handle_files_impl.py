"""Worker that globs raw paths, runs sort_files, or reloads existing table."""
from __future__ import annotations

import glob
import logging
import os
import sys

import numpy as np
from astropy.io import ascii

from .sort_table import sort_files


def _handle_files_worker(file_list, paths, tel, incl_bad=False, proc=None,
        no_redo=False, log=None):
    """Implementation for :class:`SortFilesPrimitive` / :func:`handle_files`."""
    file_table = None

    if os.path.exists(file_list):
        if no_redo:
            file_table = ascii.read(file_list, format='fixed_width')

            file_table['Target'] = file_table['Target'].astype(str)
            file_table['TargType'] = file_table['TargType'].astype(str)
            file_table['Filter'] = file_table['Filter'].astype(str)
            file_table['Amp'] = file_table['Amp'].astype(str)
            file_table['Binning'] = file_table['Binning'].astype(str)
            file_table['Exp'] = file_table['Exp'].astype(str)
            file_table['Type'] = file_table['Type'].astype(str)
            file_table['CalType'] = file_table['CalType'].astype(str)
            file_table['Time'] = file_table['Time'].astype(np.float64)

            return file_table
        elif os.path.exists(file_list):
            os.remove(file_list)

    files = glob.glob(os.path.join(paths['raw'], tel.raw_format(proc)))+\
            glob.glob(os.path.join(paths['data'], tel.raw_format(proc)))+\
            glob.glob(os.path.join(paths['bad'], tel.raw_format(proc)))

    if log:
        log.info('Sorting files and creating file lists.')
    if len(files) != 0:
        if log:
            log.info(f'{len(files)} files found.')
        file_table = sort_files(files, file_list, tel, paths, incl_bad=incl_bad,
            log=log)
    else:
        if log:
            log.critical('No files found, please check data path and rerun.')
        logging.shutdown()
        sys.exit(-1)

    if file_table is None or len(file_table) == 0:
        if log:
            log.critical('No good files found, please check files and rerun.')
        logging.shutdown()
        sys.exit(-1)

    return file_table
