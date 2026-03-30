"""Master bias frame creation from pipeline file tables."""
from __future__ import annotations

import logging
import os
import sys
import time

import numpy as np

from potpyri.primitives.base_primitive import BasePrimitive

__all__ = ["BiasPrimitive", "do_bias", "bias_calibration_worker"]


def bias_calibration_worker(file_table, tel, paths, nmin_images, log):
    """Core bias master creation (used by :class:`BiasPrimitive`)."""
    kwds = tel.filetype_keywords
    bias_match = tel.match_type_keywords(kwds['BIAS'], file_table)
    bias_table = file_table[bias_match]

    bias_num = 0
    for cal_type in np.unique(bias_table['CalType']):
        mask = bias_table['CalType'] == cal_type
        cal_table = bias_table[mask]

        if len(cal_table) < nmin_images:
            if log:
                log.info('No bias images were provided for this setup.')
            else:
                print('No bias images were provided for this setup.')
            continue
        if log:
            log.info(f'Generating bias image with {len(cal_table)} images.')
        else:
            print(f'Generating bias image with {len(cal_table)} images.')

        bias_num += 1
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        bias_name = tel.get_mbias_name(paths, amp, binn)

        if os.path.exists(bias_name):
            if log:
                log.info(f'Master bias {bias_name} exists.')
        else:
            if log:
                log.info('Master bias is being created...')
            t1 = time.time()
            if log:
                log.info('Processing bias files.')
            tel.create_bias(cal_table['File'], amp, binn, paths, log=log)
            t2 = time.time()

            if log:
                log.info(f'Master bias creation completed in {t2-t1} sec')

    if bias_num == 0:
        if log:
            log.critical('No bias present, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)


class BiasPrimitive(BasePrimitive):
    """Build master bias frames from a file table (see :func:`do_bias`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'nmin_images': 3,
    }

    def _perform(self):
        bias_calibration_worker(
            self.input, self.tel, self.paths, self.nmin_images, self.log)
        return {'output': None}


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
        Master bias FITS written to paths; exits if no bias and instrument requires it.
    """
    if not tel.bias:
        return None
    return BiasPrimitive(nmin_images=nmin_images).apply(
        input=file_table, tel=tel, paths=paths, log=log)
