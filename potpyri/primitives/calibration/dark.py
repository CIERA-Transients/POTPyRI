"""Master dark frame creation from pipeline file tables."""
from __future__ import annotations

import logging
import os
import time

import numpy as np

from potpyri.primitives.base_primitive import BasePrimitive

__all__ = ["DarkPrimitive", "do_dark", "dark_calibration_worker"]


def dark_calibration_worker(file_table, tel, paths, nmin_images, log):
    """Core dark master creation (used by :class:`DarkPrimitive`)."""
    kwds = tel.filetype_keywords
    dark_match = tel.match_type_keywords(kwds['DARK'], file_table)
    dark_table = file_table[dark_match]

    for cal_type in np.unique(dark_table['CalType']):
        mask = dark_table['CalType'] == cal_type
        cal_table = dark_table[mask]

        if len(cal_table) < nmin_images:
            if log:
                log.info('No dark images were provided for this setup.')
            else:
                print('No dark images were provided for this setup.')
            continue
        if log:
            log.info(f'Generating dark image with {len(cal_table)} images.')
        else:
            print(f'Generating dark image with {len(cal_table)} images.')

        exp = cal_table['Exp'][0]
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        dark_name = tel.get_mdark_name(paths, amp, binn)

        if os.path.exists(dark_name):
            if log:
                log.info(f'Master dark {dark_name} exists.')
        else:
            if log:
                log.info('Master dark is being created...')
            mbias = None
            if tel.bias:
                if log:
                    log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(paths, amp, binn)
                except Exception:
                    if log:
                        log.error(
                            f'No master bias found for this configuration, '
                            f'skipping master dark creation for exposure {exp}, '
                            f'{amp} amps and {binn} binning.')
                    continue

            t1 = time.time()
            tel.create_dark(cal_table['File'], amp, binn, paths, mbias=mbias, log=log)
            t2 = time.time()
            if log:
                log.info(f'Master dark creation completed in {t2-t1} sec.')


class DarkPrimitive(BasePrimitive):
    """Build master dark frames from a file table (see :func:`do_dark`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'nmin_images': 3,
    }

    def _perform(self):
        dark_calibration_worker(
            self.input, self.tel, self.paths, self.nmin_images, self.log)
        return {'output': None}


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
        Master dark FITS written to paths.
    """
    if not tel.dark:
        return None
    return DarkPrimitive(nmin_images=nmin_images).apply(
        input=file_table, tel=tel, paths=paths, log=log)
