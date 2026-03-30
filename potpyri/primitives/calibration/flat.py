"""Master flat frame creation from pipeline file tables."""
from __future__ import annotations

import logging
import os
import time

import numpy as np

from potpyri.primitives.base_primitive import BasePrimitive

__all__ = ["FlatPrimitive", "do_flat", "flat_calibration_worker"]


def flat_calibration_worker(file_table, tel, paths, nmin_images, log):
    """Core flat master creation (used by :class:`FlatPrimitive`)."""
    kwds = tel.filetype_keywords
    flat_match = tel.match_type_keywords(kwds['FLAT'], file_table)
    flat_table = file_table[flat_match]

    if tel.flat and len(flat_table) == 0:
        paths['cal'] = paths['caldb']
        return

    for cal_type in np.unique(flat_table['CalType']):
        mask = flat_table['CalType'] == cal_type
        cal_table = flat_table[mask]

        if len(cal_table) < nmin_images:
            if log:
                log.info('No flat images were provided for this setup.')
            else:
                print('No flat images were provided for this setup.')
            continue
        if log:
            log.info(f'Generating flat image with {len(cal_table)} images.')
        else:
            print(f'Generating flat image with {len(cal_table)} images.')

        fil = cal_table['Filter'][0]
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        is_science = np.any([f == 'SCIENCE' for f in cal_table['Type']])

        flat_name = tel.get_mflat_name(paths, fil, amp, binn)

        if os.path.exists(flat_name):
            if log:
                log.info(f'Master flat {flat_name} exists.')
        else:
            if log:
                log.info('Master flat is being created...')
            mbias = None
            mdark = None
            if tel.bias:
                if log:
                    log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(paths, amp, binn)
                except Exception:
                    if log:
                        log.error(
                            f'No master bias found for this configuration, '
                            f'skipping master flat creation for filter {fil}, '
                            f'{amp} amps, {binn} binning.')
                    continue

            if tel.dark:
                if log:
                    log.info('Loading master dark.')
                try:
                    mdark = tel.load_dark(paths, amp, binn)
                except Exception:
                    if log:
                        log.error(
                            f'No master dark found for this configuration, '
                            f'skipping master flat creation for filter {fil}, '
                            f'{amp} amps, {binn} binning.')
                    continue

            t1 = time.time()
            tel.create_flat(
                cal_table['File'], fil, amp, binn, paths,
                mbias=mbias, mdark=mdark, is_science=is_science, log=log)
            t2 = time.time()
            if log:
                log.info(f'Master flat creation completed in {t2-t1} sec')


class FlatPrimitive(BasePrimitive):
    """Build master flat frames from a file table (see :func:`do_flat`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'nmin_images': 3,
    }

    def _perform(self):
        flat_calibration_worker(
            self.input, self.tel, self.paths, self.nmin_images, self.log)
        return {'output': None}


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
        Master flat FITS written to paths.
    """
    if not tel.flat:
        return None
    return FlatPrimitive(nmin_images=nmin_images).apply(
        input=file_table, tel=tel, paths=paths, log=log)
