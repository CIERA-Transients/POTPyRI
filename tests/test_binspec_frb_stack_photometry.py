"""Integration test: photometry + zeropoint on FRB BINOSPEC stack.

Motivation: photutils 3.x ``DAOStarFinder`` uses ``x_centroid``/``y_centroid``
column names; POTPyRI historically expected ``xcentroid``/``ycentroid``, which
raised ``KeyError: 'xcentroid'`` inside ``extract_aperture_stats``.

The stack must be provided locally (large file, not in git). Set::

    export POTPYRI_FRB_STACK_TEST=/absolute/path/to/FRB20260213A.r.ut260322.2.11.stk.fits

or use the default ``$HOME/FRB20260213A.r.ut260322.2.11.stk.fits``. The test skips
if that path is missing.

Network (Vizier / Pan-STARRS1) is required for the zeropoint step; unreachable
catalog service skips the test similarly to ``test_absphot``.
"""
from __future__ import annotations

import os
import shutil

import numpy as np
import pytest
import requests
from astropy.io import fits

from potpyri.instruments import instrument_getter
from potpyri.primitives import absphot, photometry
from potpyri.utils import logger, options

_DEFAULT_FRB_STACK = os.path.join(
    os.path.expanduser('~'), 'FRB20260213A.r.ut260322.2.11.stk.fits'
)


def _frb_stack_path() -> str:
    return os.environ.get('POTPYRI_FRB_STACK_TEST', _DEFAULT_FRB_STACK)


def _strip_optional_hdu(hdul):
    for key in ('APPPHOT', 'PSFPHOT', 'PSFSTARS', 'RESIDUAL', 'PSF'):
        if key in hdul:
            del hdul[key]


def _sanitize_error(hdul):
    err = hdul['ERROR'].data
    med = np.nanmedian(err)
    err[np.isnan(err)] = med
    err[err < 0.0] = med
    err[np.isinf(err)] = np.max(hdul['SCI'].data)


def _prepare_stack_copy(src: str, dest: str) -> None:
    shutil.copy2(src, dest)
    with fits.open(dest, mode='update') as hdul:
        _strip_optional_hdu(hdul)
        _sanitize_error(hdul)
        hdul.flush()


@pytest.mark.integration
def test_binspec_frb_stack_photometry_and_absphot(tmp_path):
    """BINOSPEC FRB stack completes photloop and absphot.find_zeropoint."""
    stack_src = _frb_stack_path()
    if not os.path.isfile(stack_src):
        pytest.skip(
            f'Missing FRB stack at {stack_src!r}; copy the FITS locally or set '
            'POTPYRI_FRB_STACK_TEST to its absolute path.'
        )

    stack = str(tmp_path / 'frb_stack.fits')
    _prepare_stack_copy(stack_src, stack)

    tel = instrument_getter('BINOSPEC')
    paths = options.add_paths(str(tmp_path), 'files.txt', tel)
    log = logger.get_log(paths['log'])

    try:
        try:
            photometry.photloop(
                stack,
                phot_sn_min=3.0,
                phot_sn_max=20.0,
                fwhm_init=5.0,
                log=log,
            )
        except photometry.PhotometryError as exc:
            pytest.fail(f'photometry.photloop failed: {exc}')

        with fits.open(stack) as hdul:
            names = [h.name for h in hdul]
            assert 'APPPHOT' in names, f'expected APPPHOT HDU, got {names!r}'
            assert 'PSFPHOT' in names
            assert 'PSFSTARS' in names
            assert 'PSF' in names
            n_app = len(hdul['APPPHOT'].data)
            assert n_app > 100, f'expected substantial catalog, got NOBJECT={n_app}'

        try:
            ok = absphot.find_zeropoint(stack, tel, log=log)
        except (requests.exceptions.ConnectionError, OSError) as exc:
            pytest.skip(f'Catalog / network unavailable for zeropoint: {exc}')

        assert ok is True, 'absphot.find_zeropoint returned False (see log)'

        with fits.open(stack) as hdul:
            for hdr in (hdul['PRIMARY'].header, hdul['SCI'].header):
                assert 'ZPTMAG' in hdr, 'ZPTMAG not written to PRIMARY/SCI'
                assert np.isfinite(hdr['ZPTMAG'])
                assert hdr.get('ZPTNSTAR', 0) >= 5
    finally:
        log.close()
