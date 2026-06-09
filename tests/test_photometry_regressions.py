"""Regression tests for photometry failures seen in production pipeline runs.

Covers stack-trace scenarios where:
- ``DAOStarFinder`` returned ``None`` → ``'NoneType' object is not subscriptable``
- photutils 3.x centroid columns → ``KeyError: 'xcentroid'`` in ``extract_aperture_stats``
- ``photloop`` swallowed errors → stack had no ``APPPHOT`` → ``absphot`` ``KeyError``
"""
from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Column, Table

from potpyri.instruments import instrument_getter
from potpyri.primitives import absphot, photometry


def _write_stack_fits(path, shape=(64, 64), add_star=False):
    """Minimal stack FITS (SCI, MASK, ERROR) for photometry unit tests."""
    data = np.ones(shape, dtype=np.float32) * 50.0
    if add_star:
        cy, cx = shape[0] // 2, shape[1] // 2
        y, x = np.ogrid[: shape[0], : shape[1]]
        data += 2000.0 * np.exp(
            -((x - cx) ** 2 + (y - cy) ** 2) / (2 * 2.0**2)
        ).astype(np.float32)
    mask = np.zeros(shape, dtype=np.uint8)
    err = np.ones(shape, dtype=np.float32) * 5.0
    hdr = fits.Header()
    hdr['FILTER'] = 'r'
    hdr['EXTNAME'] = 'SCI'
    sci = fits.PrimaryHDU(data=data, header=hdr)
    sci.name = 'SCI'
    hdu = fits.HDUList(
        [sci, fits.ImageHDU(data=mask, name='MASK'), fits.ImageHDU(data=err, name='ERROR')]
    )
    hdu.writeto(path, overwrite=True)
    return path


def test_get_star_catalog_daofind_none_raises_photometry_error(monkeypatch):
    """DAOStarFinder returning None must raise PhotometryError, not TypeError."""

    class _FinderReturnsNone:
        def __init__(self, *args, **kwargs):
            pass

        def __call__(self, *args, **kwargs):
            return None

    monkeypatch.setattr(photometry, 'DAOStarFinder', _FinderReturnsNone)
    img = np.ones((32, 32), dtype=float) * 50.0
    mask = np.zeros((32, 32), dtype=bool)
    err = np.ones((32, 32), dtype=float)

    with pytest.raises(photometry.PhotometryError, match='None'):
        photometry.get_star_catalog(img, mask, err, fwhm_init=4.0, threshold=1.0)


def test_extract_aperture_stats_row_access_photutils3_centroids():
    """photutils 3 DeprecatedColumnQTable rows must not KeyError on xcentroid."""
    stars = Table(
        {
            'x_centroid': Column([16.0]),
            'y_centroid': Column([16.0]),
            'peak': Column([100.0]),
            'flux': Column([50.0]),
        },
    )
    stars.deprecation_map = {'xcentroid': 'x_centroid', 'ycentroid': 'y_centroid'}
    stars = photometry._normalize_daofind_catalog(stars)

    img = np.ones((32, 32), dtype=float) * 50.0
    img[16, 16] = 500.0
    mask = np.zeros((32, 32), dtype=bool)
    err = np.ones((32, 32), dtype=float) * 5.0

    stats = photometry.extract_aperture_stats(
        img, mask, err, stars, aperture_radius=5.0, log=None,
    )
    assert len(stats) == 1
    assert stats['flux_best'][0] > 0


def test_photloop_all_attempts_fail_raises_and_leaves_no_apphot(tmp_path, monkeypatch):
    """Exhausted photloop attempts must raise PhotometryError and not write APPPHOT."""
    stack = _write_stack_fits(tmp_path / 'stack.fits')

    def _always_fail(*args, **kwargs):
        raise photometry.PhotometryError('simulated photometry failure')

    monkeypatch.setattr(photometry, 'do_phot', _always_fail)

    with pytest.raises(photometry.PhotometryError, match='no successful run'):
        photometry.photloop(
            str(stack),
            phot_sn_min=3.0,
            phot_sn_max=20.0,
            fwhm_init=5.0,
            log=None,
        )

    with fits.open(stack) as hdul:
        assert 'APPPHOT' not in [h.name for h in hdul]


def test_failed_photometry_then_absphot_does_not_keyerror_appphot(tmp_path, monkeypatch):
    """Regression: trace ended with KeyError APPPHOT when photometry never succeeded."""
    stack = _write_stack_fits(tmp_path / 'stack.fits')

    monkeypatch.setattr(
        photometry,
        'do_phot',
        lambda *a, **k: (_ for _ in ()).throw(
            photometry.PhotometryError('simulated failure'),
        ),
    )

    with pytest.raises(photometry.PhotometryError):
        photometry.photloop(
            str(stack), phot_sn_min=3.0, phot_sn_max=5.0, fwhm_init=5.0, log=None,
        )

    tel = instrument_getter('LRIS')
    ok = absphot.find_zeropoint(str(stack), tel, log=None)
    assert ok is False

    with fits.open(stack) as hdul:
        assert 'APPPHOT' not in [h.name for h in hdul]


def test_photloop_succeeds_writes_appphot(tmp_path, monkeypatch):
    """Successful photloop must append APPPHOT (guards silent no-op regressions)."""
    stack = _write_stack_fits(tmp_path / 'stack.fits', add_star=True)

    def _fake_do_phot(img_file, star_param=None, log=None, **kwargs):
        with fits.open(img_file, mode='update') as hdul:
            tbl = Table({'flux': [1.0]})
            ext = fits.BinTableHDU(tbl)
            ext.name = 'APPPHOT'
            hdul.append(ext)
            hdul.flush()

    monkeypatch.setattr(photometry, 'do_phot', _fake_do_phot)

    photometry.photloop(
        str(stack), phot_sn_min=3.0, phot_sn_max=5.0, fwhm_init=5.0, log=None,
    )

    with fits.open(stack) as hdul:
        assert 'APPPHOT' in [h.name for h in hdul]
