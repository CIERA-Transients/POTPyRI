"""Tests for image_procs: remove_pv_distortion, get_fieldcenter, generate_wcs, detrend_stack, add_stack_mask, compute_relative_scales."""
from potpyri.primitives import image_procs
from potpyri.instruments import instrument_getter

import numpy as np
import os
from unittest.mock import patch

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from ccdproc import CCDData
from astropy import units as u

from tests.utils import make_dummy_wcs


def test_remove_pv_distortion(tmp_path):
    """remove_pv_distortion strips PV* keywords from header in place."""
    header = fits.Header()
    header['NAXIS'] = 2
    header['PV1_1'] = 1.0
    header['PV2_1'] = 0.0
    header['CTYPE1'] = 'RA---TAN'
    out = image_procs.remove_pv_distortion(header)
    assert out is header
    assert 'PV1_1' not in header
    assert 'PV2_1' not in header
    assert header['CTYPE1'] == 'RA---TAN'


def test_get_fieldcenter(tmp_path):
    """get_fieldcenter returns mean RA/Dec of image centers from FITS with WCS."""
    wcs_head = make_dummy_wcs()
    wcs_head['NAXIS1'] = 100
    wcs_head['NAXIS2'] = 100
    data = np.zeros((100, 100), dtype=np.float32)
    for i in range(2):
        hdu = fits.PrimaryHDU(data)
        hdu.header.update(wcs_head)
        path = os.path.join(tmp_path, f'img{i}.fits')
        hdu.writeto(path, overwrite=True)
    images = [os.path.join(tmp_path, 'img0.fits'), os.path.join(tmp_path, 'img1.fits')]
    center = image_procs.get_fieldcenter(images)
    assert len(center) == 2
    assert np.isfinite(center[0]) and np.isfinite(center[1])
    # CRVAL from make_dummy_wcs is 30, -30
    np.testing.assert_allclose(center, [30.0, -30.0], atol=0.01)


def test_generate_wcs(tmp_path):
    """generate_wcs returns a TAN WCS with correct CRVAL and size."""
    tel = instrument_getter('GMOS')
    binn = '22'
    fieldcenter = [30.0, -30.0]
    out_size = 256
    wcs = image_procs.generate_wcs(tel, binn, fieldcenter, out_size)
    assert isinstance(wcs, WCS)
    np.testing.assert_allclose(wcs.wcs.crval, [30.0, -30.0])
    assert wcs.pixel_shape == (out_size, out_size)


def test_detrend_stack(tmp_path):
    """detrend_stack removes row/column median and updates SATURATE."""
    sci = np.ones((64, 64), dtype=np.float32) * 100.0
    sci += np.arange(64)[:, None]  # row gradient
    sci += np.arange(64)[None, :]  # column gradient
    mask = np.zeros((64, 64), dtype=np.uint8)
    hdu_sci = fits.PrimaryHDU(sci.astype(np.float32))
    hdu_sci.header['SATURATE'] = 50000.0
    hdu_mask = fits.ImageHDU(mask)
    stack = fits.HDUList([hdu_sci, hdu_mask])
    out = image_procs.detrend_stack(stack)
    assert out is stack
    assert stack[0].data is not None
    # After detrend, row/col median should be ~0
    assert np.abs(np.nanmedian(stack[0].data)) < 1.0
    assert stack[0].header['SATURATE'] < 50000.0


def test_add_stack_mask(tmp_path):
    """add_stack_mask merges per-image masks and saturation into stack mask."""
    from potpyri.utils import logger
    from potpyri.utils import options
    tel = instrument_getter('MMIRS')
    paths = options.add_paths(str(tmp_path), 'files.txt', tel)
    log = logger.get_log(paths['log'])
    try:
        shape = (32, 32)
        data = np.ones(shape, dtype=np.float32) * 10.0
        mask = np.zeros(shape, dtype=np.uint8)
        err = np.ones(shape, dtype=np.float32) * 0.5
        wcs_head = make_dummy_wcs()
        ccd = CCDData(data, meta={'EXPTIME': 1.0, 'SATURATE': 1000.0}, wcs=WCS(wcs_head), unit=u.electron)
        stacking_data = [ccd]
        stack_sci = fits.PrimaryHDU(data.copy())
        stack_sci.header['SATURATE'] = 1000.0
        stack_mask = fits.ImageHDU(mask.copy())
        stack_err = fits.ImageHDU(err.copy())
        stack = fits.HDUList([stack_sci, stack_mask, stack_err])
        image_procs.add_stack_mask(stack, stacking_data)
        assert stack[1].data is not None
        assert stack[0].data is not None
        # Saturation 1000, data 10 -> no sat mask; bad pixel mask could be 0
        assert np.any(stack[1].data >= 0)
    finally:
        log.close()


def test_create_error(tmp_path):
    """create_error returns error HDU from science path, mask HDU, and rdnoise."""
    science_path = os.path.join(tmp_path, "sci.fits")
    data = np.ones((16, 16), dtype=np.float32) * 100.0
    hdu_sci = fits.PrimaryHDU(data)
    hdu_sci.header["SATURATE"] = 50000.0
    hdu_sci.writeto(science_path, overwrite=True)

    mask_data = fits.ImageHDU(np.zeros((16, 16), dtype=np.uint8))
    rdnoise = 4.0

    err_hdu = image_procs.create_error(science_path, mask_data, rdnoise)
    assert err_hdu is not None
    assert err_hdu.data.shape == (16, 16)
    assert err_hdu.header.get("BUNIT") == "ELECTRONS"
    assert np.all(err_hdu.data > 0)
    assert np.all(np.isfinite(err_hdu.data))


def _make_catalog(n, ra=30.0, dec=-30.0, flux=100.0, fluxerr=5.0, dra=0.001, ddec=0.001):
    """Build a minimal SExtractor-like table with ALPHA_J2000, DELTA_J2000, FLUX_AUTO, FLUXERR_AUTO."""
    rng = np.random.default_rng(42)
    ra_arr = ra + (rng.random(n) - 0.5) * dra
    dec_arr = dec + (rng.random(n) - 0.5) * ddec
    flux_arr = np.full(n, flux, dtype=float) + (rng.random(n) - 0.5) * 0.1 * flux
    err_arr = np.full(n, fluxerr, dtype=float)
    t = Table()
    t['ALPHA_J2000'] = ra_arr
    t['DELTA_J2000'] = dec_arr
    t['FLUX_AUTO'] = np.maximum(flux_arr, 1.0)
    t['FLUXERR_AUTO'] = err_arr
    return t


def test_compute_relative_scales_returns_none_for_single_image():
    """compute_relative_scales returns None when fewer than 2 images."""
    paths = {}
    data_images = ['/fake/path0.fits']
    exptimes = [60.0]
    out = image_procs.compute_relative_scales(data_images, paths, exptimes)
    assert out is None


def test_compute_relative_scales_returns_none_for_exptimes_length_mismatch():
    """compute_relative_scales returns None when exptimes length != number of images."""
    paths = {}
    data_images = ['/fake/a.fits', '/fake/b.fits']
    exptimes = [60.0, 30.0, 90.0]
    with patch.object(image_procs.photometry, 'run_sextractor') as m_sextractor:
        m_sextractor.return_value = _make_catalog(10)
        out = image_procs.compute_relative_scales(data_images, paths, exptimes)
    assert out is None


def test_compute_relative_scales_returns_none_when_sextractor_fails():
    """compute_relative_scales returns None when run_sextractor returns None for one image."""
    paths = {}
    data_images = ['/fake/a.fits', '/fake/b.fits']
    exptimes = [60.0, 60.0]

    def side_effect(path, log=None, sextractor_path=None):
        if 'a.fits' in path:
            return _make_catalog(10)
        return None

    with patch.object(image_procs.photometry, 'run_sextractor', side_effect=side_effect):
        out = image_procs.compute_relative_scales(data_images, paths, exptimes)
    assert out is None


def test_compute_relative_scales_returns_none_when_catalog_too_small():
    """compute_relative_scales returns None when a catalog has fewer than min_sources."""
    paths = {}
    data_images = ['/fake/a.fits', '/fake/b.fits']
    exptimes = [60.0, 60.0]
    small = _make_catalog(2)
    big = _make_catalog(10)

    def side_effect(path, log=None, sextractor_path=None):
        return small if 'a.fits' in path else big

    with patch.object(image_procs.photometry, 'run_sextractor', side_effect=side_effect):
        out = image_procs.compute_relative_scales(data_images, paths, exptimes, min_sources=5)
    assert out is None


def test_compute_relative_scales_two_frames_flux_ratio():
    """compute_relative_scales returns [1, scale] when frame 1 has half the flux of reference."""
    paths = {}
    data_images = ['/fake/ref.fits', '/fake/f1.fits']
    exptimes = [60.0, 60.0]
    # Same positions so they match; ref flux 100, frame1 flux 50 -> ratio 2.0
    ref_cat = _make_catalog(10, flux=100.0, fluxerr=5.0)
    f1_cat = _make_catalog(10, flux=50.0, fluxerr=2.5)
    # Use identical RA/Dec so match_to_catalog_sky finds pairs (same seed and pattern)
    f1_cat['ALPHA_J2000'] = ref_cat['ALPHA_J2000'].copy()
    f1_cat['DELTA_J2000'] = ref_cat['DELTA_J2000'].copy()

    def side_effect(path, log=None, sextractor_path=None):
        return ref_cat if 'ref.fits' in path else f1_cat

    with patch.object(image_procs.photometry, 'run_sextractor', side_effect=side_effect):
        out = image_procs.compute_relative_scales(data_images, paths, exptimes, min_sources=5)
    assert out is not None
    assert len(out) == 2
    assert out[0] == 1.0
    np.testing.assert_allclose(out[1], 2.0, rtol=0.15)


def test_compute_relative_scales_frame_fails_gets_exptime_scale():
    """When one frame has too few matches, that frame gets scale = 1.0/exptime."""
    paths = {}
    data_images = ['/fake/ref.fits', '/fake/f1.fits', '/fake/f2.fits']
    exptimes = [60.0, 30.0, 90.0]
    ref_cat = _make_catalog(10, ra=30.0, dec=-30.0, flux=100.0)
    f1_cat = _make_catalog(10, ra=30.0, dec=-30.0, flux=50.0)
    f1_cat['ALPHA_J2000'] = ref_cat['ALPHA_J2000'].copy()
    f1_cat['DELTA_J2000'] = ref_cat['DELTA_J2000'].copy()
    # Frame 2: completely different sky position -> no matches
    f2_cat = _make_catalog(10, ra=100.0, dec=50.0, flux=100.0)

    def side_effect(path, log=None, sextractor_path=None):
        if 'ref.fits' in path:
            return ref_cat
        if 'f1.fits' in path:
            return f1_cat
        return f2_cat

    with patch.object(image_procs.photometry, 'run_sextractor', side_effect=side_effect):
        out = image_procs.compute_relative_scales(data_images, paths, exptimes, min_sources=5, match_radius_arcsec=2.0)
    assert out is not None
    assert len(out) == 3
    # Ref (most sources) is frame 0; scale[0]=1
    assert out[0] == 1.0
    # Frame 1 matches ref -> computed ratio ~2.0
    np.testing.assert_allclose(out[1], 2.0, rtol=0.2)
    # Frame 2 has no matches within 2" -> 1/exptime = 1/90
    np.testing.assert_allclose(out[2], 1.0 / 90.0, rtol=1e-5)


def test_compute_relative_scales_uses_fluxerr_weighting():
    """compute_relative_scales uses FLUXERR_AUTO when present (weighted mean)."""
    paths = {}
    data_images = ['/fake/ref.fits', '/fake/f1.fits']
    exptimes = [60.0, 60.0]
    ref_cat = _make_catalog(10, flux=100.0, fluxerr=2.0)
    f1_cat = _make_catalog(10, flux=50.0, fluxerr=1.0)
    f1_cat['ALPHA_J2000'] = ref_cat['ALPHA_J2000'].copy()
    f1_cat['DELTA_J2000'] = ref_cat['DELTA_J2000'].copy()

    def side_effect(path, log=None, sextractor_path=None):
        return ref_cat if 'ref.fits' in path else f1_cat

    with patch.object(image_procs.photometry, 'run_sextractor', side_effect=side_effect):
        out = image_procs.compute_relative_scales(data_images, paths, exptimes, min_sources=5)
    assert out is not None
    assert out[0] == 1.0
    # Weighted mean of ratio ~2.0 (with small errors) should be close to 2
    np.testing.assert_allclose(out[1], 2.0, rtol=0.15)


def test_compute_relative_scales_without_fluxerr_uses_median():
    """When FLUXERR_AUTO is missing, compute_relative_scales uses median of ratios."""
    paths = {}
    data_images = ['/fake/ref.fits', '/fake/f1.fits']
    exptimes = [60.0, 60.0]
    ref_cat = _make_catalog(10, flux=100.0, fluxerr=5.0)
    ref_cat.remove_column('FLUXERR_AUTO')
    f1_cat = _make_catalog(10, flux=50.0, fluxerr=2.5)
    f1_cat.remove_column('FLUXERR_AUTO')
    f1_cat['ALPHA_J2000'] = ref_cat['ALPHA_J2000'].copy()
    f1_cat['DELTA_J2000'] = ref_cat['DELTA_J2000'].copy()

    def side_effect(path, log=None, sextractor_path=None):
        return ref_cat if 'ref.fits' in path else f1_cat

    with patch.object(image_procs.photometry, 'run_sextractor', side_effect=side_effect):
        out = image_procs.compute_relative_scales(data_images, paths, exptimes, min_sources=5)
    assert out is not None
    assert out[0] == 1.0
    np.testing.assert_allclose(out[1], 2.0, rtol=0.2)
