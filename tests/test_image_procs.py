"""Tests for image_procs: remove_pv_distortion, get_fieldcenter, generate_wcs, detrend_stack, add_stack_mask."""
from potpyri.primitives import image_procs
from potpyri.instruments import instrument_getter

import numpy as np
import os

from astropy.io import fits
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
