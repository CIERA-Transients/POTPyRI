"""Unit tests for instrument base class methods and instrument_getter."""
import os

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import CCDData
from astropy import units as u

from potpyri.instruments import instrument_getter
from potpyri.instruments.GMOS import GMOS
from potpyri.instruments.instrument import (
    Instrument,
    _sanitize_calibration_header,
    _read_calibration_ccd,
)
from potpyri.utils import options
from potpyri.utils import logger


def test_instrument_getter_unsupported_raises():
    """instrument_getter raises when instrument is not supported and log is None."""
    with pytest.raises(Exception, match="not supported"):
        instrument_getter("UNKNOWN_INSTRUMENT", log=None)


def test_instrument_getter_unsupported_with_log(tmp_path):
    """instrument_getter with log calls log.error and returns None for unsupported name."""
    from potpyri.instruments import __init__ as instruments_init
    tel = instrument_getter("GMOS")
    paths = options.add_paths(str(tmp_path), "files.txt", tel)
    log = logger.get_log(paths["log"])
    try:
        # When log is provided, getter still raises (code path: log.error then no return)
        # Actually re-reading the code: if log, it calls log.error and then falls through
        # and tel stays None, so it returns None. So we get None, not an exception.
        result = instrument_getter("UNSUPPORTED", log=log)
        assert result is None
    finally:
        log.close()


def test_match_type_keywords():
    """Instrument.match_type_keywords returns mask for Type column."""
    tel = GMOS()
    file_table = Table({"Type": ["SCIENCE", "BIAS", "SCIENCE", "FLAT"]})
    mask = tel.match_type_keywords("SCIENCE,FLAT", file_table)
    assert np.array_equal(mask, [True, False, True, True])


def test_needs_sky_subtraction():
    """GMOS needs_sky_subtraction True for z-band, False otherwise."""
    tel = GMOS()
    assert tel.needs_sky_subtraction("z") is True
    assert tel.needs_sky_subtraction("r") is False
    assert tel.needs_sky_subtraction("Z") is True


def test_get_pixscale():
    """get_pixscale returns instrument pixel scale."""
    tel = GMOS()
    assert tel.get_pixscale() == 0.0803


def test_get_rdnoise_get_gain():
    """get_rdnoise and get_gain use header if present else default."""
    tel = GMOS()
    hdr = fits.Header()
    assert tel.get_rdnoise(hdr) == 4.14
    assert tel.get_gain(hdr) == 1.63
    hdr["RDNOISE"] = 5.0
    hdr["GAIN"] = 2.0
    assert tel.get_rdnoise(hdr) == 5.0
    assert tel.get_gain(hdr) == 2.0


def test_get_target_get_filter_get_exptime():
    """Header getters for target, filter, exptime."""
    tel = GMOS()
    hdr = fits.Header({"OBJECT": "  NGC1234  ", "FILTER2": " r_G0326 ", "EXPTIME": 60.0})
    assert tel.get_target(hdr) == "NGC1234"
    assert tel.get_filter(hdr) == "r"
    assert tel.get_exptime(hdr) == 60.0


def test_get_ampl_get_binning():
    """GMOS get_ampl and get_binning from header."""
    tel = GMOS()
    hdr = fits.Header({"NCCDS": "1", "CCDSUM": "2 2"})
    assert tel.get_ampl(hdr) == "4"
    hdr["NCCDS"] = "2"
    assert tel.get_ampl(hdr) == "12"
    assert tel.get_binning(hdr) == "22"


def test_get_out_size():
    """get_out_size scales out_size by binning."""
    tel = GMOS()
    hdr = fits.Header({"CCDSUM": "2 2"})
    # GMOS out_size 3200, binn 2 -> 1600
    assert tel.get_out_size(hdr) == 1600


def test_get_time_get_number():
    """get_time and get_number from DATE-OBS and TIME-OBS."""
    tel = GMOS()
    hdr = fits.Header({
        "DATE-OBS": "2024-06-18",
        "TIME-OBS": "12:00:00",
        "NCCDS": "1",
        "CCDSUM": "2 2",
    })
    t = tel.get_time(hdr)
    assert t > 0 and np.isfinite(t)
    n = tel.get_number(hdr)
    assert isinstance(n, (int, np.integer))


def test_get_instrument_name():
    """get_instrument_name returns lowercase name."""
    tel = GMOS()
    hdr = fits.Header()
    assert tel.get_instrument_name(hdr) == "gmos"


def test_get_catalog():
    """GMOS get_catalog returns SkyMapper for dec < -30, PS1 otherwise."""
    tel = GMOS()
    hdr_n = fits.Header({"RA": 180.0, "DEC": 0.0})
    hdr_s = fits.Header({"RA": 180.0, "DEC": -35.0})
    assert tel.get_catalog(hdr_n) == "PS1"
    assert tel.get_catalog(hdr_s) == "SkyMapper"


def test_format_datasec():
    """format_datasec converts section string with binning."""
    tel = Instrument()
    out = tel.format_datasec("[1055:3024,217:3911]", binning=2)
    assert "[527:1512,108:1955]" == out or out == "[527:1512,109:1956]"


def test_raw_format():
    """Base raw_format and GMOS raw_format."""
    base = Instrument()
    assert base.raw_format(True) == "sci_img_*.fits"
    assert base.raw_format(False) == "sci_img*[!proc].fits"
    gmos = GMOS()
    assert gmos.raw_format("dragons") == "*.fits"
    assert gmos.raw_format("other") == "*.fits.bz2"


def test_get_stk_name_get_sci_name_get_bkg_name():
    """Naming helpers return paths under red_path."""
    tel = GMOS()
    hdr = fits.Header({
        "OBJECT": "Target",
        "FILTER2": "r",
        "DATE-OBS": "2024-06-18",
        "TIME-OBS": "12:00:00",
        "NCCDS": "1",
        "CCDSUM": "2 2",
    })
    red_path = "/data/red"
    stk = tel.get_stk_name(hdr, red_path)
    assert "Target" in stk and "r" in stk and "stk.fits" in stk and red_path in stk
    sci = tel.get_sci_name(hdr, red_path)
    assert "Target" in sci and ".fits" in sci and red_path in sci
    bkg = tel.get_bkg_name(hdr, red_path)
    assert "_bkg.fits" in bkg and red_path in bkg


def test_get_mbias_name_get_mdark_name_get_mflat_name_get_msky_name():
    """Calibration filename helpers."""
    tel = GMOS()
    paths = {"cal": "/data/red/cals"}
    assert "mbias_1_22.fits" in tel.get_mbias_name(paths, "1", "22")
    assert "mdark_1_22.fits" in tel.get_mdark_name(paths, "1", "22")
    assert "mflat_r_1_22.fits" in tel.get_mflat_name(paths, "r", "1", "22")
    assert "msky_r_1_22.fits" in tel.get_msky_name(paths, "r", "1", "22")


def test_base_instrument_get_ampl_get_binning_missing_keyword():
    """Base Instrument get_ampl/get_binning return default when keyword missing."""
    base = Instrument()
    hdr = fits.Header()
    assert base.get_ampl(hdr) == "2"
    assert base.get_binning(hdr) == "CCDSUM"


def test_sanitize_calibration_header():
    """_sanitize_calibration_header removes WCS/coord keywords so saved cals don't trigger InvalidTransformError."""
    h = fits.Header()
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["CRVAL1"] = 0.0
    h["CRVAL2"] = -100.0
    h["RA"] = "00:00:00"
    h["DEC"] = "-100:00:00"
    h["PV1_1"] = 1.0
    h["EXPTIME"] = 60.0
    h["VER"] = "1.0"
    _sanitize_calibration_header(h)
    assert "CTYPE1" not in h
    assert "CRVAL2" not in h
    assert "RA" not in h
    assert "DEC" not in h
    assert "PV1_1" not in h
    assert h["EXPTIME"] == 60.0
    assert h["VER"] == "1.0"


def test_read_calibration_ccd(tmp_path):
    """_read_calibration_ccd loads calibration FITS without parsing WCS (avoids ill-conditioned header errors)."""
    path = tmp_path / "cal.fits"
    hdu = fits.PrimaryHDU(np.zeros((10, 10), dtype=np.float32))
    hdu.header["CRVAL2"] = -100.0  # invalid; would raise if WCS were parsed
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.writeto(path, overwrite=True)
    ccd = _read_calibration_ccd(str(path), u.electron, hdu_index=0)
    assert ccd.wcs is None
    assert ccd.unit == u.electron
    assert ccd.data.shape == (10, 10)


def test_sky_subtraction_units():
    """Scaled sky (normalized * med * electron) has same unit as science so subtract is valid."""
    # Normalized master sky (dimensionless) * (med * u.electron) -> electron; then science - sky is valid
    sky = CCDData(np.ones((5, 5)), unit=u.dimensionless_unscaled)
    frame = CCDData(np.ones((5, 5)) * 100.0, unit=u.electron)
    med = 50.0
    science_unit = frame.unit if frame.unit is not None else u.electron
    frame_sky = sky.multiply(med * science_unit, propagate_uncertainties=True, handle_meta="first_found")
    result = frame.subtract(frame_sky, propagate_uncertainties=True, handle_meta="first_found")
    assert result.unit == u.electron
    np.testing.assert_allclose(result.data, 50.0)


def test_get_staticmask_filename(tmp_path):
    """get_staticmask_filename returns [path] when mask exists, [None] otherwise."""
    tel = GMOS()
    hdr = fits.Header({"CCDSUM": "2 2"})
    # Path is paths['code']/../data/staticmasks/{instname}.{binn}.staticmask.fits.fz
    code_dir = tmp_path / "code"
    code_dir.mkdir()
    data_dir = tmp_path / "data" / "staticmasks"
    data_dir.mkdir(parents=True)
    paths = {"code": str(code_dir)}
    # When file does not exist
    out = tel.get_staticmask_filename(hdr, paths)
    assert out == [None]
    # When file exists
    mask_file = data_dir / "gmos.22.staticmask.fits.fz"
    mask_file.touch()
    out = tel.get_staticmask_filename(hdr, paths)
    assert out[0] is not None
    assert "gmos.22.staticmask.fits.fz" in out[0]


def test_load_bias_load_dark_load_flat_load_sky(tmp_path):
    """load_bias, load_dark, load_flat, load_sky load from temp FITS via _read_calibration_ccd."""
    tel = GMOS()
    cal_dir = tmp_path / "cals"
    cal_dir.mkdir()
    paths = {"cal": str(cal_dir)}

    # Bias: mbias_1_22.fits
    bias_path = cal_dir / "mbias_1_22.fits"
    fits.PrimaryHDU(np.zeros((10, 10), dtype=np.float32)).writeto(bias_path, overwrite=True)
    mbias = tel.load_bias(paths, "1", "22")
    assert mbias.unit == u.electron
    assert mbias.data.shape == (10, 10)

    # Dark: mdark_1_22.fits
    dark_path = cal_dir / "mdark_1_22.fits"
    fits.PrimaryHDU(np.zeros((10, 10), dtype=np.float32)).writeto(dark_path, overwrite=True)
    mdark = tel.load_dark(paths, "1", "22")
    assert mdark.unit == u.electron

    # Flat: mflat_r_1_22.fits
    flat_path = cal_dir / "mflat_r_1_22.fits"
    fits.PrimaryHDU(np.ones((10, 10), dtype=np.float32)).writeto(flat_path, overwrite=True)
    mflat = tel.load_flat(paths, "r", "1", "22")
    assert mflat.unit == u.dimensionless_unscaled

    # Sky: msky_r_1_22.fits
    sky_path = cal_dir / "msky_r_1_22.fits"
    fits.PrimaryHDU(np.ones((10, 10), dtype=np.float32)).writeto(sky_path, overwrite=True)
    msky = tel.load_sky(paths, "r", "1", "22")
    assert msky.unit == u.dimensionless_unscaled


def test_load_bias_raises_when_missing(tmp_path):
    """load_bias raises when no bias file exists at expected path."""
    tel = GMOS()
    paths = {"cal": str(tmp_path)}
    with pytest.raises(Exception, match="Could not find bias"):
        tel.load_bias(paths, "1", "22")


def test_expand_mask():
    """expand_mask combines NaN, inf, zero pixels with optional input mask."""
    tel = Instrument()
    data = np.ones((8, 8), dtype=float)
    data[0, 0] = np.nan
    data[1, 1] = np.inf
    data[2, 2] = 0.0
    ccd = CCDData(data, unit=u.electron)
    out = tel.expand_mask(ccd, input_mask=None)
    assert out[0, 0]
    assert out[1, 1]
    assert out[2, 2]
    assert not out[3, 3]
    # With input mask
    extra = np.zeros((8, 8), dtype=bool)
    extra[4, 4] = True
    out2 = tel.expand_mask(ccd, input_mask=extra)
    assert out2[4, 4]


def test_base_get_catalog():
    """Base Instrument get_catalog returns catalog_zp."""
    base = Instrument()
    hdr = fits.Header()
    assert base.get_catalog(hdr) == "PS1"


def test_format_datasec_binning_one():
    """format_datasec with binning=1 preserves integer bounds."""
    base = Instrument()
    out = base.format_datasec("[100:200,50:150]", binning=1)
    assert "[100:200,50:150]" == out


def test_create_sky_masks_bright_sources(tmp_path):
    """create_sky uses iterative sigma clipping and masks bright sources so sky estimate is robust."""
    from astropy.stats import sigma_clipped_stats

    tel = GMOS()
    cal_dir = tmp_path / "cal"
    cal_dir.mkdir()
    sky_dir = tmp_path / "sky"
    sky_dir.mkdir()
    paths = {"cal": str(cal_dir)}

    # Two small sky frames: constant background + one frame with a bright source
    sky_val = 1000.0
    shape = (12, 12)
    for i in range(2):
        data = np.full(shape, sky_val, dtype=np.float32)
        if i == 1:
            data[5, 5] = 8000.0  # bright source
            data[6, 5] = 6000.0
        ccd = CCDData(data, unit=u.electron)
        path = sky_dir / f"sky{i}.fits"
        ccd.write(str(path), overwrite=True)

    sky_list = [str(sky_dir / "sky0.fits"), str(sky_dir / "sky1.fits")]
    tel.create_sky(sky_list, "r", "1", "22", paths, log=None,
        sky_sigma_upper=3.0, sky_sigma_lower=4.0, sky_maxiters=5,
        sky_n_sigma_high=3.0, sky_n_sigma_low=4.0,
        msky_n_sigma_high=5.0, msky_n_sigma_low=5.0, msky_maxiters=5)

    msky_path = cal_dir / "msky_r_1_22.fits"
    assert msky_path.exists()
    with fits.open(msky_path) as hdu:
        msky_data = hdu[0].data
    # Normalized sky: masked (bright) pixels set to 1.0, rest ~1.0
    assert msky_data.shape == shape
    # At (5,5) and (6,5) we had bright flux; after masking they should be 1.0
    assert msky_data[5, 5] == 1.0
    assert msky_data[6, 5] == 1.0
    # Median of combined sky should be close to 1.0 (normalized)
    med = np.nanmedian(msky_data)
    assert 0.95 < med < 1.05


def test_create_sky_accepts_default_and_custom_sigma_params(tmp_path):
    """create_sky runs with default params and with custom sigma/maxiters."""
    tel = GMOS()
    cal_dir = tmp_path / "cal"
    cal_dir.mkdir()
    sky_dir = tmp_path / "sky"
    sky_dir.mkdir()
    paths = {"cal": str(cal_dir)}
    data = np.full((8, 8), 500.0, dtype=np.float32)
    ccd = CCDData(data, unit=u.electron)
    sky_path = sky_dir / "sky.fits"
    ccd.write(str(sky_path), overwrite=True)
    sky_list = [str(sky_path)]

    # Default params
    tel.create_sky(sky_list, "r", "1", "22", paths, log=None)
    assert (cal_dir / "msky_r_1_22.fits").exists()

    # Custom params (tight sigma)
    (cal_dir / "msky_r_1_22.fits").unlink()
    tel.create_sky(sky_list, "r", "1", "22", paths, log=None,
        sky_n_sigma_high=2.0, msky_n_sigma_high=4.0)
    assert (cal_dir / "msky_r_1_22.fits").exists()
