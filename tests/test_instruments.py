"""Unit tests for instrument base class methods and instrument_getter."""
import os

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from potpyri.instruments import instrument_getter
from potpyri.instruments.GMOS import GMOS
from potpyri.instruments.instrument import Instrument
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
