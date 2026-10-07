"""Unit tests for instrument base class methods and instrument_getter."""
import os

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import CCDData
from astropy import units as u

from potpyri.instruments import (
    UnknownInstrumentError,
    __all__ as INSTRUMENT_NAMES,
    instrument_getter,
    resolve_instrument_name,
    supported_instruments,
)
from potpyri.instruments.GMOS import GMOS
from potpyri.instruments.instrument import (
    Instrument,
    _sanitize_calibration_header,
    _read_calibration_ccd,
    fix_deprecated_wcs_header_cards,
)
from potpyri.utils import options
from potpyri.utils import logger


def test_fix_deprecated_wcs_header_cards_radecsys_and_mjd():
    """fix_deprecated_wcs_header_cards migrates RADECSYS and DATE/TIME from MJD-OBS."""
    h = fits.Header()
    h['RADECSYS'] = 'FK5 '
    h['MJD-OBS'] = 60480.0
    fix_deprecated_wcs_header_cards(h)
    assert 'RADECSYS' not in h
    assert h['RADESYS'] == 'FK5'
    assert 'DATE-OBS' in h and 'TIME-OBS' in h
    assert 'MJD-OBS' not in h


def test_ensure_exptime_keyword_fills_missing_for_mosfire():
    """MOSFIRE SCI headers lack ELAPTIME; ensure_exptime_keyword fills from get_exptime."""
    tel = instrument_getter('MOSFIRE')
    h = fits.Header()
    h['TRUITIME'] = 1.45
    h['COADDONE'] = 8
    assert 'ELAPTIME' not in h
    tel.ensure_exptime_keyword(h)
    assert h['ELAPTIME'] == pytest.approx(1.45 * 8)
    # Existing value is left alone
    h['ELAPTIME'] = 99.0
    tel.ensure_exptime_keyword(h)
    assert h['ELAPTIME'] == 99.0


def test_resolve_instrument_name_aliases():
    """resolve_instrument_name maps aliases and canonical names."""
    assert resolve_instrument_name('gmos') == 'GMOS'
    assert resolve_instrument_name('BINO') == 'BINOSPEC'
    assert resolve_instrument_name('MMIR') == 'MMIRS'
    assert resolve_instrument_name('BINOSPEC') == 'BINOSPEC'


def test_supported_instruments_matches_all():
    assert set(supported_instruments()) == set(INSTRUMENT_NAMES)


def test_instrument_getter_unsupported_raises():
    """instrument_getter raises UnknownInstrumentError when log is None."""
    with pytest.raises(UnknownInstrumentError, match="not supported"):
        instrument_getter("UNKNOWN_INSTRUMENT", log=None)


def test_instrument_getter_invalid_name_type():
    with pytest.raises(TypeError):
        instrument_getter("", log=None)


def test_instrument_getter_unsupported_with_log(tmp_path):
    """instrument_getter with log calls log.error and returns None for unsupported name."""
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


def test_mosfire_get_number_fallbacks():
    """MOSFIRE get_number prefers FRAMENO, then FRAMENUM, then DATAFILE digits."""
    tel = instrument_getter('MOSFIRE')

    assert tel.get_number(fits.Header({'FRAMENO': 42})) == '00042'
    assert tel.get_number(fits.Header({'FRAMENUM': 180})) == '00180'
    # FRAMENO wins over FRAMENUM when both exist
    assert tel.get_number(fits.Header({'FRAMENO': 7, 'FRAMENUM': 180})) == '00007'
    assert tel.get_number(fits.Header({'DATAFILE': 'm260724_0180'})) == '00180'
    assert tel.get_number(
        fits.Header({'ORGFILE': '/data/raw/m260724_0233.fits'})
    ) == '00233'
    with pytest.raises(KeyError):
        tel.get_number(fits.Header({'OBJECT': 'FRB'}))


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
    """raw_format and --proc fits override across instruments."""
    base = Instrument()
    assert base.raw_format(True) == "sci_img_*.fits"
    assert base.raw_format(False) == "sci_img*[!proc].fits"
    assert base.raw_format("fits") == "*.fits"
    gmos = GMOS()
    assert gmos.raw_format("dragons") == "*.fits"
    assert gmos.raw_format("other") == "*.fits.bz2"
    assert gmos.raw_format("fits") == "*.fits"
    f2 = instrument_getter("F2")
    assert f2.raw_format(True) == "*.fits.bz2"
    assert f2.raw_format("fits") == "*.fits"
    assert f2.raw_format(".fits") == "*.fits"
    assert f2.raw_format("uncompressed") == "*.fits"
    mosfire = instrument_getter("MOSFIRE")
    assert mosfire.raw_format(None) == "*.fits.gz"
    assert mosfire.raw_format("fits") == "*.fits"
    lris = instrument_getter("LRIS")
    assert lris.raw_format("archive") == "*.fits*"
    assert lris.raw_format("fits") == "*.fits"


def test_discover_raw_files_and_message(tmp_path):
    """discover_raw_files reports search paths and glob matches."""
    raw_dir = tmp_path / 'raw'
    data_dir = tmp_path / 'data'
    bad_dir = tmp_path / 'bad'
    for d in (raw_dir, data_dir, bad_dir):
        d.mkdir()
    (raw_dir / 'science.fits').write_bytes(b'\x00')
    (data_dir / 'other.fits.bz2').write_bytes(b'\x00')

    tel = instrument_getter("F2")
    paths = {'raw': str(raw_dir), 'data': str(data_dir), 'bad': str(bad_dir)}

    info = tel.discover_raw_files(paths, proc='fits')
    assert info['pattern'] == '*.fits'
    assert len(info['files']) == 1
    assert info['files'][0].endswith('science.fits')

    msg = tel.format_raw_discovery_message(paths, proc='fits')
    assert 'Glob pattern' in msg
    assert str(raw_dir) in msg
    assert 'science.fits' not in msg or '1 file(s)' in msg
    assert 'Search directories' in msg

    info_bz2 = tel.discover_raw_files(paths, proc=True)
    assert info_bz2['pattern'] == '*.fits.bz2'
    assert len(info_bz2['files']) == 1
    assert info_bz2['files'][0].endswith('other.fits.bz2')


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


def test_load_flat_falls_back_to_caldb(tmp_path):
    """load_flat uses packaged caldb when the night cal is missing."""
    tel = GMOS()
    night = tmp_path / "cals"
    caldb = tmp_path / "caldb"
    night.mkdir()
    caldb.mkdir()
    fits.PrimaryHDU(np.ones((8, 8), dtype=np.float32)).writeto(
        caldb / "mflat_r_1_22.fits", overwrite=True,
    )
    paths = {"cal": str(night), "caldb": str(caldb)}
    mflat = tel.load_flat(paths, "r", "1", "22")
    assert mflat.data.shape == (8, 8)


def test_do_bias_no_frames_does_not_exit(tmp_path, capsys):
    """A night with no bias frames must not abort the whole pipeline."""
    from potpyri.primitives import calibration

    class _Tel:
        bias = True
        filetype_keywords = {'BIAS': 'BIAS', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

        def get_mbias_name(self, paths, amp, binn):
            return os.path.join(paths['cal'], f'mbias_{amp}_{binn}.fits')

        def _find_cal_file(self, paths, filename):
            return (None, None)

    table = Table({
        'Type': ['SCIENCE'],
        'CalType': ['R_1R_11'],
        'Amp': ['1R'],
        'Binning': ['11'],
        'File': ['d'],
    })
    calibration.do_bias(table, _Tel(), {'cal': str(tmp_path)}, nmin_images=3, log=None)
    out = capsys.readouterr().out
    assert 'No usable master bias' in out
    assert 'No master bias for science setup amp=1R' in out


def test_do_bias_reports_missing_science_setup(tmp_path, capsys):
    """do_bias warns when science amp/bin has no master and does not abort."""
    from potpyri.primitives import calibration

    class _Tel:
        bias = True
        filetype_keywords = {'BIAS': 'BIAS', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

        def get_mbias_name(self, paths, amp, binn):
            return os.path.join(paths['cal'], f'mbias_{amp}_{binn}.fits')

        def create_bias(self, files, amp, binn, paths, log=None):
            with open(self.get_mbias_name(paths, amp, binn), 'w') as fh:
                fh.write('bias')

        def _find_cal_file(self, paths, filename):
            path = os.path.join(paths['cal'], filename)
            return (path, 0) if os.path.exists(path) else (None, None)

    cal = tmp_path / 'cals'
    cal.mkdir()
    table = Table({
        'Type': ['BIAS', 'BIAS', 'BIAS', 'SCIENCE'],
        'CalType': ['4B_11', '4B_11', '4B_11', 'R_1R_11'],
        'Amp': ['4B', '4B', '4B', '1R'],
        'Binning': ['11', '11', '11', '11'],
        'File': ['a', 'b', 'c', 'd'],
    })
    calibration.do_bias(table, _Tel(), {'cal': str(cal)}, nmin_images=3, log=None)
    out = capsys.readouterr().out
    assert 'No master bias for science setup amp=1R' in out
    assert (cal / 'mbias_4B_11.fits').exists()


def test_do_dark_no_frames_does_not_exit(tmp_path, capsys):
    """A required dark with no frames must not abort the whole pipeline."""
    from potpyri.primitives import calibration

    class _Tel:
        dark = True
        filetype_keywords = {'DARK': 'DARK', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

        def get_mdark_name(self, paths, amp, binn):
            return os.path.join(paths['cal'], f'mdark_{amp}_{binn}.fits')

        def _find_cal_file(self, paths, filename):
            return (None, None)

    table = Table({
        'Type': ['SCIENCE'],
        'CalType': ['R_1R_11'],
        'Amp': ['1R'],
        'Binning': ['11'],
        'File': ['d'],
    })
    calibration.do_dark(table, _Tel(), {'cal': str(tmp_path)}, nmin_images=3, log=None)
    out = capsys.readouterr().out
    assert 'No usable master dark' in out
    assert 'No master dark for science setup amp=1R' in out


def test_do_dark_reports_missing_science_setup(tmp_path, capsys):
    """do_dark warns when science amp/bin has no master and does not abort."""
    from potpyri.primitives import calibration

    class _Tel:
        dark = True
        bias = False
        filetype_keywords = {'DARK': 'DARK', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

        def get_mdark_name(self, paths, amp, binn):
            return os.path.join(paths['cal'], f'mdark_{amp}_{binn}.fits')

        def create_dark(self, files, amp, binn, paths, mbias=None, log=None):
            with open(self.get_mdark_name(paths, amp, binn), 'w') as fh:
                fh.write('dark')

        def _find_cal_file(self, paths, filename):
            path = os.path.join(paths['cal'], filename)
            return (path, 0) if os.path.exists(path) else (None, None)

    cal = tmp_path / 'cals'
    cal.mkdir()
    table = Table({
        'Type': ['DARK', 'DARK', 'DARK', 'SCIENCE'],
        'CalType': ['4B_11', '4B_11', '4B_11', 'R_1R_11'],
        'Amp': ['4B', '4B', '4B', '1R'],
        'Binning': ['11', '11', '11', '11'],
        'Exp': [100.0, 100.0, 100.0, 30.0],
        'File': ['a', 'b', 'c', 'd'],
    })
    calibration.do_dark(table, _Tel(), {'cal': str(cal)}, nmin_images=3, log=None)
    out = capsys.readouterr().out
    assert 'No master dark for science setup amp=1R' in out
    assert (cal / 'mdark_4B_11.fits').exists()


def test_do_dark_skipped_when_not_required(tmp_path, capsys):
    """Instruments that do not require darks emit no missing-dark errors."""
    from potpyri.primitives import calibration

    class _Tel:
        dark = False
        filetype_keywords = {'DARK': 'DARK', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

    table = Table({
        'Type': ['SCIENCE'],
        'CalType': ['R_1R_11'],
        'Amp': ['1R'],
        'Binning': ['11'],
        'File': ['d'],
    })
    calibration.do_dark(table, _Tel(), {'cal': str(tmp_path)}, nmin_images=3, log=None)
    assert capsys.readouterr().out == ''


def test_do_flat_reports_missing_science_setup(tmp_path, capsys):
    """do_flat warns when science filter/amp/bin has no master and does not abort."""
    from potpyri.primitives import calibration

    class _Tel:
        flat = True
        bias = False
        dark = False
        filetype_keywords = {'FLAT': 'FLAT', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

        def get_mflat_name(self, paths, fil, amp, binn):
            return os.path.join(paths['cal'], f'mflat_{fil}_{amp}_{binn}.fits')

        def create_flat(self, files, fil, amp, binn, paths, mbias=None,
                mdark=None, is_science=False, log=None):
            with open(self.get_mflat_name(paths, fil, amp, binn), 'w') as fh:
                fh.write('flat')

        def _find_cal_file(self, paths, filename):
            path = os.path.join(paths['cal'], filename)
            return (path, 0) if os.path.exists(path) else (None, None)

    cal = tmp_path / 'cals'
    cal.mkdir()
    table = Table({
        'Type': ['FLAT', 'FLAT', 'FLAT', 'SCIENCE'],
        'CalType': ['g_4B_11', 'g_4B_11', 'g_4B_11', 'R_1R_11'],
        'Filter': ['g', 'g', 'g', 'R'],
        'Amp': ['4B', '4B', '4B', '1R'],
        'Binning': ['11', '11', '11', '11'],
        'File': ['a', 'b', 'c', 'd'],
    })
    calibration.do_flat(table, _Tel(), {'cal': str(cal)}, nmin_images=3, log=None)
    out = capsys.readouterr().out
    assert 'No master flat for science setup filter=R, amp=1R' in out
    assert (cal / 'mflat_g_4B_11.fits').exists()


def test_do_flat_caldb_satisfies_science_setup(tmp_path, capsys):
    """A packaged caldb flat counts as present; no missing-flat error."""
    from potpyri.primitives import calibration

    cal = tmp_path / 'cals'
    caldb = tmp_path / 'caldb'
    cal.mkdir()
    caldb.mkdir()
    (caldb / 'mflat_R_1R_11.fits').write_text('flat')

    class _Tel:
        flat = True
        filetype_keywords = {'FLAT': 'FLAT', 'SCIENCE': 'SCIENCE'}

        def match_type_keywords(self, kwds, table):
            return table['Type'] == kwds

        def get_mflat_name(self, paths, fil, amp, binn):
            return os.path.join(paths['cal'], f'mflat_{fil}_{amp}_{binn}.fits')

        def _find_cal_file(self, paths, filename):
            for key in ('cal', 'caldb'):
                path = os.path.join(paths[key], filename)
                if os.path.exists(path):
                    return (path, 0)
            return (None, None)

    table = Table({
        'Type': ['SCIENCE'],
        'CalType': ['R_1R_11'],
        'Filter': ['R'],
        'Amp': ['1R'],
        'Binning': ['11'],
        'File': ['d'],
    })
    calibration.do_flat(
        table, _Tel(), {'cal': str(cal), 'caldb': str(caldb)},
        nmin_images=3, log=None,
    )
    out = capsys.readouterr().out
    assert 'No night master flat was built' in out
    assert 'No master flat for science setup' not in out


def test_do_flat_empty_does_not_change_cal_path(capsys):
    """No flat frames must not retarget paths['cal'] to the packaged caldb."""
    from potpyri.primitives import calibration

    tel = instrument_getter("LRIS")
    paths = {"cal": "/night/cals", "caldb": "/pkg/caldb"}
    table = Table({
        "Type": ["SCIENCE"],
        "CalType": ["x"],
        "Filter": ["R"],
        "Amp": ["1R"],
        "Binning": ["11"],
    })
    calibration.do_flat(table, tel, paths)
    assert paths["cal"] == "/night/cals"
    out = capsys.readouterr().out
    assert 'No master flat for science setup filter=R, amp=1R' in out


def test_mask_flat_sources_excludes_stars_from_flat_level():
    """Stars above the twilight-flat level are masked and do not change the median."""
    tel = instrument_getter("LRIS")
    data = np.full((64, 64), 10000.0)
    data[28:36, 28:36] = 40000.0
    unmasked_mean = float(np.mean(data))
    frame = CCDData(data.copy(), unit=u.electron)
    tel.mask_flat_sources(frame)
    assert np.any(np.isnan(frame.data))
    assert np.all(np.isnan(frame.data[30:34, 30:34]))
    after = float(np.nanmedian(frame.data))
    assert abs(after - 10000.0) < 1.0
    assert after < unmasked_mean


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


def test_set_zeropoint_catalog_overrides_gmos_latitude_switch():
    """set_zeropoint_catalog forces DECALS even for southern GMOS fields."""
    tel = instrument_getter("GMOS")
    tel.set_zeropoint_catalog("decals")
    assert tel.catalog_zp == "DECALS"
    hdr_s = fits.Header()
    hdr_s["RA"] = 0.0
    hdr_s["DEC"] = -40.0
    assert tel.get_catalog(hdr_s) == "DECALS"


def test_format_datasec_binning_one():
    """format_datasec with binning=1 preserves integer bounds."""
    base = Instrument()
    out = base.format_datasec("[100:200,50:150]", binning=1)
    assert "[100:200,50:150]" == out


def test_create_sky_masks_bright_sources(tmp_path):
    """create_sky uses iterative sigma clipping and masks bright sources so sky estimate is robust."""
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
