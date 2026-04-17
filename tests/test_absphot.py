"""Tests for zeropoint (find_zeropoint, Vizier/PS1) and ZeropointFitter methods."""
import os

import numpy as np
import pytest
import requests
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.coordinates import SkyCoord

from potpyri.utils import options, logger
from potpyri.primitives import zeropoint
from potpyri.instruments import instrument_getter

from tests.utils import download_gdrive_file


@pytest.mark.integration
def test_absphot(tmp_path):
    """Run find_zeropoint on LRIS stack and check ZPTMAG, MAGSYS, ZPTCAT, ZPTMUCER in header."""
    instrument = 'LRIS'
    file_list_name = 'files.txt'

    # Raw science file (download to tmp_path so log is writable)
    file_path = download_gdrive_file('LRISr/red/R155_host.R.ut240629.1R.11.stk.fits.fz', output_dir=str(tmp_path), use_cached=True)
    cat_path = download_gdrive_file('LRISr/red/test_catalog.dat', output_dir=str(tmp_path), use_cached=True)
    input_catalog = ascii.read(cat_path)

    data_path, basefile = os.path.split(file_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    hdu = fits.open(file_path)
    hdu.close()

    try:
        zeropoint.find_zeropoint(file_path, tel, log=log)
    except (requests.exceptions.ConnectionError, OSError) as e:
        pytest.skip(f"VizieR unreachable (zeropoint needs PS1): {e}")
    finally:
        log.close()

    hdu = fits.open(file_path)
    header = hdu['PRIMARY'].header
    hdu.close()

    assert header['ZPTCAT']=='PS1'
    assert header['FILTER']=='R'
    assert header['MAGSYS'] == zeropoint.DEFAULT_MAGSYS
    assert np.abs(27.576871-header['ZPTMAG'])<0.01
    assert header['ZPTMUCER']<0.01


def test_get_zeropoint():
    """get_zeropoint returns (zpt, zpterr) consistent with mag = zpt - 2.5*log10(flux)."""
    cal = zeropoint.ZeropointFitter()
    flux = np.array([1000.0, 2000.0, 500.0])
    fluxerr = np.array([10.0, 14.0, 7.0])
    mag = np.array([25.0, 24.25, 25.75])
    magerr = np.array([0.01, 0.01, 0.01])
    zpt, zpterr = cal.get_zeropoint(flux, fluxerr, mag, magerr)
    assert np.isfinite(zpt) and np.isfinite(zpterr)
    mag_derived = zpt - 2.5 * np.log10(flux)
    np.testing.assert_allclose(mag, mag_derived, atol=0.1)


def test_zpt_iteration():
    """zpt_iteration returns (zpt, zpterr, master_mask) with mask length = input length."""
    cal = zeropoint.ZeropointFitter(iterations=2, sigma=3.0)
    np.random.seed(42)
    n = 20
    flux = np.random.uniform(500, 3000, n).astype(float)
    fluxerr = flux * 0.01
    zpt_true = 26.0
    mag = zpt_true - 2.5 * np.log10(flux) + np.random.normal(0, 0.02, n)
    magerr = np.full(n, 0.02)
    zpt, zpterr, master_mask = cal.zpt_iteration(flux, fluxerr, mag, magerr, log=None)
    assert np.isfinite(zpt) and np.isfinite(zpterr)
    assert len(master_mask) == n
    assert np.sum(master_mask) <= n


def test_convert_filter_name():
    """convert_filter_name maps instrument filter names to catalog names."""
    cal = zeropoint.ZeropointFitter()
    assert cal.convert_filter_name('gG0301') == 'g'
    assert cal.convert_filter_name('rG0303') == 'r'
    assert cal.convert_filter_name('RG850') == 'z'
    assert cal.convert_filter_name('Y') == 'J'
    assert cal.convert_filter_name('V') == 'g'
    assert cal.convert_filter_name('I') == 'i'
    # K, Ks, Kspec all map to 2MASS K-band
    assert cal.convert_filter_name('K') == 'K'
    assert cal.convert_filter_name('Ks') == 'K'
    assert cal.convert_filter_name('Kspec') == 'K'


def test_get_minmag():
    """get_minmag returns bright limit by filter."""
    cal = zeropoint.ZeropointFitter()
    assert cal.get_minmag('J') == 15.5
    assert cal.get_minmag('K') == 13.0
    assert cal.get_minmag('Y') == 15.0
    assert cal.get_minmag('g') == 16.0


def test_Y_band():
    """GF18 eq. (6): VISTA Y_V = J + C_Y*(J-Ks)_2MASS with C_Y=0.46 (Vega-like Y)."""
    cal = zeropoint.ZeropointFitter()
    J, J_err = 15.0, 0.02
    Ks, Ks_err = 14.0, 0.02
    Yv, Y_err = cal.y_band_from_jk(J, J_err, Ks, Ks_err)
    expected_Yv = J + zeropoint.GF18_TWOMASS_TO_VISTA_Y_SLOPE * (J - Ks)
    np.testing.assert_allclose(Yv, expected_Yv)
    assert np.isfinite(Y_err) and Y_err > 0
    Y2, Y2_err = cal.Y_band(J, J_err, Ks, Ks_err)
    np.testing.assert_allclose(Y2, Yv)
    np.testing.assert_allclose(Y2_err, Y_err)


def test_twomass_y_includes_gf18_ab_offset():
    """GF18 Appendix D (D3): Y_AB = Y_V + 0.600 for twomass_y_mag_to_ab."""
    cal = zeropoint.ZeropointFitter()
    J, J_err = 15.0, 0.02
    Ks, Ks_err = 14.0, 0.02
    Yv, _ = cal.twomass_j_ks_to_vista_y_vega(J, J_err, Ks, Ks_err)
    Yab, _ = cal.twomass_y_mag_to_ab(J, J_err, Ks, Ks_err)
    np.testing.assert_allclose(Yab, Yv + zeropoint.GF18_VISTA_Y_VEGA_TO_AB)


def test_twomass_vega_to_ab():
    """twomass_vega_to_ab applies standard 2MASS Vega -> AB offsets."""
    m = 10.0
    assert zeropoint.ZeropointFitter.twomass_vega_to_ab(m, 'J') == m + (4.56 - 3.65)
    assert zeropoint.ZeropointFitter.twomass_vega_to_ab(m, 'H') == m + (4.71 - 3.32)
    assert zeropoint.ZeropointFitter.twomass_vega_to_ab(m, 'K') == m + (5.14 - 3.29)
    assert zeropoint.ZeropointFitter.twomass_vega_to_ab(m, 'Ks') == m + (5.14 - 3.29)


def test_magsys_get_set():
    """get_magsys / set_magsys and DEFAULT_MAGSYS."""
    cal = zeropoint.ZeropointFitter()
    assert cal.get_magsys() == zeropoint.DEFAULT_MAGSYS == 'ABMAG'
    cal.set_magsys('VEGAMAG')
    assert cal.get_magsys() == 'VEGAMAG'
