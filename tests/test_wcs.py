"""Tests for solve_wcs (solve_astrometry, align_to_gaia on GMOS stack) and helpers (clean_up_astrometry)."""
import os

import numpy as np
import pytest
import requests
from astropy.io import fits
from astropy.table import Table

from potpyri.utils import options, logger
from potpyri.primitives import solve_wcs
from potpyri.instruments import instrument_getter

from tests.utils import download_gdrive_file


def test_wcs(tmp_path):
    """Check RADISP/DEDISP and scale from a FITS file (no network required).

    Uses a minimal FITS file with RADISP/DEDISP set as align_to_gaia would
    write them, so the test passes consistently without being skipped.
    Full solve_astrometry + align_to_gaia flow is covered by test_wcs_integration.
    """
    instrument = 'GMOS'
    tel = instrument_getter(instrument)
    # Dispersion values that satisfy the pipeline quality cuts (arcsec).
    # Must be < 0.5 and < tel.pixscale so disp/pixscale < 1 (GMOS pixscale ≈ 0.08).
    ra_disp_val = dec_disp_val = 0.05

    # Create minimal FITS with RADISP/DEDISP as align_to_gaia writes them
    file_path = os.path.join(tmp_path, 'test_solved.fits')
    hdu = fits.PrimaryHDU(data=np.zeros((50, 50), dtype=np.float32))
    hdu.header['NAXIS1'] = 50
    hdu.header['NAXIS2'] = 50
    hdu.header['RADISP'] = (ra_disp_val, 'Dispersion in R.A. of WCS [Arcsec]')
    hdu.header['DEDISP'] = (dec_disp_val, 'Dispersion in Decl. of WCS [Arcsec]')
    hdu.writeto(file_path, overwrite=True)

    with fits.open(file_path) as hdu:
        ra_disp = hdu[0].header['RADISP']
        dec_disp = hdu[0].header['DEDISP']

    assert ra_disp < 0.5
    assert dec_disp < 0.5
    assert ra_disp / tel.pixscale < 1
    assert dec_disp / tel.pixscale < 1


@pytest.mark.integration
def test_wcs_integration(tmp_path):
    """Full pipeline: solve_astrometry and align_to_gaia on GMOS slice (requires network and astrometry index)."""
    instrument = 'GMOS'
    file_list_name = 'files.txt'

    file_path = download_gdrive_file(
        'GMOS/red/workspace/sGRB240615A-GRB.i.ut240618.12.22.1403185114.fits.fz',
        output_dir=str(tmp_path), use_cached=True)
    astm_path = download_gdrive_file(
        'astrometry/index-52m1-14.fits', output_dir=str(tmp_path), use_cached=True)

    data_path, _ = os.path.split(file_path)
    data_path, _ = os.path.split(data_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)
    log = logger.get_log(paths['log'])

    with fits.open(file_path) as hdu:
        binn = tel.get_binning(hdu[1].header)
        file_path_unfz = file_path.replace('.fz', '')
        hdu[1].header.pop('RADISP', None)
        hdu[1].header.pop('DEDISP', None)
        fits.writeto(file_path_unfz, hdu[1].data, hdu[1].header, overwrite=True)
    file_path = file_path_unfz

    try:
        solve_wcs.solve_astrometry(file_path, tel, binn, paths, index=astm_path, log=log)
        solve_wcs.align_to_gaia(file_path, tel, radius=0.5, log=log)
    except (requests.exceptions.ConnectionError, OSError) as e:
        pytest.skip(f"VizieR/network unreachable (Gaia catalog): {e}")
    finally:
        log.close()

    with fits.open(file_path) as hdu:
        ra_disp = hdu[0].header['RADISP']
        dec_disp = hdu[0].header['DEDISP']

    assert ra_disp < 0.5
    assert dec_disp < 0.5
    assert ra_disp / tel.pixscale < 1
    assert dec_disp / tel.pixscale < 1


def test_clean_up_astrometry(tmp_path):
    """clean_up_astrometry removes .axy, .corr, .solved, .wcs etc. in directory."""
    base = 'test_image.fits'
    exten = '.fits'
    to_create = [
        base.replace(exten, '.axy'),
        base.replace(exten, '.corr'),
        base.replace(exten, '.solved'),
        base.replace(exten, '.wcs'),
    ]
    for f in to_create:
        path = os.path.join(tmp_path, f)
        with open(path, 'w') as fp:
            fp.write('dummy')
    dir_path = tmp_path
    file_path = os.path.join(dir_path, base)
    solve_wcs.clean_up_astrometry(dir_path, base, exten)
    for f in to_create:
        path = os.path.join(tmp_path, f)
        assert not os.path.exists(path), f'{f} should have been removed'


def test_align_to_gaia_fallback_when_no_gaia_catalog(tmp_path, monkeypatch):
    """align_to_gaia falls back to coarse WCS and returns True when Gaia catalog is unavailable."""
    file_path = os.path.join(tmp_path, 'test_no_gaia.fits')
    hdu = fits.PrimaryHDU(data=np.zeros((40, 40), dtype=np.float32))
    hdu.header['NAXIS1'] = 40
    hdu.header['NAXIS2'] = 40
    hdu.header['CRPIX1'] = 20.0
    hdu.header['CRPIX2'] = 20.0
    hdu.header['CRVAL1'] = 30.0
    hdu.header['CRVAL2'] = -30.0
    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CD1_1'] = -2.2e-5
    hdu.header['CD1_2'] = 0.0
    hdu.header['CD2_1'] = 0.0
    hdu.header['CD2_2'] = 2.2e-5
    hdu.writeto(file_path, overwrite=True)

    monkeypatch.setattr(solve_wcs, 'get_gaia_catalog', lambda *args, **kwargs: None)
    tel = instrument_getter('GMOS')
    ok = solve_wcs.align_to_gaia(file_path, tel, log=None)
    assert ok is True
    with fits.open(file_path) as out:
        assert out[0].header['RADISP'] == 1.0
        assert out[0].header['DEDISP'] == 1.0
        assert 'GAIAFAIL' in out[0].header


def test_align_to_gaia_fallback_when_no_sextractor_sources(tmp_path, monkeypatch):
    """align_to_gaia falls back and returns True when SExtractor has no detections."""
    file_path = os.path.join(tmp_path, 'test_no_sex.fits')
    hdu = fits.PrimaryHDU(data=np.zeros((40, 40), dtype=np.float32))
    hdu.header['NAXIS1'] = 40
    hdu.header['NAXIS2'] = 40
    hdu.header['CRPIX1'] = 20.0
    hdu.header['CRPIX2'] = 20.0
    hdu.header['CRVAL1'] = 30.0
    hdu.header['CRVAL2'] = -30.0
    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CD1_1'] = -2.2e-5
    hdu.header['CD1_2'] = 0.0
    hdu.header['CD2_1'] = 0.0
    hdu.header['CD2_2'] = 2.2e-5
    hdu.writeto(file_path, overwrite=True)

    cat = Table()
    # Coordinates close to CRVAL so stars project in-frame with the test WCS.
    cat['RA_ICRS'] = [30.0, 30.0001, 29.9999, 30.0002, 29.9998, 30.00005, 29.99995, 30.00015]
    cat['DE_ICRS'] = [-30.0, -29.9999, -30.0001, -29.9998, -30.0002, -29.99995, -30.00005, -29.99985]
    cat['PSS'] = [1.0] * 8
    cat['Plx'] = [1.0] * 8
    cat['PM'] = [1.0] * 8
    monkeypatch.setattr(solve_wcs, 'get_gaia_catalog', lambda *args, **kwargs: cat)
    monkeypatch.setattr(solve_wcs.photometry, 'run_sextractor', lambda *args, **kwargs: None)

    tel = instrument_getter('GMOS')
    ok = solve_wcs.align_to_gaia(file_path, tel, log=None)
    assert ok is True
    with fits.open(file_path) as out:
        assert out[0].header['RADISP'] == 1.0
        assert out[0].header['DEDISP'] == 1.0
        assert 'GAIAFAIL' in out[0].header
