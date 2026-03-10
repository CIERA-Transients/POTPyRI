"""Tests for photometry.photloop (PSF and aperture photometry on stacks) and helpers (create_conv, create_params, extract_*, get_star_catalog, extract_fwhm_from_epsf)."""
import warnings

from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import photometry
from potpyri.instruments import instrument_getter

import os
import numpy as np

from astropy.io import fits
from astropy.table import Table
from photutils.psf import EPSFModel

import pytest
from tests.utils import download_gdrive_file


@pytest.mark.integration
def test_photometry(tmp_path):
    """Run photloop on GMOS stack; assert APPPHOT, PSFPHOT, PSFSTARS, PSF extensions and catalog size."""
    instrument = 'GMOS'
    file_list_name = 'files.txt'

    # Raw science file (download to tmp_path so log and overwrites are writable)
    file_path = download_gdrive_file('GMOS/red/sGRB240615A-GRB.i.ut240618.12.22.stk.fits.fz', output_dir=str(tmp_path), use_cached=True)

    # Strip out extensions added by photometry loop if they exist
    hdu = fits.open(file_path)
    for key in ['APPPHOT','PSFPHOT','PSFSTARS','RESIDUAL','PSF']:
        if key in [h.name for h in hdu]:
            del hdu[key]

    # Sanitize error array
    mask = np.isnan(hdu['ERROR'].data)
    hdu['ERROR'].data[mask] = np.nanmedian(hdu['ERROR'].data)
    mask = hdu['ERROR'].data < 0.0
    hdu['ERROR'].data[mask] = np.nanmedian(hdu['ERROR'].data)
    maxval = np.max(hdu['SCI'].data)
    mask = np.isinf(hdu['ERROR'].data)
    hdu['ERROR'].data[mask] = maxval

    hdu.writeto(file_path, overwrite=True)

    data_path, basefile = os.path.split(file_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    phot_sn_min = 20.0
    phot_sn_max = 100.0
    fwhm_init = 5.0

    try:
        photometry.photloop(file_path, phot_sn_min=phot_sn_min,
            fwhm_init=fwhm_init, phot_sn_max=phot_sn_max, log=log)
    finally:
        log.close()

    with fits.open(file_path) as hdu:
        assert np.any(['PSFPHOT' in h.name for h in hdu])
        assert np.any(['PSFSTARS' in h.name for h in hdu])
        assert np.any(['APPPHOT' in h.name for h in hdu])
        assert np.any(['PSF' in h.name for h in hdu])
        assert len(hdu['APPPHOT'].data) > 100


def test_create_conv(tmp_path):
    """create_conv writes a 3x3 CONV NORM kernel file."""
    outfile = os.path.join(tmp_path, 'test.conv')
    photometry.create_conv(outfile)
    assert os.path.exists(outfile)
    with open(outfile) as f:
        content = f.read()
    assert 'CONV NORM' in content
    assert '1 2 1' in content
    assert '2 4 2' in content


def test_create_params(tmp_path):
    """create_params writes parameter file with NUMBER, X_IMAGE, FWHM_IMAGE, etc."""
    outfile = os.path.join(tmp_path, 'test.param')
    photometry.create_params(outfile)
    assert os.path.exists(outfile)
    with open(outfile) as f:
        content = f.read()
    assert 'NUMBER' in content
    assert 'X_IMAGE' in content
    assert 'FWHM_IMAGE' in content
    assert 'MAG_AUTO' in content


def test_extract_aperture_stats(tmp_path):
    """extract_aperture_stats returns table with fwhm, flux_best, Xpos, Ypos columns."""
    from photutils.aperture import CircularAperture, ApertureStats
    shape = (64, 64)
    img_data = np.ones(shape, dtype=float) * 100.0
    img_data[32, 32] = 500.0
    img_mask = np.zeros(shape, dtype=bool)
    img_error = np.ones(shape, dtype=float) * 5.0
    stars = Table()
    stars['xcentroid'] = [32.0]
    stars['ycentroid'] = [32.0]
    tel = instrument_getter('MMIRS')
    paths = options.add_paths(str(tmp_path), 'files.txt', tel)
    log = logger.get_log(paths['log'])
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='.*Units from inserted quantities.*', category=UserWarning)
            stats = photometry.extract_aperture_stats(
                img_data, img_mask, img_error, stars,
                aperture_radius=5.0, log=log)
        assert len(stats) == 1
        assert 'fwhm' in stats.colnames
        assert 'flux_best' in stats.colnames
        assert 'Xpos' in stats.colnames
        assert 'Ypos' in stats.colnames
    finally:
        log.close()


def test_get_star_catalog(tmp_path):
    """get_star_catalog finds sources and returns table with aperture stats."""
    shape = (128, 128)
    y, x = np.ogrid[:shape[0], :shape[1]]
    img_data = np.ones(shape, dtype=float) * 50.0
    # Add a clear Gaussian star so DAOStarFinder and sharpness/roundness pass
    for cy, cx in [(64, 64), (32, 96)]:
        img_data += 2000.0 * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * 9.0))
    img_mask = np.zeros(shape, dtype=bool)
    img_error = np.ones(shape, dtype=float) * 5.0
    tel = instrument_getter('MMIRS')
    paths = options.add_paths(str(tmp_path), 'files.txt', tel)
    log = logger.get_log(paths['log'])
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='.*Units from inserted quantities.*', category=UserWarning)
            stars = photometry.get_star_catalog(
                img_data, img_mask, img_error,
                fwhm_init=4.0, threshold=10.0, log=log)
        assert stars is not None
        assert len(stars) >= 1
        assert 'xcentroid' in stars.colnames or 'Xpos' in stars.colnames
        assert 'flux_best' in stars.colnames or 'flux' in stars.colnames
    finally:
        log.close()


def test_extract_fwhm_from_epsf():
    """extract_fwhm_from_epsf returns finite FWHM from EPSFModel."""
    # Build a minimal EPSF-like array (Gaussian blob)
    size = 15
    y, x = np.ogrid[-size//2:size//2+1, -size//2:size//2+1]
    sigma = 2.0
    data = np.exp(-(x*x + y*y) / (2 * sigma**2)).astype(float)
    with warnings.catch_warnings():
        try:
            from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning
            warnings.simplefilter('ignore', AstropyDeprecationWarning)
            warnings.simplefilter('ignore', AstropyUserWarning)
        except ImportError:
            warnings.simplefilter('ignore', DeprecationWarning)
            warnings.simplefilter('ignore', UserWarning)
        epsf = EPSFModel(data)
        fwhm = photometry.extract_fwhm_from_epsf(epsf, fwhm_init=3.0)
    assert np.isfinite(fwhm)
    assert fwhm > 0
