"""Integration test: PS1 vs DECaLS zeropoints on a dummy field in survey overlap.

Uses a sky position covered by both Pan-STARRS1 and DECaLS / Legacy Surveys,
downloads *r*-band point-source photometry from each, builds a synthetic
observation whose fluxes are tied to PS1 magnitudes at a known zeropoint, then
checks that recovering the zeropoint with PS1 and with DECaLS yields similar
values (within expected PS1–DECam filter differences).
"""
from __future__ import annotations

import numpy as np
import pytest
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from potpyri.instruments import instrument_getter
from potpyri.primitives import absphot
from potpyri.utils import catalogs


# Field in the PS1 footprint and Legacy Surveys / DECaLS equatorial coverage.
_OVERLAP_RA = 150.0   # deg
_OVERLAP_DEC = 10.0   # deg
_BAND = 'r'
_TRUE_ZPT = 27.0
# PS1 r vs DECam r typically agree to ~few × 0.01 mag for stars; allow margin.
_ZPT_TOL = 0.15
_MIN_MATCHES = 15
_SEARCH_WIDTH = 0.4 * u.deg


def _skip_network(exc):
    pytest.skip(f'Catalog / network unavailable: {exc}')


def _write_dummy_stack(path, ra, dec, flux, flux_err, filt=_BAND):
    """Minimal stack FITS (SCI + APPPHOT) for find_zeropoint."""
    n = len(ra)
    shape = (64, 64)
    w = WCS(naxis=2)
    w.wcs.crpix = [32, 32]
    w.wcs.crval = [float(np.median(ra)), float(np.median(dec))]
    w.wcs.cdelt = [-0.0001, 0.0001]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    sci = fits.PrimaryHDU(np.ones(shape, dtype=np.float32))
    sci.name = 'SCI'
    sci.header.update(w.to_header())
    sci.header['FILTER'] = filt
    sci.header['EXTNAME'] = 'SCI'

    phot = Table({
        'RA': np.asarray(ra, dtype=float),
        'Dec': np.asarray(dec, dtype=float),
        'flux': np.asarray(flux, dtype=float),
        'flux_err': np.asarray(flux_err, dtype=float),
        'Xpos': np.full(n, 32.0),
        'Ypos': np.full(n, 32.0),
        'mag': np.full(n, 20.0),
        'mag_err': np.full(n, 0.05),
        'SN': np.full(n, 50.0),
        'FWHM': np.full(n, 5.0),
        'sky': np.zeros(n),
    })
    app = fits.BinTableHDU(phot, name='APPPHOT')
    mask = fits.ImageHDU(np.zeros(shape, dtype=np.uint8), name='MASK')
    err = fits.ImageHDU(np.ones(shape, dtype=np.float32), name='ERROR')
    fits.HDUList([sci, mask, err, app]).writeto(path, overwrite=True)


def _download_ps1_r(center, width):
    cols = ['RAJ2000', 'DEJ2000', 'rmag', 'e_rmag', 'rKmag']
    tab = catalogs.query_vizier_region(center, width, 'II/349', cols, log=None)
    if tab is None or len(tab) == 0:
        raise RuntimeError('PS1 Vizier query returned no rows')
    tab = tab[np.isfinite(tab['rmag']) & (tab['e_rmag'] > 0)]
    tab = catalogs.apply_point_source_cut(tab, 'II/349', _BAND, mag_col='rmag')
    out = Table()
    out['ra'] = np.asarray(tab['RAJ2000'], dtype=float)
    out['dec'] = np.asarray(tab['DEJ2000'], dtype=float)
    out['mag'] = np.asarray(tab['rmag'], dtype=float)
    out['mag_err'] = np.asarray(tab['e_rmag'], dtype=float)
    return out


def _download_decals_r(center, width):
    tab = catalogs.query_decals_region(center, width, _BAND, log=None)
    if tab is None or len(tab) == 0:
        raise RuntimeError('DECaLS Data Lab query returned no rows')
    tab = catalogs.apply_point_source_cut(
        tab, catalogs.DECALS_TRACTOR_TABLE, _BAND, mag_col='mag')
    return tab


@pytest.mark.integration
def test_ps1_decals_overlap_dummy_zeropoint_r_band(tmp_path):
    """Download PS1 and DECaLS in overlap; dummy obs yields similar r-band ZP."""
    center = SkyCoord(_OVERLAP_RA, _OVERLAP_DEC, unit='deg')

    try:
        ps1 = _download_ps1_r(center, _SEARCH_WIDTH)
        decals = _download_decals_r(center, _SEARCH_WIDTH)
    except (requests.exceptions.RequestException, OSError, RuntimeError,
            ValueError, TimeoutError) as exc:
        _skip_network(exc)
    except Exception as exc:
        # Data Lab / Vizier client errors vary by version.
        _skip_network(exc)

    assert len(ps1) >= _MIN_MATCHES, f'PS1 too sparse: {len(ps1)} sources'
    assert len(decals) >= _MIN_MATCHES, f'DECaLS too sparse: {len(decals)} sources'

    # Cross-match surveys (common stars for a fair ZP comparison).
    c_ps1 = SkyCoord(ps1['ra'], ps1['dec'], unit='deg')
    c_dec = SkyCoord(decals['ra'], decals['dec'], unit='deg')
    idx, sep, _ = c_ps1.match_to_catalog_sky(c_dec)
    match = sep < 1.0 * u.arcsec
    # Prefer unsaturated, well-measured stars used for ZP (above bright limit).
    bright_ok = (ps1['mag'] > 16.0) & (ps1['mag'] < 21.0)
    match &= bright_ok
    assert np.sum(match) >= _MIN_MATCHES, (
        f'Need >= {_MIN_MATCHES} PS1–DECaLS matches, got {np.sum(match)}'
    )

    ps1_m = ps1[match]
    dec_m = decals[idx[match]]

    # Catalog magnitudes for the same stars should be similar in r.
    dmag = np.asarray(ps1_m['mag'], float) - np.asarray(dec_m['mag'], float)
    med_dmag = float(np.nanmedian(dmag))
    assert abs(med_dmag) < 0.25, (
        f'Median PS1-DECaLS r mag offset {med_dmag:.3f} too large for overlap field'
    )

    # Dummy observation: fluxes implied by PS1 mags at a known zeropoint.
    flux = 10.0 ** (-0.4 * (np.asarray(ps1_m['mag'], float) - _TRUE_ZPT))
    flux_err = flux * 0.02
    stack = tmp_path / 'dummy_overlap_r.stk.fits'
    _write_dummy_stack(
        str(stack),
        np.asarray(ps1_m['ra'], float),
        np.asarray(ps1_m['dec'], float),
        flux,
        flux_err,
        filt=_BAND,
    )

    tel = instrument_getter('GMOS')
    cal = absphot.absphot(iterations=3, sigma=3.0)

    # Recover ZP with PS1 (live download via get_catalog).
    tel.set_zeropoint_catalog('PS1')
    try:
        ok_ps1 = cal.find_zeropoint(str(stack), tel, log=None)
    except (requests.exceptions.RequestException, OSError) as exc:
        _skip_network(exc)
    assert ok_ps1 is True
    with fits.open(stack) as hdul:
        zpt_ps1 = float(hdul['PRIMARY'].header['ZPTMAG'])
        n_ps1 = int(hdul['PRIMARY'].header.get('ZPTNSTAR', 0))

    # Recover ZP with DECaLS / Legacy on the same dummy photometry.
    tel.set_zeropoint_catalog('DECALS')
    try:
        ok_dec = cal.find_zeropoint(str(stack), tel, log=None)
    except (requests.exceptions.RequestException, OSError) as exc:
        _skip_network(exc)
    assert ok_dec is True
    with fits.open(stack) as hdul:
        zpt_dec = float(hdul['PRIMARY'].header['ZPTMAG'])
        n_dec = int(hdul['PRIMARY'].header.get('ZPTNSTAR', 0))
        zpt_cat = hdul['PRIMARY'].header.get('ZPTCAT')

    assert zpt_cat == 'DECALS'
    assert n_ps1 >= 5
    assert n_dec >= 5
    assert abs(zpt_ps1 - _TRUE_ZPT) < 0.05, (
        f'PS1 ZP {zpt_ps1:.3f} should recover injected {_TRUE_ZPT}'
    )
    # Fluxes were built from PS1 mags, so DECaLS recovers
    #   ZPT_DECaLS ≈ ZPT_PS1 - median(m_PS1 - m_DECaLS)
    # for the same stars (PS1–DECam filter zero-point / color offset).
    expected_dec = zpt_ps1 - med_dmag
    assert abs(zpt_dec - expected_dec) < _ZPT_TOL, (
        f'DECaLS ZP {zpt_dec:.3f} vs expected {expected_dec:.3f} '
        f'(PS1 ZP {zpt_ps1:.3f} minus median PS1-DECaLS offset {med_dmag:.3f})'
    )
    # Direct survey ZPs should still be close for stellar r-band.
    assert abs(zpt_dec - zpt_ps1) < 0.30, (
        f'Raw PS1/DECaLS ZP difference {abs(zpt_dec - zpt_ps1):.3f} mag unexpectedly large'
    )
