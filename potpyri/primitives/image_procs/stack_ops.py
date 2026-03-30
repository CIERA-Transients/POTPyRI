"""Detrending, relative flux calibration, and median/average stacking."""
from __future__ import annotations

import logging
import sys

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from ccdproc import combine

from potpyri.primitives import photometry


def detrend_stack(stack):
    """Remove row and column median from stack science data (detrend).

    Parameters
    ----------
    stack : astropy.io.fits.HDUList
        HDU list with [0]=SCI, [1]=MASK. Modified in place.

    Returns
    -------
    astropy.io.fits.HDUList
        The same stack reference (modified in place).
    """
    data = stack[0].data
    mask = stack[1].data.astype(bool)

    mean, med, stddev = sigma_clipped_stats(data, mask=mask, axis=1,
        sigma_upper=2.5)
    data = data - med[:,None]

    row_med = np.nanmedian(med)

    mean, med, stddev = sigma_clipped_stats(data, mask=mask, axis=0,
        sigma_upper=2.5)
    data = data - med[None,:]

    col_med = np.nanmedian(med)

    # Reapply mask to data
    data[mask] = 0.0

    stack[0].data = data
    stack[0].header['SATURATE'] = stack[0].header['SATURATE'] - (row_med+col_med)

    return(stack)


def compute_relative_scales(data_images, paths, exptimes, log=None,
        match_radius_arcsec=2.0, min_sources=5, min_rel_err=0.01):
    """Compute relative flux scale factors from SExtractor catalogs and RA/Dec cross-matching.

    Runs Source Extractor on each aligned science image, cross-matches sources
    between frames in RA/Dec, and uses the inverse-variance-weighted mean of
    flux ratios (using FLUXERR_AUTO) to define scale factors so that combine()
    can stack on a common flux scale. For frames where matching fails (too few
    matches), the scale is set to 1.0/exptime for that frame.

    Parameters
    ----------
    data_images : list of str
        Paths to FITS with SCI extension (e.g. *_data.fits from image_proc).
    paths : dict
        Paths dict; paths.get('source_extractor') used for SExtractor binary.
    exptimes : array-like
        Exposure time in seconds for each image (same order as data_images).
    log : ColoredLogger, optional
        Logger for progress.
    match_radius_arcsec : float, optional
        Match radius in arcsec for cross-matching sources. Default 2.0.
    min_sources : int, optional
        Minimum matched sources per frame to compute scale. Default 5.
    min_rel_err : float, optional
        Minimum relative flux error (err/flux) used in variance of ratio to
        avoid zero variance. Default 0.01.

    Returns
    -------
    np.ndarray or None
        Relative scale factors, one per image (reference frame scale = 1.0).
        Frames that fail get scale = 1.0/exptime. None if catalogs cannot be
        built at all (caller should use exposure-time-only scaling).
    """
    nimg = len(data_images)
    if nimg < 2:
        return None
    exptimes = np.atleast_1d(np.asarray(exptimes, dtype=float))
    if len(exptimes) != nimg:
        if log:
            log.warning('Relative calibration: exptimes length does not match data_images.')
        return None

    sextractor_path = paths.get('source_extractor')
    def _col(cat, name):
        """Get column by case-insensitive name."""
        for c in cat.colnames:
            if c.upper() == name.upper():
                return c
        return None

    catalogs = []
    for p in data_images:
        cat = photometry.run_sextractor(p, log=log, sextractor_path=sextractor_path)
        if cat is None or len(cat) < min_sources:
            if log:
                log.warning(f'Relative calibration: insufficient catalog for {p}, skipping.')
            return None
        ra_col = _col(cat, 'ALPHA_J2000')
        dec_col = _col(cat, 'DELTA_J2000')
        flux_col = _col(cat, 'FLUX_AUTO')
        fluxerr_col = _col(cat, 'FLUXERR_AUTO')
        if not (ra_col and dec_col and flux_col):
            if log:
                log.warning('Relative calibration: catalog missing RA/Dec or FLUX_AUTO.')
            return None
        catalogs.append((cat, ra_col, dec_col, flux_col, fluxerr_col))

    # Reference = frame with most sources
    ref_idx = int(np.argmax([len(c[0]) for c in catalogs]))
    ref, ref_ra, ref_dec, ref_flux, ref_fluxerr = catalogs[ref_idx]
    ref_coords = SkyCoord(ref[ref_ra], ref[ref_dec], unit='deg')

    match_radius = match_radius_arcsec * u.arcsec
    relative_scales = np.ones(nimg, dtype=float)
    # Reference frame scale = 1.0; failed frames will get 1.0/exptime
    for i in range(nimg):
        if i == ref_idx:
            continue
        cat, ra_col, dec_col, flux_col, fluxerr_col = catalogs[i]
        coords = SkyCoord(cat[ra_col], cat[dec_col], unit='deg')
        # For each source in this frame, find nearest in reference
        idx_ref, d2d, _ = coords.match_to_catalog_sky(ref_coords)
        keep = d2d < match_radius
        n_match = np.sum(keep)
        if n_match < min_sources:
            if log:
                log.warning(f'Relative calibration: frame {i} has {n_match} matches (< {min_sources}), using 1/exptime.')
            relative_scales[i] = 1.0 / exptimes[i]
            continue
        flux_ref = np.array(ref[ref_flux][idx_ref[keep]], dtype=float)
        flux_i = np.array(cat[flux_col][keep], dtype=float)
        # Flux errors: use FLUXERR_AUTO if present, else fall back to unweighted
        if ref_fluxerr and fluxerr_col:
            err_ref = np.array(ref[ref_fluxerr][idx_ref[keep]], dtype=float)
            err_i = np.array(cat[fluxerr_col][keep], dtype=float)
        else:
            err_ref = err_i = None
        # Avoid zeros and negatives in flux
        valid = (flux_i > 0) & (flux_ref > 0)
        if np.sum(valid) < min_sources:
            if log:
                log.warning(f'Relative calibration: frame {i} has too few valid flux pairs, using 1/exptime.')
            relative_scales[i] = 1.0 / exptimes[i]
            continue
        ratio = flux_ref[valid] / flux_i[valid]
        if err_ref is not None and err_i is not None:
            err_ref = err_ref[valid]
            err_i = err_i[valid]
            # Relative errors with floor to avoid zero variance
            rel_err_ref = np.maximum(np.abs(err_ref) / np.maximum(flux_ref[valid], 1e-30), min_rel_err)
            rel_err_i = np.maximum(np.abs(err_i) / np.maximum(flux_i[valid], 1e-30), min_rel_err)
            # Variance of ratio R = A/B: var(R) ≈ R^2 * ( (σ_A/A)^2 + (σ_B/B)^2 )
            var_ratio = ratio ** 2 * (rel_err_ref ** 2 + rel_err_i ** 2)
            var_ratio = np.maximum(var_ratio, 1e-30)
            weight = 1.0 / var_ratio
            relative_scales[i] = float(np.average(ratio, weights=weight))
        else:
            relative_scales[i] = float(np.median(ratio))

    if log:
        log.info(f'Relative calibration: scales (ref frame {ref_idx} = 1.0): {relative_scales}')
    return relative_scales


def stack_data(stacking_data, tel, masks, errors, mem_limit=8.0e9, log=None,
        relative_scales=None):
    """Combine aligned CCDData list with exposure scaling (median or average).

    Masks are applied (masked pixels set to nan) before combination.
    Optionally uses relative flux scales from source cross-matching so that
    frames are combined on a common flux scale.

    Parameters
    ----------
    stacking_data : list of ccdproc.CCDData
        Aligned science images.
    tel : Instrument
        Instrument instance (for stack_method and get_exptime).
    masks : list of astropy.io.fits.ImageHDU
        Per-image masks.
    errors : list of astropy.io.fits.ImageHDU
        Per-image error arrays.
    mem_limit : float, optional
        Memory limit in bytes for ccdproc.combine. Default is 8e9.
    log : ColoredLogger, optional
        Logger for progress.
    relative_scales : np.ndarray, optional
        Relative flux scale factors (one per image; reference = 1.0). If
        provided, scale = relative_scales (exposure-time scaling is not applied).

    Returns
    -------
    astropy.io.fits.HDUList
        HDU list [science, mask, error] (mask/error from first image).
    """
    stack_method = tel.stack_method
    if not stack_method:
        if log:
            log.critical('Could not get stacking method for these images')
            logging.shutdown()
            sys.exit(-1)
        else:
            raise Exception('Could not get stacking method for these images')

    exptimes = []
    for i,stk in enumerate(stacking_data):
        mask = masks[i].data.astype(bool)
        stacking_data[i].data[mask] = np.nan
        exptimes.append(float(tel.get_exptime(stk.header)))

    exptimes = np.array(exptimes)
    scale = 1.0 / exptimes
    if relative_scales is not None:
        scale = np.asarray(relative_scales, dtype=float)

    if len(scale)==1:
        stack_method='average'

    sci_med = combine(stacking_data, scale=scale,
            method=stack_method, mem_limit=mem_limit)

    sci_med = sci_med.to_hdu()

    return(sci_med)


def add_stack_mask(stack, stacking_data):
    """Update stack mask from per-image masks and saturation; set masked pixels to 0.

    Parameters
    ----------
    stack : astropy.io.fits.HDUList
        HDU list with [0]=SCI, [1]=MASK, [2]=ERROR. Modified in place.
    stacking_data : list of ccdproc.CCDData
        Per-image data (headers used for SATURATE).

    Returns
    -------
    astropy.io.fits.HDUList
        The same stack reference (modified in place).
    """
    data = stack[0].data
    mask = stack[1].data
    error = stack[2].data

    # Keep track of original masked values
    bp_mask = (mask > 0) | (data==0.0) | np.isnan(data) | np.isnan(mask) |\
        np.isnan(error) | (error==0.0) | np.isinf(data)

    # Reset mask
    mask = np.zeros(mask.shape).astype(np.uint8)

    # Add bad pixels back to mask
    mask[bp_mask] = 1

    # Add saturation to mask
    sat = np.median([d.header['SATURATE'] for d in stacking_data])
    stack[0].header['SATURATE'] = sat
    sat_mask = data >= sat
    mask[sat_mask] = 4

    # Set masked values to 0
    data[mask>0] = 0.0
    error[mask>0] = 0.0

    # Make sure that nan and inf values are removed from data
    data[np.isnan(data)]=0.0
    data[np.isinf(data)]=0.0

    stack[0].data = data
    stack[1].data = mask
    stack[2].data = error

    return(stack)
