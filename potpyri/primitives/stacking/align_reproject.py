"""Align reduced images to a common WCS (astrometry.net + Gaia + reprojection)."""
from __future__ import annotations

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from ccdproc import CCDData
from ccdproc import wcs_project

from potpyri.primitives import solve_wcs

from .wcs_geom import generate_wcs, get_fieldcenter, remove_pv_distortion


def align_images(reduced_files, paths, tel, binn, use_wcs=None, fieldcenter=None,
    out_size=None, skip_gaia=False, keep_all_astro=False, log=None):
    """Solve WCS and align reduced images to a common grid.

    Runs astrometry.net and optionally Gaia alignment, then reprojects to
    a common WCS.

    Parameters
    ----------
    reduced_files : list of str
        Paths to reduced FITS files.
    paths : dict
        Paths dict from options.add_paths.
    tel : Instrument
        Instrument instance.
    binn : str
        Binning string.
    use_wcs : astropy.wcs.WCS, optional
        If provided, use this WCS instead of solving.
    fieldcenter : sequence, optional
        [ra, dec] for generating WCS.
    out_size : int, optional
        Output image size.
    skip_gaia : bool, optional
        If True, skip Gaia alignment. Default is False.
    keep_all_astro : bool, optional
        If True, keep all images regardless of astrometric dispersion.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    tuple or None
        (aligned_data, masks, errors) as lists, or None if no images aligned.
    """
    solved_images = []
    aligned_images = []
    aligned_data = []
    for file in reduced_files:
        # Coarse WCS solution using astrometry.net
        success = solve_wcs.solve_astrometry(file, tel, binn, paths, 
            shift_only=False, log=log)
        if not success: continue

        # Fine WCS solution using Gaia DR3 point sources
        if skip_gaia:
            hdu = fits.open(file)
            hdu[0].header['RADISP']=0.0
            hdu[0].header['DEDISP']=0.0
            hdu.writeto(file, overwrite=True)
        else:
            success = solve_wcs.align_to_gaia(file, tel, log=log)
            if not success: continue

        solved_images.append(file)

    # Reject images from stack where either RADISP or DEDISP is >5-sigma outlier
    if not keep_all_astro:
        if len(solved_images)>2:
            radisp = [] ; dedisp = []
            for file in solved_images:
                hdu = fits.open(file, mode='readonly')
                radisp.append(hdu[0].header['RADISP'])
                dedisp.append(hdu[0].header['DEDISP'])

            radisp = np.array(radisp) ; dedisp = np.array(dedisp)
            ramean, ramedian, rastddev = sigma_clipped_stats(radisp)
            demean, demedian, destddev = sigma_clipped_stats(dedisp)

            if log: log.info(f'Median dispersion in R.A.={ramedian}')
            if log: log.info(f'Median dispersion in Decl.={demedian}')

            mask = (radisp-ramedian <= 5 * rastddev) &\
                   (dedisp-demedian <= 5 * destddev)

            if log:
                log.info('Rejecting the following images for high astrometric dispersion:')
                if np.all(mask):
                    log.info('No images rejected')
                else:
                    for i,m in enumerate(mask):
                        if not m: log.info(solved_images[i])

            solved_images = np.array(solved_images)[mask]
    else:
        if log:
            log.info('Keeping all images (even poorly aligned).')
        solved_images = np.array(solved_images)


    if len(solved_images)==0:
        return(None, None)

    # Determine what the value of use_wcs should be
    if use_wcs is None:
        if fieldcenter is None:
            fieldcenter = get_fieldcenter(solved_images)
            use_wcs = generate_wcs(tel, binn, fieldcenter, out_size)

    for file in solved_images:
        if log: log.info(f'Reprojecting {file}...')

        # Create a reprojection of the file in a common astrometric frame
        hdu = fits.open(file)
        header = remove_pv_distortion(hdu[0].header)
        wcs = WCS(hdu[0].header)

        if use_wcs is None:
            use_wcs = wcs

        image = CCDData(hdu[0].data, meta=hdu[0].header, wcs=wcs, 
            unit=u.electron)

        if out_size:
            target_shape = (out_size, out_size)
        else:
            target_shape = None

        reprojected = wcs_project(image, target_wcs=use_wcs, order='bilinear',
            target_shape=target_shape)
        # Get rid of mask in data
        reprojected.mask = None
        outfile = file.replace('.fits','_reproj.fits')
        reprojected.write(file.replace('.fits','_reproj.fits'), overwrite=True)

        aligned_images.append(outfile)
        aligned_data.append(reprojected)

    return(aligned_images, aligned_data)
