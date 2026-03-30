"""Orchestrates alignment, masking, stacking, and optional detrend for one target."""
from __future__ import annotations

import logging
import os
import sys
import time

import numpy as np
from astropy.io import fits

from .align import align_images
from .mask_error import create_error, create_mask
from .satellites import mask_satellites
from .stack_ops import (
    add_stack_mask,
    compute_relative_scales,
    detrend_stack,
    stack_data,
)
from .wcs_geom import generate_wcs


def _image_proc_worker(image_data, tel, paths, skip_skysub=False,
    fieldcenter=None, out_size=None, satellites=True, cosmic_ray=True,
    skip_gaia=False, keep_all_astro=False, relative_calibration=False,
    log=None):
    """Full image processing: align, mask, stack, and optionally detrend.

    Orchestrates WCS solving, alignment, satellite/cosmic-ray masking,
    stacking, and optional sky subtraction. Optionally calibrates frames
    to each other using SExtractor catalogs and RA/Dec cross-matching before stacking.

    Calibrated workspace ``*.fits`` and intermediate ``*_reproj.fits`` are
    removed from disk once they are no longer needed; final per-frame products
    are ``*_data.fits`` under ``paths['work']``.

    Parameters
    ----------
    image_data : astropy.table.Table
        File table subset (e.g. one TargType) with 'File' column.
    tel : Instrument
        Instrument instance.
    paths : dict
        Paths dict from options.add_paths.
    skip_skysub : bool, optional
        If True, skip sky subtraction. Default is False.
    fieldcenter : sequence, optional
        [ra, dec] for alignment.
    out_size : int, optional
        Output image size.
    satellites : bool, optional
        If True, mask satellite trails. Default is True.
    cosmic_ray : bool, optional
        If True, run cosmic-ray rejection. Default is True.
    skip_gaia : bool, optional
        If True, skip Gaia alignment. Default is False.
    keep_all_astro : bool, optional
        If True, keep all images regardless of astrometric dispersion.
    relative_calibration : bool, optional
        If True, run SExtractor on each aligned frame, cross-match sources in
        RA/Dec, and use relative fluxes to set scale factors for the stack
        combine. Default is False.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    str or None
        Path to stacked output FITS, or None if processing failed.
    """

    wavelength = tel.wavelength

    red_path = paths['red']
    work_path = paths['work']

    # All data in the table should be for one target
    assert np.all([image_data['Target']==image_data['Target'][0]])
    assert np.all([image_data['TargType']==image_data['TargType'][0]])
    
    cal_type = image_data['TargType'][0]
    target = image_data['Target'][0]
    fil = image_data['Filter'][0]
    amp = image_data['Amp'][0]
    binn = image_data['Binning'][0]

    # Load bias frame
    if tel.bias:
        if log: log.info('Loading master bias.')
        try:
            mbias = tel.load_bias(paths, amp, binn)
        except Exception as e:
            if log: 
                log.error('No master bias found for this configuration')
                log.error(f'Skipping reduction for: {cal_type}')
                log.error(e)
            return(None)
    else:
        mbias = None

    # Load bias frame
    if tel.dark:
        if log: log.info('Loading master dark.')
        try:
            mdark = tel.load_dark(paths, amp, binn)
        except Exception as e:
            if log: 
                log.error('No master dark found for this configuration')
                log.error(f'Skipping reduction for: {cal_type}')
                log.error(e)
            return(None)
    else:
        mdark = None

    # Load flat frame
    if tel.flat:
        if log: log.info('Loading master flat.')
        try:
            mflat = tel.load_flat(paths, fil, amp, binn)
        except Exception as e:
            if log: 
                log.error('No master bias found for this configuration')
                log.error(f'Skipping reduction for: {cal_type}')
                log.error(e)
            return(None)
    else:
        mflat = None
        
    t1 = time.time()
    if log: log.info(f'Processing data for {cal_type}')

    # Bias subtraction, gain correction, flat correction, and flat fielding
    files = image_data['File']
    processed = tel.process_science(files, fil, amp, binn, paths,
        mbias=mbias, mflat=mflat, mdark=mdark, skip_skysub=skip_skysub, log=log)

    # Get filenames for output processed data
    reduced_files = [p.header['FILENAME'] for p in processed]

    t2 = time.time()
    if log: log.info(f'Data processed in {t2-t1} sec')
    if log: log.info('Aligning images.')

    # If masking satellite trails, this needs to be done before reprojection so
    # that the acstools algorithm accurately models the edges of the detector
    if satellites:
        mask_satellites(processed, reduced_files, log=log)

    # Sort files so deepest exposure is first
    exptimes = []
    for file in reduced_files:
        hdu = fits.open(file)
        exptimes.append(tel.get_exptime(hdu[0].header))
    idx = np.flip(np.argsort(exptimes))
    reduced_files = np.array(reduced_files)[idx]

    if out_size is None:
        out_size = tel.get_out_size(processed[0].header)

    if fieldcenter is not None:
        use_wcs = generate_wcs(tel, binn, fieldcenter, out_size)
    else:
        use_wcs = None

    aligned_images, aligned_data = align_images(reduced_files, paths, tel, binn,
        use_wcs=use_wcs, fieldcenter=fieldcenter, out_size=out_size,
        skip_gaia=skip_gaia, keep_all_astro=keep_all_astro, log=log)

    if aligned_images is None or aligned_data is None:
        return(None)

    # Calibrated workspace frames are no longer needed once reprojected.
    for reproj in aligned_images:
        cal = reproj.replace('_reproj.fits', '.fits')
        if os.path.isfile(cal):
            try:
                os.remove(cal)
            except OSError as e:
                if log:
                    log.warning(f'Could not remove calibrated frame {cal}: {e}')

    if log: log.info('Creating mask and error arrays.')
    masks = []
    errors = []
    data_images = []
    for stack_img in aligned_images:
        hdu = fits.open(stack_img, mode='readonly')
        hdr = hdu[0].header
        clean, mask = create_mask(stack_img,
            hdr['SATURATE'], np.mean(tel.get_rdnoise(hdr)),
            log=log, cosmic_ray=cosmic_ray, outpath=work_path)
        error = create_error(stack_img, mask, np.mean(tel.get_rdnoise(hdr)))

        masks.append(mask)
        errors.append(error)

        # Set nan values to 0 now that they are masked
        data = hdu[0].data
        data[np.isnan(data)] = 0.0

        # Create final image triplet before stacking with data, mask, and error
        hdulist = fits.HDUList()
        sci_hdu = fits.ImageHDU()
        sci_hdu.data = data
        sci_hdu.header = hdu[0].header
        msk_hdu = fits.ImageHDU()
        msk_hdu.data = mask.data
        msk_hdu.header = mask.header
        err_hdu = fits.ImageHDU()
        err_hdu.data = error.data
        err_hdu.header = error.header

        hdulist.append(sci_hdu)
        hdulist.append(msk_hdu)
        hdulist.append(err_hdu)

        hdulist[0].name='SCI'
        hdulist[1].name='MASK'
        hdulist[2].name='ERROR'

        hdulist[0].header['EXTNAME']='SCI'
        hdulist[1].header['EXTNAME']='MASK'
        hdulist[2].header['EXTNAME']='ERROR'

        filename = os.path.basename(hdulist[0].header['FILENAME']).replace(
            '.fits', '_data.fits')
        fullfilename = os.path.join(work_path, filename)

        if log:
            log.info(f'Writing out all file data: {fullfilename}')
        else:
            print(f'Writing out all file data: {fullfilename}')

        hdulist.writeto(fullfilename, overwrite=True, output_verify='silentfix')
        data_images.append(fullfilename)
        hdu.close()

        if os.path.isfile(stack_img):
            try:
                os.remove(stack_img)
            except OSError as e:
                if log:
                    log.warning(f'Could not remove reprojected frame {stack_img}: {e}')

    if log: log.info('Creating median stack.')
    if len(aligned_data)>1:
        relative_scales = None
        if relative_calibration:
            exptimes = np.array([tel.get_exptime(d.header) for d in aligned_data])
            relative_scales = compute_relative_scales(data_images, paths, exptimes, log=log)
        sci_med = stack_data(aligned_data, tel, masks, errors, log=log,
            relative_scales=relative_scales)
        sci_med = add_stack_mask(sci_med, aligned_data)

        if tel.detrend:
            if log: log.info('Detrending stack')
            sci_med = detrend_stack(sci_med)
        else:
            if log: log.info('Skipping detrending')

    else:
        sci_med = fits.open(data_images[0])

    sci_med[0].data = sci_med[0].data.astype(np.float32)
    sci_med[1].data = sci_med[1].data.astype(np.uint8)
    sci_med[2].data = sci_med[2].data.astype(np.float32)
    sci_med[0].header['EXTNAME']='SCI'
    sci_med[1].header['EXTNAME']='MASK'
    sci_med[2].header['EXTNAME']='ERROR'
    sci_med[0].header['BITPIX'] = -32
    sci_med[1].header['BITPIX'] = 8
    sci_med[2].header['BITPIX'] = -32

    # Get time parameters from aligned data
    mid_time = np.average([tel.get_time(d.header) for d in aligned_data])
    exptimes = np.array([tel.get_exptime(d.header) for d in aligned_data])
    eff_time = np.sum(exptimes**2)/np.sum(exptimes)
    total_time = np.sum(exptimes)

    # Rescale both data and error images by effective exposure times so they 
    # are in e- instead of e-/s
    sci_med[0].data = sci_med[0].data * eff_time
    sci_med[2].data = sci_med[2].data * eff_time

    # Explicitly note that data and error extensions are in ELECTRONS
    sci_med[0].header['BUNIT'] = 'ELECTRONS'
    sci_med[2].header['BUNIT'] = 'ELECTRONS'

    # Add both formats for code that requires either
    sci_med[0].header['MJD-OBS'] = (mid_time, 
        'Mid-MJD of the observation sequence.')
    sci_med[0].header['MJD'] = (mid_time, 
        'Mid-MJD of the observation sequence.')
    # Since we generated stack with median, effective exposure should be a 
    # weighted average of the input exposure times
    sci_med[0].header['EXPTIME'] = (eff_time, 
        'Effective expsoure time in seconds.')
    sci_med[0].header['EXPTOT'] = (total_time, 
        'Total exposure time in seconds')
    sci_med[0].header['GAIN'] = (len(aligned_data), 
        'Effecetive gain for stack.')

    # Calculate read noise
    rdnoise = np.mean(tel.get_rdnoise(sci_med[0].header))/np.sqrt(len(aligned_data))
    sci_med[0].header['RDNOISE'] = (rdnoise, 'Readnoise of stack.')
    
    sci_med[0].header['NFILES'] = (len(aligned_data), 
        'Number of images in stack')
    sci_med[0].header['FILTER'] = fil
    sci_med[0].header['OBSTYPE'] = 'OBJECT'
    
    # Generate stack name and write out
    stackname = tel.get_stk_name(sci_med[0].header, red_path)
    sci_med.writeto(stackname, overwrite=True, output_verify='silentfix')
    
    if log: 
        log.info(f'Stack made for {stackname}')
    else:
        print(f'Stack made for {stackname}')

    return(stackname)
