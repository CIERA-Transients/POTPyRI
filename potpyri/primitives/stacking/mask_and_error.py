"""Bad-pixel, cosmic-ray, and saturation masking; per-pixel error arrays."""
from __future__ import annotations

import copy
import time

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from ccdproc import CCDData
from ccdproc import cosmicray_lacosmic

from potpyri._version import __version__
from potpyri.primitives import photometry


def create_mask(science_data, saturation, rdnoise, sigclip=3.0,
    sigfrac=0.10, objlim=4.0, niter=6, outpath='', grow=0, cosmic_ray=True,
    fsmode='convolve', cleantype='medmask', log=None):
    """Build a bad-pixel/saturation/cosmic-ray mask using LACosmic and optional growth.

    Parameters
    ----------
    science_data : str
        Path to science FITS file.
    saturation : float
        Saturation level (ADU); pixels above this get flag 4.
    rdnoise : float
        Read noise for LACosmic (electrons).
    sigclip : float, optional
        LACosmic sigma clipping. Default is 3.0.
    sigfrac : float, optional
        LACosmic sigfrac. Default is 0.10.
    objlim : float, optional
        LACosmic objlim. Default is 4.0.
    niter : int, optional
        LACosmic iterations. Default is 6.
    outpath : str, optional
        Output path for debug files. Default is ''.
    grow : int, optional
        Pixels to grow mask (adds flag 8). Default is 0.
    cosmic_ray : bool, optional
        If True, run LACosmic. Default is True.
    fsmode : str, optional
        LACosmic fsmode. Default is 'convolve'.
    cleantype : str, optional
        LACosmic cleantype. Default is 'medmask'.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    tuple
        (cleaned_data_ndarray, mask_ImageHDU). Mask uses additive flags:
        1=bad, 2=cosmic ray, 4=saturated, 8=neighbor.
    """
    t_start = time.time()
    
    if log:
        log.info(f'Running astroscrappy on {science_data}')
    else:
        print(f'Running astroscrappy on {science_data}')

    hdu = fits.open(science_data)
    data = np.ascontiguousarray(hdu[0].data.astype(np.float32))

    # Estimate FWHM size for LACosmic
    table = photometry.run_sextractor(science_data, log=log)
    if table is not None:
        # Clip fwhm_stars by fwhm
        fwhm, meanfwhm, stdfwhm = sigma_clipped_stats(table['FWHM_IMAGE'])
        if log: 
            log.info(f'Using FWHM for cosmic rays: {fwhm}')
        else:
            print(f'Using FWHM for cosmic rays: {fwhm}')
    else:
        fwhm = 3.5

    # Astroscrappy requires added sky background, so add this value back
    # Set the sky background to some nominal value if it is too low for CR rej
    skybkg = hdu[0].header['SKYBKG']
    if log:
        log.info(f'Sky background in science frame is {skybkg}')
    else:
        print(f'Sky background in science frame is {skybkg}')

    if skybkg < 2000.0: 
        skybkg = 2000.0
        if log:
            log.info('Setting sky background to 2000.0')
        else:
            print('Setting sky background to 2000.0')

    # This needs to be done before data is modified to preserve bad pixels
    mask_bp = (data==0.0) | np.isnan(data)

    # Set data background to skybkg
    data = data + skybkg
    # Also need to adjust the saturation level by SKYBKG for saturated pixels
    saturation += skybkg
    
    if log: log.info('Masking saturated pixels.')

    mask_sat = np.zeros(data.shape).astype(bool) # create empty mask
    mask_sat = mask_sat.astype(np.uint8) #set saturated star mask type
    mask_sat[data >= saturation] = 4 #set saturated pixel flag

    if log: log.info('Cleaning and masking cosmic rays.')
    if log: log.info(f'Using sigclip={sigclip}, sigfrac={sigfrac}, objlim={objlim}')
    #clean science image of cosmic rays and create cosmic ray mask
    inmask = (mask_sat+mask_bp).astype(bool)

    scidata = CCDData(data, unit=u.electron, mask=inmask, 
        wcs=WCS(hdu[0].header), meta=hdu[0].header)

    mask_cr = np.zeros(data.shape)
    mask_cr = mask_cr.astype(np.uint8)

    psfsize = int(np.round(2.5*fwhm))
    if psfsize%2==0: psfsize+=1

    if cosmic_ray:
        newdata, mask_cr = cosmicray_lacosmic(data,
            readnoise=rdnoise, satlevel=saturation, verbose=True,
            sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=niter,
            psffwhm=fwhm, psfsize=psfsize, fsmode=fsmode,
            cleantype=cleantype)

        mask_cr = mask_cr.astype(np.uint8) #set cosmic ray mask type
    else:
        newdata = copy.copy(data)

    mask_cr[mask_cr == 1] = 2 #set cosmic ray flag

     #combine bad pixel, cosmic ray, saturated star and satellite trail masks
    mask = mask_bp+mask_sat+mask_cr

    if grow>0:
        shape = mask.shape
        grow_mask = np.zeros(shape)
        shape = mask.shape
        mask = data > 0
        xs, ys = np.where(mask.astype(bool))
        # Grow mask by grow factor
        g = int(np.ceil((grow-1)/2))
        for p in zip(xs, ys):
            grow_mask[np.max([p[0]-g,0]):np.min([p[0]+g,shape[0]-1]),
                      np.max([p[1]-g,0]):np.min([p[1]+g,shape[1]-1])]=8

        grow_mask = grow_mask.astype(np.uint8)

        mask = mask_bp+mask_sat+mask_cr+grow_mask

    # Get number of bad, cosmic-ray flagged, and saturated pixels
    nbad = np.sum(mask & 1 == 1)
    ncr = np.sum(mask & 2 == 2)
    nsat = np.sum(mask & 4 == 4)
    ngrow = np.sum(mask & 8 == 8)
    
    mask_hdu = fits.ImageHDU(mask) #create mask Primary HDU
    mask_hdu.header['VER'] = (__version__, 
        'Version of image procedures used.') #name of log file
    mask_hdu.header['USE'] = 'Complex mask using additive flags.'#header comment
    mask_hdu.header['M-BP'] = (1, 'Value of masked bad pixels.')
    mask_hdu.header['M-BPNUM'] = (nbad, 'Number of bad pixels.')
    mask_hdu.header['M-CR'] = (2, 'Value of masked cosmic ray pixels.')
    mask_hdu.header['M-CRNUM'] = (ncr, 'Number of cosmic ray pixels.')
    mask_hdu.header['SATURATE'] = (saturation, 'Level of saturation.')
    mask_hdu.header['M-SP'] = (4, 'Value of masked saturated pixels.')
    mask_hdu.header['M-SPNUM'] = (nsat, 'Number of saturated pixels.')
    mask_hdu.header['M-NE'] = (8, 'Value of masked neighbor pixels.')
    mask_hdu.header['M-NENUM'] = (ngrow, 'Number of neighboring masked pixels.')
    
    if log: 
        log.info('Mask created.')
        log.info(f'{nbad} bad, {ncr} cosmic-ray, and {nsat} saturated pixels')
        log.info(f'{ngrow} neighbor masked pixels.')
    else:
        print('Mask created.')
        print(f'{nbad} bad, {ncr} cosmic-ray, and {nsat} saturated pixels')
        print(f'{ngrow} neighbor masked pixels.')

    t_end = time.time()
    if log: 
        log.info(f'Mask creation completed in {t_end-t_start} sec')
    else:
        print(f'Mask creation completed in {t_end-t_start} sec')

    return(newdata, mask_hdu)


def create_error(science_data, mask_data, rdnoise):
    """Compute per-pixel error array from science image, mask, and read noise.

    Error = sqrt(poisson + rms^2 + rdnoise^2). Masked pixels are excluded
    from RMS estimation. science_data can be a path or HDU; mask_data is
    an HDU with mask array.

    Parameters
    ----------
    science_data : str or astropy.io.fits.HDUList
        Path to science FITS or open HDU list.
    mask_data : astropy.io.fits.PrimaryHDU or ImageHDU
        HDU containing boolean or integer mask (masked pixels excluded from rms).
    rdnoise : float
        Read noise in electrons.

    Returns
    -------
    astropy.io.fits.PrimaryHDU
        Error array in electrons (BUNIT='ELECTRONS').
    """
    with fits.open(science_data) as hdu:
        img_data = hdu[0].data.astype(np.float32)
        mask = mask_data.data.astype(bool)

        # Check for issue with all of the file being masked
        if np.all(mask):
            rms = hdu[0].header['SATURATE']
        else:
            rms = 0.5 * (
                np.percentile(img_data[~mask], 84.13)
                - np.percentile(img_data[~mask], 15.86)
            )

        poisson = img_data.copy()
        poisson[poisson < 0.] = 0.
        error = np.sqrt(poisson + rms**2 + rdnoise)

        # Sanitize error array
        mask = np.isnan(error)
        error[mask] = np.nanmedian(error)
        mask = error < 0.0
        error[mask] = np.nanmedian(error)
        maxval = np.float32(hdu[0].header['SATURATE'])
        mask = np.isinf(error)
        error[mask] = maxval

    error_hdu = fits.PrimaryHDU(error)  # create mask Primary HDU
    error_hdu.header['VER'] = (__version__, 
        'Version of image procedures used used.')
    error_hdu.header['USE'] = 'Error array for Poisson, read, and RMS noise'
    error_hdu.header['BUNIT'] = ('ELECTRONS', 'Units of the error array.')

    return(error_hdu)
