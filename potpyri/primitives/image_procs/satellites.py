"""Satellite trail detection and masking via acstools."""
import os

import acstools
import numpy as np
from astropy.io import fits


def mask_satellites(images, filenames, log=None):
    """Detect and mask satellite trails in science images using acstools.satdet.

    Parameters
    ----------
    images : list of ccdproc.CCDData
        Science images; data is modified in place (trail pixels set to nan).
    filenames : list of str
        Paths corresponding to images (used for temp files).
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
    """
    out_data = []
    for i,science_data in enumerate(images):

        if log: log.info(f'Checking for satellite trails in: {filenames[i]}')
        
        data = np.ascontiguousarray(science_data.data.astype(np.float32))

        tmphdu = fits.PrimaryHDU()
        tmpfile = filenames[i].replace('.fits','.tmp.fits')

        tmphdu.data = data
        tmphdu.writeto(tmpfile, overwrite=True)

        tmpshape = data.shape
        buf = int(np.round(np.max(tmpshape)/10.))

        results, errors = acstools.satdet.detsat(tmpfile,
            chips=[0], sigma=4, h_thresh=0.1, buf=buf, verbose=False)

        trail_coords = results[(tmpfile, 0)]

        if len(trail_coords)>0:
            n_trails = len(trail_coords)
            trail_masks = []
            for i,coord in enumerate(trail_coords):
                if log: log.info(f'Masking for satellite trail {i+1}/{n_trails}')
                try:
                    mask = acstools.satdet.make_mask(tmpfile, 0, coord, 
                        verbose=False)
                except ValueError:
                    continue
                trail_masks.append(mask)

            trail_mask = np.any(trail_masks, axis=0)
            
            # Set actual pixel values to nan so they don't contribute to image
            science_data.data[trail_mask] = np.nan

        out_data.append(science_data)

        os.remove(tmpfile)
