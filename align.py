#!/usr/bin/env python

"Function to pixel align a list of images using astroalign"
"(Beroiz, M., Cabral, J. B., & Sanchez, B. Astronomy and Computing, Volume 32, July 2020, 100384)."
"Author: Kerry Paterson"

import numpy as np
from astropy.io import fits
import astroalign as aa

def align_stars(images,ref=None,hdu=0):
    """
    Pixel align a list of images using astroalign
    
    Parameters
    ----------
    images: list of numpy arrays or list of file names
        List of images to align.

    ref: numpy array or file name, optional
        Image to use as a reference point for alginment, all images
        in the list will be aligned to this image. If none
        the first image in the list will be used.

    hdu: int
        Hdu of data if list of file names is provided. Default is 0.
    
    Returns
    -------
    algined_arrays: numpy array
        List of pixel aligned numpy arrays
    """

    #determine type of input and read data to list if file names
    list_type = type(images[0])
    if list_type is str:
        image_arrays = []
        for f in images:
            with fits.open(f) as fo:
                image_data = fo[hdu].data
            image_arrays.append(np.nan_to_num(image_data).astype(np.float64))
    else:
        image_arrays = [np.nan_to_num(x).astype(np.float64) for x in images]

    #check for a reference image and determine type
    if ref is None:
        ref_array = image_arrays.pop(0)
    else:
        ref_type = type(ref)
        if ref_type is str:
            with fits.open(ref) as fo:
                ref_data = fo[hdu].data
            ref_array = np.nan_to_num(ref_data).astype(np.float64)
        else:
            ref_array = np.nan_to_num(ref).astype(np.float64)

    #align using astroalign
    aligned_arrays = []
    aligned_arrays.append(ref_array)
    for image in image_arrays:
        try:
            aligned, _ = aa.register(image,ref_array)
            aligned_arrays.append(aligned)
        except:
            print('MaxIterError: Max iterations exceeded while trying to find acceptable transformation: image not aligned.' )

    return aligned_arrays