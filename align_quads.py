#!/usr/bin/env python

"Function to pixel align a list of images using quads."
"Author: Kerry Paterson"

__version__ = "1.2" #last updated 18/03/2021

import numpy as np
from astropy.io import fits
from skimage import transform as tf
from astropy.nddata import CCDData
import astropy.units.astrophys as u
import astropy.units as u
import importlib
import solve_wcs
import tel_params

def align_stars(images,telescope,hdu=0):
    """
    Pixel align a list of images
    
    Parameters
    ----------
    images: list of file names
        List of images to align.

    tel: telescope parameter file
        Telescope needed for sextractor.
    
    hdu: int
        Hdu of data if list of file names is provided. Default is 0.
    
    Returns
    -------
    algined_arrays: numpy array
        List of pixel aligned numpy arrays
    """

    #import telescope parameter file
    try:
        tel = importlib.import_module('tel_params.'+telescope)
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)
    
    #run sextractor
    stars_list = []
    d_list = []
    ds_list = []
    ratios_list = []
    for f in images:
        cat_name = f.replace('.fits','.cat')
        table = solve_wcs.run_sextractor(f, cat_name, tel, sex_config_dir='./Config')
        table = table[(table['FLAGS']==0)&(table['IMAFLAGS_ISO']==0)&(table['EXT_NUMBER']==tel.wcs_extension()+1)]
        table.sort('MAG_BEST')
        stars, d, ds, ratios = solve_wcs.make_quads(table['XWIN_IMAGE'],
                table['YWIN_IMAGE'], use=10)
        stars_list.append(stars)
        d_list.append(d)
        ds_list.append(ds)
        ratios_list.append(ratios)

    #match quads and calcculate shifts
    shift_x = []
    shift_y = []
    for i in range(len(stars_list)-1):
        starsx1, starsy1, starsx2, starsy2 = solve_wcs.match_quads(stars_list[0],stars_list[i+1],
                d_list[0],d_list[i+1],ds_list[0],ds_list[i+1],ratios_list[0],ratios_list[i+1],sky_coords=False)
        shift_x.append(np.median([np.array(starsx1[j])-np.array(starsx2[j]) for j in range(len(starsx1))]))
        shift_y.append(np.median([np.array(starsy1[j])-np.array(starsy2[j]) for j in range(len(starsy1))]))

    #apply shifts
    image_arrays = []
    header_arrays = []
    for i, f in enumerate(images):
        with fits.open(f) as fo:
            image_data = fo[hdu].data
            image_arrays.append(np.nan_to_num(image_data).astype(np.float64))
            header_arrays.append(fo[hdu].header)
    aligned_arrays = []
    aligned_arrays.append(CCDData(image_arrays[0],meta=header_arrays[0],unit=u.electron/u.second))
    for i in range(len(image_arrays)-1):
        tform = tf.SimilarityTransform(scale=1, rotation=0, translation=(-shift_x[i], -shift_y[i]))
        aligned_arrays.append(CCDData(tf.warp(image_arrays[i+1], tform),meta=header_arrays[i],unit=u.electron/u.second))
    
    for i in range(len(aligned_arrays)):
        aligned_arrays[i].write(images[i].replace('.fits','_aligned.fits'),overwrite=True)

    return aligned_arrays