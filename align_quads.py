#!/usr/bin/env python

"Function to pixel align a list of images using quads."
"Author: Kerry Paterson"

__version__ = "1.8" #last updated 30/07/2021

import time
import numpy as np
from astropy.io import fits
from skimage import transform as tf
from astropy.nddata import CCDData
import astropy.units.astrophys as u
import astropy.units as u
import importlib
import solve_wcs
import tel_params

def align_stars(images,telescope,hdu=0,mask=None,log=None):
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

    t_start = time.time()

    if log:
        log.info('Running align_quads version '+str(__version__))

    #import telescope parameter file
    try:
        tel = importlib.import_module('tel_params.'+telescope)
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)
    
    #check for flipped images
    crpix1, crpix2, cd11, cd12, cd21, cd22 = wcs_keyword = tel.WCS_keywords()
    for i, f in enumerate(images):
        rewrite = False
        with fits.open(f) as fo:
            header = fo[hdu].header
            data = fo[hdu].data
        if i==0:
            xsign = header[cd11]
            ysign = header[cd22]
        else:
            xflip = header[cd11]
            yflip = header[cd22]
            try:
                if xsign/xflip < 0:
                    data = np.fliplr(data)
                    rewrite = True
            except ZeroDivisionError:
                pass
            try:
                if ysign/yflip < 0:
                    data = np.flipud(data)
                    rewrite = True
            except ZeroDivisionError:
                pass
            if rewrite:
                fits.writeto(f.replace('.fits','_flipped.fits'),data,header,overwrite=True)
                images[i] = f.replace('.fits','_flipped.fits')
    
    #run sextractor
    stars_list = []
    d_list = []
    ds_list = []
    ratios_list = []
    for f in images:
        if log:
            log.info('Loading file: '+f)
            log.info('Running SExtracor.')
        cat_name = f.replace('.fits','.cat')
        table = solve_wcs.run_sextractor(f, cat_name, tel, sex_config_dir='./Config', log=log)
        table = table[(table['FLAGS']==0)&(table['EXT_NUMBER']==tel.wcs_extension()+1)]
        if mask:
            with fits.open(mask) as mask_hdu:
                stat_mask = mask_hdu[0].data
            table = table[(stat_mask[table['YWIN_IMAGE'].astype(int)-1,table['XWIN_IMAGE'].astype(int)-1]!=0)]
        table.sort('MAG_BEST')
        if log:
            log.info('Creating quads with the 10 brightest non-saturated stars.')
        stars, d, ds, ratios = solve_wcs.make_quads(table['XWIN_IMAGE'],
                table['YWIN_IMAGE'], use=10)
        stars_list.append(stars)
        d_list.append(d)
        ds_list.append(ds)
        ratios_list.append(ratios)

    #match quads and calculate shifts
    if log:
        log.info('Matching quads and calculating shifts.')
        log.info('Using '+images[0]+' as the reference image.')
    shift_x = []
    shift_y = []
    aligned_images = []
    aligned_images.append(images[0])
    for i in range(len(stars_list)-1):
        try:
            starsx1, starsy1, starsx2, starsy2 = solve_wcs.match_quads(stars_list[0],stars_list[i+1],
                    d_list[0],d_list[i+1],ds_list[0],ds_list[i+1],ratios_list[0],ratios_list[i+1],sky_coords=False)
            if log:
                log.info('Found '+str(len(starsx1))+' unique star matches for image '+images[i+1])
        except:
            if log:
                log.error(images[i+1]+' failed to alignment, removing from stack.')
            else:
                print(images[i+1]+' failed to alignment, removing from stack.')
            continue
        x_shift = np.median([np.array(starsx1[j])-np.array(starsx2[j]) for j in range(len(starsx1))])
        y_shift = np.median([np.array(starsy1[j])-np.array(starsy2[j]) for j in range(len(starsy1))])
        shift_x.append(x_shift)
        shift_y.append(y_shift)
        if log:
            log.info('Shifts: (x,y) = %.3f, %.3f'%(x_shift,y_shift))
        aligned_images.append(images[i+1])

    #apply shifts
    if log:
        log.info('Applying shifts.')
    image_arrays = []
    header_arrays = []
    for i, f in enumerate(aligned_images):
        with fits.open(f) as fo:
            image_data = fo[hdu].data
            image_arrays.append(np.nan_to_num(image_data).astype(np.float64))
            header_arrays.append(fo[hdu].header)
    aligned_arrays = []
    aligned_arrays.append(CCDData(image_arrays[0],meta=header_arrays[0],unit=u.electron/u.second))
    for i in range(len(image_arrays)-1):
        tform = tf.EuclideanTransform(rotation=0, translation=(-shift_x[i], -shift_y[i]))
        aligned_arrays.append(CCDData(tf.warp(image_arrays[i+1], tform),meta=header_arrays[i],unit=u.electron/u.second))

    if log:
        log.info('Writing out aligned images.')

    for i in range(len(aligned_arrays)):
        aligned_images[i] = aligned_images[i].replace('.fits','_aligned.fits')
        aligned_arrays[i].write(aligned_images[i],overwrite=True)

    t_end = time.time()
    if log:
        log.info('Align_quads ran in '+str(t_end-t_start)+' sec')
    else:
        print('Align_quads ran in '+str(t_end-t_start)+' sec')

    return aligned_images, aligned_arrays 