#!/usr/bin/env python

"Function to assess image quality for stacking."
"Author: Kerry Paterson"

__version__ = "1.0" #last updated 20/04/2021

import time
import numpy as np
from astropy.io import fits, ascii
import importlib
import tel_params

def quality_check(aligned_images, aligned_data, telescope, log):

    t_start = time.time()

    log.info('Running quality_check version '+str(__version__))

    #import telescope parameter file
    try:
        tel = importlib.import_module('tel_params.'+telescope)
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)

    param_file = './Config/params'
    params = []
    with open(param_file) as f:
        for line in f:
            params.append(line.split()[0].strip())
    
    fwhm = []
    elong = []
    for f in aligned_images:
        log.info('Loading catalog for: '+f)
        cat = f.replace('.fits','.cat')
        table = ascii.read(cat, names=params, comment='#')
        table = table[(table['FLAGS']==0)&(table['EXT_NUMBER']==tel.wcs_extension()+1)]
        fwhm_image = np.median(table['FWHM_IMAGE'])*tel.pixscale()
        elong_image = np.median(table['ELONGATION'])
        fwhm.append(fwhm_image)
        elong.append(elong_image)
        log.info('FWHM of image in arcsec: '+str(fwhm_image))
        log.info('Elongation of image: '+str(elong_image))
    
    log.info('Calculating median parameters.')
    fwhm_med = np.median(fwhm)
    fwhm_std = np.std(fwhm)
    elong_med = np.median(elong)
    elong_std = np.std(elong)
    log.info('Median FWHM of images in arcsec = %.3f +/- %.3f'%(fwhm_med,fwhm_std))
    log.info('Median Elongation of images = %.3f +/- %.3f'%(elong_med,elong_std))

    log.info('Using 3 sigma cuts.')
    stacking_images = []
    stacking_arrays = []
    for i,f in enumerate(aligned_images):
        if (fwhm[i]>fwhm_med+3*fwhm_std) or (elong[i]>elong_med+3*elong_std):
            log.info('Removing '+f+' from stack due to bad quality.')
        else:
            stacking_images.append(f)
            stacking_arrays.append(aligned_data[i])

    log.info('Calculating mid-time of images.')
    time_list = []
    for image in stacking_arrays:
        time_list.append(tel.time_format(image.header))

    sorted_times = sorted(time_list)
    mid_time = sorted_times[0]+(sorted_times[-1]-sorted_times[0])/2

    t_end = time.time()
    log.info('Quailty check completed in '+str(t_end-t_start)+' sec')

    return stacking_arrays, mid_time