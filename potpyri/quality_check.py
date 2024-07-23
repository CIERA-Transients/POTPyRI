#!/usr/bin/env python

"Function to assess image quality for stacking."
"Author: Kerry Paterson"

__version__ = "1.7" #last updated 11/01/2022

import time
import numpy as np
from astropy.io import fits, ascii
import importlib
import solve_wcs
import params
import sys
import os

def quality_check(aligned_images, aligned_data, telescope, log,
    nsource_align=100, nsource_align_min=20):

    t_start = time.time()

    log.info('Running quality_check version '+str(__version__))

    #import telescope parameter file
    try:
        tel = importlib.import_module(f'params.{telescope}')
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)   

    fwhm = []
    elong = []
    stars = []
    i = 0
    for f in aligned_images:
        log.info('Loading catalog for: '+f)
        cat = f.replace('.fits','.cat')
        config_dir = os.path.abspath(os.path.join('potpyri','config'))
        table = solve_wcs.run_sextractor(f, cat, tel, 
            sex_config_dir=config_dir, log=log)
        if len(table)==0:
            log.info('No sources found by SExtractor.')
            fwhm.append(0)
            elong.append(0)
            stars.append(0)
        elif len(table)<nsource_align_min:
            log.info('Not enough sources found by SExtractor.')
            fwhm.append(0)
            elong.append(0)
            stars.append(0)
        else:
            i += 1
            if i==1:
                table = table[(table['FLAGS']==0)&(table['EXT_NUMBER']==tel.wcs_extension()+1)&(table['MAGERR_BEST']!=99)]
                table.sort('MAG_BEST')
                nsource = np.min([len(table), nsource_align])
                x_pos = table['XWIN_IMAGE'][0:nsource]
                y_pos = table['YWIN_IMAGE'][0:nsource]
                match = nsource
            else:
                matches = []
                nsource = np.min([len(table), nsource_align])
                for j in range(nsource):
                    x_pos_search = table['XWIN_IMAGE'][(table['XWIN_IMAGE']>x_pos[j]-5)&(table['XWIN_IMAGE']<x_pos[j]+5)&(table['YWIN_IMAGE']>y_pos[j]-5)&(table['YWIN_IMAGE']<y_pos[j]+5)]
                    y_pos_search = table['YWIN_IMAGE'][(table['XWIN_IMAGE']>x_pos[j]-5)&(table['XWIN_IMAGE']<x_pos[j]+5)&(table['YWIN_IMAGE']>y_pos[j]-5)&(table['YWIN_IMAGE']<y_pos[j]+5)]
                    if len(x_pos_search)!=0 and len(y_pos_search)!=0:
                        matches.append(x_pos_search[0])
                match = len(matches)
            fwhm_image = np.median(table['FWHM_IMAGE'])*tel.pixscale()
            elong_image = np.median(table['ELONGATION'])
            fwhm.append(fwhm_image)
            elong.append(elong_image)
            stars.append(match)
            log.info('FWHM of image in arcsec: '+str(fwhm_image))
            log.info('Elongation of image: '+str(elong_image))
    
    log.info('Calculating median parameters.')
    fwhm_med = np.median(fwhm)
    fwhm_std = np.std(fwhm)
    log.info('Median FWHM of images in arcsec = %.3f +/- %.3f'%(fwhm_med,fwhm_std))
    elong_med = np.median(elong)
    elong_std = np.std(elong)
    log.info('Median Elongation of images = %.3f +/- %.3f'%(elong_med,elong_std))

    log.info('Using 3 sigma cuts to calculate the number of bad images.')
    fwhm = np.array(fwhm)
    elong = np.array(elong)
    fwhm_3sigma_cut = fwhm[fwhm>fwhm_med+3*fwhm_std]
    elong_3sigma_cut = elong[elong>elong_med+3*elong_std]
    n_fwhm = len(fwhm_3sigma_cut)
    n_elong = len(elong_3sigma_cut)
    log.info('Number of images with FWHM > 3sigma = %d'%n_fwhm)
    log.info('Number of images with elongation > 3sigma = %d'%n_elong)
    ten_percent = int(len(aligned_images)*0.1)
    if n_fwhm > ten_percent or n_elong > ten_percent:
        log.info('More than 10 percent of the images have bad quality.')
        log.info('Restricting cuts to 10 percent.')
        fwhm_sorted = sorted(fwhm_3sigma_cut,reverse=True)
        fwhm_cut = fwhm_sorted[ten_percent-1]
        elong_sorted = sorted(elong_3sigma_cut,reverse=True)
        elong_cut = elong_sorted[ten_percent-1]
        log.info('Using new cuts of %.2f sigma for FWHM and %.2f sigma for elongation'
            %((fwhm_cut-fwhm_med)/fwhm_std,(elong_cut-elong_med)/elong_std))
    else:
        fwhm_cut = fwhm_med+3*fwhm_std
        elong_cut = elong_med+3*elong_std
    stacking_images = []
    stacking_arrays = []
    log.info('Checking alignment of images.')
    for i,f in enumerate(aligned_images):
        if (fwhm[i]>fwhm_cut) or (elong[i]>elong_cut):
            log.info('Removing '+f+' from stack due to bad quality.')
        elif stars[i] < 5:
            log.info('Removing '+f+' from stack due to bad alignment.')
        else:
            stacking_images.append(f)
            stacking_arrays.append(aligned_data[i])

    log.info('Calculating mid-time and total exposure time of stack.')
    time_list = []
    total_time = 0
    for image in stacking_arrays:
        time_list.append(tel.time_format(image.header))
        total_time += tel.exptime(image.header)

    sorted_times = sorted(time_list)
    if len(sorted_times) == 1:
        mid_time = sorted_times[0]
    else:
        mid_time = sorted_times[0]+(sorted_times[-1]-sorted_times[0])/2

    t_end = time.time()
    log.info('Quailty check completed in '+str(t_end-t_start)+' sec')

    return stacking_arrays, mid_time, total_time
