#!/usr/bin/env python

"Function to create mask."
"Author: Kerry Paterson"

__version__ = "1.1" #last updated 24/05/2021

import os
import time
import numpy as np
from astropy.io import fits
from acstools.satdet import detsat, make_mask
import astroscrappy
import importlib
import tel_params

def create_mask(science_file,red,suffix,static_mask,source_mask,saturation,binning,rdnoise,sigclip,sigfrac,objlim,log):

    t_start = time.time()
    
    log.info('Running create_mask version: '+str(__version__))

    with fits.open(static_mask) as hdr:
        mask_bp = -~-hdr[0].data

    data = np.ascontiguousarray(red.data.astype(np.float32))
    
    log.info('Masking saturated pixels.')
    mask_sat = np.zeros((np.shape(data)[0],np.shape(data)[1])).astype(np.bool) #create empty mask
    astroscrappy.update_mask(data,mask_sat,saturation,True) #add saturated stars to mask
    mask_sat = mask_sat.astype(np.uint8) #set saturated star mask type
    mask_sat[data >= saturation] = 4 #set saturated pixel flag
    mask_sat[mask_sat == 1] = 8 #set connected to saturated pixel flag

    log.info('Cleaning and masking cosmic rays.')
    log.info('Using sigclip = %d, sigfrac = %.1f, objlim = %d'%(sigclip,sigfrac,objlim))
    mask_cr, clean = astroscrappy.detect_cosmics(data,inmask=(mask_bp+mask_sat+source_mask).astype(np.bool),sigclip=sigclip,sigfrac=sigfrac,readnoise=rdnoise,satlevel=saturation,objlim=objlim) #clean science image of cosmic rays and create cosmic ray mask
    mask_cr = mask_cr.astype(np.uint8) #set cosmic ray mask type
    mask_cr[mask_cr == 1] = 2 #set cosmic ray flag
    
    log.info('Maksing satellite trails.')
    satellite_fitting = False
    red_path = os.path.dirname(science_file.replace('/raw/','/red/'))
    binned_data = clean.reshape(int(np.shape(clean)[0]/binning[0]),binning[0],int(np.shape(clean)[1]/binning[1]),binning[1]).sum(3).sum(1) #bin data
    for j in range(3):
        fits.PrimaryHDU(binned_data).writeto(red_path+'/binned_mask.fits',overwrite=True) #write binned data to tmp file
        results, errors = detsat(red_path+'/binned_mask.fits', chips=[0], n_processes=1, buf=40, sigma=3, h_thresh=0.2) #detect sateliite trails
        trail_coords = results[(red_path+'/binned_mask.fits',0)] #create satellite trail if found
        if len(trail_coords) > 0: #continue if sateliite trail found
            trail_segment = trail_coords[0]
            try:
                mask_binned = make_mask(red_path+'/binned_mask.fits', 0, trail_segment, sublen=5, pad=0, sigma=5, subwidth=5000).astype(np.uint8) #create satellite trail mask
            except ValueError:
                log.error('Satellite trail could not be fitted for file '+science_file+' and is not included in the mask. ') #if error occurs, add comment
                break
            satellite_fitting = True
            binned_data[mask_binned == 1] = np.median(binned_data)
            try:
                open_old_mask = fits.open(red_path+'/old_mask.fits')
                old_mask = open_old_mask[0].data
                open_old_mask.close()
                mask_binned = old_mask+mask_binned
            except IOError:
                pass
            fits.writeto(red_path+'/old_mask.fits',mask_binned,overwrite=True)
        else:
            break
    os.remove(red_path+'/binned_mask.fits')
    if os.path.exists(red_path+'/old_mask.fits'):
        os.remove(red_path+'/old_mask.fits')
    if satellite_fitting == True:
        mask_sate = np.kron(mask_binned, np.ones((binning[0],binning[1]))).astype(np.uint8) #unbin mask
        mask_sate[mask_sate == 1] = 16 #set satellite trail flag
    else: #if no satellite trails are found, create empty mask
        mask_sate = (np.zeros([np.shape(clean)[0],np.shape(clean)[1]])).astype(np.uint8)

    mask = mask_bp+mask_cr+mask_sat+mask_sate #combine bad pixel, cosmic ray, saturated star and satellite trail masks
    mask_hdu = fits.PrimaryHDU(mask) #create mask Primary HDU
    mask_hdu.header['VER'] = (__version__, 'Version of create mask used.') #name of log file
    mask_hdu.header['USE'] = ('Complex mask using additive flags.', 'e.g. 6 = 2 + 4') #header comment
    mask_hdu.header['M-BP'] = (1, 'Value of masked bad pixels.')
    mask_hdu.header['M-BPNUM'] = (np.sum(mask & 1 == 1), 'Number of bad pixels.')
    mask_hdu.header['M-CR'] = (2, 'Value of masked cosmic ray pixels.')
    mask_hdu.header['M-CRNUM'] = (np.sum(mask & 2 == 2), 'Number of cosmic ray pixels.')
    mask_hdu.header['SATURATE'] = (saturation, 'Level of saturation.')
    mask_hdu.header['M-SP'] = (4, 'Value of masked saturated pixels.')
    mask_hdu.header['M-SPNUM'] = (np.sum(mask & 4 == 4), 'Number of saturated pixels.')
    mask_hdu.header['M-CSP'] = (8, 'Value of masked saturated-connected pixels.')
    mask_hdu.header['M-CSPNUM'] = (np.sum(mask & 8 == 8), 'Number of saturated-connected pixels.')
    mask_hdu.header['M-SAT'] = (16, 'Value of masked satellite trail pixels.')
    mask_hdu.header['M-SATNUM'] = (np.sum(mask & 16 == 16), 'Number of satellite trail pixels.')
    mask_hdu.writeto(science_file.replace('/raw/','/red/').replace('.fits',suffix),overwrite=True) #write mask to file
    log.info('Mask created.')

    t_end = time.time()
    log.info('Mask creation completed in '+str(t_end-t_start)+' sec')

    return clean, mask