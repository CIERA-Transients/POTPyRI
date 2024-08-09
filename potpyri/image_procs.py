import os 
import ccdproc
import time
import importlib
import numpy as np
import astroscrappy
import logging

import astropy.units as u
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
from astropy.nddata import CCDData
from astropy.wcs import WCS

from ccdproc import wcs_project
from ccdproc import CCDData
from ccdproc import combine
from ccdproc import cosmicray_lacosmic

from photutils.background import Background2D
from photutils.background import MeanBackground
from photutils.segmentation import detect_threshold, detect_sources

# Internal dependencies
import solve_wcs

# Edited on 2024-07-29
__version__ = "1.3"

def align_images(reduced_files, paths, tel, use_wcs=None, log=None):

    aligned_images = []
    aligned_data = []
    for file in reduced_files:
        success = solve_wcs.solve_astrometry(paths['work'], file, 
            tel, log=log)
        if success:
            hdu = fits.open(file)
            wcs = WCS(hdu[0].header)
            image = CCDData(hdu[0].data, meta=hdu[0].header,
                wcs=wcs, unit=u.electron)
            
            if not use_wcs:
                use_wcs = wcs

            reprojected = wcs_project(image, use_wcs, order='bilinear')
            outfile = file.replace('.fits','_reproj.fits')
            reprojected.write(file.replace('.fits','_reproj.fits'),
                overwrite=True)

            aligned_images.append(outfile)
            aligned_data.append(reprojected)

    return(aligned_images, aligned_data)

def image_proc(image_data, tel, paths, proc=None, log=None):

    wavelength = tel.wavelength()

    red_path = paths['red']
    cal_path = paths['cal']
    work_path = paths['work']

    # All data in the table should be for one target
    assert np.all([image_data['Target']==image_data['Target'][0]])
    assert np.all([image_data['CalType']==image_data['CalType'][0]])
    
    cal_type = image_data['CalType'][0]
    target = image_data['Target'][0]
    fil = image_data['Filter'][0]
    amp = image_data['Amp'][0]
    binn = image_data['Binning'][0]

    # Load bias frame
    if tel.bias():
        if log: log.info('Loading master bias.')
        try:
            mbias = tel.load_bias(cal_path, amp, binn)
        except:
            if log: log.error(f'''No master bias found for this configuration, 
                skipping reduction for: {cal_type}''')
            return(None)
    else:
        mbias = None

    # Load flat frame
    if log: log.info('Loading master flat.')
    master_flat = tel.get_mflat_name(cal_path, fil, amp, binn)
    if not os.path.exists(master_flat):
        if log: log.error(f'''No master flat present for filter {fil}, skipping 
            data reduction for {tar}. Check data before rerunning.''')
        return(None)
    flat_data = tel.load_flat(master_flat)
        
    t1 = time.time()
    if log: log.info(f'Processing data for {cal_type}')

    # Bias subtraction, gain correction, flat correction, and flat fielding
    files = image_data['File']
    staticmask = tel.static_mask(paths)
    processed = tel.process_science(files, fil, amp, binn, work_path,
        mbias=mbias, mflat=flat_data, proc=proc, staticmask=staticmask, log=log)
        
    t2 = time.time()
    if log: log.info(f'Data processed in {t2-t1} sec')

    # Background subtraction
    # IR image?
    if wavelength=='NIR':
        processed=sky_background_subtraction(processed, image_data, log=log)
            
    # Write reduced data to file and stack images with the same pointing
    if log: log.info('Writing out reduced data.')
    
    reduced_files = [os.path.join(work_path, tel.get_base_science_name(sci))
        for sci in image_data['File']]
    
    for j,process_data in enumerate(processed):
        process_data.header['FILENAME'] = os.path.basename(reduced_files[j])
        process_data.write(reduced_files[j], overwrite=True)

    if log: log.info('Aligning images.')

    # Sort files so deepest exposure is first
    exptimes = []
    for file in reduced_files:
        hdu = fits.open(file)
        exptimes.append(tel.exptime(hdu[0].header))
    idx = np.flip(np.argsort(exptimes))
    reduced_files = np.array(reduced_files)[idx]

    aligned_images, aligned_data = align_images(reduced_files, paths, tel, 
        log=log)

    if log: log.info('Creating mask and error arrays.')
    masks = []
    errors = []
    for stack_img in aligned_images:
        hdu = fits.open(stack_img, mode='readonly')
        clean, mask = create_mask(stack_img,
            hdu[0].header['SATURATE'], tel.rdnoise(hdu[0].header),
            log=log, outpath=work_path)
        error = create_error(stack_img, mask, tel.rdnoise(hdu[0].header))

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

        filename = hdulist[0].header['FILENAME'].replace('.fits',
            '_data.fits')
        fullfilename = os.path.join(work_path, filename)

        if log: 
            log.info(f'Writing out all file data: {fullfilename}')
        else:
            print(f'Writing out all file data: {fullfilename}')
        hdulist.writeto(fullfilename, overwrite=True, output_verify='silentfix')
    
    if log: log.info('Creating median stack.') 
    sci_med = stack_data(aligned_data, tel, masks, errors, log=log)
    sci_med = sci_med.to_hdu()

    sci_med = add_stack_mask(sci_med, aligned_data)

    if tel.detrend(sci_med[0].header):
        sci_med = detrend_stack(sci_med)

    sci_med[0].data = sci_med[0].data.astype(np.float32)
    sci_med[1].data = sci_med[1].data.astype(np.uint8)
    sci_med[2].data = sci_med[2].data.astype(np.float32)
    sci_med[0].header['BITPIX'] = -32
    sci_med[1].header['BITPIX'] = 8
    sci_med[2].header['BITPIX'] = -32

    # Get time parameters from aligned data
    mid_time = np.average([tel.time_format(d.header) for d in aligned_data])
    exptimes = np.array([tel.exptime(d.header) for d in aligned_data])
    eff_time = np.sum(exptimes**2)/np.sum(exptimes)
    total_time = np.sum(exptimes)

    sci_med[0].header['MJD-OBS'] = (mid_time, 
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
    rdnoise = tel.rdnoise(sci_med[0].header)/np.sqrt(len(aligned_data))
    sci_med[0].header['RDNOISE'] = (rdnoise, 'Readnoise of stack.')
    
    sci_med[0].header['NFILES'] = (len(aligned_data), 
        'Number of images in stack')
    sci_med[0].header['FILTER'] = fil
    sci_med[0].header['OBSTYPE'] = 'OBJECT'

    sci_med = tel.edit_stack_headers(sci_med)
    
    # Generate stack name and write out
    stackname = tel.stacked_image(image_data[0], red_path)
    sci_med.writeto(stackname, overwrite=True, output_verify='silentfix')
                
    if log: 
        log.info(f'Median stack made for {stackname}')
    else:
        print(f'Median stack made for {stackname}')

    return(stackname)

def sky_background_subtraction(processed, image_data, log=None):

    t1 = time.time()
    if log: 
        log.info('NIR data, creating NIR sky maps.')

    for j,n in enumerate(processed):
        
        if log: log.info(f'Creating map for file {j+1}/{len(processed)}')

        time_diff = sorted([(abs(image_data['Time']-n2),k) 
            for k,n2 in enumerate(image_data['Time'])])

        sky_list = [files[k] for _,k in time_diff[0:5]]
        sky_data = [processed[k] for _,k in time_diff[0:5]]
        sky_mask = [masks[k] for _,k in time_diff[0:5]]
        sky_masked_data = []
                
        for k in range(len(sky_data)):
            # Make a background array from sky_data
            bkg = Background2D(sky_data[k], (20, 20), 
                filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), 
                bkg_estimator=MeanBackground(), mask=sky_mask[k], 
                exclude_percentile=80)

            # Do masking on sky
            masked = np.array(sky_data[k])
            masked[sky_mask[k]] = bkg.background[sky_mask[k]]
            sky_masked_data.append(CCDData(masked, unit=u.electron/u.second))

        # Generate sky output file and edit headers
        sky_hdu = fits.PrimaryHDU()
        sky_hdu.header['FILE'] = (os.path.basename(image_data['File'][j]), 
            'NIR sky flat for file.')
                
        for k,m in enumerate(sky_list):
            sky_hdu.header['FILE'+str(k+1)] = (os.path.basename(str(m)), 
                'Name of file used in creation of sky.')
                
        # Image combination procedure
        sky = combine(sky_masked_data,method='median',
            sigma_clip=True,sigma_clip_func=np.ma.median,mask=sky_mask)
                
        # Apply header
        sky.header = sky_hdu.header

        # Write out file to 
        baseskyfile = os.path.basename(image_data['File'][j])
        baseskyfile = baseskyfile.replace('.fits','_sky.fits')
        baseskyfile = baseskyfile.replace('.gz','')
        baseskyfile = baseskyfile.replace('.bz2','')
        sky.write(os.path.join(work_path, baseskyfile), overwrite=True)
                
        # Process files by subtracting sky data
        processed[j] = n.subtract(sky, propagate_uncertainties=True,
            handle_meta='first_found')

    t2 = time.time()
    if log: log.info(f'Sky maps complete and subtracted in {t2-t1} sec')

    return(processed)

def detrend_stack(stack):

    data = stack[0].data
    mask = stack[1].data.astype(bool)

    mean, med, stddev = sigma_clipped_stats(data, mask=mask, axis=1)
    med[np.isnan(med)]=0.
    data = data - med[:,None]

    row_med = np.nanmedian(med)

    mean, med, stddev = sigma_clipped_stats(data, mask=mask, axis=0)
    data = data - med[None,:]

    col_med = np.nanmedian(med)

    # Reapply mask to data
    data[mask] = 0.0

    stack[0].data = data
    stack[0].header['SATURATE'] = stack[0].header['SATURATE'] - (row_med+col_med)

    return(stack)

def stack_data(stacking_data, tel, masks, errors, log=None):

    stack_method = tel.stack_method(stacking_data[0].header)
    if not stack_method:
        if log:
            log.critical('Could not get stacking method for these images')
            logging.shutdown()
            sys.exit(-1)
        else:
            raise Exception('Could not get stacking method for these images')

    weights = []
    for i,error in enumerate(errors):
        ivar = 1./error.data**2
        mask = masks[i]
        ivar[mask.data.astype(bool)]=0.
        weights.append(ivar)
    weights=np.array(weights)

    new_data = []
    for i,stk in enumerate(stacking_data):
        mask = masks[i].data.astype(bool)
        stacking_data[i].data[mask] = np.nan

    all_data = np.array([s.data for s in stacking_data])

    sci_med = combine(stacking_data, weights=weights, 
        method=stack_method)

    return(sci_med)

def add_stack_mask(stack, stacking_data, sat_grow=3):

    data = stack[0].data
    mask = stack[1].data
    error = stack[2].data

    # Keep track of original masked values
    bp_mask = (mask > 0) | (data==0.0) | np.isnan(data) | np.isnan(mask) |\
        np.isnan(error) | (error==0.0)

    # Reset mask
    mask = np.zeros(mask.shape).astype(np.uint8)

    # Add bad pixels back to mask
    mask[(data==0.0) | bp_mask] = 1

    # Add saturation to mask
    sat = np.min([d.header['SATURATE'] for d in stacking_data])
    stack[0].header['SATURATE'] = sat
    sat_mask = data >= sat
    mask[sat_mask] = 4

    data[mask>0] = 0.0
    error[mask>0] = 0.0

    stack[0].data = data
    stack[1].data = mask
    stack[2].data = error

    return(stack)

def create_mask(science_data, saturation, rdnoise, sigclip=3.5, 
    sigfrac=0.2, objlim=4.5, niter=6, outpath='', log=None):

    t_start = time.time()
    
    if log:
        log.info(f'Running astroscrappy on {science_data}')
    else:
        print(f'Running astroscrappy on {science_data}')

    hdu = fits.open(science_data)
    data = np.ascontiguousarray(hdu[0].data.astype(np.float32))

    # Astroscrappy requires added sky background, so add this value back
    data = data + hdu[0].header['SKYBKG']
    
    if log: log.info('Masking saturated pixels.')

    mask_bp = (data==0.0) | np.isnan(data)

    mask_sat = np.zeros(data.shape).astype(bool) # create empty mask
    mask_sat = mask_sat.astype(np.uint8) #set saturated star mask type
    mask_sat[data >= saturation] = 4 #set saturated pixel flag

    if log: log.info('Cleaning and masking cosmic rays.')
    if log: log.info(f'Using sigclip={sigclip}, sigfrac={sigfrac}, objlim={objlim}')
    #clean science image of cosmic rays and create cosmic ray mask
    inmask = (mask_sat+mask_bp).astype(bool)

    scidata = CCDData(data, unit=u.electron, mask=inmask, 
        wcs=WCS(hdu[0].header), meta=hdu[0].header)

    newdata, mask_cr = cosmicray_lacosmic(data,
        readnoise=rdnoise, satlevel=saturation, verbose=True,
        sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=niter)

    mask_cr = mask_cr.astype(np.uint8) #set cosmic ray mask type
    mask_cr[mask_cr == 1] = 2 #set cosmic ray flag

     #combine bad pixel, cosmic ray, saturated star and satellite trail masks
    mask = mask_bp+mask_sat+mask_cr

    # Get number of bad, cosmic-ray flagged, and saturated pixels
    nbad = np.sum(mask & 1 == 1)
    ncr = np.sum(mask & 2 == 2)
    nsat = np.sum(mask & 4 == 4)
    
    mask_hdu = fits.PrimaryHDU(mask) #create mask Primary HDU
    mask_hdu.header['VER'] = (__version__, 
        'Version of image procedures used.') #name of log file
    mask_hdu.header['USE'] = 'Complex mask using additive flags.' #header comment
    mask_hdu.header['M-BP'] = (1, 'Value of masked bad pixels.')
    mask_hdu.header['M-BPNUM'] = (nbad, 'Number of bad pixels.')
    mask_hdu.header['M-CR'] = (2, 'Value of masked cosmic ray pixels.')
    mask_hdu.header['M-CRNUM'] = (ncr, 'Number of cosmic ray pixels.')
    mask_hdu.header['SATURATE'] = (saturation, 'Level of saturation.')
    mask_hdu.header['M-SP'] = (4, 'Value of masked saturated pixels.')
    mask_hdu.header['M-SPNUM'] = (nsat, 'Number of saturated pixels.')
    
    if log: 
        log.info('Mask created.')
        log.info(f'{nbad} bad, {ncr} cosmic-ray, and {nsat} saturated pixels')
    else:
        print('Mask created.')
        print(f'{nbad} bad, {ncr} cosmic-ray, and {nsat} saturated pixels')

    t_end = time.time()
    if log: 
        log.info(f'Mask creation completed in {t_end-t_start} sec')
    else:
        print(f'Mask creation completed in {t_end-t_start} sec')

    return newdata, mask_hdu


def create_error(science_data, mask_data, rdnoise):

    hdu = fits.open(science_data)
    img_data = hdu[0].data.astype(np.float32)
    mask = mask_data.data.astype(bool)

    rms = 0.5 * (
        np.percentile(img_data[~mask], 84.13)
        - np.percentile(img_data[~mask], 15.86)
    )

    poisson = img_data[img_data < 0.]=0.
    error = np.sqrt(img_data + rms**2 + rdnoise)

    error_hdu = fits.PrimaryHDU(error) #create mask Primary HDU
    error_hdu.header['VER'] = (__version__, 
        'Version of image procedures used used.')
    error_hdu.header['USE'] = 'Error array for Poisson, read, and RMS noise'
    error_hdu.header['BUNIT'] = ('ELECTRONS', 'Units of the error array.')

    return(error_hdu)


if __name__=="__main__":

    t = ascii.read('/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/file_list.txt',
        format='fixed_width')
    mask = t['Target']=='R155_host'
    t = t[mask]

    t.sort('File')

    global tel
    tel = importlib.import_module('params.LRIS')

    paths={}
    paths['red']='/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/red'
    paths['work']='/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/red/workspace'
    paths['cal']='/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/red/cals'

    image_proc(t, tel, paths, proc=None, log=None)

