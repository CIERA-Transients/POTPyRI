import os 
import ccdproc
import time
import importlib
import numpy as np
import astroscrappy

import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.stats import SigmaClip
from astropy.nddata import CCDData

from photutils.background import Background2D
from photutils.background import MeanBackground
from photutils.segmentation import detect_threshold, detect_sources

# Internal dependencies
import align_quads
import quality_check

__version__ = "1.2"

def image_proc(image_data, tel, paths, proc=None, log=None):

    wavelength = tel.wavelength()

    red_path = paths['red']
    cal_path = paths['cal']
    work_path = paths['work']

    # All data in the table should be for one target
    assert np.all([image_data['Target']==image_data['Target'][0]])
    assert np.all([image_data['CalType']==image_data['CalType'][0]])

    stack = tel.stacked_image(image_data[0], red_path)
    
    cal_type = image_data['CalType'][0]
    target = image_data['Target'][0]
    fil = image_data['Filter'][0]
    amp = image_data['Amp'][0]
    binn = image_data['Binning'][0]

    # Load bias frame
    if tel.bias():
        if log: log.info('Loading master bias.')
        try:
            mbias = tel.load_bias(cal_path,amp,binn)
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
    masks = []
    processed = tel.process_science(files, fil, amp, binn, work_path,
        mbias=mbias, mflat=flat_data, proc=proc, log=log)
        
    t2 = time.time()
    if log: log.info(f'Data processed in {t2-t1} sec')

    # Background subtraction
    # IR image?
    if wavelength=='NIR':
        processed=background_subtraction(processed, masks, image_data, log=log)

    # Optical image?
    if wavelength=='OPT' and tel.fringe_correction(fil):
        processed=fringe_correction(stack, processed, masks, work_path, log=log)
            
    # Write reduced data to file and stack images with the same pointing
    if log: log.info('Writing out reduced data.')
    
    reduced_files = [os.path.join(work_path, 
        os.path.basename(sci).replace('.fits','_red.fits').replace('.gz','').replace('.bz2',''))
        for sci in image_data['File']]
    
    for j,process_data in enumerate(processed):
        process_data.header['FILENAME'] = os.path.basename(reduced_files[j])
        process_data.write(reduced_files[j], overwrite=True)

    if log: log.info('Aligning images.')
    static_mask = tel.static_mask(proc)
    aligned_images, aligned_data = align_quads.align_stars(reduced_files,
        tel, hdu=tel.wcs_extension(), mask=static_mask, log=log)
    
    if log: log.info('Checking qualty of images.')
    stacking_data, mid_time, total_time = quality_check.quality_check(
        aligned_images, aligned_data, tel, log=log)

    if log: log.info('Creating mask and error arrays.')
    masks = []
    errors = []
    for stack_img in stacking_data:
        clean, mask = create_mask(stack_img, static_mask,
            stack_img.header['SATURATE'], tel.rdnoise(stack_img.header),
            log=log, outpath=work_path)
        error = create_error(stack_img, mask, tel.rdnoise(stack_img.header))

        masks.append(mask)
        errors.append(error)

        # Create final image triplet before stacking with data, mask, and error
        hdulist = fits.HDUList()
        sci_hdu = fits.ImageHDU()
        sci_hdu.data = stack_img.data
        sci_hdu.header = stack_img.header
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

        filename = hdulist[0].header['FILENAME'].replace('_red.fits',
            '_data.fits')
        fullfilename = os.path.join(work_path, filename)

        if log: 
            log.info(f'Writing out all file data: {fullfilename}')
        else:
            print(f'Writing out all file data: {fullfilename}')
        hdulist.writeto(fullfilename, overwrite=True, output_verify='silentfix')
    
    if log: log.info('Creating median stack.') 
    sci_med = stack_data(stacking_data, masks, errors)
    sci_med = sci_med.to_hdu()

    sci_med = add_stack_mask(sci_med, stacking_data)

    sci_med[0].data = sci_med[0].data.astype(np.float32)
    sci_med[1].data = sci_med[1].data.astype(np.uint8)
    sci_med[2].data = sci_med[2].data.astype(np.float32)
    sci_med[0].header['BITPIX'] = -32
    sci_med[1].header['BITPIX'] = 8
    sci_med[2].header['BITPIX'] = -32

    sci_med[0].header['MJD-OBS'] = (mid_time, 
        'Mid-MJD of the observation sequence calculated using DATE-OBS.')
    sci_med[0].header['EXPTIME'] = (total_time, 
        'Effective expsoure tiime for the stack in seconds.')
    sci_med[0].header['EXPTOT'] = (total_time, 
        'Total exposure time of stack in seconds')
    sci_med[0].header['GAIN'] = (len(stacking_data), 
        'Effecetive gain for stack.')

    # Calculate read noise
    rdnoise = tel.rdnoise(sci_med[0].header)/np.sqrt(len(stacking_data))
    sci_med[0].header['RDNOISE'] = (rdnoise, 'Readnoise of stack.')
    
    sci_med[0].header['NFILES'] = (len(stacking_data), 
        'Number of images in stack')
    sci_med[0].header['FILTER'] = fil
    sci_med[0].header['OBSTYPE'] = 'OBJECT'

    sci_med = tel.edit_stack_headers(sci_med)
    
    sci_med.writeto(stack, overwrite=True, output_verify='silentfix')
                
    if log: 
        log.info(f'Median stack made for {stack}')
    else:
        print(f'Median stack made for {stack}')

    return(stack)

def background_subtraction(processed, masks, image_data, log=None):

    t1 = time.time()
    if log: log.info('NIR data, creating NIR sky maps.')
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
        sky = ccdproc.combine(sky_masked_data,method='median',
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

def fringe_correction(stack, processed, masks, work_path, log=None):

    t1 = time.time()
    dimen = len(stack)
    for m in range(dimen):
        fringe_data = []
        if dimen == 1:
            suffix = ['.fits']
            for k,n in enumerate(processed):
                bkg = Background2D(n, (20, 20), filter_size=(3, 3),
                    sigma_clip=SigmaClip(sigma=3), 
                    bkg_estimator=MeanBackground(), 
                    mask=masks[k], exclude_percentile=80)
                masked = np.array(n)
                masked[masks[k]] = bkg.background[masks[k]]
                fringe_data.append(CCDData(masked,unit=u.electron/u.second))
            
            fringe_map = ccdproc.combine(fringe_data,method='median',
                sigma_clip=True,sigma_clip_func=np.ma.median,mask=masks)
            fringe_map.write(os.path.join(work_path, 
                'fringe_map_'+fil+'_'+amp+'_'+binn+suffix[m]),overwrite=True)
            for j,n in enumerate(processed):
                processed[j] = n.subtract(fringe_map,
                    propagate_uncertainties=True, handle_meta='first_found')
        else:
            suffix = [s.replace('_red','') for s in tel.suffix()]
            for k,n in enumerate(processed[m]):
                bkg = Background2D(n, (20, 20), filter_size=(3, 3),
                    sigma_clip=SigmaClip(sigma=3), 
                    bkg_estimator=MeanBackground(), 
                    mask=masks[m][k], exclude_percentile=80)
                masked = np.array(n)
                masked[masks[m][k]] = bkg.background[masks[m][k]]
                fringe_data.append(CCDData(masked,unit=u.electron/u.second))
            
            fringe_map = ccdproc.combine(fringe_data,method='median',
                sigma_clip=True,sigma_clip_func=np.ma.median,mask=masks)
            fringe_map.write(os.path.join(work_path, 
                'fringe_map_'+fil+'_'+amp+'_'+binn+suffix[m]),overwrite=True)
            for j,n in enumerate(processed[m]):
                processed[m][j] = n.subtract(fringe_map,
                    propagate_uncertainties=True,handle_meta='first_found')
    
    t2 = time.time()
    if log: log.info(f'Fringe correction complete and subtracted in {t2-t1} sec')

    return(processed)

def stack_data(stacking_data, masks, errors):

    weights = []
    for i,error in enumerate(errors):
        ivar = 1./error.data**2
        mask = masks[i]
        ivar[mask.data.astype(bool)]=0.
        weights.append(ivar)
    weights=np.array(weights)

    sci_med = ccdproc.combine(stacking_data, weights=weights, 
        method='median', sigma_clip=True, sigma_clip_func=np.ma.median)

    return(sci_med)

def add_stack_mask(stack, stacking_data, sat_grow=3):

    data = stack[0].data
    mask = stack[1].data

    # Keep track of original masked values
    bp_mask = mask > 0

    # Reset mask
    mask = np.zeros(mask.shape).astype(np.uint8)

    # Add bad pixels back to mask
    mask[(data==0.0) | bp_mask] = 1

    # Add saturation to mask
    sat = np.min([d.header['SATURATE'] for d in stacking_data])
    stack[0].header['SATURATE'] = sat
    sat_mask = data >= sat
    mask[sat_mask] = 4

    stack[0].data = data
    stack[1].data = mask

    return(stack)

def create_mask(science_data, static_mask, saturation, rdnoise, sigclip=2.5, 
    sigfrac=0.1, objlim=2.5, niter=6, outpath='', log=None):

    t_start = time.time()
    
    if log: log.info(f'Running create_mask version: {__version__}')

    data = np.ascontiguousarray(science_data.data.astype(np.float32))

    # Astroscrappy requires added sky background, so add this value back
    data = data + science_data.header['SKYBKG']
    
    if log: log.info('Masking saturated pixels.')

    mask_bp = data==0.0

    mask_sat = np.zeros(data.shape).astype(bool) # create empty mask
    astroscrappy.update_mask(data, mask_sat, saturation, True) #add saturated stars to mask
    mask_sat = mask_sat.astype(np.uint8) #set saturated star mask type
    mask_sat[data >= saturation] = 4 #set saturated pixel flag
    mask_sat[mask_sat == 1] = 8 #set connected to saturated pixel flag

    if log: log.info('Cleaning and masking cosmic rays.')
    if log: log.info(f'Using sigclip={sigclip}, sigfrac={sigfrac}, objlim={objlim}')
    #clean science image of cosmic rays and create cosmic ray mask
    inmask = (mask_sat+mask_bp).astype(bool)
    mask_cr, clean = astroscrappy.detect_cosmics(data,
        inmask=inmask, sigclip=sigclip, sigfrac=sigfrac, readnoise=rdnoise, 
        satlevel=saturation, objlim=objlim, niter=niter)

    mask_cr = mask_cr.astype(np.uint8) #set cosmic ray mask type
    mask_cr[mask_cr == 1] = 2 #set cosmic ray flag

     #combine bad pixel, cosmic ray, saturated star and satellite trail masks
    mask = mask_bp+mask_sat+mask_cr
    
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
    
    if log: log.info('Mask created.')

    t_end = time.time()
    if log: log.info(f'Mask creation completed in {t_end-t_start} sec')

    return clean, mask_hdu


def create_error(science_data, mask_data, rdnoise):

    img_data = science_data.data.astype(np.float32)
    mask = mask_data.data.astype(bool)

    rms = 0.5 * (
        np.percentile(img_data[~mask], 84.13)
        - np.percentile(img_data[~mask], 15.86)
    )

    poisson = img_data[img_data < 0.]=0.
    error = np.sqrt(img_data + rms**2 + rdnoise)

    error_hdu = fits.PrimaryHDU(error) #create mask Primary HDU
    error_hdu.header['VER'] = (__version__, 'Version of create mask used.') #name of log file
    error_hdu.header['USE'] = ('Error array based on Poisson, read, and RMS noise', '') #header comment
    error_hdu.header['BUNIT'] = ('ELECTRONS', 'Units of the error array.')

    return(error_hdu)


if __name__=="__main__":

    t = Table.read('/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/file_list.txt',
        format='ascii.ecsv')
    mask = t['Target']=='R155_host'
    t = t[mask]

    t.sort('File')

    global tel
    tel = importlib.import_module('params.LRIS')

    red_path='/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/red'

    image_proc(t, tel, red_path, proc=None, log=None)

