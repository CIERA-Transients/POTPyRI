import os 
import ccdproc
import time
import importlib
import numpy as np
import astroscrappy
import logging
import sys

# Needed for satellite trail detection - dependency can be removed later if 
# we turn off satellite masking or use some other algorithm (e.g., the Rubin one)
import acstools

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
from utilities import util

# Edited on 2024-07-29
__version__ = "1.3"

def remove_pv_distortion(header):

    done = False
    while not done:
        bad_key = False
        for key in header.keys():
            if key.startswith('PV'):
                if key in header.keys():
                    bad_key = True
                    del header[key]

        if not bad_key: done = True

    return(header)

def get_fieldcenter(images):

    ras=[] ; des = []
    for file in images:
        hdu = fits.open(file)
        w = WCS(hdu[0].header)

        data_shape = hdu[0].data.shape
        center_pix = (data_shape[0]/2., data_shape[1]/2.)

        coord = w.pixel_to_world(*center_pix)

        ras.append(coord.ra.degree)
        des.append(coord.dec.degree)

    mean_ra = np.mean(ras)
    mean_de = np.mean(des)

    return([mean_ra, mean_de])

def generate_wcs(tel, binn, fieldcenter, out_size):

    w = {'NAXES':2, 'NAXIS1': out_size, 'NAXIS2': out_size,
        'EQUINOX': 2000.0, 'LONPOLE': 180.0, 'LATPOLE': 0.0}

    pixscale = tel.pixscale(None)
    cdelt = pixscale * int(str(binn)[0])/3600.0

    w['CDELT1'] = -1.0 * cdelt
    w['CDELT2'] = cdelt

    w['CTYPE1']='RA---TAN'
    w['CTYPE2']='DEC--TAN'
    w['CUNIT1']='deg'
    w['CUNIT2']='deg'

    w['CRPIX1']=float(out_size)/2 + 0.5
    w['CRPIX2']=float(out_size)/2 + 0.5

    coord = util.parse_coord(fieldcenter[0], fieldcenter[1])

    w['CRVAL1']=coord.ra.degree
    w['CRVAL2']=coord.dec.degree

    w = WCS(w)

    return(w)

def align_images(reduced_files, paths, tel, binn, use_wcs=None, fieldcenter=None,
    out_size=None, log=None):

    solved_images = []
    aligned_images = []
    aligned_data = []
    for file in reduced_files:
        # Coarse WCS solution using astrometry.net
        success = solve_wcs.solve_astrometry(file, tel, shift_only=False, 
            log=log)
        if not success: continue

        # Fine WCS solution using Gaia DR3 point sources
        success = solve_wcs.align_to_gaia(file, tel, log=log)
        if not success: continue

        solved_images.append(file)

    # Determine what the value of use_wcs should be
    if use_wcs is None:
        if fieldcenter is None:
            fieldcenter = get_fieldcenter(solved_images)
            use_wcs = generate_wcs(tel, binn, fieldcenter, out_size)

    for file in solved_images:
        if log: log.info(f'Reprojecting {file}...')

        # Create a reprojection of the file in a common astrometric frame
        hdu = fits.open(file)
        header = remove_pv_distortion(hdu[0].header)
        wcs = WCS(hdu[0].header)

        if use_wcs is None:
            use_wcs = wcs

        image = CCDData(hdu[0].data, meta=hdu[0].header, wcs=wcs, 
            unit=u.electron)

        if out_size:
            target_shape = (out_size, out_size)
        else:
            target_shape = None

        reprojected = wcs_project(image, target_wcs=use_wcs, order='bilinear',
            target_shape=target_shape)
        # Get rid of mask in data
        reprojected.mask = None
        outfile = file.replace('.fits','_reproj.fits')
        reprojected.write(file.replace('.fits','_reproj.fits'), overwrite=True)

        aligned_images.append(outfile)
        aligned_data.append(reprojected)

    return(aligned_images, aligned_data)

def image_proc(image_data, tel, paths, proc=None, skip_skysub=False, 
    fieldcenter=None, out_size=None, satellites=True, log=None):

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

    # Load bias frame
    if tel.dark():
        if log: log.info('Loading master dark.')
        try:
            mdark = tel.load_dark(cal_path, amp, binn)
        except:
            if log: log.error(f'''No master dark found for this configuration, 
                skipping reduction for: {cal_type}''')
            return(None)
    else:
        mdark = None

    # Load flat frame
    if tel.flat():
        if log: log.info('Loading master flat.')
        master_flat = tel.get_mflat_name(cal_path, fil, amp, binn)
        if not os.path.exists(master_flat):
            if log: log.error(f'''No master flat present for filter {fil}, skipping 
                data reduction for {tar}. Check data before rerunning.''')
            return(None)
        mflat = tel.load_flat(master_flat)
    else:
        mflat = None
        
    t1 = time.time()
    if log: log.info(f'Processing data for {cal_type}')

    # Bias subtraction, gain correction, flat correction, and flat fielding
    files = image_data['File']
    staticmask = tel.static_mask(paths)
    processed = tel.process_science(files, fil, amp, binn, work_path,
        mbias=mbias, mflat=mflat, mdark=mdark, proc=proc, 
        staticmask=staticmask, skip_skysub=skip_skysub, log=log)

    # If masking satellite trails, this needs to be done before reprojection so
    # that the acstools algorithm accurately models the edges of the detector
    if satellites:
        mask_satellites(processed, files, log=log)
        
    t2 = time.time()
    if log: log.info(f'Data processed in {t2-t1} sec')

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

    if out_size is None:
        out_size = tel.out_size()

    if fieldcenter is not None:
        use_wcs = generate_wcs(tel, binn, fieldcenter, out_size)
    else:
        use_wcs = None

    aligned_images, aligned_data = align_images(reduced_files, paths, tel, binn,
        use_wcs=use_wcs, fieldcenter=fieldcenter, out_size=out_size, log=log)

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
        log.info(f'Stack made for {stackname}')
    else:
        print(f'Stack made for {stackname}')

    return(stackname)

def detrend_stack(stack):

    data = stack[0].data
    mask = stack[1].data.astype(bool)

    mean, med, stddev = sigma_clipped_stats(data, mask=mask, axis=1)
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
    exptimes = []
    for i,stk in enumerate(stacking_data):
        mask = masks[i].data.astype(bool)
        stacking_data[i].data[mask] = np.nan
        exptimes.append(float(tel.exptime(stk.header)))

    exptimes = np.array(exptimes)
    scale = 1./exptimes

    all_data = np.array([s.data for s in stacking_data])

    sci_med = combine(stacking_data, weights=weights, scale=scale,
        method=stack_method, mem_limit=64e9)

    return(sci_med)

def add_stack_mask(stack, stacking_data):

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

def mask_satellites(images, filenames, log=None):

    out_data = []
    for i,science_data in enumerate(images):

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

def create_mask(science_data, saturation, rdnoise, sigclip=3.5, 
    sigfrac=0.2, objlim=4.5, niter=6, outpath='', grow=0, satellites=True,
    log=None):

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
    
    mask_hdu = fits.PrimaryHDU(mask) #create mask Primary HDU
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

    return newdata, mask_hdu


def create_error(science_data, mask_data, rdnoise):

    hdu = fits.open(science_data)
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

    poisson = img_data[img_data < 0.]=0.
    error = np.sqrt(img_data + rms**2 + rdnoise)

    error_hdu = fits.PrimaryHDU(error) #create mask Primary HDU
    error_hdu.header['VER'] = (__version__, 
        'Version of image procedures used used.')
    error_hdu.header['USE'] = 'Error array for Poisson, read, and RMS noise'
    error_hdu.header['BUNIT'] = ('ELECTRONS', 'Units of the error array.')

    return(error_hdu)


if __name__=="__main__":

    t = ascii.read('/Users/ckilpatrick/Dropbox/Data/POTPyRI/test/Binospec/file_list.txt',
        format='fixed_width')
    mask = t['Target']=='GRB240809A_r_ep1'
    t = t[mask]

    t.sort('File')

    global tel
    tel = importlib.import_module('params.BINOSPEC')

    paths={}
    paths['red']='/Users/ckilpatrick/Dropbox/Data/POTPyRI/test/Binospec/red'
    paths['work']='/Users/ckilpatrick/Dropbox/Data/POTPyRI/test/Binospec/red/workspace'
    paths['cal']='/Users/ckilpatrick/Dropbox/Data/POTPyRI/test/Binospec/red/cals'
    paths['code']=os.path.dirname(sys.argv[0])

    image_proc(t, tel, paths, fieldcenter=['15:50:08.0920','-2:16:08.617'],
        out_size=5000, proc=None, log=None)

