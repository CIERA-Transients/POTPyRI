"Functions for calibrating, masking, creating error images, and stacking images."
"Authors: Charlie Kilpatrick"

# Last updated on 09/29/2024
__version__ = "2.1"

import os 
import time
import importlib
import numpy as np
import logging
import sys
import copy

# Needed for satellite trail detection - dependency can be removed later if 
# we turn off satellite masking or use some other algorithm (e.g., the Rubin one)
import acstools

import astropy.units as u
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS

from ccdproc import wcs_project
from ccdproc import CCDData
from ccdproc import combine
from ccdproc import cosmicray_lacosmic

# Internal dependencies
from potpyri.stages import solve_wcs
from potpyri.utils import utilities

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
        center_pix = (data_shape[1]/2., data_shape[0]/2.)

        coord = w.pixel_to_world(*center_pix)

        ras.append(coord.ra.degree)
        des.append(coord.dec.degree)

    mean_ra = np.mean(ras)
    mean_de = np.mean(des)

    return([mean_ra, mean_de])

def generate_wcs(tel, binn, fieldcenter, out_size):

    w = {'NAXES':2, 'NAXIS1': out_size, 'NAXIS2': out_size,
        'EQUINOX': 2000.0, 'LONPOLE': 180.0, 'LATPOLE': 0.0}

    pixscale = tel.get_pixscale()
    cdelt = pixscale * int(str(binn)[0])/3600.0

    w['CDELT1'] = -1.0 * cdelt
    w['CDELT2'] = cdelt

    w['CTYPE1']='RA---TAN'
    w['CTYPE2']='DEC--TAN'
    w['CUNIT1']='deg'
    w['CUNIT2']='deg'

    w['CRPIX1']=float(out_size)/2 + 0.5
    w['CRPIX2']=float(out_size)/2 + 0.5

    coord = utilities.parse_coord(fieldcenter[0], fieldcenter[1])

    w['CRVAL1']=coord.ra.degree
    w['CRVAL2']=coord.dec.degree

    w = WCS(w)

    return(w)

def align_images(reduced_files, paths, tel, binn, use_wcs=None, fieldcenter=None,
    out_size=None, skip_gaia=False, log=None):

    solved_images = []
    aligned_images = []
    aligned_data = []
    for file in reduced_files:
        # Coarse WCS solution using astrometry.net
        success = solve_wcs.solve_astrometry(file, tel, binn, paths, 
            shift_only=False, log=log)
        if not success: continue

        # Fine WCS solution using Gaia DR3 point sources
        if skip_gaia:
            hdu = fits.open(file)
            hdu[0].header['RADISP']=0.0
            hdu[0].header['DEDISP']=0.0
            hdu.writeto(file, overwrite=True)
        else:
            success = solve_wcs.align_to_gaia(file, tel, log=log)
            if not success: continue

        solved_images.append(file)

    # Reject images from stack where either RADISP or DEDISP is >5-sigma outlier
    if len(solved_images)>2:
        radisp = [] ; dedisp = []
        for file in solved_images:
            hdu = fits.open(file, mode='readonly')
            radisp.append(hdu[0].header['RADISP'])
            dedisp.append(hdu[0].header['DEDISP'])

        radisp = np.array(radisp) ; dedisp = np.array(dedisp)
        ramean, ramedian, rastddev = sigma_clipped_stats(radisp)
        demean, demedian, destddev = sigma_clipped_stats(dedisp)

        if log: log.info(f'Median dispersion in R.A.={ramedian}')
        if log: log.info(f'Median dispersion in Decl.={demedian}')

        mask = (radisp-ramedian <= 5 * rastddev) &\
               (dedisp-demedian <= 5 * destddev)

        if log:
            log.info('Rejecting the following images for high astrometric dispersion:')
            for i,m in enumerate(mask):
                if not m: log.info(solved_images[i])

        solved_images = np.array(solved_images)[mask]


    if len(solved_images)==0:
        return(None, None)

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

def image_proc(image_data, tel, paths, skip_skysub=False, 
    fieldcenter=None, out_size=None, satellites=True, cosmic_ray=True,
    skip_gaia=False, log=None):

    wavelength = tel.wavelength

    red_path = paths['red']
    work_path = paths['work']

    # All data in the table should be for one target
    assert np.all([image_data['Target']==image_data['Target'][0]])
    assert np.all([image_data['TargType']==image_data['TargType'][0]])
    
    cal_type = image_data['TargType'][0]
    target = image_data['Target'][0]
    fil = image_data['Filter'][0]
    amp = image_data['Amp'][0]
    binn = image_data['Binning'][0]

    # Load bias frame
    if tel.bias:
        if log: log.info('Loading master bias.')
        try:
            mbias = tel.load_bias(paths, amp, binn)
        except Exception as e:
            if log: 
                log.error('No master bias found for this configuration')
                log.error(f'Skipping reduction for: {cal_type}')
                log.error(e)
            return(None)
    else:
        mbias = None

    # Load bias frame
    if tel.dark:
        if log: log.info('Loading master dark.')
        try:
            mdark = tel.load_dark(paths, amp, binn)
        except Exception as e:
            if log: 
                log.error('No master dark found for this configuration')
                log.error(f'Skipping reduction for: {cal_type}')
                log.error(e)
            return(None)
    else:
        mdark = None

    # Load flat frame
    if tel.flat:
        if log: log.info('Loading master flat.')
        try:
            mflat = tel.load_flat(paths, fil, amp, binn)
        except Exception as e:
            if log: 
                log.error('No master bias found for this configuration')
                log.error(f'Skipping reduction for: {cal_type}')
                log.error(e)
            return(None)
    else:
        mflat = None
        
    t1 = time.time()
    if log: log.info(f'Processing data for {cal_type}')

    # Bias subtraction, gain correction, flat correction, and flat fielding
    files = image_data['File']
    processed = tel.process_science(files, fil, amp, binn, paths,
        mbias=mbias, mflat=mflat, mdark=mdark, skip_skysub=skip_skysub, log=log)
    
    # Get filenames for output processed data
    reduced_files = [p.header['FILENAME'] for p in processed]

    t2 = time.time()
    if log: log.info(f'Data processed in {t2-t1} sec')
    if log: log.info('Aligning images.')

    # If masking satellite trails, this needs to be done before reprojection so
    # that the acstools algorithm accurately models the edges of the detector
    if satellites:
        mask_satellites(processed, reduced_files, log=log)

    # Sort files so deepest exposure is first
    exptimes = []
    for file in reduced_files:
        hdu = fits.open(file)
        exptimes.append(tel.get_exptime(hdu[0].header))
    idx = np.flip(np.argsort(exptimes))
    reduced_files = np.array(reduced_files)[idx]

    if out_size is None:
        out_size = tel.get_out_size(processed[0].header)

    if fieldcenter is not None:
        use_wcs = generate_wcs(tel, binn, fieldcenter, out_size)
    else:
        use_wcs = None

    aligned_images, aligned_data = align_images(reduced_files, paths, tel, binn,
        use_wcs=use_wcs, fieldcenter=fieldcenter, out_size=out_size, 
        skip_gaia=skip_gaia, log=log)

    if aligned_images is None or aligned_data is None:
        return(None)

    if log: log.info('Creating mask and error arrays.')
    masks = []
    errors = []
    data_images = []
    for stack_img in aligned_images:
        hdu = fits.open(stack_img, mode='readonly')
        hdr = hdu[0].header
        clean, mask = create_mask(stack_img,
            hdr['SATURATE'], np.mean(tel.get_rdnoise(hdr)),
            log=log, cosmic_ray=cosmic_ray, outpath=work_path)
        error = create_error(stack_img, mask, np.mean(tel.get_rdnoise(hdr)))

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

        hdulist[0].header['EXTNAME']='SCI'
        hdulist[1].header['EXTNAME']='MASK'
        hdulist[2].header['EXTNAME']='ERROR'

        filename = hdulist[0].header['FILENAME'].replace('.fits',
            '_data.fits')
        fullfilename = os.path.join(work_path, filename)

        if log: 
            log.info(f'Writing out all file data: {fullfilename}')
        else:
            print(f'Writing out all file data: {fullfilename}')

        hdulist.writeto(fullfilename, overwrite=True, output_verify='silentfix')
        data_images.append(fullfilename)
    
    if log: log.info('Creating median stack.')
    if len(aligned_data)>1:
        sci_med = stack_data(aligned_data, tel, masks, errors, log=log)
        sci_med = add_stack_mask(sci_med, aligned_data)

        if tel.detrend:
            sci_med = detrend_stack(sci_med)

    else:
        sci_med = fits.open(data_images[0])

    sci_med[0].data = sci_med[0].data.astype(np.float32)
    sci_med[1].data = sci_med[1].data.astype(np.uint8)
    sci_med[2].data = sci_med[2].data.astype(np.float32)
    sci_med[0].header['EXTNAME']='SCI'
    sci_med[1].header['EXTNAME']='MASK'
    sci_med[2].header['EXTNAME']='ERROR'
    sci_med[0].header['BITPIX'] = -32
    sci_med[1].header['BITPIX'] = 8
    sci_med[2].header['BITPIX'] = -32

    # Get time parameters from aligned data
    mid_time = np.average([tel.get_time(d.header) for d in aligned_data])
    exptimes = np.array([tel.get_exptime(d.header) for d in aligned_data])
    eff_time = np.sum(exptimes**2)/np.sum(exptimes)
    total_time = np.sum(exptimes)

    # Rescale both data and error images by effective exposure times so they 
    # are in e- instead of e-/s
    sci_med[0].data = sci_med[0].data * eff_time
    sci_med[2].data = sci_med[2].data * eff_time

    # Explicitly note that data and error extensions are in ELECTRONS
    sci_med[0].header['BUNIT'] = 'ELECTRONS'
    sci_med[2].header['BUNIT'] = 'ELECTRONS'

    # Add both formats for code that requires either
    sci_med[0].header['MJD-OBS'] = (mid_time, 
        'Mid-MJD of the observation sequence.')
    sci_med[0].header['MJD'] = (mid_time, 
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
    rdnoise = np.mean(tel.get_rdnoise(sci_med[0].header))/np.sqrt(len(aligned_data))
    sci_med[0].header['RDNOISE'] = (rdnoise, 'Readnoise of stack.')
    
    sci_med[0].header['NFILES'] = (len(aligned_data), 
        'Number of images in stack')
    sci_med[0].header['FILTER'] = fil
    sci_med[0].header['OBSTYPE'] = 'OBJECT'
    
    # Generate stack name and write out
    stackname = tel.get_stk_name(sci_med[0].header, red_path)
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

def stack_data(stacking_data, tel, masks, errors, mem_limit=8.0e9, log=None):

    stack_method = tel.stack_method
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
        exptimes.append(float(tel.get_exptime(stk.header)))

    exptimes = np.array(exptimes)
    scale = 1./exptimes

    all_data = np.array([s.data for s in stacking_data])

    if len(scale)==1:
        stack_method='average'

    # Determine if intermediate stacks are needed
    # Get size of individual frame, account for data, noise, and mask
    exdata = stacking_data[0].data
    size = 64.0/8.0 * exdata.shape[0]*exdata.shape[1] * 2 + exdata.shape[0]*exdata.shape[1]

    if log:
        log.info(f'Image size is {size}')
    else:
        print(f'Image size is {size}')

    # Get maximum size of chunk
    max_chunk = int(np.floor(mem_limit / 4.0 / (size)))
    if max_chunk==0: max_chunk=1

    if len(stacking_data)<=max_chunk:
        sci_med = combine(stacking_data, weights=weights, scale=scale,
            method=stack_method, mem_limit=mem_limit)
    else:
        nimgs = len(stacking_data) * 1.0
        nchunks = int(np.ceil(nimgs/(1.0*max_chunk)))

        if log:
            log.info(f'Splitting stacking into {nchunks} chunks')
        else:
            print(f'Splitting stacking into {nchunks} chunks')

        inter_med = []
        for i in np.arange(nchunks):
            if log:
                log.info(f'Stacking chunk {i+1}/{nchunks}')
            else:
                print(f'Stacking chunk {i+1}/{nchunks}')
            start = i*max_chunk
            end = (i+1)*max_chunk
            if end>len(stacking_data): end=None

            sci_med = combine(stacking_data[start:end], 
                weights=weights[start:end], scale=scale[start:end],
                method=stack_method, mem_limit=mem_limit)

            inter_med.append(sci_med)

        sci_med = combine(inter_med, method=stack_method)

    sci_med = sci_med.to_hdu()

    return(sci_med)

def add_stack_mask(stack, stacking_data):

    data = stack[0].data
    mask = stack[1].data
    error = stack[2].data

    # Keep track of original masked values
    bp_mask = (mask > 0) | (data==0.0) | np.isnan(data) | np.isnan(mask) |\
        np.isnan(error) | (error==0.0) | np.isinf(data)

    # Reset mask
    mask = np.zeros(mask.shape).astype(np.uint8)

    # Add bad pixels back to mask
    mask[bp_mask] = 1

    # Add saturation to mask
    sat = np.min([d.header['SATURATE'] for d in stacking_data])
    stack[0].header['SATURATE'] = sat
    sat_mask = data >= sat
    mask[sat_mask] = 4

    # Set masked values to 0
    data[mask>0] = 0.0
    error[mask>0] = 0.0

    # Make sure that nan and inf values are removed from data
    data[np.isnan(data)]=0.0
    data[np.isinf(data)]=0.0

    stack[0].data = data
    stack[1].data = mask
    stack[2].data = error

    return(stack)

def mask_satellites(images, filenames, log=None):

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

def create_mask(science_data, saturation, rdnoise, sigclip=3.5, 
    sigfrac=0.2, objlim=4.5, niter=6, outpath='', grow=0, satellites=True,
    cosmic_ray=True, log=None):

    t_start = time.time()
    
    if log:
        log.info(f'Running astroscrappy on {science_data}')
    else:
        print(f'Running astroscrappy on {science_data}')

    hdu = fits.open(science_data)
    data = np.ascontiguousarray(hdu[0].data.astype(np.float32))

    # Astroscrappy requires added sky background, so add this value back
    # Set the sky background to some nominal value if it is too low for CR rej
    skybkg = hdu[0].header['SKYBKG']
    if log:
        log.info(f'Sky background in science frame is {skybkg}')
    else:
        print(f'Sky background in science frame is {skybkg}')

    if skybkg < 2000.0: 
        skybkg = 2000.0
        if log:
            log.info('Setting sky background to 2000.0')
        else:
            print('Setting sky background to 2000.0')

    # Set data background to skybkg
    data = data + skybkg
    # Also need to adjust the saturation level by SKYBKG for saturated pixels
    saturation += skybkg
    
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

    mask_cr = np.zeros(data.shape)
    mask_cr = mask_cr.astype(np.uint8)

    if cosmic_ray:
        newdata, mask_cr = cosmicray_lacosmic(data,
            readnoise=rdnoise, satlevel=saturation, verbose=True,
            sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=niter)

        mask_cr = mask_cr.astype(np.uint8) #set cosmic ray mask type
    else:
        newdata = copy.copy(data)

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

    return(newdata, mask_hdu)


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

    poisson = img_data
    poisson[poisson<0.]=0.
    error = np.sqrt(poisson + rms**2 + rdnoise)

    error_hdu = fits.PrimaryHDU(error) #create mask Primary HDU
    error_hdu.header['VER'] = (__version__, 
        'Version of image procedures used used.')
    error_hdu.header['USE'] = 'Error array for Poisson, read, and RMS noise'
    error_hdu.header['BUNIT'] = ('ELECTRONS', 'Units of the error array.')

    return(error_hdu)


if __name__=="__main__":

    t = ascii.read('/Users/ckilpatrick/Downloads/FRB20250316A/2025.0323/file_list.txt',
        format='fixed_width')
    mask = t['TargType']=='FRB20250316A_r_2_11'
    t = t[mask]

    t.sort('File')

    global tel
    import importlib
    module = importlib.import_module('instruments.BINOSPEC')
    tel = getattr(module, 'BINOSPEC')()

    paths={}
    paths['red']='/Users/ckilpatrick/Downloads/FRB20250316A/2025.032/red'
    paths['work']='/Users/ckilpatrick/Downloads/FRB20250316A/2025.032/workspace'
    paths['cal']='/Users/ckilpatrick/Downloads/FRB20250316A/2025.032/red/cals'
    paths['code']=os.path.dirname(sys.argv[0])

    import glob
    files = glob.glob('/Users/ckilpatrick/Downloads/FRB20250316A/2025.0323/red/workspace/*_data.fits')
    aligned_data=[]
    masks=[]
    errors=[]
    for file in files:
        hdu = fits.open(file)
        wcs = WCS(hdu[0].header)
        image = CCDData(hdu[0].data, meta=hdu[0].header, wcs=wcs, 
            unit=u.electron)
        aligned_data.append(image)
        masks.append(hdu['MASK'])
        errors.append(hdu['ERROR'])

    sci_med = stack_data(aligned_data, tel, masks, errors, log=None)

    #image_proc(t, tel, paths, fieldcenter=None,
    #    out_size=4200, cosmic_ray=False, log=None)

