#!/usr/bin/env python

"Automatic generalized pipeline for imaging reduction. Creates median coadded images for each target."
"Individual images should be checked and removed from sci path."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "1.0" #last updated 20/04/2021

import sys
import numpy as np
import os
import datetime
import time
import astropy
import astropy.units.astrophys as u
import astropy.units as u
from astropy.io import fits
import ccdproc
import glob
import argparse
import logging
from astropy.nddata import CCDData
from photutils import Background2D, MeanBackground
from astropy.stats import SigmaClip
import importlib
import Sort_files
import align_quads
import solve_wcs
import quality_check

def main_pipeline(telescope,data_path,cal_path=None,target=None,skip_red=None,proc=None,use_dome_flats=None,phot=None):
    #start time
    t_start = time.time()
    #import telescope parameter file
    global tel
    try:
        tel = importlib.import_module('tel_params.'+telescope)
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)

    raw_path = data_path+'/raw/' #path containing the raw data
    if not os.path.exists(raw_path): #create reduced file path if it doesn't exist
        os.makedirs(raw_path)
    bad_path = data_path+'/bad/' #path containing the raw data
    if not os.path.exists(bad_path): #create reduced file path if it doesn't exist
        os.makedirs(bad_path)
    spec_path = data_path+'/spec/' #path containing the raw data
    if not os.path.exists(spec_path): #create reduced file path if it doesn't exist
        os.makedirs(spec_path)
    red_path = data_path+'/red/' #path to write the reduced files
    if not os.path.exists(red_path): #create reduced file path if it doesn't exist
        os.makedirs(red_path)

    if cal_path is not None:
        cal_path = cal_path
    else:
        cal_path = tel.cal_path() #os.getenv("HOME")+'/Pipelines/MMIRS_calib/'
    if cal_path:             
        flat_path = cal_path
    else:
        flat_path = red_path


    wavelength = tel.wavelength()

    log_file_name = red_path+telescope+'_log_'+datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')+'.log' #create log file name
    log = logging.getLogger(log_file_name) #create logger
    log.setLevel(logging.INFO) #set level of logger
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s") #set format of logger
    logging.Formatter.converter = time.gmtime #convert time in logger to UCT
    filehandler = logging.FileHandler(log_file_name, 'w+') #create log file
    filehandler.setFormatter(formatter) #add format to log file
    log.addHandler(filehandler) #link log file to logger
    streamhandler = logging.StreamHandler() #add log stream to logger
    streamhandler.setFormatter(formatter) #add format to log stream
    log.addHandler(streamhandler) #link logger to log stream

    log.info('Running main pipeline version '+str(__version__))
    log.info('Running telescope paramater file version '+str(tel.__version__))

    if os.path.exists(data_path+'/file_list.txt'):
        log.info('Previous file list exists, loading lists.')
        cal_list, sci_list, sky_list, time_list = Sort_files.load_files(data_path+'/file_list.txt', telescope,log)
    else:
        log.info('Sorting files and creating file lists.')
        files = sorted(glob.glob(data_path+tel.raw_format(proc)))
        if len(files) != 0:
            log.info(str(len(files))+' files found.')
            cal_list, sci_list, sky_list, time_list = Sort_files.sort_files(files,telescope,data_path,log)
        else:
            log.critical('No files found, please check data path and rerun.')
            logging.shutdown()
            sys.exit(-1)
    
    if tel.bias():
        if len(cal_list['BIAS']) != 0:
            process_bias = True
            if skip_red:
                log.info('User input to skip reduction.')
                if os.path.exists(red_path+'mbias.fits'):
                    log.info('Found previous master bias, loading.')
                    mbias = CCDData.read(red_path+'mbias.fits', unit=u.electron)
                    process_bias = False
                else:
                    log.error('No master bias found, creating master bias.')
            if process_bias:
                t1 = time.time()
                log.info('Processing bias files.')
                mbias = tel.create_bias(cal_list,red_path,log)
                t2 = time.time()
                log.info('Master bias creation completed in '+str(t2-t1)+' sec')
        else:
            log.critical('No bias files present, check data before rerunning.')
            logging.shutdown()
            sys.exit(-1)
    else:
        mbias = None

    if tel.dark():
        for cal in cal_list:
            if 'DARK' in cal:
                process_dark = True
                if skip_red:
                    log.info('User input to skip reduction.')
                    if os.path.exists(red_path+cal+'.fits'):
                        log.info('Found previous master dark, loading.')
                        process_dark = False
                    else:
                        log.error('No master dark found, creating master dark.')
                if process_dark:
                    t1 = time.time()
                    tel.create_dark(cal_list,cal,mbias,red_path,log)
                    t2 = time.time()
                    log.info('Master dark creation completed in '+str(t2-t1)+' sec')

    if tel.flat():
        for cal in cal_list:
            if 'FLAT' in cal:
                process_flat = True
                fil = cal.split('_')[-1]
                if skip_red:
                    log.info('User input to skip reduction.')
                    master_flat = tel.flat_name(flat_path, fil)
                    if np.all([os.path.exists(mf) for mf in master_flat]):
                        log.info('Found previous master flat for filter '+fil)
                        process_flat = False
                    else:
                        log.info('No master flat found for filter '+fil+', creating master flat.')
                if process_flat:
                    if wavelength=='OPT':
                        t1 = time.time()
                        tel.create_flat(cal_list[cal],fil,red_path,mbias=mbias,log=log)
                        t2 = time.time()
                        log.info('Master flat creation completed in '+str(t2-t1)+' sec')
                    elif wavelength=='NIR':
                        if use_dome_flats: #use dome flats instead of sky flats for NIR
                            log.info('User input to use dome flats to create master flat')
                            flat_type = 'dome'
                            t1 = time.time()
                            tel.create_flat(cal_list[cal],fil,red_path,mbias=mbias,log=log)
                            t2 = time.time()
                            log.info('Master flat creation completed in '+str(t2-t1)+' sec')
        if wavelength=='NIR' and not use_dome_flats: #default to use science files for master flat creation
            for fil_list in sky_list:
                process_flat = True
                if skip_red:
                    log.info('User input to skip reduction.')
                    master_flat = tel.flat_name(flat_path, fil_list)
                    if np.all([os.path.exists(mf) for mf in master_flat]):
                        log.info('Found previous master flat for filter '+fil_list)
                        process_flat = False
                    else:
                        log.info('No master flat found for filter '+fil_list+', creating master flat.')
                if process_flat:              
                    flat_type = 'sky'
                    log.info('Using science files to create master flat')
                    t1 = time.time()
                    tel.create_flat(sky_list[fil_list],fil_list,red_path,mbias=mbias,log=log)
                    t2 = time.time()
                    log.info('Master flat creation completed in '+str(t2-t1)+' sec')
    
    if len(sci_list) == 0:
        log.critical('No science files to process, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)     
    log.info('User input target for reduction: '+str(target))
    for tar in sci_list:
        stack = tel.stacked_image(tar,red_path)
        target = tar.split('_')[0]
        fil = tar.split('_')[-1]
        if target is not None:
            if target not in tar:
                continue
            else:
                log.info('Matching target found: '+tar)
        process_data = True
        if skip_red:
            log.info('User input to skip reduction.')
            if tel.run_wcs():
                final_stack = [st.replace('.fits','_wcs.fits') for st in stack]
            else:
                final_stack = stack
            if np.all([os.path.exists(st) for st in stack]):
                process_data = False
            else:
                log.error('Missing stacks, processing data.')   
        if process_data:
            log.info('Loading master flat.')
            master_flat = tel.flat_name(flat_path, fil)
            if not np.all([os.path.exists(mf) for mf in master_flat]):
                log.error('No master flat present for filter '+fil+', skipping data reduction for '+tar+'. Check data before rerunning')
                continue
            flat_data = tel.load_flat(master_flat)
            t1 = time.time()
            log.info('Processing data for '+str(tar))
            processed, masks = tel.process_science(sci_list[tar],fil,red_path,mbias=mbias,mflat=flat_data,proc=proc,log=log)
            t2 = time.time()
            log.info('Data processed in '+str(t2-t1)+' sec')
            if wavelength=='NIR':
                t1 = time.time()
                log.info('NIR data, creating NIR sky maps.')
                for j,n in enumerate(processed):
                    time_diff = sorted([(abs(time_list[tar][j]-n2),k) for k,n2 in enumerate(time_list[tar])])
                    sky_list = [sci_list[tar][k] for _,k in time_diff[0:3]]
                    sky_data = [processed[k] for _,k in time_diff[0:3]]
                    sky_mask = [masks[k] for _,k in time_diff[0:3]]
                    sky_masked_data = []
                    for k in range(3): 
                        bkg = Background2D(sky_data[k], (20, 20), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=sky_mask[k], exclude_percentile=80)
                        masked = np.array(sky_data[k])
                        masked[sky_mask[k]] = bkg.background[sky_mask[k]]
                        sky_masked_data.append(CCDData(masked,unit=u.electron/u.second))
                    sky_hdu = fits.PrimaryHDU()
                    sky_hdu.header['FILE'] = (os.path.basename(sci_list[tar][j]), 'NIR sky flat for file.')
                    for k,m in enumerate(sky_list):
                        sky_hdu.header['FILE'+str(k+1)] = (os.path.basename(str(m)), 'Name of file used in creation of sky.')
                    sky = ccdproc.combine(sky_masked_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median,mask=sky_mask)
                    sky.header = sky_hdu.header
                    sky.write(red_path+os.path.basename(sci_list[tar][j]).replace('.fits','_sky.fits'),overwrite=True)
                    processed[j] = n.subtract(sky,propagate_uncertainties=True,handle_meta='first_found')
                t2 = time.time()
                log.info('Sky maps complete and subtracted in '+str(t2-t1)+' sec')
            log.info('Writing out reduced data.')
            dimen = len(stack)
            if dimen == 1:
                suffix = ['_red.fits']
            else:
                log.info('Multiple extensions to stack.')
                suffix = tel.suffix()
            mask = tel.static_mask(proc)
            for k in range(dimen):
                red_list = [red_path+os.path.basename(sci).replace('.fits',suffix[k]) for sci in sci_list[tar]]
                if dimen == 1:
                    for j,process_data in enumerate(processed):
                        process_data.write(red_list[j],overwrite=True)
                else:
                    for j,process_data in enumerate(processed[k]):
                        process_data.write(red_list[j],overwrite=True)
                log.info('Aligning images.')
                aligned_images, aligned_data = align_quads.align_stars(red_list,telescope,hdu=tel.wcs_extension(),mask=mask[k],log=log)
                log.info('Checking qualty of images.')
                stacking_data, mid_time = quality_check.quality_check(aligned_images, aligned_data, telescope, log)
                log.info('Creating median stack.')
                sci_med = ccdproc.combine(stacking_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
                sci_med.header['MJD-OBS'] = (mid_time, 'Mid-MJD of the observation sequence calculated using DATE-OBS.')
                sci_med.header['EXPTIME'] = (1, 'Effective expsoure tiime for the stack in seconds.')
                sci_med.header['GAIN'] = (1, 'Effecetive gain for stack.')
                sci_med.header['RDNOISE'] = (tel.rdnoise(sci_med.header)/np.sqrt(len(aligned_images)), 'Readnoise of stack.')
                sci_med.header['NFILES'] = (len(aligned_images), 'Number of images in stack')
                sci_med.write(stack[k],overwrite=True)
                log.info('Median stack made for '+stack[k])
                if tel.run_wcs():
                    log.info('Solving WCS')
                    wcs_error = solve_wcs.solve_wcs(stack[k],telescope,log=log)   
                    log.info(wcs_error)
                    stack[k] = stack[k].replace('.fits','_wcs.fits')
        #do auto phot
        #if phot: do manual
    t_end = time.time()
    log.info('Pipeline finshed.')
    log.info('Total runtime: '+str(t_end-t_start)+' sec')


def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('telescope', default=None, help='Path of data, in default format from KOA.') #path to MOSFIRE data: required
    params.add_argument('data_path', default=None, help='Path of data, in default format from KOA.') #path to MOSFIRE data: required
    params.add_argument('--use_dome_flats', type=str, default=None, help='Use dome flats for flat reduction.') #use dome flat instead of sci images to create master flat
    params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.') #
    params.add_argument('--target', type=str, default=None, help='Option to only reduce this target.') #
    params.add_argument('--proc', type=str, default=None, help='If working with the _proc data from MMT.')
    params.add_argument('--cal_path', type=str, default=None, help='Use dome flats for flat reduction.') #use dome flat instead of sci images to create master flat
    params.add_argument('--phot', type=str, default=None, help='Option to use IRAF to perform photometry.') #must have pyraf install and access to IRAF to use
    args = params.parse_args()
    
    main_pipeline(args.telescope,args.data_path,args.cal_path,target=args.target,skip_red=args.skip_red,proc=args.proc,use_dome_flats=args.use_dome_flats,phot=args.phot)

if __name__ == "__main__":
    main()
    