#!/usr/bin/env python

"Automatic generalized pipeline for imaging reduction. Creates median coadded images for each target."
"Individual images should be checked and removed from sci path."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "1.0" #last updated 18/03/2021

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

def main_pipeline(telescope,data_path,cal_path=None,target=None,skip_red=None,proc=None):
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

    if os.path.exists(data_path+'/file_list.txt'):
        cal_list, sci_list, sky_list, time_list = Sort_files.load_files(data_path+'/file_list.txt', telescope)
    else:
        files = sorted(glob.glob(data_path+tel.raw_format(proc)))
        if len(files) != 0:
            cal_list, sci_list, sky_list, time_list = Sort_files.sort_files(files,telescope,data_path)
        else:
            log.critical('No files found, please check data path and rerun.')
            logging.shutdown()
            sys.exit(-1)
    
    if tel.bias():
        if len(cal_list['BIAS']) != 0:
            process_bias = True
            if skip_red == 'True' or skip_red == 'yes':
                if os.path.exists(red_path+'mbias.fits'):
                    mbias = CCDData.read(red_path+'mbias.fits', unit=u.electron)
                    process_bias = False
                else:
                    log.error('No master bias found, creating master bias.')
            if process_bias:
                processed = []
                for bias in cal_list['BIAS']:
                    try:
                        raw = CCDData.read(bias, hdu=1, unit='adu')
                    except astropy.io.registry.IORegistryError:
                        log.error('File '+bias+' not recognized.')
                    processed.append(ccdproc.ccd_process(raw, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron))
                mbias = ccdproc.combine(processed,method='median')
                mbias.write(red_path+'mdark.fits',overwrite=True)
        else:
            log.critical('No bias present, check data before rerunning.')
            logging.shutdown()
            sys.exit(-1)
    else:
        mbias = None

    if tel.dark():
        if len(cal_list['DARK']) != 0:
            process_dark = True
            if skip_red == 'True' or skip_red == 'yes':
                if os.path.exists(red_path+'mdark.fits'):
                    mdark = CCDData.read(red_path+'mdark.fits', unit=u.electron)
                    process_dark = False
                else:
                    log.error('No master dark found, creating master dark.')
            if process_dark:
                processed = []
                for dark in cal_list['DARK']:
                    try:
                        raw = CCDData.read(dark, hdu=1, unit='adu')
                    except astropy.io.registry.IORegistryError:
                        log.error('File '+dark+' not recognized.')
                    processed.append(ccdproc.ccd_process(raw, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron, master_bias=mbias))
                mdark = ccdproc.combine(processed,method='median')
                mdark.write(red_path+'mdark.fits',overwrite=True)
        else:
            log.critical('No darks present, check data before rerunning.')
            logging.shutdown()
            sys.exit(-1)
    else:
        mdark = None

    if tel.flat():
            for cal in cal_list:
                if 'FLAT' in cal:
                    if wavelength=='OPT':
                        for cal in cal_list:
                            if 'FLAT' in cal:
                                tel.create_flat(cal_list[cal])
                    elif wavelength=='NIR':
                        fil = cal.split('_')[-1]
                        if use_dome_flats == 'yes' or use_dome_flats == 'True': #use dome flats instead of sky flats for NIR
                            flat_type = 'dome'
                            log.info('User set option to use dome flats to create master flat')
                            log.info(str(len(cal_list[cal]))+' dome flats found for filter '+fil)
                            tel.create_flat(cal_list[cal])
                        else: #default to use science files for master flat creation
                            flat_type = 'sky'
                            log.info('Using science files to create master flat')
                            tel.create_flat(sky_list[fil])

    if len(sci_list) == 0:
        log.critical('No science files to process, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)     
    log.info('User specified a target for reduction: '+str(target))
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
        if skip_red == 'True' or skip_red == 'yes':
            if np.all([os.path.exists(st) for st in stack]):
                # with fits.open(red_path+tar+'.fits',mode='update') as hdr:
                #     sci_med = hdr[0].data
                #     wcs_object = wcs.WCS(hdr[0].header)
                # processed = [x.replace(sci_path,red_path).replace('.fits','_red.fits') for x in sci_list[tar]]
                process_data = False
            else:
                log.error('No reduced image found, processing data.')      
        if process_data:             
            master_flat = tel.flat_name(cal_path, fil)
            if not np.all([os.path.exists(mf) for mf in master_flat]):
                log.error('No master flat present for filter '+fil+', skipping data reduction for '+tar+'. Check data before rerunning')
                continue
            flat_data = tel.load_flat(master_flat)
            log.info('Processing data for '+str(tar))
            processed, masks = tel.process_science(sci_list[tar],fil,cal_path,mdark=mdark,mbias=mbias,mflat=flat_data,proc=proc)
            log.info('Data processed.')
            if wavelength=='NIR':
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
            dimen = len(stack)
            if dimen == 1:
                red_list = [red_path+os.path.basename(sci).replace('.fits','_red.fits') for sci in sci_list[tar]]
                for j,process_data in enumerate(processed):
                    process_data.write(red_list[j],overwrite=True)
                log.info('Aligning images.')
                aligned = align_quads.align_stars(red_list,telescope,hdu=tel.wcs_extension(),mask=tel.static_mask(proc))
                log.info('Images aligned, creating median stack.')
                sci_med = ccdproc.combine(aligned,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
                sci_med.header['RDNOISE'] = tel.rdnoise(sci_med.header)/np.sqrt(len(aligned))
                sci_med.header['NFILES'] = len(aligned)
                sci_med.write(stack[0],overwrite=True)
                log.info('Median stack made for '+str(tar))
            else:
                log.info('Multiple extensions to stack.')
                suffix = tel.suffix()
                mask = tel.static_mask(proc)
                for k in range(dimen):
                    red_list = [red_path+os.path.basename(sci).replace('.fits',suffix[k]) for sci in sci_list[tar]]
                    for j,process_data in enumerate(processed[k]):
                        process_data.write(red_list[j],overwrite=True)
                    log.info('Aligning images.')
                    aligned = align_quads.align_stars(red_list,telescope,hdu=tel.wcs_extension(),mask=mask[k])
                    log.info('Images aligned, creating median stack.')
                    sci_med = ccdproc.combine(aligned,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
                    sci_med.header['RDNOISE'] = tel.rdnoise(sci_med.header)/np.sqrt(len(aligned))
                    sci_med.header['NFILES'] = len(aligned)
                    sci_med.write(stack[k],overwrite=True)
                    log.info('Median stack made for '+str(tar))                    


def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('telescope', default=None, help='Path of data, in default format from KOA.') #path to MOSFIRE data: required
    params.add_argument('data_path', default=None, help='Path of data, in default format from KOA.') #path to MOSFIRE data: required
    params.add_argument('--use_dome_flats', type=str, default=None, help='Use dome flats for flat reduction.') #use dome flat instead of sci images to create master flat
    params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.') #
    params.add_argument('--target', type=str, default=None, help='Option to only reduce this target.') #
    params.add_argument('--proc', type=str, default=None, help='If working with the _proc data from MMT.')
    params.add_argument('--cal_path', type=str, default=None, help='Use dome flats for flat reduction.') #use dome flat instead of sci images to create master flat
    params.add_argument('--wcs', type=str, default=None, help='Option to use IRAF to apply WCS.') #must have pyraf install and access to IRAF to use
    params.add_argument('--phot', type=str, default=None, help='Option to use IRAF to perform photometry.') #must have pyraf install and access to IRAF to use
    args = params.parse_args()
    
    main_pipeline(args.telescope,args.data_path,args.cal_path,target=args.target,skip_red=args.skip_red,proc=args.proc)

if __name__ == "__main__":
    main()
    