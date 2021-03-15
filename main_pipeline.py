#!/usr/bin/env python

"Automatic generalized pipeline for imaging reduction. Creates median coadded images for each target."
"Individual images should be checked and removed from sci path."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "1.0" #last updated 25/02/2021

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
import importlib
import Sort_files
import align_quads
import solve_wcs

def main_pipeline(telescope,data_path,cal_path=None,target=None,skip_red=None):
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
        cal_list, sci_list, sky_list, time_list = Sort_files.load_files(data_path+'/file_list.txt')
    else:
        files = glob.glob(data_path+tel.raw_format())
        if len(files) != 0:
            cal_list, sci_list, sky_list, time_list = Sort_files.sort_files(files,telescope,data_path)
        else:
            log.critical('No files found, please check data path and rerun.')
            logging.shutdown()
            sys.exit(-1)

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
                    processed.append(ccdproc.ccd_process(raw, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron))
                mdark = ccdproc.combine(processed,method='median')
                mdark.write(red_path+'mdark.fits',overwrite=True)
        else:
            con = input('No darks present. Continue without dark subtraction? (True or False) ')
            if con == 'True':
                log.error('No darks present, continuing without dark subtraction.')
                mdark = None
            else:
                log.critical('No darks present, check data before rerunning.')
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
            con = input('No bias present. Continue without bias subtraction? (True or False) ')
            if con == 'True':
                log.error('No bias present, continuing without bias subtraction.')
                mbias = None
            else:
                log.critical('No bias present, check data before rerunning.')
                logging.shutdown()
                sys.exit(-1)

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
    log.info('User specified a target for reduction: '+target)
    for tar in sci_list:
        target = tar.split('_')[0]
        fil = tar.split('_')[-1]
        if target is not None:
            if target not in tar:
                continue
            else:
                log.info('Matching target found: '+tar)
        process_data = True
        if skip_red == 'True' or skip_red == 'yes':
            if os.path.exists(red_path+tar+'.fits'):
                with fits.open(red_path+tar+'.fits',mode='update') as hdr:
                    sci_med = hdr[0].data
                    wcs_object = wcs.WCS(hdr[0].header)
                processed = [x.replace(sci_path,red_path).replace('.fits','.fits') for x in sci_list[tar]]
                process_data = False
            else:
                log.error('No reduced image found, processing data.')        
        if process_data:             
            fil = target.split('_')[-1]
            master_flat = tel.flat_name(fil)
            if not os.path.exists(master_flat):
                log.error('No master flat present for filter '+fil+', skipping data reduction for '+tar+'. Check data before rerunning')
                continue
            flat_data = tel.load_flat(master_flat)
            processed, masks = tel.process_science(sci_list[tar])
            if wavelength=='NIR':
                for j,n in enumerate(processed):
                    time_diff = sorted([(abs(time_list[tar][j]-n2),k) for k,n2 in enumerate(time_list[target])])
                    sky_list = [time_list[tar][k] for _,k in time_diff[1:10]]
                    sky_data = [processed[k] for _,k in time_diff[1:10]]
                    sky_mask = [masks[k] for _,k in time_diff[1:10]]
                    sky_hdu = fits.PrimaryHDU()
                    sky_hdu.header['FILE'] = (time_list[tar][j], 'NIR sky flat for file.')
                    for k,m in enumerate(sky_list):
                        sky_hdu.header['FILE'+str(k+1)] = (os.path.basename(m), 'Name of file used in creation of sky.')
                    sky = ccdproc.combine(sky_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median,mask=sky_mask)
                    sky.header = sky_hdu.header
                    sky.write(red_path+os.path.basename(time_list[tar][j]).replace('.fits','_sky.fits'),overwrite=True)
                    processed[j] = n.subtract(sky,propagate_uncertainties=True,handle_meta='first_found')
            for process_data in processed:
                process_data.write(red_path+os.path.basename(sci_list[tar][j]).replace('.fits','_red.fits'),overwrite=True)
            red_list = [sci.replace('.fits','_red.fits') for sci in sci_list[tar]]
            aligned = align_quads.align_stars(red_list,telescope,hdu=tel.wcs_extension())
            sci_med = ccdproc.combine(aligned,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
            sci_med.header['RDNOISE'] = sci_med.header['RDNOISE']/len(aligned)
            sci_med.header['NFILES'] = len(aligned)
            sci_med.write(red_path+tar+'.fits',overwrite=True)


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
    
    main_pipeline(args.telescope,args.data_path,args.cal_path,target=args.target,skip_red=args.skip_red)

if __name__ == "__main__":
    main()
    