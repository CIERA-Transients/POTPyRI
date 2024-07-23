#!/usr/bin/env python

"Automatic generalized pipeline for imaging reduction. Creates median coadded images for each target."
"Individual images should be checked and removed from sci path."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

# Refactor on 2024-07-19
__version__ = "1.22"

import sys
import numpy as np
import os
import datetime
import time
import shutil
import astropy
import astropy.units.astrophys as u
import astropy.units as u
from astropy.io import fits
import ccdproc
import glob
import argparse
import logging
import astropy.wcs as wcs
from astropy.nddata import CCDData
from photutils.background import Background2D, MeanBackground
from astropy.stats import SigmaClip
from astropy.coordinates import SkyCoord
import importlib
import Sort_files
import align_quads
import solve_wcs
import quality_check
import psf
import absphot
import cal_procs

from custom_logger import ColoredLogger
from astroquery.astrometry_net import AstrometryNet

# Update to new photutils methods and deprecate previous Find_target_phot when possible
import Find_target_phot as tp

from utilities import extinction
from utilities import util
from colorama import init, Fore, Back, Style
init()

def main_pipeline(telescope:str,data_path:str,
                  cal_path:str=None,
                  input_target:list=None,
                  skip_red:bool=None,
                  proc:str=None,
                  use_dome_flats:bool=None,
                  phot:bool=None,
                  reset:bool=None,
                  use_anet:bool=None,
                  no_verify:bool=None)->None:

    # start time
    t_start = time.time()
    # import telescope parameter file
    global tel
    try:
        tel = importlib.import_module('params.'+telescope)
    except ImportError:
        telfile=f'params.{telescope}'
        print(f'''No such telescope file {telfile}, please check that you have 
            entered the correct name or this telescope is available.''')
        sys.exit(-1)

    raw_path = os.path.join(data_path, 'raw') #path containing the raw data
    bad_path = os.path.join(data_path, 'bad') #path containing the bad (unused) data
    red_path = os.path.join(data_path, 'red') #path to write the reduced files
    for dirname in [raw_path, bad_path, red_path]:
        if not os.path.exists(dirname): os.makedirs(dirname)

    # Get path to code directory
    abspath = os.path.abspath(sys.argv[0])
    code_dir = os.path.split(abspath)[0]

    # Copy config directory to data_path
    config_dir = os.path.join(code_dir, 'config')
    if not os.path.exists(os.path.join(data_path, 'config')):
        shutil.copytree(config_dir, os.path.join(data_path, 'config'))


    if cal_path is not None:
        cal_path = cal_path
    else:
        cal_path = tel.cal_path()
    if cal_path:
        flat_path = cal_path
    else:
        flat_path = red_path

    wavelength = tel.wavelength()

    # Generate logfile
    datestr = datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')
    base_logname = f'{telescope}_log_{datestr}.log'
    log_file_name = os.path.join(red_path, base_logname)
    log = ColoredLogger(log_file_name)

    log.info(f'Running main pipeline version {__version__}')
    log.info(f'Running telescope paramater file version {tel.__version__}')

    file_list = os.path.join(data_path, 'file_list.txt')

    if reset is not None:
        if os.path.exists(file_list):
            os.remove(file_list)
        if reset=='all':
            files = glob.glob(os.path.join(raw_path,'*'))+\
                    glob.glob(os.path.join(bad_path,'*'))
        if reset=='raw':
            files = glob.glob(os.path.join(raw_path,'*'))
        for f in files:
            shutil.move(f,data_path)

    # CDK - added editing for headers with LRIS raw data files (*[b,r]*.fits)
    if ((telescope=='LRIS' and proc and str(proc)=='raw') or
        (telescope=='DEIMOS' and proc and str(proc)=='raw')):
        log.info('Edit raw headers')
        tel.edit_raw_headers(data_path)

    if os.path.exists(file_list):
        log.info('Previous file list exists, loading lists.')
        cal_list, sci_list, sky_list, time_list = Sort_files.load_files(file_list, telescope,log)
    else:
        log.info('Sorting files and creating file lists.')
        # CDK - updated this to a os.path.join so it doesn't require trailing /
        glob_str = os.path.join(data_path,tel.raw_format(proc))
        log.info(f'Looking for files with format: {glob_str}')
        files = sorted(glob.glob(glob_str))
        if len(files) != 0:
            log.info(str(len(files))+' files found.')
            cal_list, sci_list, sky_list, time_list = Sort_files.sort_files(files,telescope,data_path,log)
            if not no_verify:
                print('''Please check the created file list and make any 
                    nessecary edits (e.g. remove bad files mentioned by the 
                    observing log etc.)''')
                con = input(Back.RED+'Did you edit the file list (yes or no)? '+Style.RESET_ALL)
            else:
                con='yes'
            if con=='yes':
                log.info('Edits made to the file list by the user, reloading file list.')
                cal_list, sci_list, sky_list, time_list = Sort_files.load_files(file_list, telescope,log)
        else:
            log.critical('No files found, please check data path and rerun.')
            logging.shutdown()
            sys.exit(-1)

    # CDK - added editing for headers with LRIS raw data files (*[b,r]*.fits)
    if ((telescope=='LRIS' and proc and str(proc)=='raw') or
        (telescope=='DEIMOS' and proc and str(proc)=='raw')):
        tel.edit_raw_headers(raw_path)

    # Master bias creation
    if tel.bias():
        cal_procs.do_bias(cal_list, tel, red_path, skip_red=skip_red, log=log)

    # Master dark creation
    if tel.dark():
        cal_procs.do_dark(cal_list, tel, red_path, skip_red=skip_red, log=log)
    
    # Master flat creation
    if tel.flat():
        for cal in cal_list:
            if 'FLAT' in cal:
                process_flat = True
                fil = cal.split('_')[1]
                amp = cal.split('_')[2]
                binn = cal.split('_')[3]
                if skip_red:
                    log.info('User input to skip reduction.')
                    master_flat = tel.get_mflat_name(flat_path, fil, amp, binn)
                    if os.path.exists(master_flat):
                        log.info('Found previous master flat for filter '+fil+', '+amp+' amps and '+binn+' binning.')
                        process_flat = False
                    else:
                        log.info('No master flat found for filter '+fil+', '+amp+' amps and '+binn+' binning, creating master flat.')
                if process_flat:
                    if tel.bias():
                        log.info('Loading master bias.')
                        try:
                            mbias = tel.load_bias(red_path,amp,binn)
                        except:
                            log.error('No master bias found for this configuration, skipping master flat creation for filter '+fil+', '+amp+' amps and '+binn+' binning.')
                            continue
                    else:
                        mbias = None
                    if wavelength=='OPT':
                        t1 = time.time()
                        tel.create_flat(cal_list[cal],fil,amp,binn,red_path,mbias=mbias,log=log)
                        t2 = time.time()
                        log.info('Master flat creation completed in '+str(t2-t1)+' sec')
                    elif wavelength=='NIR':
                        if use_dome_flats: #use dome flats instead of sky flats for NIR
                            log.info('User input to use dome flats to create master flat')
                            flat_type = 'dome'
                            t1 = time.time()
                            tel.create_flat(cal_list[cal],fil,amp,binn,red_path,mbias=mbias,log=log)
                            t2 = time.time()
                            log.info('Master flat creation completed in '+str(t2-t1)+' sec')
        if wavelength=='NIR' and not use_dome_flats: #default to use science files for master flat creation
            for fil_list in sky_list:
                fil = fil_list.split('_')[0]
                amp = fil_list.split('_')[1]
                binn = fil_list.split('_')[2]
                process_flat = True
                if skip_red:
                    log.info('User input to skip reduction.')
                    master_flat = tel.get_mflat_name(flat_path, fil, amp, binn)
                    if os.path.exists(master_flat):
                        log.info('Found previous master flat for filter '+fil+', '+amp+' amps and '+binn+' binning.')
                        process_flat = False
                    else:
                        log.info('No master flat found for filter '+fil+', '+amp+' amps and '+binn+' binning, creating master flat.')
                if process_flat:
                    if tel.bias():
                        log.info('Loading master bias.')
                        try:
                            mbias = tel.load_bias(red_path,amp,binn)
                        except:
                            log.error('No master bias found for this configuration, skipping master flat creation for filter '+fil+', '+amp+' amps and '+binn+' binning.')
                            continue
                    else:
                        mbias = None
                    flat_type = 'sky'
                    log.info('Using science files to create master flat')
                    t1 = time.time()
                    tel.create_flat(sky_list[fil_list],fil,amp,binn,red_path,mbias=mbias,log=log)
                    t2 = time.time()
                    log.info('Master flat creation completed in '+str(t2-t1)+' sec')

    # Science files
    ##################
    # Basic processing
    # Does it even need to run?
    if len(sci_list) == 0:
        log.critical('No science files to process, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)
    
    # Begin processing
    log.info('User input target for reduction: '+str(input_target))
    for tar in sci_list:
        stack = tel.stacked_image(tar,red_path)
        target = tar.split('_')[0]
        fil = tar.split('_')[-3]
        amp = tar.split('_')[-2]
        binn = tar.split('_')[-1]
        if input_target is not None:
            if input_target not in tar:
                continue
            else:
                log.info('Matching target found: '+tar)
        if tel.run_wcs():
            final_stack = [st.replace('.fits','_wcs.fits') for st in stack]
        else:
            final_stack = stack
        process_data = True
        # Skip reduction if requested
        if skip_red:
            log.info('User input to skip reduction.')
            if np.all([os.path.exists(st) for st in final_stack]):
                process_data = False
            else:
                log.error('Missing stacks, processing data.')
        # Proceed otherwise
        if process_data:

            # Load bias frame
            if tel.bias():
                log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(red_path,amp,binn)
                except:
                    log.error('No master bias found for this configuration, skipping reduction for: '+tar)
                    continue
            else:
                mbias = None

            # Load flat frame
            log.info('Loading master flat.')
            master_flat = tel.get_mflat_name(flat_path, fil, amp, binn)
            if not os.path.exists(master_flat):
                log.error('No master flat present for filter '+fil+', skipping data reduction for '+tar+'. Check data before rerunning.')
                continue
            flat_data = tel.load_flat(master_flat)
            t1 = time.time()
            log.info('Processing data for '+str(tar))

            # Bias subtraction, gain correction and flat correction. and flat fielding
            processed, masks = tel.process_science(sci_list[tar],fil,amp,binn,red_path,mbias=mbias,mflat=flat_data,proc=proc,log=log)
            t2 = time.time()
            log.info('Data processed in '+str(t2-t1)+' sec')

            # Background subtraction
            # IR image?
            if wavelength=='NIR':
                t1 = time.time()
                log.info('NIR data, creating NIR sky maps.')
                for j,n in enumerate(processed):
                    print('Creating map for file %d/%d'%(j+1,len(processed)))
                    time_diff = sorted([(abs(time_list[tar][j]-n2),k) for k,n2 in enumerate(time_list[tar])])
                    sky_list = [sci_list[tar][k] for _,k in time_diff[0:5]]
                    sky_data = [processed[k] for _,k in time_diff[0:5]]
                    sky_mask = [masks[k] for _,k in time_diff[0:5]]
                    sky_masked_data = []
                    for k in range(len(sky_data)):
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

                    sky.write(os.path.join(red_path, os.path.basename(sci_list[tar][j]).replace('.fits','_sky.fits').replace('.gz','').replace('.bz2','')),overwrite=True)
                    processed[j] = n.subtract(sky,propagate_uncertainties=True,handle_meta='first_found')
                t2 = time.time()
                log.info('Sky maps complete and subtracted in '+str(t2-t1)+' sec')
            # Optical image?
            if wavelength=='OPT':
                t1 = time.time()
                # Fringe correction needed?
                if tel.fringe_correction(fil):
                    dimen = len(stack)
                    for m in range(dimen):
                        fringe_data = []
                        if dimen == 1:
                            suffix = ['.fits']
                            for k,n in enumerate(processed):
                                bkg = Background2D(n, (20, 20), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=masks[k], exclude_percentile=80)
                                masked = np.array(n)
                                masked[masks[k]] = bkg.background[masks[k]]
                                fringe_data.append(CCDData(masked,unit=u.electron/u.second))
                            fringe_map = ccdproc.combine(fringe_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median,mask=masks)
                            fringe_map.write(os.path.join(red_path, 'fringe_map_'+fil+'_'+amp+'_'+binn+suffix[m]),overwrite=True)
                            for j,n in enumerate(processed):
                                processed[j] = n.subtract(fringe_map,propagate_uncertainties=True,handle_meta='first_found')
                        else:
                            suffix = [s.replace('_red','') for s in tel.suffix()]
                            for k,n in enumerate(processed[m]):
                                bkg = Background2D(n, (20, 20), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=masks[m][k], exclude_percentile=80)
                                masked = np.array(n)
                                masked[masks[m][k]] = bkg.background[masks[m][k]]
                                fringe_data.append(CCDData(masked,unit=u.electron/u.second))
                            fringe_map = ccdproc.combine(fringe_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median,mask=masks)
                            fringe_map.write(os.path.join(red_path, 'fringe_map_'+fil+'_'+amp+'_'+binn+suffix[m]),overwrite=True)
                            for j,n in enumerate(processed[m]):
                                processed[m][j] = n.subtract(fringe_map,propagate_uncertainties=True,handle_meta='first_found')
                    t2 = time.time()
                    log.info('Fringe correction complete and subtracted in '+str(t2-t1)+' sec')
            
            # Write reduced data to file and stack images with the same pointing
            log.info('Writing out reduced data.')
            dimen = len(stack)
            if dimen == 1:
                suffix = ['_red.fits']
            else:
                log.info('Multiple extensions to stack.')
                suffix = tel.suffix()
            mask = tel.static_mask(proc)
            for k in range(dimen):
                red_list = [os.path.join(red_path, os.path.basename(sci).replace('.fits',suffix[k]).replace('.gz','').replace('.bz2','')) for sci in sci_list[tar]]
                if dimen == 1:
                    for j,process_data in enumerate(processed):
                        process_data.write(red_list[j],overwrite=True)
                else:
                    for j,process_data in enumerate(processed[k]):
                        process_data.write(red_list[j],overwrite=True)
                log.info('Aligning images.')
                aligned_images, aligned_data = align_quads.align_stars(red_list,telescope,hdu=tel.wcs_extension(),mask=mask[k],log=log)
                log.info('Checking qualty of images.')
                stacking_data, mid_time, total_time = quality_check.quality_check(aligned_images, aligned_data, telescope, log)
                log.info('Creating median stack.')
                sci_med = ccdproc.combine(stacking_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)

                sci_med.data = sci_med.data.astype('float32')
                sci_med.header['BITPIX'] = -32
                sci_med.header['MJD-OBS'] = (mid_time, 'Mid-MJD of the observation sequence calculated using DATE-OBS.')
                sci_med.header['EXPTIME'] = (1, 'Effective expsoure tiime for the stack in seconds.')
                sci_med.header['EXPTOT'] = (total_time, 'Total exposure time of stack in seconds')
                sci_med.header['GAIN'] = (len(stacking_data), 'Effecetive gain for stack.')
                sci_med.header['RDNOISE'] = (tel.rdnoise(sci_med.header)/np.sqrt(len(stacking_data)), 'Readnoise of stack.')
                sci_med.header['NFILES'] = (len(stacking_data), 'Number of images in stack')
                sci_med.header['FILTER'] = fil
                sci_med.header['OBSTYPE'] = 'OBJECT'
                sci_med.write(stack[k],overwrite=True)
                log.info('Median stack made for '+stack[k])
                print('Please check the quality of the stack made.')

                # Auto-WCS solution
                if tel.run_wcs():

                    log.info('Solving WCS.')

                    basefile = os.path.basename(stack[k])
                    dirname = os.path.dirname(stack[k])
                    solve_wcs.solve_astrometry(dirname, basefile, tel,
                            replace=True, log=log)

                if tel.run_phot():
                    log.info('Running psf photometry.')
                    try:
                        epsf, fwhm = psf.do_phot(stack[k], log=log)
                        log.info('FWHM = %2.4f"'%(fwhm*tel.pixscale()))
                        log.info('Calculating zeropint.')
                        zp_catalogs = tel.catalog_zp()
                        zp_cal = absphot.absphot()
                        for zp_cat in zp_catalogs:
                            catfile = stack[k].replace('.fits','.pcmp')
                            zp, zp_err = zp_cal.find_zeropoint(catfile, fil, 
                                zp_cat, plot=True, log=log)
                            if zp:
                                break
                    except Exception as e:
                        log.error('PSF photometry failed due to: '+str(e))
        
    t_end = time.time()
    log.info('Pipeline finshed.')
    log.info('Total runtime: '+str(t_end-t_start)+' sec')


def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('instrument', 
        default=None, 
        choices=['Binospec','DEIMOS','GMOS','LRIS','MMIRS','MOSFIRE','TEST'],
        help='''Name of instrument (must be in params folder) of data to 
        reduce. Required to run pipeline. Use TEST to run the pipeline test.''')
    params.add_argument('data_path', 
        default=None, 
        help='''Path of data to reduce. See manual for specific details. 
        Required to run pipeline.''')
    params.add_argument('--use_dome_flats', 
        type=str, 
        default=None, 
        help='''Use dome flats instead of science images to create master flat 
        for NIR. Optional parameter. Default is False.''')
    params.add_argument('--skip_red', 
        type=str, 
        default=None, 
        help='''Option to skip reduction i.e. the creation of master files used 
        for calibration or the stacked image for science targets. See manual 
        for specific details. Optional parameter.''')
    params.add_argument('--target', 
        type=str, 
        default=None, 
        help='''Option to only reduce a specific target. String used here must 
        be contained within the target name in file headers. Optional 
        parameter.''')
    params.add_argument('--proc', 
        type=str, 
        default=True, 
        help='''Option to use the _proc data from MMT. Optional parameter. 
        Default is False.''')
    params.add_argument('--cal_path', 
        type=str, 
        default=None, 
        help='''Full path of the folder containing calibrations if different 
        from the default. Optional parameter.''')
    params.add_argument('--phot', 
        type=str, 
        default=None, 
        help='''Option to perform aperture photometry on stacked image. See 
        manual for specific details. Optional parameter.''')
    params.add_argument('--reset', 
        type=str, 
        default=None, 
        help='''Option to reset data files sorted previously by pipeline. 
        Optional parameter.''')
    params.add_argument('--use-anet',
        default=False,
        action='store_true',
        help='''Use astrometry.net in initial WCS solver.''')
    params.add_argument('--no-verify',
        default=False,
        action='store_true',
        help='''Do not prompt the user for input.''')
    args = params.parse_args()

    main_pipeline(args.instrument,args.data_path,args.cal_path,
        input_target=args.target,skip_red=args.skip_red,proc=args.proc,
        use_dome_flats=args.use_dome_flats,phot=args.phot,reset=args.reset,
        use_anet=args.use_anet,no_verify=args.no_verify)

if __name__ == "__main__":
    main()
