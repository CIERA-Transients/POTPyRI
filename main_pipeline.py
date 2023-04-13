#!/usr/bin/env python

"Automatic generalized pipeline for imaging reduction. Creates median coadded images for each target."
"Individual images should be checked and removed from sci path."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "1.21" #last updated 24/02/2023

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

# Update to new photutils methods and deprecate previous Find_target_phot when possible
try:
    import Find_target_phot as tp
except ImportError:
    import Find_target_phot_new as tp

import extinction
from utilities.util import *
from colorama import init, Fore, Back, Style
init()

def main_pipeline(telescope,data_path,cal_path=None,input_target=None,skip_red=None,proc=None,use_dome_flats=None,phot=None,reset=None):
    #start time
    t_start = time.time()
    #import telescope parameter file
    global tel
    try:
        tel = importlib.import_module('tel_params.'+telescope)
    except ImportError:
        telfile='tel_params.'+telescope
        print(f'No such telescope file {telfile}, please check that you have entered the'+\
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

    # Get path to code directory
    abspath = os.path.abspath(sys.argv[0])
    code_dir = os.path.split(abspath)[0]

    # Copy config directory to data_path
    config_dir = os.path.join(code_dir, 'Config')
    if not os.path.exists(os.path.join(data_path, 'Config')):
        shutil.copytree(config_dir, os.path.join(data_path, 'Config'))


    if cal_path is not None:
        cal_path = cal_path
    else:
        cal_path = tel.cal_path()
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

    if reset is not None:
        if os.path.exists(data_path+'/file_list.txt'):
            os.remove(data_path+'/file_list.txt')
        if reset=='all':
            files = glob.glob(os.path.join(raw_path,'*'))+\
                    glob.glob(os.path.join(bad_path,'*'))+\
                    glob.glob(os.path.join(spec_path,'*'))
        if reset=='raw':
            files = glob.glob(os.path.join(raw_path,'*'))
        for f in files:
            shutil.move(f,data_path)

    # CDK - added editing for headers with LRIS raw data files (*[b,r]*.fits)
    if ((telescope=='LRIS' and proc and str(proc)=='raw') or
        (telescope=='DEIMOS' and proc and str(proc)=='raw')):
        print('Edit raw headers')
        tel.edit_raw_headers(data_path)

    if os.path.exists(data_path+'/file_list.txt'):
        log.info('Previous file list exists, loading lists.')
        cal_list, sci_list, sky_list, time_list = Sort_files.load_files(data_path+'/file_list.txt', telescope,log)
    else:
        log.info('Sorting files and creating file lists.')
        # CDK - updated this to a os.path.join so it doesn't require trailing /
        files = sorted(glob.glob(os.path.join(data_path,tel.raw_format(proc))))
        if len(files) != 0:
            log.info(str(len(files))+' files found.')
            cal_list, sci_list, sky_list, time_list = Sort_files.sort_files(files,telescope,data_path,log)
            print('Please check the created file list and make any nessecary edits (e.g. remove bad files mentioned by the observing log etc.)')
            con = input(Back.GREEN+'Did you edit the file list (yes or no)? '+Style.RESET_ALL)
            if con=='yes':
                log.info('Edits made to the file list by the user, reloading file list.')
                cal_list, sci_list, sky_list, time_list = Sort_files.load_files(data_path+'/file_list.txt', telescope,log)
        else:
            log.critical('No files found, please check data path and rerun.')
            logging.shutdown()
            sys.exit(-1)

    # CDK - added editing for headers with LRIS raw data files (*[b,r]*.fits)
    if ((telescope=='LRIS' and proc and str(proc)=='raw') or
        (telescope=='DEIMOS' and proc and str(proc)=='raw')):
        tel.edit_raw_headers(raw_path)

    if tel.bias():
        bias_num = 0
        for cal in cal_list:
            if 'BIAS' in cal:
                bias_num += 1
                process_bias = True
                amp = cal.split('_')[1]
                binn = cal.split('_')[2]
                if skip_red:
                    log.info('User input to skip reduction.')
                    if os.path.exists(red_path+'mbias_'+amp+'_'+binn+'.fits'):
                        log.info('Found previous master bias.')
                        process_bias = False
                    else:
                        log.error('No master bias found, creating master bias.')
                if process_bias:
                    t1 = time.time()
                    log.info('Processing bias files.')
                    tel.create_bias(cal_list,cal,red_path,log)
                    t2 = time.time()
                    log.info('Master bias creation completed in '+str(t2-t1)+' sec')
        if bias_num==0:
            log.critical('No bias files present, check data before rerunning.')
            logging.shutdown()
            sys.exit(-1)

    if tel.dark():
        for cal in cal_list:
            if 'DARK' in cal:
                process_dark = True
                exp = cal.split('_')[1]
                amp = cal.split('_')[2]
                binn = cal.split('_')[3]
                if skip_red:
                    log.info('User input to skip reduction.')
                    if os.path.exists(red_path+cal+'.fits'):
                        log.info('Found previous master dark.')
                        process_dark = False
                    else:
                        log.error('No master dark found, creating master dark.')
                if process_dark:
                    if tel.bias():
                        log.info('Loading master bias.')
                        try:
                            mbias = tel.load_bias(red_path,amp,binn)
                        except:
                            log.error('No master bias found for this configuration, skipping master dark creation for exposure '+exp+', '+amp+' amps and '+binn+' binning.')
                            continue
                    else:
                        mbias = None
                    t1 = time.time()
                    tel.create_dark(cal_list,cal,mbias,red_path,log)
                    t2 = time.time()
                    log.info('Master dark creation completed in '+str(t2-t1)+' sec')

    if tel.flat():
        for cal in cal_list:
            if 'FLAT' in cal:
                process_flat = True
                fil = cal.split('_')[1]
                amp = cal.split('_')[2]
                binn = cal.split('_')[3]
                if skip_red:
                    log.info('User input to skip reduction.')
                    master_flat = tel.flat_name(flat_path, fil, amp, binn)
                    if np.all([os.path.exists(mf) for mf in master_flat]):
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
                    master_flat = tel.flat_name(flat_path, fil, amp, binn)
                    if np.all([os.path.exists(mf) for mf in master_flat]):
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

    if len(sci_list) == 0:
        log.critical('No science files to process, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)
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
        if skip_red:
            log.info('User input to skip reduction.')
            if np.all([os.path.exists(st) for st in final_stack]):
                process_data = False
            else:
                log.error('Missing stacks, processing data.')
        if process_data:
            if tel.bias():
                log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(red_path,amp,binn)
                except:
                    log.error('No master bias found for this configuration, skipping reduction for: '+tar)
                    continue
            else:
                mbias = None
            log.info('Loading master flat.')
            master_flat = tel.flat_name(flat_path, fil, amp, binn)
            if not np.all([os.path.exists(mf) for mf in master_flat]):
                log.error('No master flat present for filter '+fil+', skipping data reduction for '+tar+'. Check data before rerunning.')
                continue
            flat_data = tel.load_flat(master_flat)
            t1 = time.time()
            log.info('Processing data for '+str(tar))
            processed, masks = tel.process_science(sci_list[tar],fil,amp,binn,red_path,mbias=mbias,mflat=flat_data,proc=proc,log=log)
            t2 = time.time()
            log.info('Data processed in '+str(t2-t1)+' sec')
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
                    sky.write(red_path+os.path.basename(sci_list[tar][j]).replace('.fits','_sky.fits').replace('.gz','').replace('.bz2',''),overwrite=True)
                    processed[j] = n.subtract(sky,propagate_uncertainties=True,handle_meta='first_found')
                t2 = time.time()
                log.info('Sky maps complete and subtracted in '+str(t2-t1)+' sec')
            if wavelength=='OPT':
                t1 = time.time()
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
                            fringe_map.write(red_path+'fringe_map_'+fil+'_'+amp+'_'+binn+suffix[m],overwrite=True)
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
                            fringe_map.write(red_path+'fringe_map_'+fil+'_'+amp+'_'+binn+suffix[m],overwrite=True)
                            for j,n in enumerate(processed[m]):
                                processed[m][j] = n.subtract(fringe_map,propagate_uncertainties=True,handle_meta='first_found')
                    t2 = time.time()
                    log.info('Fringe correction complete and subtracted in '+str(t2-t1)+' sec')
            log.info('Writing out reduced data.')
            dimen = len(stack)
            if dimen == 1:
                suffix = ['_red.fits']
            else:
                log.info('Multiple extensions to stack.')
                suffix = tel.suffix()
            mask = tel.static_mask(proc)
            for k in range(dimen):
                red_list = [red_path+os.path.basename(sci).replace('.fits',suffix[k]).replace('.gz','').replace('.bz2','') for sci in sci_list[tar]]
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
                sci_med.header['MJD-OBS'] = (mid_time, 'Mid-MJD of the observation sequence calculated using DATE-OBS.')
                sci_med.header['EXPTIME'] = (1, 'Effective expsoure tiime for the stack in seconds.')
                sci_med.header['EXPTOT'] = (total_time, 'Total exposure time of stack in seconds')
                sci_med.header['GAIN'] = (len(stacking_data), 'Effecetive gain for stack.')
                sci_med.header['RDNOISE'] = (tel.rdnoise(sci_med.header)/np.sqrt(len(stacking_data)), 'Readnoise of stack.')
                sci_med.header['NFILES'] = (len(stacking_data), 'Number of images in stack')
                sci_med.write(stack[k],overwrite=True)
                log.info('Median stack made for '+stack[k])
                print('Please check the quality of the stack made.')
                con = input(Back.GREEN+'Would you like to continue with the reduction (solving the WCS) for this stack (yes or no)? '+Style.RESET_ALL)
                if con=='no':
                    continue
                if tel.run_wcs():
                    log.info('Solving WCS.')
                    try:
                        wcs_error = solve_wcs.solve_wcs(stack[k],telescope,log=log)
                        log.info(wcs_error)
                        stack[k] = stack[k].replace('.fits','_wcs.fits')
                    except:
                        log.error('Automatic WCS solution failed.')
                    redo_wcs = input(Back.GREEN+'Please review the WCS plots and errors. Do you wish to manually redo wcs (yes or no)?  '+Style.RESET_ALL)
                    while redo_wcs=='yes':
                        log.info('Manually redoing WCS.')
                        cat = input(Fore.RED+'Please enter the name of the catalog you wish to use. The options are "gaia" (GAIA DR3), "sdssdr12" (SDSSDR12), "2mass" (2MASS). '+Style.RESET_ALL)
                        if cat == 'gaia':
                            cat_header, cat_stars = import_catalog(stack[k].replace('_wcs','').replace('.fits','_wcs.gaia'))
                            cat_stars_ra = cat_stars['ra']
                            cat_stars_dec = cat_stars['dec']
                        else:
                            with fits.open(stack[k].replace('_wcs','')) as hdr:
                                header = hdr[0].header
                                data = hdr[0].data
                            ra_center, dec_center = (wcs.WCS(header)).all_pix2world(np.shape(data)[1]/2,np.shape(data)[0]/2,1)
                            coord = SkyCoord(ra_center, dec_center,unit=u.deg)
                            cat_stars = search_catalogs(coord, [cat], search_radius=10*u.arcmin, outfile=stack[k].replace('_wcs','').replace('.fits','_wcs.'+cat))
                            cat_stars_ra = cat_stars[cat+'_'+viziercat[cat]['columns'][0]]
                            cat_stars_dec = cat_stars[cat+'_'+viziercat[cat]['columns'][1]]
                        if len(cat_stars) == 0:
                            log.info('No stars found in this catalog.')
                            continue
                        else:
                            wcs_error = solve_wcs.man_wcs(telescope, stack[k].replace('_wcs',''), cat, cat_stars_ra, cat_stars_dec)
                            log.info(wcs_error)
                            redo_wcs = input(Back.GREEN+'Please review the WCS plots and errors. Do you wish to manually redo wcs (yes or no)?  '+Style.RESET_ALL)
                    try:
                        if wcs_error!=0:
                            stack[k] = stack[k].replace('_wcs','').replace('.fits','_wcs.fits')
                    except UnboundLocalError:
                        log.error('WCS solution could not be calculated. Skipping the rest of the reduction for this target.')
                        continue
                if tel.run_phot():
                    con = input(Back.GREEN+'Based on the WCS solution, would you like to continue with the PSF photometry for this stack (yes or no)? '+Style.RESET_ALL)
                    if con=='no':
                        continue
                    log.info('Running psf photometry.')
                    try:
                        epsf, fwhm = psf.do_phot(stack[k])
                        log.info('FWHM = %2.4f"'%(fwhm*tel.pixscale()))
                        con = input(Back.GREEN+'Please check the PSF determined. Would you like to continue with the reduction (calculating zero point) for this stack (yes or no)? '+Style.RESET_ALL)
                        if con=='yes':
                            log.info('Calculating zeropint.')
                            zp_catalogs = tel.catalog_zp()
                            zp_cal = absphot.absphot()
                            for zp_cat in zp_catalogs:
                                zp, zp_err = zp_cal.find_zeropoint(stack[k].replace('.fits','.pcmp'), fil, zp_cat, plot=True, log=log)
                                if zp:
                                    break
                    except Exception as e:
                        log.error('PSF photometry failed due to: '+str(e))
        if phot:
            log.info('User input to perform manual aperture photometry.')
            log.info('List of final stacks: '+str(final_stack))
            k = int(input(Back.GREEN+'Index (starting from 0) of the stack you want to perform aperture photometry on? '+Style.RESET_ALL))
            log.info('Performing aperture photometry on '+final_stack[k])
            if not os.path.exists(final_stack[k].replace('.fits','.pcmp')):
                log.info('Running psf photometry.')
                epsf, fwhm = psf.do_phot(final_stack[k], log=log)
                log.info('Calculating zeropint.')
                zp_catalogs = tel.catalog_zp()
                zp_cal = absphot.absphot()
                for zp_cat in zp_catalogs:
                    zp, zp_err = zp_cal.find_zeropoint(final_stack[k].replace('.fits','.pcmp'), fil, zp_cat, plot=True, log=log)
                    if zp:
                        break
            log.info('Loading FWHM from psf photometry.')
            header, table = import_catalog(final_stack[k].replace('.fits','.pcmp'))
            fwhm = header['FWHM']
            log.info('FWHM = %2.4f pixels'%fwhm)
            enter_zp = input(Back.GREEN+'Enter user zeropiont instead of loading from psf photometry ("yes" or "no")? '+Style.RESET_ALL)
            if enter_zp == 'yes':
                zp = float(input(Fore.RED+'Please enter zeropoint in AB mag: '+Style.RESET_ALL))
                zp_err = float(input(Fore.RED+'Please enter zeropoint error in AB mag: '+Style.RESET_ALL))
                log.info('User entered zpt = %2.4f +/- %2.4f AB mag'%(zp,zp_err))
            else:
                try:
                    zp = header['ZPTMAG']
                    zp_err = header['ZPTMUCER']
                    log.info('zpt = %2.4f +/- %2.4f AB mag'%(zp,zp_err))
                except:
                    log.info('No zeropint found.')
                    zp = float(input(Fore.RED+'Please enter zeropoint in AB mag: '+Style.RESET_ALL))
                    zp_err = float(input(Fore.RED+'Please enter zeropoint error in AB mag: '+Style.RESET_ALL))
                    log.info('User entered zpt = %2.4f +/- %2.4f AB mag'%(zp,zp_err))
            with fits.open(final_stack[k]) as hdr:
                gain_eff = hdr[0].header['GAIN']
            pos = input(Back.GREEN+'Would you like to enter the RA and Dec ("wcs") or x and y ("xy") position of the target? '+Style.RESET_ALL)
            if pos == 'wcs':
                ra = float(input(Fore.RED+'Enter RA in degrees: '+Style.RESET_ALL))
                dec = float(input(Fore.RED+'Enter Dec in degrees: '+Style.RESET_ALL))
                tp.find_target_phot(final_stack[k], fil, fwhm, zp, zp_err, tel.pixscale(), show_phot=True, ra=ra, dec=dec, log=log, gain=gain_eff)
            elif pos == 'xy':
                x = float(input(Fore.RED+'Enter x position in pixel: '+Style.RESET_ALL))
                y = float(input(Fore.RED+'Enter y position in pixel: '+Style.RESET_ALL))
                tp.find_target_phot(final_stack[k], fil, fwhm, zp, zp_err, tel.pixscale(), show_phot=True, x=x, y=y, log=log, gain=gain_eff)
                ra, dec = (wcs.WCS(header)).all_pix2world(x,y,1)
            log.info('Calculating extinction correction from Schlafly & Finkbeiner (2011).')
            coords = SkyCoord(ra,dec,unit='deg')
            ext_cor = extinction.calculate_mag_extinction(coords,fil)
            log.info('Galactic extinction = %2.4f mag'%ext_cor)
    t_end = time.time()
    log.info('Pipeline finshed.')
    log.info('Total runtime: '+str(t_end-t_start)+' sec')


def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('telescope', default=None, help='Name of instrument (must be in tel_params folder) of data to reduce. Required to run pipeline. Use TEST to run the pipeline test.')
    params.add_argument('data_path', default=None, help='Path of data to reduce. See manual for specific details. Required to run pipeline.')
    params.add_argument('--use_dome_flats', type=str, default=None, help='Use dome flats instead of science images to create master flat for NIR. Optional parameter. Default is False.')
    params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction i.e. the creation of master files used for calibration or the stacked image for science targets. See manual for specific details. Optional parameter.')
    params.add_argument('--target', type=str, default=None, help='Option to only reduce a specific target. String used here must be contained within the target name in file headers. Optional parameter.')
    params.add_argument('--proc', type=str, default=True, help='Option to use the _proc data from MMT. Optional parameter. Default is False')
    params.add_argument('--cal_path', type=str, default=None, help='Full path of the folder containing calibrations if different from the default. Optional parameter.')
    params.add_argument('--phot', type=str, default=None, help='Option to perform aperture photometry on stacked image. See manual for specific details. Optional parameter.')
    params.add_argument('--reset', type=str, default=None, help='Option to reset data files sorted previously by pipeline. Optional parameter.')
    args = params.parse_args()

    main_pipeline(args.telescope,args.data_path,args.cal_path,input_target=args.target,skip_red=args.skip_red,proc=args.proc,use_dome_flats=args.use_dome_flats,phot=args.phot,reset=args.reset)

if __name__ == "__main__":
    main()
