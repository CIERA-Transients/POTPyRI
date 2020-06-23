#!/usr/bin/env python

"Python pipeline to reduce BINOSPEC (imaging) data."
"Assumes single object per file."
"Pipeline is based on the IDL pipeline by Igor Chilingarian"
"(Kansky J., et al., 2019, ascl.soft, ascl:1905.004"
"https://bitbucket.org/chil_sai/binospec/src/master/)."
"Author: Kerry Paterson."
"This project was funded by the National Science Foundation: AST-1814782 and AST-1909358"
"If you use this code for your work, please consider citing."

import os
import sys
import datetime
import time
import glob
import argparse
import logging
import shutil
import ccdproc
import numpy as np
import astropy.wcs as wcs
import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData
from photutils import MeanBackground
from astropy.stats import SigmaClip
from astropy.modeling import models
from Sort_files import sort_files_bino
from Find_zero_point import find_zero_point
from Find_target_phot import find_target_phot
from align import align_stars
from astropy.utils.exceptions import AstropyWarning
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

params = argparse.ArgumentParser(description='Path of data.')
params.add_argument('data_path', default=None, help='Path to data, assumes it is in the default format from KOA.')
params.add_argument('--proc', type=str, default=None, help='If working with the _proc data from MMT.')
params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.')
params.add_argument('--phot', type=str, default=None, help='Option to perform photometry.')
args = params.parse_args()

gain = [1.085, 1.04649118, 1.04159151, 0.97505369, 1.028, 1.16341855, 1.04742053, 1.0447564]

raw_path = args.data_path+'/raw/' #path containing the raw data
if not os.path.exists(raw_path): #create reduced file path if it doesn't exist
    os.makedirs(raw_path)
if not os.path.exists(raw_path+'bad/'): #create reduced file path if it doesn't exist
    os.makedirs(raw_path+'bad/')
red_path = args.data_path+'/red/' #path to write the reduced files
if not os.path.exists(red_path): #create reduced file path if it doesn't exist
    os.makedirs(red_path)

log_file_name = red_path+'BINOSPEC_log_'+datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')+'.log' #create log file name
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

if args.proc == 'True' or args.proc == 'yes':
    files = glob.glob(args.data_path+'sci_img_*proc.fits')
else:
    files_all = glob.glob(args.data_path+'sci_img_*.fits')
    files = [x for x in files_all if 'proc' not in x]

if len(files) != 0:
    cal_list, sci_list = sort_files_bino(files,raw_path)
    for i,f in enumerate(cal_list):
        shutil.move(f,raw_path)
        cal_list[i] = f.replace(args.data_path,raw_path)
    for i,f in enumerate(sci_list):
        shutil.move(f,raw_path)
        sci_list[i] = f.replace(args.data_path,raw_path)
else:
    if args.proc == 'True' or args.proc == 'yes':
        files = glob.glob(raw_path+'sci_img_*proc.fits')
    else:
        files_all = glob.glob(raw_path+'sci_img_*.fits')
        files = [x for x in files_all if 'proc' not in x]
    if len(files) != 0:
        cal_list, sci_list = sort_files_bino(files,raw_path)
    else:
        log.critical('No files found, please check data path and rerun.')
        logging.shutdown()
        sys.exit(-1)

if len(sci_list) != 0:
    process_data = True
    with fits.open(sci_list[0]) as hdr:
        header = hdr[1].header
    target = header['OBJECT']
    fil = header['FILTER'].split('_')[0]
    if args.skip_red == 'True' or args.skip_red == 'yes':
        if os.path.exists(red_path+target+'_left.fits') and os.path.exists(red_path+target+'_right.fits'):
            with fits.open(red_path+target+'_left.fits',mode='update') as hdr_left:
                sci_left = hdr_left[0].data
            with fits.open(red_path+target+'_right.fits',mode='update') as hdr_right:
                sci_right = hdr_right[0].data
            process_data = False
        else:
            log.error('No reduced image found, processing data.')
    if process_data:
        flat_left = CCDData.read('./BINOSPEC_calib/mflat_'+fil+'_left.fits', hdu=1, unit=u.electron)
        flat_right = CCDData.read('./BINOSPEC_calib/mflat_'+fil+'_right.fits', hdu=1, unit=u.electron)
        left_list = []
        right_list = []
        for j,sci in enumerate(sci_list):
            bkg = MeanBackground(SigmaClip(sigma=3.))
            if args.proc == 'True' or args.proc == 'yes':
                left = CCDData.read(sci, hdu=1, unit=u.electron)
                left = ccdproc.flat_correct(left,flat_left)
                bkg_value = bkg.calc_background(left.data)
                left.subtract(bkg_value*u.electron,propagate_uncertainties=True,handle_meta='first_found')
                left.divide(left.header['EXPTIME'],propagate_uncertainties=True,handle_meta='first_found')
                left.write(red_path+os.path.basename(sci).replace('.fits','_left.fits'),overwrite=True)
                left_list.append(left)
                right = CCDData.read(sci, hdu=2, unit=u.electron)
                right = ccdproc.flat_correct(right,flat_right)
                bkg_value = bkg.calc_background(left.data)
                right.subtract(bkg_value*u.electron,propagate_uncertainties=True,handle_meta='first_found')
                right.divide(right.header['EXPTIME'],propagate_uncertainties=True,handle_meta='first_found')
                right.write(red_path+os.path.basename(sci).replace('.fits','_right.fits'),overwrite=True)
                right_list.append(right)
            else:
                with fits.open(sci) as hdr:
                    header_left = hdr[1].header
                    header_right = hdr[6].header
                data_list = []
                for i in range(8):
                    data = ccdproc.CCDData.read(sci,hdu=i+1,unit=u.adu)
                    red = ccdproc.ccd_process(data, oscan=data[:,0:50], oscan_model=models.Chebyshev1D(3), trim='[1200:2098,210:2056]', gain=gain[i]*u.electron/u.adu, readnoise=4*u.electron)
                    data_list.append(np.asarray(red).astype(np.float32))
                top_left = np.concatenate([data_list[0],np.fliplr(data_list[1])],axis=1)
                bot_left = np.flipud(np.concatenate([data_list[3],np.fliplr(data_list[2])],axis=1))
                left = CCDData(np.concatenate([top_left,bot_left]),unit=u.electron,header=header_left,wcs=wcs.WCS(header_left))
                left = ccdproc.flat_correct(left,flat_left[209:3903,1199:2997])
                bkg_value = bkg.calc_background(left.data)
                left.subtract(bkg_value*u.electron,propagate_uncertainties=True,handle_meta='first_found')
                left.divide(left.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
                left.header['DATASEC'] = '[1:1798,1:3694]'
                left.write(red_path+os.path.basename(sci).replace('.fits','_left.fits'),overwrite=True)
                left_list.append(left)
                top_right = np.concatenate([data_list[6],np.fliplr(data_list[7])],axis=1)
                bot_right = np.flipud(np.concatenate([data_list[5],np.fliplr(data_list[4])],axis=1))
                right = CCDData(np.concatenate([top_right,bot_right]),unit=u.electron,header=header_right,wcs=wcs.WCS(header_right))
                right = ccdproc.flat_correct(right,flat_right[209:3903,1199:2997])
                bkg_value = bkg.calc_background(right.data)
                right.subtract(bkg_value*u.electron,propagate_uncertainties=True,handle_meta='first_found')
                right.divide(right.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
                right.header['DATASEC'] = '[1:1798,1:3694]'
                right.write(red_path+os.path.basename(sci).replace('.fits','_right.fits'),overwrite=True)
                right_list.append(right)
        align_left = align_stars(left_list)
        com_left = [CCDData(x,unit=u.electron/u.second,header=left_list[0].header,wcs=left_list[0].wcs) for x in align_left] 
        sci_left = ccdproc.combine(com_left,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
        sci_left.write(red_path+target+'_left.fits',overwrite=True)
        with fits.open(red_path+target+'_left.fits',mode='update') as hdr_left:
            hdr_left[0].header['RDNOISE'] = 4.0/np.sqrt(len(align_left))
            hdr_left[0].header['NFILES'] = len(align_left)
        align_right = align_stars(right_list)
        com_right = [CCDData(x,unit=u.electron/u.second,header=right_list[0].header,wcs=right_list[0].wcs) for x in align_right]       
        sci_right = ccdproc.combine(com_right,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
        sci_right.write(red_path+target+'_right.fits',overwrite=True)
        with fits.open(red_path+target+'_right.fits',mode='update') as hdr_right:
            hdr_right[0].header['RDNOISE'] = 4.0/np.sqrt(len(align_right))
            hdr_right[0].header['NFILES'] = len(align_right)
    if args.phot == 'True' or args.phot == 'yes':
        side = input('Perform photometry on left or right image? (type left or right) ')
        if side == 'left':
            stack = red_path+target+'_left.fits'
        elif side == 'right':
            stack = red_path+target+'_right.fits'
        zp, zp_err, fwhm = find_zero_point(stack, fil, 0, 'SDSS', log=log)
        method = input('Use RA and Dec or pixel coords to find source? (type radec or xy) ')
        if method == 'radec':
            ra = input('Enter RA of source in degrees: ')
            dec = input('Enter Dec of source in degrees: ')
            mag, mag_err = find_target_phot(stack, fil, fwhm=fwhm, zp=zp, zp_err=zp_err, show_phot=True, log=log, log2=log, ra=ra, dec=dec)
        elif method == 'xy':
            x = input('Enter x position of source in pixels: ')
            y = input('Enter y position of source in pixels: ')
            mag, mag_err = find_target_phot(stack, fil, fwhm=fwhm, zp=zp, zp_err=zp_err, show_phot=True, log=log, log2=log, x=x, y=y)