#!/usr/bin/env python

"Python pipeline to reduce LRIS (imaging) data."
"Author: Kerry Paterson."
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

import re
import sys
import numpy as np
import os
import datetime
import time
import astropy
import astropy.units.astrophys as u
from astropy.io import fits
import ccdproc
import glob
import argparse
import subprocess
import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, ZScaleInterval
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.io import ascii
import astropy.wcs as wcs
import shutil
from astropy.nddata import CCDData
from astropy.modeling import models
from photutils import make_source_mask, MeanBackground, StdBackgroundRMS, CircularAperture, CircularAnnulus, aperture_photometry, Background2D, MedianBackground, DAOStarFinder
from ccdproc import wcs_project
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.utils import calc_total_error
from astroquery.irsa import Irsa
from align import align_stars

params = argparse.ArgumentParser(description='Path of data.')
params.add_argument('data_path', default=None, help='Path of data, in default format from Gemini.') #path to GMOS data: required
params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.') #
params.add_argument('--target', type=str, default=None, help='Option to only reduce this target.') #
params.add_argument('--individual', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--wcs', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--phot', type=str, default=None, help='Option to perform photometry.') #must have pyraf install and access to IRAF to use
args = params.parse_args()

def sort_files(files): #sort the calibration files: 
    cal_list = {}
    cal_list.update({'BIAS4':[]})
    cal_list.update({'BIAS8':[]})
    sci_list = {}
    for f in files:
        with fits.open(f) as file_open:
            hdr = file_open[0].header
            target = hdr['OBJECT'].replace(' ','')
            fil = hdr['DWFILNAM']
            if hdr['AMPMODE']=='SINGLE:B':
                amps = '4'
            elif hdr['AMPMODE']=='DUAL:A+B':
                amps = '8'
            if hdr['MOSMODE']=='Direct':
                if hdr['OBSTYPE']=='Bias' or target=='Bias':
                    cal_list['BIAS'+amps].append(f)
                elif hdr['OBSTYPE']=='DmFlat' or target=='Domeflat':
                    try:
                        cal_list[fil+amps]
                    except KeyError:
                        cal_list.update({fil+amps:[]})
                    cal_list[fil+amps].append(f)
                else:#if hdr['OBSTYPE']=='Object': #IntFlat because of single amp mode?
                    try:
                        sci_list[target+'_'+fil+amps]
                    except KeyError:
                        sci_list.update({target+'_'+fil+amps:[]})
                    sci_list[target+'_'+fil+amps].append(f)
                # else: #other things?
                #     pass
            else:
                pass #spec data
    return cal_list, sci_list

raw_path = args.data_path+'/raw/' #path containing the raw data
bad_path = raw_path+'bad/'
cal_path = raw_path+'cal/' #path containing the calibration files
sci_path = raw_path+'sci/' #path containing the scie files
red_path = args.data_path+'/red/' #path to write the reduced files
if not os.path.exists(raw_path): #create reduced file path if it doesn't exist
    os.makedirs(raw_path)
if not os.path.exists(bad_path): #create reduced file path if it doesn't exist
    os.makedirs(bad_path)
if not os.path.exists(cal_path): #create reduced file path if it doesn't exist
    os.makedirs(cal_path)
if not os.path.exists(sci_path): #create reduced file path if it doesn't exist
    os.makedirs(sci_path)
if not os.path.exists(red_path): #create reduced file path if it doesn't exist
    os.makedirs(red_path)

log_file_name = red_path+'DEIMOS_log_'+datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')+'.log' #create log file name
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

gain = 1.2
rdnoise = 2.5

files = glob.glob(args.data_path+'/*.fits')
if len(files) != 0:
    cal_list, sci_list = sort_files(files)
    for cal in cal_list:
        for i,f in enumerate(cal_list[cal]):
            shutil.move(f,cal_path)
            cal_list[cal][i] = f.replace(args.data_path,cal_path)
    for sci in sci_list:
        for i,f in enumerate(sci_list[sci]):
            shutil.move(f,sci_path)
            sci_list[sci][i] = f.replace(args.data_path,sci_path)
else:
    files = glob.glob(cal_path+'*.fits')+glob.glob(sci_path+'*.fits')
    cal_list, sci_list = sort_files(files)

if len(cal_list['BIAS4'])+len(cal_list['BIAS8']) == 0:
    log.critical('No bias data present, exiting.')
    logging.shutdown()
    sys.exit(-1)

def create_bias(cal_list):
    for i, bias_file in enumerate(cal_list['BIAS4']):
        with fits.open(bias_file) as hdr:
            header = hdr[0].header
        bias_raw = [CCDData.read(bias_file, hdu=x+1, unit='adu') for x in range(4)]
        bias_processed = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim='[13:2060,1:2601]', gain=gain*u.electron/u.adu, readnoise=rdnoise*u.electron) for j,x in enumerate(bias_raw)]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in bias_processed: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias_file.replace(cal_path,red_path),overwrite=True)
        cal_list['BIAS4'][i] = bias_file.replace(cal_path,red_path)
    if len(cal_list['BIAS4']) != 0:
        mbias4 = [ccdproc.combine(cal_list['BIAS4'],hdu=x+1,unit=u.electron) for x in range(4)]
        mbias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in mbias4: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        mbias_hdu.writeto(red_path+'mbias4.fits',overwrite=True)
        for bias in cal_list['BIAS4']: os.remove(bias)
    else:
        mbias4 = None
    for i, bias_file in enumerate(cal_list['BIAS8']):
        with fits.open(bias_file) as hdr:
            header = hdr[0].header
        bias_raw = [CCDData.read(bias_file, hdu=x+1, unit='adu') for x in range(8)]
        bias_processed = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim='[13:1036,1:2601]', gain=gain*u.electron/u.adu, readnoise=rdnoise*u.electron) for j,x in enumerate(bias_raw)]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in bias_processed: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias_file.replace(cal_path,red_path),overwrite=True)
        cal_list['BIAS8'][i] = bias_file.replace(cal_path,red_path)
    if len(cal_list['BIAS8']) != 0:
        mbias8 = [ccdproc.combine(cal_list['BIAS8'],hdu=x+1,unit=u.electron) for x in range(8)]
        mbias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in mbias8: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        mbias_hdu.writeto(red_path+'mbias8.fits',overwrite=True)
        for bias in cal_list['BIAS8']: os.remove(bias)
    else:
        mbias8 = None
    return mbias4, mbias8

if args.skip_red == 'True' or args.skip_red == 'yes':
    process_bias = False
    mbias4, mbias8 = None, None
else:
    process_bias = True
if process_bias:
    mbias4, mbias8 = create_bias(cal_list)

def load_bias(amps,mbias4,mbias8):
    if amps=='4':
        mbias = mbias4
    elif amps=='8':
        mbias = mbias8
    if mbias is None:
        if os.path.exists(red_path+'mbias'+amps+'.fits'):
            mbias = [CCDData.read(red_path+'mbias'+amps+'.fits', hdu=x+1, unit=u.electron) for x in range(int(amps))]
        else:
            log.error('No master bias found, creating master bias.')
            mbias4, mbias8 = mcreate_bias(cal_list)
            if amps=='4':
                mbias = mbias4
            elif amps=='8':
                mbias = mbias8
            if mbias is None:
                log.critical('No bias data present for this configuration, skipping.')
    return mbias        

for cal in cal_list:
    if 'BIAS' not in cal:
        process_flat = True
        if args.skip_red == 'True' or args.skip_red == 'yes':
            if os.path.exists(red_path+'mflat_'+cal+'.fits'):
                process_flat = False
            else:
                log.error('No master flat found, creating master flat.')
        if process_flat:
            mbias = load_bias(cal[-1],mbias4,mbias8)
            if mbias is None:
                break
            scale = []
            for i, flat_file in enumerate(cal_list[cal]):
                with fits.open(flat_file) as hdr:
                    header = hdr[0].header
                    data = hdr[int(int(cal[-1])/2)].data
                scale.append(1/np.median(data[1200:2000,800:1600]))
                flat_raw = [CCDData.read(flat_file, hdu=x+1, unit='adu') for x in range(int(cal[-1]))]
                flat_processed = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=gain*u.electron/u.adu, readnoise=rdnoise*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(flat_raw)]
                flat_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
                for x in flat_processed: flat_hdu.append(fits.ImageHDU(x.data,header=x.header))
                flat_hdu.writeto(flat_file.replace(cal_path,red_path),overwrite=True)
                cal_list[cal][i] = flat_file.replace(cal_path,red_path)
            mflat = [ccdproc.combine(cal_list[cal],hdu=x+1,unit=u.electron,method='median',scale=scale,sigma_clip=True) for x in range(int(cal[-1]))]
            mflat_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
            for x in mflat: mflat_hdu.append(fits.ImageHDU(x.data,header=x.header))
            mflat_hdu.writeto(red_path+'mflat_'+cal+'.fits',overwrite=True)
            for flat in cal_list[cal]: os.remove(flat)

for sci in sci_list:
    if args.target is not None:
        if args.target not in sci:
            continue
    process_data = True
    if args.skip_red == 'True' or args.skip_red == 'yes':
        if os.path.exists(red_path+sci+'.fits'):
            with fits.open(red_path+sci+'.fits',mode='update') as hdr:
                wcs_object = wcs.WCS(hdr[0].header)
                sci_med = hdr[0].data
            processed = [x.replace(sci_path,red_path) for x in sci_list[sci]]
            process_data = False
        else:
            log.error('No reduced image found for '+red_path+sci+'.fits, processing data.')
    if process_data:
        fil = sci.split('_')[-1]
        mbias = load_bias(cal[-1],mbias4,mbias8)
        if mbias is None:
            continue
        if os.path.exists(red_path+'mflat_'+cal+'.fits'):
            mflat = [CCDData.read(red_path+'mflat_'+fil+'.fits', hdu=x+1, unit=u.electron) for x in range(int(fil[-1]))]
        else:
            log.critical('No bias data present for this configuration, skipping.')
            continue
        reprojected = []
        for i, sci_file in enumerate(sci_list[sci]):
            with fits.open(sci_file) as hdr:
                header = hdr[0].header
            sci_raw = [CCDData.read(sci_file, hdu=x+1, unit='adu') for x in range(int(fil[-1]))]
            sci_processed = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=gain*u.electron/u.adu, readnoise=rdnoise*u.electron, master_bias=mbias[k], master_flat=mflat[k], gain_corrected=True) for k,x in enumerate(sci_raw)]
            bkg = MeanBackground(SigmaClip(sigma=3.))
            bkg_value = [bkg.calc_background(x) for x in sci_processed]
            sci_final = [x.subtract(bkg_value[k]*u.electron,propagate_uncertainties=True,handle_meta='first_found').divide(header['ELAPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found') for k,x in enumerate(sci_processed)]
            if fil[-1] == '4':
                sci_full = CCDData(np.concatenate(sci_final,axis=1),header=header,unit=u.electron/u.second)
            elif fil[-1] == '8':
                sci_full = CCDData(np.concatenate([sci_final[1],np.fliplr(sci_final[0]),sci_final[3],np.fliplr(sci_final[2]),sci_final[5],np.fliplr(sci_final[4]),sci_final[7],np.fliplr(sci_final[6])],axis=1),header=header,unit=u.electron/u.second)
            sci_full = ccdproc.trim_image(sci_full[700:2400,4100:6100])#[700:2400,300:7300])
            sci_full.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(sci_full.data)[1],np.shape(sci_full.data)[0]))
            reprojected.append(sci_full)
            reprojected[i].write(sci_file.replace(sci_path,red_path),overwrite=True)
        align = align_stars(reprojected)
        com = [CCDData(x,unit=u.electron/u.second) for x in align]
        sci_med = ccdproc.combine(reprojected,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
        sci_med.header = reprojected[0].header
        ra = sci_med.header['RA'].split(':')
        dec = sci_med.header['DEC'].split(':')
        coords = SkyCoord(ra[0]+'h'+ra[1]+'m'+ra[2]+'s',dec[0]+'d'+dec[1]+'m'+dec[2]+'s',frame='icrs')
        sci_med.header['RADECSYS'] = 'ICRS'
        sci_med.header['CUNIT1'] = 'deg'
        sci_med.header['CUNIT2'] = 'deg'
        sci_med.header['CTYPE1'] = 'RA---TAN'
        sci_med.header['CTYPE2'] = 'DEC--TAN'
        sci_med.header['CRPIX1'] = 5000
        sci_med.header['CRPIX2'] = 1000
        sci_med.header['CRVAL1'] = coords.ra.deg
        sci_med.header['CRVAL2'] = coords.dec.deg
        sci_med.header['CD1_2'] = -0.1185/3600*np.cos(np.pi/180.*(sci_med.header['ROTPOSN']+90))
        sci_med.header['CD1_1'] = 0.1185/3600*np.sin(np.pi/180.*(sci_med.header['ROTPOSN']+90))
        sci_med.header['CD2_2'] = -0.1185/3600*np.sin(np.pi/180.*(sci_med.header['ROTPOSN']+90))
        sci_med.header['CD2_1'] = 0.1185/3600*np.cos(np.pi/180.*(sci_med.header['ROTPOSN']+90))
        sci_med.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(sci_med)[1],np.shape(sci_med)[0]))
        sci_med.header['RDNOISE'] = rdnoise/np.sqrt(len(align))
        sci_med.header['NFILES'] = len(align)
            # for k, n in enumerate(sci_list[sci]):
            #     hdr[0].header['FILE'+str(k+1)] = (os.path.basename(n), 'Name of file used in median.')
        sci_med.wcs = wcs.WCS(sci_med.header)
        sci_med.write(red_path+sci+'.fits',overwrite=True)
