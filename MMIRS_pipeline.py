#!/usr/bin/env python

"Python pipeline to reduce MMIRS (imaging) data."
"Assumes single object per folder."
"Pipeline uses the IDL pipeline by Igor Chilingarian as a reference"
"(Chilingarian I., et al., 2015, PASP, 127, 406 -"
"https://bitbucket.org/chil_sai/mmirs-pipeline/src/master/)."
"Default is to use calibration files from the MMRIS IDL pipeline."
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
import random
from astropy.nddata import CCDData
from astropy.modeling import models
from photutils import make_source_mask, MeanBackground, StdBackgroundRMS, CircularAperture, CircularAnnulus, aperture_photometry, Background2D, MedianBackground, DAOStarFinder
from ccdproc import wcs_project
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.utils import calc_total_error
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier

params = argparse.ArgumentParser(description='Path of data.')
params.add_argument('data_path', default=None, help='Path of data, in default format from Gemini.') #path to GMOS data: required
params.add_argument('--cal_path', type=str, default=None, help='Use dome flats for flat reduction.') #use dome flat instead of sci images to create master flat
params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.') #
params.add_argument('--individual', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--wcs', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--phot', type=str, default=None, help='Option to perform photometry.') #must have pyraf install and access to IRAF to use
args = params.parse_args()

def sort_files(files):
    dark_list = []
    sci_list = []
    time_list = []
    for f in files:
        if 'dark' in f:
            dark_list.append(f)
        else:
            try:
                with fits.open(f) as hdr:
                    header = hdr[1].header
                if header['OBSMODE'] == 'imaging' and header['APTYPE'] == 'open':
                    sci_list.append(f)
                    time_list.append(datetime.datetime.strptime(header['DATE-OBS'],'%Y-%m-%dT%H:%M:%S'))
                else:
                    shutil.move(f,raw_path+'bad/')
            except AttributeError:
                shutil.move(f,raw_path+'bad/')
    return dark_list, sci_list, time_list

def onclick(event,source_list):
    if event.dblclick:
        source_list.append((event.xdata, event.ydata))
        log.info('You selected x, y position: %d, %d' %(event.xdata, event.ydata)) 

def onpick(event,source_list,gaia_list,std_list):
    label = event.artist.get_label()
    log.info(event.artist.get_label())
    if 'Source' in label:
        source_list.append([np.float(label.split('x = ')[1].split(',')[0]),np.float(label.split('y = ')[1])])
    elif 'Gaia' in label:
        gaia_list.append([np.float(label.split('RA = ')[1].split(',')[0]),np.float(label.split('Dec = ')[1])])
    else:
        std_list.append(np.float(label.split('mag = ')[1]))

def AB_conversion(fil):
    if fil == 'K':
        cor = 1.85
    elif fil == 'J':
        cor = 0.91
    elif fil == 'H':
        cor = 1.39
    elif fil == 'Y':
        cor = 0.634
    else:
        cor = None
    return cor

raw_path = args.data_path+'/raw/' #path containing the raw data
if not os.path.exists(raw_path): #create reduced file path if it doesn't exist
    os.makedirs(raw_path)
if not os.path.exists(raw_path+'bad/'): #create reduced file path if it doesn't exist
    os.makedirs(raw_path+'bad/')
red_path = args.data_path+'/red/' #path to write the reduced files
if not os.path.exists(red_path): #create reduced file path if it doesn't exist
    os.makedirs(red_path)

if args.cal_path is not None:
    cal_path = args.cal_path
else:
    cal_path = os.getenv("HOME")+'/Pipelines/MMIRS_calib/'

log_file_name = red_path+'MMIRS_log_'+datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')+'.log' #create log file name
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

files = glob.glob(args.data_path+'*.fits')
if len(files) != 0:
    dark_list, sci_list, time_list = sort_files(files)
    for i,f in enumerate(dark_list):
        shutil.move(f,raw_path)
        dark_list[i] = f.replace(args.data_path,raw_path)
    for i,f in enumerate(sci_list):
        shutil.move(f,raw_path)
        sci_list[i] = f.replace(args.data_path,raw_path)
else:
    files = glob.glob(raw_path+'*.fits')
    if len(files) != 0:
        dark_list, sci_list, time_list = sort_files(files)
    else:
        log.critical('No files found, please check data path and rerun.')
        logging.shutdown()
        sys.exit(-1)

if len(dark_list) != 0:
    process_dark = True
    if args.skip_red == 'True' or args.skip_red == 'yes':
        if os.path.exists(red_path+'mdark.fits'):
            mdark = CCDData.read(red_path+'mdark.fits', unit=u.electron)
            process_bias = False
        else:
            log.error('No master dark found, creating master dark.')
    if process_dark:
        processed = []
        for dark in dark_list:
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

if len(sci_list) != 0:
    with fits.open(sci_list[0]) as hdr:
        header = hdr[1].header
    target = header['OBJECT']
    fil = header['FILTER']
    process_data = True
    if args.skip_red == 'True' or args.skip_red == 'yes':
        if os.path.exists(red_path+target+'.fits'):
            with fits.open(red_path+target+'.fits',mode='update') as hdr:
                sci_med = hdr[0].data
                wcs_object = wcs.WCS(hdr[0].header)
            process_data = False
        else:
            log.error('No reduced image found, processing data.')
    if process_data:
        wcs_object = wcs.WCS(header)
        if os.path.exists(cal_path+'/pixflat_'+fil+'.fits.gz'):
            flat = cal_path+'/pixflat_'+fil+'.fits.gz'
        else:
            flat = input('Could not find flat '+cal_path+detector+'/pixflat_'+fil+'.fits.gz, please enter full path of file or hit enter to exit. ')
            if not flat:
                log.critical('No master flat present, check data before rerunning.')
                logging.shutdown()
                sys.exit(-1)
        with fits.open(flat) as hdr:
            mflat = hdr[1].data
            mflat[np.isnan(mflat)] = np.nanmedian(mflat)
            mflat = CCDData(mflat,unit=u.electron)
        masks = []
        processed = []
        reprojected = []
        for sci in sci_list:
            raw = CCDData.read(sci,hdu=1,unit=u.adu)
            red = ccdproc.subtract_overscan(raw, overscan=raw[0:4,:], overscan_axis=0, model=models.Chebyshev1D(3))
            red = ccdproc.ccd_process(red, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron)
            red = ccdproc.subtract_dark(red, mdark, exposure_time='EXPTIME', exposure_unit=u.second)
            red = ccdproc.flat_correct(red, mflat)
            processed_data = ccdproc.ccd_process(red, trim=raw.header['DATASEC'])
            sigma_clip = SigmaClip(sigma=3)
            bkg_estimator = MedianBackground()
            _, median, std = sigma_clipped_stats(processed_data, sigma=3.0) 
            m = np.ma.masked_greater(processed_data,median+3*std)
            masks.append(m.mask)
            bkg = Background2D(processed_data, (30, 30), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=m.mask, exclude_percentile=80)
            processed.append(processed_data.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found').divide(red.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found'))
        for i,n in enumerate(processed):
            time_diff = sorted([(abs((time_list[i]-n2).total_seconds()),k) for k,n2 in enumerate(time_list)])
            sky_data = [processed[k] for _,k in time_diff[1:10]]
            sky_mask = [masks[k] for _,k in time_diff[1:10]]
            sky = ccdproc.combine(sky_data,method='median',sigma_clip=True,sigma_clip_func=np.ma.median,mask=sky_mask)
            sky.write(red_path+os.path.basename(sci_list[i]).replace('.fits','_sky.fits'),overwrite=True)
            final = n.subtract(sky,propagate_uncertainties=True,handle_meta='first_found')
            reprojected.append(wcs_project(final,wcs_object))
            reprojected[i].write(red_path+os.path.basename(sci_list[i]),overwrite=True)
        sci_med = ccdproc.combine(reprojected,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
        header = sci_med.header.extend(sci_med.wcs.to_header())
        fits.writeto(red_path+target+'.fits',sci_med.data,sci_med.header,overwrite=True)
        with fits.open(red_path+target+'.fits',mode='update') as hdr:
            hdr[0].header['RDNOISE'] = sci_med.header['RDNOISE']/len(sci_list)
            hdr[0].header['NFILES'] = len(sci_list)
    if args.wcs == 'True' or args.wcs == 'yes':
        fig, ax = plt.subplots(figsize=(7,7))
        ax = plt.subplot(projection=wcs_object)  
        gaia = Irsa.query_region(SkyCoord(hdr[0].header['CRVAL1']*u.deg, hdr[0].header['CRVAL2']*u.deg,frame='fk5'), catalog="gaia_dr2_source", spatial="Cone",radius=4*u.arcmin)
        if len(gaia) == 0:
            log.info('No GAIA stars found within 4 arcmin for starlist.')
            plt.close()
        else:
            ax.imshow(sci_med, cmap='gray', norm=ImageNormalize(sci_med, interval=ZScaleInterval()))
            _, median, std = sigma_clipped_stats(sci_med, sigma=3.0)  
            daofind = DAOStarFinder(fwhm=7.0, threshold=5.*std)  
            sources = daofind(np.asarray(sci_med)) 
            for l,m in enumerate(gaia['source_id']):
                x, y = (wcs.WCS(hdr[0].header)).all_world2pix(gaia['ra'][l],gaia['dec'][l],1)
                ax.add_patch(patches.Circle((x,y),radius=4,edgecolor='g',alpha=0.5,facecolor='none',linewidth=2, label='Gaia star: RA = %f, Dec = %f'%(gaia['ra'][l],gaia['dec'][l]), picker=True))
            for i in range(len(sources)):
                ax.add_patch(patches.Circle((sources['xcentroid'][i],sources['ycentroid'][i]),radius=3,edgecolor='b',alpha=0.5,facecolor='none',linewidth=2,label='Source star x = %f, y = %f'%(sources['xcentroid'][i],sources['ycentroid'][i]), picker=True))
            ax.set_title('Target: '+target)
            source_star = []
            gaia_star = []
            fig.canvas.mpl_connect('pick_event', lambda event: onpick(event,source_star,gaia_star,[]))
            print('Displaying interactive plot to select star to center WCS solution.')
            print('Select 1 star near the center/near your target.')
            print('Default WCS is generally accurate enough but may be off-set.')
            print('First select star (blue) and corresponding Gaia match (green).')
            print('A message will confirm the selection of the star.')
            print('Note: star selection is turned off in zoom/pan mode.')
            print('Close figure when finished.')
            plt.show()
            with fits.open(red_path+target+'.fits',mode='update') as hdr:
                hdr[0].header['CRVAL1'] = gaia_star[0][0]
                hdr[0].header['CRVAL2'] = gaia_star[0][1]
                hdr[0].header['CRPIX1'] = source_star[0][0]
                hdr[0].header['CRPIX2'] = source_star[0][1]
            diff = []
            gaia_cat = SkyCoord(gaia['ra'], gaia['dec'],frame='fk5')
            source_ra, source_dec = (wcs.WCS(hdr[0].header)).all_pix2world(sources['xcentroid'][(sources['xcentroid']>500)&(sources['xcentroid']<1500)&(sources['ycentroid']>500)&(sources['ycentroid']<1500)],sources['ycentroid'][(sources['xcentroid']>500)&(sources['xcentroid']<1500)&(sources['ycentroid']>500)&(sources['ycentroid']<1500)],1)
            source_cat = SkyCoord(source_ra*u.deg,source_dec*u.deg,frame='fk5')
            idx, d2, d3 = gaia_cat.match_to_catalog_sky(source_cat)
            for i,n in enumerate(idx):
                if d2[i].deg*3600 < 1.5:
                    diff.append(d2[i].value*3600)
            log.info('Error on astrometry in central [500:1500,500:1500] of the image is %.3f arcsec.'%np.median(diff))
            if args.individual == 'True' or args.individual == 'yes':
                for sci in sci_list:
                    with fits.open(red_path+os.path.basename(sci),mode='update') as hdr_ind:
                        header = hdr_ind[0].header
                        header['CRVAL1'] = gaia_star[0][0]
                        header['CRVAL2'] = gaia_star[0][1]
                        header['CRPIX1'] = source_star[0][0]
                        header['CRPIX2'] = source_star[0][1]
    if args.phot == 'True' or args.phot == 'yes':
        zp_cal = True
        if fil == 'Y':
            twomass_filter = 'Ymag'
            twomass = Vizier.query_region(SkyCoord(hdr[0].header['CRVAL1']*u.deg, hdr[0].header['CRVAL2']*u.deg,frame='fk5'), catalog="II/319",radius=4*u.arcmin)[0]
        else:
            twomass_filter = fil.lower()+'_m'
            twomass = Irsa.query_region(SkyCoord(hdr[0].header['CRVAL1']*u.deg, hdr[0].header['CRVAL2']*u.deg,frame='fk5'), catalog="fp_psc", spatial="Cone",radius=4*u.arcmin)
        if len(twomass) != 0:
            pass
        else:
            log.info('No 2MASS stars available in this field.')
            zp_cal = False
        try:
            check = twomass[twomass_filter]
        except KeyError:
            log.info('No 2MASS stars available in this filter.')
            zp_cal = False
        record_phot = True
        while record_phot:
            fig, ax = plt.subplots(figsize=(7,7))
            ax = plt.subplot(projection=wcs_object)  
            ax.imshow(sci_med, cmap='gray', norm=ImageNormalize(sci_med, interval=ZScaleInterval()))
            if zp_cal:
                for i in range(len(twomass)):
                    if twomass_filter == 'Ymag':
                        x, y = (wcs.WCS(header)).all_world2pix(twomass['RAJ2000'][i],twomass['DEJ2000'][i],1)
                    else:
                        x, y = (wcs.WCS(header)).all_world2pix(twomass['ra'][i],twomass['dec'][i],1)
                    ax.add_patch(patches.Circle((x,y),radius=4,edgecolor='r',alpha=0.5,facecolor='none',linewidth=2,label='2MASS star mag = %f'%(twomass[twomass_filter][i]), picker=True))
                std_mag = []
                fig.canvas.mpl_connect('pick_event', lambda event: onpick(event,positions,[],std_mag))
            else:
                fig.canvas.mpl_connect('pick_event', lambda event: onpick(event,positions,[],[]))
            if args.wcs == 'True' or args.wcs == 'yes':
                pass
            else:
                _, median, std = sigma_clipped_stats(sci_med, sigma=3.0)  
                daofind = DAOStarFinder(fwhm=7.0, threshold=5.*std)  
                sources = daofind(np.asarray(sci_med)) 
            for i in range(len(sources)):
                ax.add_patch(patches.Circle((sources['xcentroid'][i],sources['ycentroid'][i]),radius=3,edgecolor='b',facecolor='none',linewidth=2,label='Source star x = %f, y = %f'%(sources['xcentroid'][i],sources['ycentroid'][i]), picker=True))
            ax.set_title('Target: '+target)
            positions = []
            fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event,positions))
            print('Displaying interactive plot to select star on which you wish to perform photometry on.')
            print('Select target first, followed by SDSS stars for zero point calculation.')
            print('Either select the source detection or double click on a star to record the x, y position of star.')
            print('Message: "You selected x, y position." will confirm selection of star.')
            print('Select SDSS object to record magnitude after recording the x, y position.')
            print('Note: star selection is turned off in zoom/pan mode.')
            print('Close figure when finished.')
            plt.show()
            if zp_cal:
                if len(positions)-len(std_mag) != 1:
                    record_phot = input('Error occurred while recording positions and standard magnitudes. Redo ? (type True or False) ')
                    if record_phot == 'False':
                        skip_phot = True
                else:
                    record_phot = False
                    skip_phot = False
            else:
                record_phot = False
                skip_phot = False
        if skip_phot:
            logging.shutdown()
            sys.exit(-1)
        sigma_clip = SigmaClip(sigma=3)
        bkg_estimator = MedianBackground()
        bkg = Background2D(sci_med, (120, 120), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        error = calc_total_error(sci_med, bkg.background_rms, 1)
        fwhm = np.float(input('FWHM of stars? '))
        log.info('Using FWHM of '+str(fwhm))
        apertures = CircularAperture(positions, r=fwhm*2.5)
        annulus_apertures = CircularAnnulus(positions, r_in=fwhm*3, r_out=fwhm*4)
        annulus_masks = annulus_apertures.to_mask(method='center')
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(sci_med)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        phot = aperture_photometry(sci_med, apertures, error=error)
        bkg_aper = bkg_median*apertures.area
        mag_inst = -2.5*np.log10((np.asarray(phot['aperture_sum'])-bkg_aper))
        mag_inst_err = 2.5/np.log(10.)*np.asarray(phot['aperture_sum_err'])/np.asarray(phot['aperture_sum'])
        log.info('Calculating magnitudes for '+target)
        for i,m in enumerate(mag_inst):
            if i==0:
                log.info('Instrumental magnitude of target = %.3f +/- %.3f.'%(mag_inst[i],mag_inst_err[i]))
            else:
                log.info('Instrumental magnitude of standard = %.3f +/- %.3f.'%(mag_inst[i],mag_inst_err[i]))
        if zp_cal:
            log.info('Calculating zero point from %i stars with sigma clipping and a maximum error of 0.1.'%len(std_mag))
            for i in np.flip(np.linspace(1,3,num=10)): 
                _, zp, zp_err = sigma_clipped_stats(std_mag-mag_inst[1:len(mag_inst)],sigma=i) 
                log.info('zp, zp_err, sigma used for clipping: %.3f, %.3f, %.3f'%(zp, zp_err, i))
                if zp_err < 0.1: 
                    break 
                else: 
                    pass
            log.info('zp = %.3f +/- %.3f.'%(zp,zp_err))
            mag = mag_inst+zp
            mag_err = np.sqrt(mag_inst_err**2+zp_err**2)
            for i,m in enumerate(mag):
                if i==0:
                    log.info('Magnitude of target = %.3f +/- %.3f AB mag.'%(mag[i]+AB_conversion(fil),mag_err[i]))
                else:
                    log.info('Magnitude of standard = %.3f +/- %.3f AB mag.'%(mag[i]+AB_conversion(fil),mag_err[i]))
        if args.individual == 'True' or args.individual == 'yes':
            log.info('User selected to extract light curve.')
            lightcurve = {}
            lightcurve_err = {}
            mjd = []
            for sci in sci_list:
                with fits.open(red_path+os.path.basename(sci)) as hdr:
                    start_time = hdr[0].header['DATE-OBS']
                    mid_time = datetime.datetime.strptime(start_time,'%Y-%m-%dT%H:%M:%S')+datetime.timedelta(seconds=hdr[0].header['EXPTIME']/2) 
                    mjd.append(astropy.time.Time(mid_time).jd)
                    sci_data = hdr[0].data
                bkg = Background2D(sci_data, (120, 120), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
                error = calc_total_error(sci_data, bkg.background_rms, 1)
                apertures = CircularAperture(positions, r=fwhm*2.5)
                annulus_apertures = CircularAnnulus(positions, r_in=fwhm*3, r_out=fwhm*4)
                annulus_masks = annulus_apertures.to_mask(method='center')
                bkg_median = []
                for mask in annulus_masks:
                    annulus_data = mask.multiply(sci_data)
                    annulus_data_1d = annulus_data[mask.data > 0]
                    _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)
                    bkg_median.append(median_sigclip)
                bkg_median = np.array(bkg_median)
                phot = aperture_photometry(sci_data, apertures, error=error)
                bkg_aper = bkg_median*apertures.area
                mag = -2.5*np.log10((np.asarray(phot['aperture_sum'])-bkg_aper))
                mag_err = 2.5/np.log(10.)*np.asarray(phot['aperture_sum_err'])/np.asarray(phot['aperture_sum'])
                if zp_cal:
                    mag+=zp
                    mag_err = np.sqrt(mag_err**2+zp_err**2)
                for j,m in enumerate(mag):
                    if i==0:
                        lightcurve.update({j:[]})
                        lightcurve_err.update({j:[]})
                    lightcurve[j].append(m)
                    lightcurve_err[j].append(mag_err[j])
            for j in range(len(lightcurve[0])):
                std_stars = []
                for k in range(len(lightcurve)-2):
                    std_stars.append(lightcurve[k+1][j])
                cor = np.nanmedian(std_stars)
                for k in range(len(lightcurve)):
                    lightcurve[k][j] = lightcurve[k][j]-cor+np.median(std_mag)
            for star in lightcurve:
                if star==0:
                    np.savetxt(red_path+target+'_lightcurve_target.txt',np.c_[mjd,lightcurve[star],lightcurve_err[star]])
                else:
                    np.savetxt(red_path+target+'_lightcurve_standard'+str(star)+'.txt',np.c_[mjd,lightcurve[star],lightcurve_err[star]])

