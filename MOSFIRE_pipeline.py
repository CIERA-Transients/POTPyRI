#!/usr/bin/env python

"Automatic pipeline for imaging with MOSFIRE. Creates median coadded images for each target."
"Individual images should be checked and removed from sci path."
"Author: Kerry Paterson"
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
from astroquery.vizier import Vizier
from astroquery.irsa import Irsa

params = argparse.ArgumentParser(description='Path of data.')
params.add_argument('data_path', default=None, help='Path of data, in default format from KOA.') #path to MOSFIRE data: required
params.add_argument('--use_dome_flats', type=str, default=None, help='Use dome flats for flat reduction.') #use dome flat instead of sci images to create master flat
params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.') #
params.add_argument('--target', type=str, default=None, help='Option to only reduce this target.') #
params.add_argument('--individual', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--wcs', type=str, default=None, help='Option to use IRAF to apply WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--phot', type=str, default=None, help='Option to use IRAF to perform photometry.') #must have pyraf install and access to IRAF to use
args = params.parse_args()

def sort_files(files): #sort the calibration files: 
    cal_list = {}
    sci_list = {}
    sky_list = {}
    time_list = {}
    for f in files:
        with fits.open(f) as file_open:
            hdr = file_open[0].header
            if hdr['MASKNAME']=='OPEN': #only collect imaging files **errors with the mask can result in this being MIRA#only collect dome flats for imaging
                if 'FLAT' in hdr['TARGNAME']:
                    fil = hdr['FILTER'] #sort by filter
                    try:
                        cal_list[fil]
                    except KeyError:
                        cal_list.update({fil:[]})
                    cal_list[fil].append(f)
                else:
                    target = hdr['TARGNAME'].replace(' ','') #sort by target
                    fil = hdr['FILTER'] #and filter
                    try:
                        sci_list[target+'_'+fil]
                    except KeyError:
                        sci_list.update({target+'_'+fil:[]})
                    sci_list[target+'_'+fil].append(f)
                    try:
                        time_list[target+'_'+fil]
                    except KeyError:
                        time_list.update({target+'_'+fil:[]})
                    time_list[target+'_'+fil].append(hdr['MJD-OBS'])
                    try:
                        sky_list[fil]
                    except KeyError:
                        sky_list.update({fil:[]})
                    sky_list[fil].append(f)
            else:
                try:
                    shutil.move(f,sci_path) #move spec data here for now
                except:
                    pass
    return cal_list, sci_list, sky_list, time_list

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

def get_twomass_filter(mosfire_fil):
    if mosfire_fil == 'Ks':
        fil = 'k_m'
    elif mosfire_fil == 'J':
        fil = 'j_m'
    elif mosfire_fil == 'H':
        fil = 'h_m'
    elif mosfire_fil == 'Y':
        fil = 'y_m'
    else:
        fil = None
    return fil

def AB_conversion(fil):
    if fil == 'Ks':
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

readnoise = {1:21,4:10.8,8:7.7,16:5.8,32:4.2,64:3.5,128:3.0}

raw_path = args.data_path+'/raw/' #path containing the raw data
cal_path = raw_path+'cal/' #path containing the calibration files
if not os.path.exists(cal_path): #create reduced file path if it doesn't exist
    os.makedirs(cal_path)
sci_path = raw_path+'sci/' #path containing the scie files
if not os.path.exists(sci_path): #create reduced file path if it doesn't exist
    os.makedirs(sci_path)
red_path = args.data_path+'/red/' #path to write the reduced files
if not os.path.exists(red_path): #create reduced file path if it doesn't exist
    os.makedirs(red_path)

log_file_name = red_path+'MOSFIRE_log_'+datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')+'.log' #create log file name
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

log.info('Reducing data from '+raw_path)
files = glob.glob(args.data_path+'/*.fits.gz')
if len(files) != 0:
    cal_list, sci_list, sky_list, time_list = sort_files(files)
    for cal in cal_list:
        for i,f in enumerate(cal_list[cal]):
            shutil.move(f,cal_path)
            cal_list[cal][i] = f.replace(args.data_path,cal_path)
    for sci in sci_list:
        for i,f in enumerate(sci_list[sci]):
            shutil.move(f,sci_path)
            sci_list[sci][i] = f.replace(args.data_path,sci_path)
    for sky in sky_list:
        for i,f in enumerate(sky_list[sky]):
            sky_list[sky][i] = f.replace(args.data_path,sci_path)
else:
    files = glob.glob(cal_path+'*.fits.gz')+glob.glob(sci_path+'*.fits')
    cal_list, sci_list, sky_list, time_list = sort_files(files)

if len(sci_list) == 0: #check if there are images to reduce, otherwise exit #to do:add cal only
    log.critical('No science images found, check data_path and mask name.')
    logging.shutdown()
    sys.exit(-1)
else:
    log.info(str(len(sci_list))+' science targets found')

filter_list = []
for target in sci_list:
    fil = target.split('_')[1]
    if args.skip_red == 'True' or args.skip_red == 'yes':
        if os.path.exists(red_path+'master_flat_'+fil+'.fits'):
            process_flat = False
        else:
            process_flat = True
    else:
        process_flat = True
    if process_flat:
        if fil not in filter_list:
            if args.use_dome_flats == 'yes' or args.use_dome_flats == 'True': #use dome flats if user used use_dome_flats arguement
                flat_type = 'dome'
                log.info('User set option to use dome flats to create master flat')
                log.info(str(len(cal_list[fil]))+' dome flats found for filter '+fil)
                flat_list = cal_list[fil]
            else: #default to use science files for master flat creation
                flat_type = 'sky'
                log.info('Using science files to create master flat')
                flat_list = sky_list[fil]
            flat_scale = lambda arr: 1.0/np.ma.median(arr[224:1824,224:1824]) #scaling function for normalization
            if len(flat_list) > 31:
                random.seed(4)
                flat_list = random.sample(flat_list,k=31)
            image_list = []
            for image in flat_list:
                sci_raw = CCDData.read(image, unit=u.adu)
                image_list.append(ccdproc.ccd_process(sci_raw, gain=sci_raw.header['SYSGAIN']*u.electron/u.adu, readnoise=readnoise[sci_raw.header['NUMREADS']]*u.electron))
                mask = [make_source_mask(x, nsigma=1, npixels=5) for x in image_list] #####make into list
            weights = np.ones([len(image_list),2048,2048])
            for i,n in enumerate(image_list):
                weights[i]*=np.median(n[224:1824,224:1824]) #flats weighted by median in flat
            master_flat = ccdproc.combine(image_list,method='median',weights=weights,scale=flat_scale,mask=mask,sigma_clip=True,sigma_clip_func=np.ma.median) #combine flats
            flat_hdu = fits.PrimaryHDU(master_flat) #create primary HDU with master flat data
            flat_hdu.header['TYPE'] = (flat_type, 'Type of observations used to create master flat.') #add info to header and log
            log.info('Files (including weight) used for master flat creation:')
            for i,n in enumerate(flat_list):
                flat_hdu.header['FLAT'+str(i+1)] = (os.path.basename(n), 'Name of flat file used in creation of master flat.')
                log.info(os.path.basename(n)+' ('+str(np.median(weights[i]))+')')
            flat_hdu.header['FLATMED'] = (np.nanmedian(master_flat), 'Median level of master flat.')
            flat_hdu.header['FLATSTD'] = (np.nanstd(master_flat), 'Std of master flat.')
            flat_hdu.header['LOG'] = (log_file_name, 'Name of the log file.')
            log.info('Median level of master flat: '+str(np.nanmedian(master_flat)))
            log.info('Std of master flat: '+str(np.nanstd(master_flat)))
            flat_hdu.writeto(red_path+'master_flat_'+fil+'.fits',overwrite=True) #write master flat
            filter_list.append(fil)

for target in sci_list:
    if args.target is not None:
        if args.target not in target:
            continue
        else:
            log.info('User specified a target for reduction: '+args.target)
            log.info('Matching target found: '+target)
    process_data = True
    if args.skip_red == 'True' or args.skip_red == 'yes':
        if os.path.exists(red_path+target+'.fits'):
            with fits.open(red_path+target+'.fits',mode='update') as hdr:
                sci_med = hdr[0].data
                wcs_object = wcs.WCS(hdr[0].header)
            processed = [x.replace(sci_path,red_path).replace('.fits','.fits') for x in sci_list[target]]
            process_data = False
        else:
            log.error('No reduced image found, processing data.')
    if process_data:
        fil = target.split('_')[-1]
        master_flat = red_path+'master_flat_'+fil+'.fits'
        flat_data = CCDData.read(master_flat,unit=u.electron)
        log.info('Subtracting sky from science files and applying flat field')
        new_sci_list = []
        for j,n in enumerate(sci_list[target]):
            with fits.open(n) as hdr:
                header = hdr[0].header
                data = hdr[0].data
            data[np.isnan(data)] = np.nanmedian(data)
            sci_data = CCDData(data,meta=header,unit=u.adu)
            red = ccdproc.ccd_process(sci_data, gain=header['SYSGAIN']*u.electron/u.adu, readnoise=readnoise[header['NUMREADS']]*u.electron, gain_corrected=True)
            flat = np.asarray(ccdproc.flat_correct(red,flat_data,add_keyword={'FLATTYPE': flat_data.header['TYPE'],
                                        'FLATNAME': red_path+'master_flat_'+flat_data.header['TYPE']+'_'+fil+'.fits',
                                        'LOG': log_file_name})) #flat field
            flat[np.isnan(flat)] = np.nanmedian(flat)
            bkg = MeanBackground(SigmaClip(sigma=3.))
            bkg_value = bkg.calc_background(flat)
            red = CCDData(flat,unit=u.electron,header=red.header,wcs=red.wcs).subtract(bkg_value*np.ones(np.shape(flat))*u.electron,propagate_uncertainties=True,handle_meta='first_found')
            red = red.divide(red.header['TRUITIME']*red.header['COADDONE']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            red.write(n.replace(sci_path,red_path).replace('.fits','_red.fits'),overwrite=True)
            new_sci_list.append(n.replace(sci_path,red_path).replace('.fits','_red.fits'))
        for j,n in enumerate(new_sci_list):
            time_diff = sorted([(abs(time_list[target][j]-n2),k) for k,n2 in enumerate(time_list[target])])
            sky_list = [new_sci_list[k] for _,k in time_diff[1:10]]
            sky_data = []
            sky_hdu = fits.PrimaryHDU()
            for k,m in enumerate(sky_list):
                sky_hdu.header['FILE'+str(k+1)] = (os.path.basename(m), 'Name of file used in creation of sky.')
                sky_data.append(CCDData.read(m,unit=u.electron/u.second))
            mask = [make_source_mask(x, nsigma=1, npixels=5) for x in sky_data]
            sky_hdu.data = ccdproc.combine(sky_data,method='median',mask=mask,sigma_clip=True,sigma_clip_func=np.ma.median)
            sky_hdu.writeto(n.replace('_red.fits','_sky.fits'),overwrite=True)
        reprojected = []
        for j,n in enumerate(new_sci_list):
            sci_data = CCDData.read(n,unit=u.electron/u.second)
            sky = CCDData.read(n.replace('_red.fits','_sky.fits'),unit=u.electron/u.second)
            sci_final = sci_data.subtract(sky,propagate_uncertainties=True,handle_meta='first_found')
            if j == 0: wcs_object = sci_data.wcs
            sci_final.wcs = sci_data.wcs
            reprojected.append(wcs_project(sci_final,wcs_object))
            reprojected[j].write(n,overwrite=True)
        sci_med = ccdproc.combine(reprojected,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
        sci_med.write(red_path+target+'.fits',overwrite=True)
        with fits.open(red_path+target+'.fits',mode='update') as hdr:
            hdr[0].header['RDNOISE'] = readnoise[hdr[0].header['NUMREADS']]/len(new_sci_list)
            hdr[0].header['NFILES'] = len(new_sci_list)
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
                for sci in new_sci_list:
                    with fits.open(red_path+os.path.basename(sci),mode='update') as hdr_ind:
                        header = hdr_ind[0].header
                        header['CRVAL1'] = gaia_star[0][0]
                        header['CRVAL2'] = gaia_star[0][1]
                        header['CRPIX1'] = source_star[0][0]
                        header['CRPIX2'] = source_star[0][1]
    if args.phot == 'True' or args.phot == 'yes':
        zp_cal = True
        twomass_filter = get_twomass_filter(fil)
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
                    x, y = (wcs.WCS(hdr[0].header)).all_world2pix(twomass['ra'][i],twomass['dec'][i],1)
                    ax.add_patch(patches.Circle((x,y),radius=4,edgecolor='r',alpha=0.5,facecolor='none',linewidth=2,label='2MASS star mag = %f'%(twomass[twomass_filter][i]), picker=True))
                std_mag = []
                fig.canvas.mpl_connect('pick_event', lambda event: onpick(event,positions,[],std_mag))
            else:
                fig.canvas.mpl_connect('pick_event', lambda event: onpick(event,positions,[],[]))
            if args.wcs == 'True' or args.wcs == 'yes':
                pass
            else:
                _, median, std = sigma_clipped_stats(sci_med, sigma=3.0)  
                daofind = DAOStarFinder(fwhm=7.0, threshold=3.*std)  
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
                        record_phot = False
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
            for sci in new_sci_list:
                with fits.open(red_path+os.path.basename(sci)) as hdr:
                    start_time = hdr[0].header['DATE-OBS']
                    mid_time = datetime.datetime.strptime(start_time,'%Y-%m-%dT%H:%M:%S')+datetime.timedelta(seconds=hdr[0].header['DARKTIME']/2) 
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