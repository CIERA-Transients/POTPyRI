#!/usr/bin/env python

"Automatic pipeline for imaging with GMOS."
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
from astropy.nddata import CCDData
from astropy.modeling import models
from photutils import make_source_mask, MeanBackground, StdBackgroundRMS, CircularAperture, CircularAnnulus, aperture_photometry, Background2D, MedianBackground, DAOStarFinder
from ccdproc import wcs_project
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.utils import calc_total_error
from astroquery.vizier import Vizier
from astroquery.irsa import Irsa
import warnings
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning

params = argparse.ArgumentParser(description='Path of data.')
params.add_argument('data_path', default=None, help='Path of data, in default format from Gemini.') #path to GMOS data: required
params.add_argument('--skip_red', type=str, default=None, help='Option to skip reduction.') #
params.add_argument('--target', type=str, default=None, help='Option to only reduce this target.') #
params.add_argument('--individual', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--wcs', type=str, default=None, help='Option to recenter WCS.') #must have pyraf install and access to IRAF to use
params.add_argument('--phot', type=str, default=None, help='Option to perform photometry.') #must have pyraf install and access to IRAF to use
params.add_argument('--verbose', type=str, default=None, help='Print warnings.') #must have pyraf install and access to IRAF to use
args = params.parse_args()

if not args.verbose:
    warnings.simplefilter('ignore',category=(AstropyWarning, AstropyDeprecationWarning))

def sort_files(files): #sort the calibration files:
    cal_list = {}
    cal_list.update({'BIAS':[]})
    sci_list = {}
    for f in files:
        with fits.open(f) as file_open:
            hdr = file_open[0].header
            if hdr['RAWGEMQA'] == 'USABLE' or hdr['RAWGEMQA'] == 'PASS' or hdr['RAWGEMQA'] == 'UNKNOWN':
                if hdr['GRATING'] == 'MIRROR':
                    if hdr['OBSCLASS'] == 'dayCal':
                        if hdr['OBSTYPE'] == 'BIAS':
                            cal_list['BIAS'].append(f)
                        elif hdr['OBSTYPE'] == 'OBJECT':
                            fil = hdr['FILTER2'].replace('_','')
                            try:
                                cal_list[fil]
                            except KeyError:
                                cal_list.update({fil:[]})
                            cal_list[fil].append(f)
                        else:
                            shutil.move(f,bad_path) #other options, here for now
                    elif hdr['OBSCLASS'] == 'science':
                        target = hdr['OBJECT'].replace(' ','').replace('_','')
                        fil = hdr['FILTER2'].replace('_','')
                        try:
                            sci_list[target+'_'+fil]
                        except KeyError:
                            sci_list.update({target+'_'+fil:[]})
                        sci_list[target+'_'+fil].append(f)
                    else:
                        shutil.move(f,bad_path) #other options, here for now
                else:
                    try:
                        shutil.move(f,sci_path)
                    except:
                        pass
            else:
                shutil.move(f,bad_path)
    return cal_list, sci_list

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

log_file_name = red_path+'GMOS_log_'+datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')+'.log' #create log file name
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

files = glob.glob(args.data_path+'/*.fits.bz2')
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
    files = glob.glob(cal_path+'*.fits.bz2')+glob.glob(sci_path+'*.fits.bz2')
    cal_list, sci_list = sort_files(files)

if len(cal_list['BIAS']) == 0:
    log.critical('No bias data present, exiting.')
    logging.shutdown()
    sys.exit(-1)

process_bias = True
if args.skip_red == 'True' or args.skip_red == 'yes':
    if os.path.exists(red_path+'mbias.fits'):
        mbias = [CCDData.read(red_path+'mbias.fits', hdu=x+1, unit=u.electron) for x in range(12)]
        process_bias = False
    else:
        log.error('No bias master bias found, creatging master bias.')
if process_bias:
    for i, bias_file in enumerate(cal_list['BIAS']):
        with fits.open(bias_file) as hdr:
            header = hdr[0].header
        bias_raw = [CCDData.read(bias_file, hdu=x+1, unit='adu') for x in range(12)]
        bias_processed = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, readnoise=x.header['RDNOISE']*u.electron) for x in bias_raw]
        for x in bias_processed: x.header['DATASEC'] = ('[1:256,1:2112]', 'updated DATASEC')
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in bias_processed: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias_file.replace(cal_path,red_path).replace('.fits.bz2','.fits'),overwrite=True)
        cal_list['BIAS'][i] = bias_file.replace(cal_path,red_path).replace('.fits.bz2','.fits')
    mbias = [ccdproc.combine(cal_list['BIAS'],hdu=x+1,unit=u.electron) for x in range(12)]
    mbias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
    for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
    mbias_hdu.writeto(red_path+'mbias.fits',overwrite=True)
    for bias in cal_list['BIAS']: os.remove(bias)

for cal in cal_list:
    if cal != 'BIAS':
        process_flat = True
        if args.skip_red == 'True' or args.skip_red == 'yes':
            if os.path.exists(red_path+'mflat_'+cal+'.fits'):
                process_flat = False
            else:
                log.error('No master flat found, creatging master flat.')
        if process_flat:
            scale = []
            for i, flat_file in enumerate(cal_list[cal]):
                with fits.open(flat_file) as hdr:
                    header = hdr[0].header
                    hdr3 = hdr[3].header
                    data3 = hdr[3].data
                xbin, ybin = int(hdr3['CCDSUM'][0]), int(hdr3['CCDSUM'][-1])
                scale.append(1/np.median(data3[int(100/ybin):int(4500/ybin),int(100/xbin):int(1800/xbin)]))
                flat_raw = [CCDData.read(flat_file, hdu=x+1, unit='adu') for x in range(12)]
                flat_processed = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, readnoise=x.header['RDNOISE']*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(flat_raw)]
                for x in flat_processed: x.header['DATASEC'] = ('[1:256,1:2112]', 'updated DATASEC')
                flat_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
                for x in flat_processed: flat_hdu.append(fits.ImageHDU(x.data,header=x.header))
                flat_hdu.writeto(flat_file.replace(cal_path,red_path).replace('.fits.bz2','.fits'),overwrite=True)
                cal_list[cal][i] = flat_file.replace(cal_path,red_path).replace('.fits.bz2','.fits')
            mflat = [ccdproc.combine(cal_list[cal],hdu=x+1,unit=u.electron,method='median',scale=scale,sigma_clip=True) for x in range(12)]
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
                sci_med = hdr[0].data
                wcs_object = wcs.WCS(hdr[0].header)
            processed = [x.replace(sci_path,red_path).replace('.fits.bz2','.fits') for x in sci_list[sci]]
            process_data = False
        else:
            log.error('No reduced image found for '+red_path+sci+'.fits, processing data.')
    if process_data:
        fil = sci.split('_')[-1]
        mflat = [CCDData.read(red_path+'mflat_'+fil+'.fits', hdu=x+1, unit=u.electron) for x in range(12)]
        wcs_objects = []
        for i, sci_file in enumerate(sci_list[sci]):
            with fits.open(sci_file) as hdr:
                header = hdr[0].header
                if header['NCCDS'] == 1: #at the moment this assumes that all data from the same target either has 1 or 3 NCCDS
                    exts = 4
                else:
                    exts = 0
                header_wcs = hdr[7-exts].header
            wcs_objects.append(wcs.WCS(header_wcs))
            sci_raw = [CCDData.read(sci_file, hdu=x+1, unit='adu') for x in range(12-2*exts)]
            sci_processed = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, readnoise=x.header['RDNOISE']*u.electron, master_bias=mbias[k+exts], master_flat=mflat[k+exts], gain_corrected=True) for k,x in enumerate(sci_raw)]
            for x in sci_processed: x.header['DATASEC'] = ('[1:256,1:2112]', 'updated DATASEC')
            header.extend(sci_raw[0].wcs.to_header())
            bkg = MeanBackground(SigmaClip(sigma=3.))
            bkg_value = [bkg.calc_background(x) for x in sci_processed]
            sci_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
            for j,x in enumerate(sci_processed): sci_hdu.append(fits.ImageHDU(x.subtract(bkg_value[j]*np.ones(np.shape(x))*u.electron,propagate_uncertainties=True,handle_meta='first_found').divide(header['ELAPSED']*u.second,propagate_uncertainties=True,handle_meta='first_found').data,header=x.header))
            sci_hdu.writeto(sci_file.replace(sci_path,red_path).replace('.fits.bz2','.fits'),overwrite=True)
        processed = [x.replace(sci_path,red_path).replace('.fits.bz2','.fits') for x in sci_list[sci]]
        fringe = []
        for i in range(12-2*exts):
            processed_data = [CCDData.read(x,hdu=i+1,unit=u.electron/u.second) for k,x in enumerate(processed)]
            mask = [make_source_mask(x, nsigma=3, npixels=5) for x in processed_data]
            fringe.append(ccdproc.combine(processed_data,method='median',mask=mask,sigma_clip=True,sigma_clip_func=np.ma.median)) #minmax_clip=True,minmax_clip_max=np.median(sky_level)+4*np.std(sky_level)
        fri_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in fringe: fri_hdu.append(fits.ImageHDU(x.data,header=x.header))
        fri_hdu.writeto(red_path+'fringe_'+sci+'.fits',overwrite=True)
        reprojected = []
        for i,sci_ind in enumerate(processed):
            with fits.open(sci_ind) as hdr:
                header = hdr[0].header
            sci_exts = [CCDData.read(sci_ind, hdu=x+1, unit=u.electron/u.second) for x in range(12-2*exts)]
            sci_final = [x.subtract(fringe[j],propagate_uncertainties=True,handle_meta='first_found') for j,x in enumerate(sci_exts)]
            header += sci_final[0].header
            header.update(wcs_objects[i].to_header())
            if header['INSTRUME'] == 'GMOS-S':
                header['CRPIX1'] = header['CRPIX1']+(np.shape(sci_final[0])[1]*(4-exts))+40-(10*exts)
                if exts == 0:
                    sci_full = CCDData(np.concatenate(sci_final[2:4]+[np.ones([np.shape(sci_final[3])[0],40])*-999]+sci_final[4:8]+[np.ones([np.shape(sci_final[7])[0],40])*-999]+sci_final[8:10],axis=1),header=header,wcs=wcs.WCS(header),unit=u.electron)
                elif exts == 4:
                    sci_full = CCDData(np.concatenate(sci_final,axis=1),header=header,wcs=wcs.WCS(header),unit=u.electron)
            elif header['INSTRUME'] == 'GMOS-N': #check
                header['CRPIX1'] = header['CRPIX1']+(np.shape(sci_final[0])[1]*(4-exts))+40-(10*exts)
                if exts == 0:
                    sci_full = CCDData(np.concatenate(sci_final[2:4]+[np.ones([np.shape(sci_final[3])[0],30])*-999]+sci_final[4:8]+[np.ones([np.shape(sci_final[7])[0],30])*-999]+sci_final[8:10],axis=1),header=header,wcs=wcs.WCS(header),unit=u.electron)
                elif exts == 4:
                    sci_full = CCDData(np.concatenate(sci_final,axis=1),header=header,wcs=wcs.WCS(header),unit=u.electron)
            sci_full.mask = np.zeros(np.shape(sci_full))
            sci_full.mask[sci_full==-999] = 1
            sci_full.mask = sci_full.mask.astype(np.bool)
            sci_full.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(sci_full)[1],np.shape(sci_full)[0]))
            if i == 0: wcs_object = wcs.WCS(header)
            wcs_objects[i] = wcs.WCS(header)
            reprojected.append(wcs_project(sci_full,wcs_objects[0]))
            reprojected[i].write(processed[i],overwrite=True)
        sci_med = ccdproc.combine(reprojected,method='median',sigma_clip=True,sigma_clip_func=np.ma.median)
        sci_med.write(red_path+sci+'.fits',overwrite=True)
        with fits.open(red_path+sci+'.fits',mode='update') as hdr:
            hdr[0].header['RDNOISE'] = header['RDNOISE']/len(sci_list[sci])
            hdr[0].header['NFILES'] = len(sci_list[sci])
            for k, n in enumerate(sci_list[sci]):
                hdr[0].header['FILE'+str(k+1)] = (os.path.basename(n), 'Name of file used in median.')
    if args.wcs == 'True' or args.wcs == 'yes':
        fig, ax = plt.subplots(figsize=(7,7))
        ax = plt.subplot(projection=wcs_object)
        gaia = Irsa.query_region(SkyCoord(hdr[0].header['CRVAL1']*u.deg, hdr[0].header['CRVAL2']*u.deg,frame='fk5'), catalog="gaia_dr2_source", spatial="Cone",radius=3*u.arcmin)
        if len(gaia) == 0:
            log.info('No GAIA stars found within 3 arcmin for starlist.')
            plt.close()
        else:
            ax.imshow(sci_med, cmap='gray', norm=ImageNormalize(sci_med, interval=ZScaleInterval()))
            _, median, std = sigma_clipped_stats(sci_med, sigma=3.0)
            daofind = DAOStarFinder(fwhm=7.0, threshold=5.*std)
            sources = daofind(np.asarray(sci_med))
            for l,m in enumerate(gaia['source_id']):
                x, y = (wcs.WCS(hdr[0].header)).all_world2pix(gaia['ra'][l],gaia['dec'][l],1)
                ax.add_patch(patches.Circle((x,y),radius=3,edgecolor='g',alpha=0.5,facecolor='none',linewidth=2, label='Gaia star: RA = %f, Dec = %f'%(gaia['ra'][l],gaia['dec'][l]), picker=True))
            for i in range(len(sources)):
                ax.add_patch(patches.Circle((sources['xcentroid'][i],sources['ycentroid'][i]),radius=3,edgecolor='b',alpha=0.5,facecolor='none',linewidth=2,label='Source star x = %f, y = %f'%(sources['xcentroid'][i],sources['ycentroid'][i]), picker=True))
            ax.set_title('Target: '+sci)
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
            with fits.open(red_path+sci+'.fits',mode='update') as hdr:
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
                for i,sci_ind in enumerate(processed):
                    with fits.open(sci_ind,mode='update') as hdr_ind:
                        header = hdr_ind[0].header
                        header['CRVAL1'] = gaia_star[0][0]
                        header['CRVAL2'] = gaia_star[0][1]
                        header['CRPIX1'] = source_star[0][0]
                        header['CRPIX2'] = source_star[0][1]
    if args.phot == 'True' or args.phot == 'yes':
        fil = sci.split('_')[-1]
        SDSS_filter = fil[0]
        Vizier.ROW_LIMIT = -1
        zp_cal = True
        try:
            sdss = Vizier.query_region(SkyCoord(hdr[0].header['CRVAL1']*u.deg, hdr[0].header['CRVAL2']*u.deg, frame='fk5'),radius=4*u.arcmin,catalog='V/139')[0]
            sdss = sdss[(sdss['mode']==1)&(sdss['q_mode']=='+')]
            if len(sdss) != 0:
                pass
            else:
                log.info('No SDSS stars available in this field.')
                zp_cal = False
            try:
                check = sdss[SDSS_filter+'mag']
            except KeyError:
                log.info('No SDSS stars available in this filter.')
                zp_cal = False
        except IndexError:
            log.info('No SDSS stars available in this field.')
            zp_cal = False
        record_phot = True
        while record_phot:
            fig, ax = plt.subplots(figsize=(7,7))
            ax = plt.subplot(projection=wcs_object)
            ax.imshow(sci_med, cmap='gray', norm=ImageNormalize(sci_med, interval=ZScaleInterval()))
            if zp_cal:
                for i in range(len(sdss)):
                    x, y = (wcs.WCS(hdr[0].header)).all_world2pix(sdss['RA_ICRS'][i],sdss['DE_ICRS'][i],1)
                    ax.add_patch(patches.Circle((x,y),radius=4,edgecolor='r',alpha=0.5,facecolor='none',linewidth=2,label='SDSS star mag = %f'%(sdss[SDSS_filter+'mag'][i]), picker=True))
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
            ax.set_title('Target: '+sci)
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
            continue
        sci_med = np.asarray(sci_med)
        dum, median, std = sigma_clipped_stats(sci_med, sigma=3.0)
        masked = np.ma.masked_outside(sci_med,median-3*std,median+3*std)
        sigma_clip = SigmaClip(sigma=3)
        bkg_estimator = MedianBackground()
        bkg = Background2D(sci_med, (60, 60), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=masked.mask, exclude_percentile=80)
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
        log.info('Calculating magnitudes for '+sci)
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
                    log.info('Magnitude of target = %.3f +/- %.3f AB mag.'%(mag[i],mag_err[i]))
                else:
                    log.info('Magnitude of standard = %.3f +/- %.3f AB mag.'%(mag[i],mag_err[i]))
        if args.individual == 'True' or args.individual == 'yes':
            log.info('User selected to extract light curve.')
            lightcurve = {}
            lightcurve_err = {}
            mjd = []
            for i,sci_ind in enumerate(processed):
                with fits.open(sci_ind) as hdr:
                    start_time = hdr[0].header['DATE']+hdr[0].header['UTSTART']
                    mid_time = datetime.datetime.strptime(start_time,'%Y-%m-%d%H:%M:%S.%f')+datetime.timedelta(seconds=hdr[0].header['ELAPSED']/2)
                    mjd.append(astropy.time.Time(mid_time).jd)
                    sci_data = hdr[0].data
                bkg = Background2D(sci_data, (133, 132), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
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
                    np.savetxt(red_path+sci+'_lightcurve_target.txt',np.c_[mjd,lightcurve[star],lightcurve_err[star]])
                else:
                    np.savetxt(red_path+sci+'_lightcurve_standard'+str(star)+'.txt',np.c_[mjd,lightcurve[star],lightcurve_err[star]])



