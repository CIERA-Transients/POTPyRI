#parameter file for BINOSPEC/MMT
import os
import astropy
import datetime
import numpy as np
from photutils import make_source_mask, Background2D, MeanBackground
from astropy.stats import SigmaClip
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData
import astropy.units.astrophys as u
import astropy.units as u
import astropy.wcs as wcs
import ccdproc
from astropy.modeling import models
import create_mask
from utilities import util

__version__ = 1.4 #last edited 24/08/2021

def static_mask(proc):
    if proc:
        return ['','']#'./staticmasks/bino_proc_left.trim.staticmask.fits','./staticmasks/bino_proc_right.trim.staticmask.fits']
    else:
        return ['','']#'./staticmasks/bino_left.staticmask.fits','./staticmasks/bino_right.staticmask.fits']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.24

def saturation(hdr):
    return 65000 #defualt hdr['DATAMAX']*hdr['GAIN']

def WCS_keywords_old(): #WCS keywords
    return ['PC1_1','PC1_2','PC2_1','PC2_2','WCSNAMEA','CUNIT1A','CUNIT2A','CTYPE1A','CTYPE2A','CRPIX1A','CRPIX2A','CRVAL1A','CRVAL2A','CD1_1A','CD1_2A','CD2_1A','CD2_2A']
    
def WCS_keywords(): #WCS keywords
    return ['CRPIX1','CRPIX2','PC1_1','PC1_2','PC2_1','PC2_2']

def cal_path():
    return str(os.getenv("PIPELINE_HOME"))+'/Imaging_pipelines/BINOSPEC_calib/'

def raw_format(proc):
    if proc:
        return 'sci_img_*proc.fits'
    else:
        return 'sci_img*[!proc].fits'

def dark():
    return False

def bias():
    return False

def flat():
    return False

def raw_header_ext():
    return 1

def science_keyword():
    return ['MASK','SCRN']

def science_files():
    return ['imaging','stowed']

def flat_keyword():
    return ['MASK','SCRN']

def flat_files():
    return ['imaging','deployed']

def bias_keyword():
    return []

def bias_files():
    return []

def dark_keyword():
    return []

def dark_files():
    return []

def spec_keyword():
    return ['MASK']

def spec_files():
    return ['spectroscopy']

def bad_keyword():
    return ['MASK']

def bad_files():
    return ['mira']

def target_keyword():
    return 'OBJECT'

def filter_keyword(hdr):
    return hdr['FILTER'].replace(' ','').split('_')[0]

def amp_keyword(hdr):
    return '1'

def bin_keyword(hdr):
    return hdr['CCDSUM'].replace(' ','')

def time_format(hdr):
    return hdr['MJD']

def wavelength():
    return 'OPT'

def flat_name(cpath,fil,amp,binn):
    return [cpath+'/mflat_'+fil+'_left.fits',cpath+'/mflat_'+fil+'_right.fits']

def load_flat(flat):
    mflat = []
    for f in flat:
        mflat.append(CCDData.read(f,hdu=1,unit=u.electron))
    return mflat

def gain():
    return [1.085, 1.04649118, 1.04159151, 0.97505369, 1.028, 1.16341855, 1.04742053, 1.0447564]

def process_science(sci_list,fil,amp,binn,red_path,mbias=None,mflat=None,proc=None,log=None):
    masks = []
    processed = []
    flat_left = mflat[0]
    flat_right = mflat[1]
    left_list = []
    right_list = []
    left_mask = []
    right_mask = []
    if proc:  
        for j,sci in enumerate(sci_list):
            log.info('Loading file: '+sci)
            log.info('Applying flat correction and trimming left image.')
            left = CCDData.read(sci, hdu=1, unit=u.electron)
            left = ccdproc.flat_correct(left,flat_left)
            left = ccdproc.ccd_process(left, trim=left.header['DATASEC'])
            log.info('Left image proccessed and trimmed.')
            log.info('Cleaning cosmic rays and creating mask.')
            mask = make_source_mask(left, nsigma=3, npixels=5)
            left_mask.append(mask)
            # clean, com_mask = create_mask.create_mask(sci,left,'_mask_left.fits',static_mask(proc)[0],mask,saturation(left.header),binning(proc,'left'),rdnoise(left.header),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
            # left.data = clean
            log.info('Calculating 2D background.')
            bkg = Background2D(left, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            log.info('Median background for left iamge: '+str(np.median(bkg.background)))
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_left.fits'),np.array(bkg.background),overwrite=True)
            left = left.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            log.info('Exposure time of left image is '+str(left.header['EXPTIME']))
            left = left.divide(left.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            log.info('Background subtracted and image divided by exposure time.')
            left.header['DATASEC'] = '[1:'+str(np.shape(left)[1])+',1:'+str(np.shape(left)[0])+']'
            left_list.append(left)
            log.info('Applying flat correction and trimming right image.')
            right = CCDData.read(sci, hdu=2, unit=u.electron)
            right = ccdproc.flat_correct(right,flat_right)
            right = ccdproc.ccd_process(right, trim=right.header['DATASEC'])
            log.info('Right image proccessed and trimmed.')
            log.info('Cleaning cosmic rays and creating mask.')
            mask = make_source_mask(right, nsigma=3, npixels=5)
            right_mask.append(mask)
            # clean, com_mask = create_mask.create_mask(sci,right,'_mask_right.fits',static_mask(proc)[1],mask,saturation(right.header),binning(proc,'right'),rdnoise(right.header),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
            # right.data = clean
            log.info('Calculating 2D background.')
            bkg = Background2D(right, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            log.info('Median background for right image : '+str(np.median(bkg.background)))
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_right.fits'),np.array(bkg.background),overwrite=True)
            right = right.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            log.info('Exposure time of right image is '+str(right.header['EXPTIME']))
            right = right.divide(right.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            log.info('Background subtracted and image divided by exposure time.')
            right.header['DATASEC'] = '[1:'+str(np.shape(right)[1])+',1:'+str(np.shape(right)[0])+']'
            right_list.append(right)
    else:
        for j,sci in enumerate(sci_list):
            log.info('Loading file: '+sci)
            log.info('Applying gain correction, overscan correction, flat correction, and trimming image.')
            with fits.open(sci) as hdr:
                header_left = hdr[1].header
                header_right = hdr[6].header
            data_list = []
            for i in range(8):
                data = ccdproc.CCDData.read(sci,hdu=i+1,unit=u.adu)
                red = ccdproc.ccd_process(data, oscan=data[:,0:50], oscan_model=models.Chebyshev1D(3), trim='[1200:2098,210:2056]', gain=gain()[i]*u.electron/u.adu, readnoise=4*u.electron)
                data_list.append(np.asarray(red).astype(np.float32))
            top_left = np.concatenate([data_list[0],np.fliplr(data_list[1])],axis=1)
            bot_left = np.flipud(np.concatenate([data_list[3],np.fliplr(data_list[2])],axis=1))
            left = CCDData(np.concatenate([top_left,bot_left]),unit=u.electron,header=header_left)
            left = ccdproc.flat_correct(left,flat_left[209:3903,1149:2947])
            log.info('Left image proccessed and trimmed.')
            log.info('Cleaning cosmic rays and creating mask.')
            mask = make_source_mask(left, nsigma=3, npixels=5)
            left_mask.append(mask)
            # clean, com_mask = create_mask.create_mask(sci,left,static_mask(proc)[0],mask,saturation(left.header),binning(proc,'left'),rdnoise(left.header),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
            # processed_data.data = clean
            log.info('Calculating 2D background.')
            bkg = Background2D(left, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            log.info('Median background for left image : '+str(np.median(bkg.background)))
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_left.fits'),bkg.background,overwrite=True)
            left = left.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            log.info('Exposure time of left image is '+str(left.header['EXPTIME']))
            left = left.divide(left.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            log.info('Background subtracted and image divided by exposure time.')
            left.header['DATASEC'] = '[1:1798,1:3694]'
            left.header['RADECSYS'] = 'ICRS'
            left.header['CUNIT1'] = 'deg'
            left.header['CUNIT2'] = 'deg'
            left.header['CTYPE1'] = 'RA---TAN'
            left.header['CTYPE2'] = 'DEC--TAN'
            left.header['CRPIX1'] = 2301
            left.header['CRPIX2'] = 1846
            coord = util.parse_coord(left.header['RA'],left.header['DEC'])
            left.header['CRVAL1'] = coord.ra.deg
            left.header['CRVAL2'] = coord.dec.deg
            left.header['PC1_1'] = -pixscale()/3600*np.sin(np.pi/180.*(left.header['POSANG']+90))
            left.header['PC1_2'] = pixscale()/3600*np.cos(np.pi/180.*(left.header['POSANG']+90))
            left.header['PC2_1'] = -pixscale()/3600*np.cos(np.pi/180.*(left.header['POSANG']+90))
            left.header['PC2_2'] = pixscale()/3600*np.sin(np.pi/180.*(left.header['POSANG']+90))
            left.write(sci.replace('/raw/','/red/').replace('.fits','_left.fits'),overwrite=True)
            left_list.append(left)
            top_right = np.concatenate([data_list[6],np.fliplr(data_list[7])],axis=1)
            bot_right = np.flipud(np.concatenate([data_list[5],np.fliplr(data_list[4])],axis=1))
            right = CCDData(np.concatenate([top_right,bot_right]),unit=u.electron,header=header_right)
            right = ccdproc.flat_correct(right,flat_right[209:3903,1149:2947])
            log.info('Right image proccessed and trimmed.')
            log.info('Cleaning cosmic rays and creating mask.')
            mask = make_source_mask(right, nsigma=3, npixels=5)
            right_mask.append(mask)
            # clean, com_mask = create_mask.create_mask(sci,right,static_mask(proc)[1],mask,saturation(right.header),binning(proc,'right'),rdnoise(right.header),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
            # processed_data.data = clean
            log.info('Calculating 2D background.')
            bkg = Background2D(right, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            log.info('Median background for right image : '+str(np.median(bkg.background)))
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_right.fits'),bkg.background,overwrite=True)
            right = right.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            log.info('Exposure time of right image is '+str(right.header['EXPTIME']))
            right = right.divide(right.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            log.info('Background subtracted and image divided by exposure time.')
            right.header['DATASEC'] = '[1:1798,1:3694]'
            right.header['RADECSYS'] = 'ICRS'
            right.header['CUNIT1'] = 'deg'
            right.header['CUNIT2'] = 'deg'
            right.header['CTYPE1'] = 'RA---TAN'
            right.header['CTYPE2'] = 'DEC--TAN'
            right.header['CRPIX1'] = -504
            right.header['CRPIX2'] = 1845
            coord = util.parse_coord(right.header['RA'],right.header['DEC'])
            right.header['CRVAL1'] = coord.ra.deg
            right.header['CRVAL2'] = coord.dec.deg
            right.header['PC1_1'] = -pixscale()/3600*np.sin(np.pi/180.*(right.header['POSANG']+90))
            right.header['PC1_2'] = pixscale()/3600*np.cos(np.pi/180.*(right.header['POSANG']+90))
            right.header['PC2_1'] = -pixscale()/3600*np.cos(np.pi/180.*(right.header['POSANG']+90))
            right.header['PC2_2'] = pixscale()/3600*np.sin(np.pi/180.*(right.header['POSANG']+90))
            right.write(sci.replace('/raw/','/red/').replace('.fits','_right.fits'),overwrite=True)
            right_list.append(right)
    return [left_list,right_list], [left_mask,right_mask]

def stacked_image(tar,red_path):
    return [red_path+tar+'_left.fits',red_path+tar+'_right.fits']

def suffix():
    return ['_red_left.fits','_red_right.fits']

def rdnoise(header):
    return 4.0

def binning(proc,side):
    if proc:
        if side=='left':
            return [4,5]
        elif side=='right':
            return [4,7]
    else:
        if side=='left':
            return [4,4]
        elif side=='right':
            return [4,4]

def cr_clean_sigclip():
    return 50

def cr_clean_sigcfrac():
    return 0.1

def cr_clean_objlim():
    return 100

def run_phot():
    return True

def catalog_zp():
    return ['SDSS','PS1']

def exptime(hdr):
    return hdr['EXPTIME']

def fringe_correction(fil):
    if fil == 'z':
        return True
    else:
        return False

def trim(f):
    return False