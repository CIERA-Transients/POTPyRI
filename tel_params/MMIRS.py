#parameter file for MMIRS/MMT
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
import ccdproc
from astropy.modeling import models
import create_mask

__version__ = 1.8 #last edited 09/11/2021

def static_mask(proc):
    return ['./staticmasks/MMIRS.staticmask.fits']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.202

def saturation(hdr):
    return hdr['DATAMAX']*hdr['GAIN']

def WCS_keywords_old(): #WCS keywords
    return ['PC1_1','PC1_2','PC2_1','PC2_2','PV2_1','PV2_2','PV2_3','PV2_4','PV2_5']
    
def WCS_keywords(): #WCS keywords
    return ['CRPIX1','CRPIX2','PC1_1','PC1_2','PC2_1','PC2_2']

def cal_path():
    return None

def raw_format(proc):
    return '*.fits'

def dark():
    return True

def bias():
    return False

def flat():
    return True

def raw_header_ext():
    return 1

def science_keyword():
    return ['OBSMODE','APTYPE']

def science_files():
    return ['imaging','open']

def flat_keyword():
    return []

def flat_files():
    return []

def bias_keyword():
    return []

def bias_files():
    return []

def dark_keyword():
    return ['OBJECT']

def dark_files():
    return ['ark']

def spec_keyword():
    return ['OBSMODE','APERTURE']

def spec_files():
    return ['spectral','open']

def bad_keyword():
    return ['MOSID']

def bad_files():
    return ['closed']

def target_keyword():
    return 'OBJECT'

def filter_keyword(hdr):
    return hdr['FILTER'].replace(' ','').split('_')[0]

def amp_keyword(hdr):
    return str(hdr['NAMPS'])

def bin_keyword(hdr):
    return '11'

def time_format(hdr):
    return Time(hdr['DATE-OBS']).mjd

def wavelength():
    return 'NIR'

def flat_name(flatpath,fil,amp,binn):
    return [flatpath+'mflat_'+fil+'_'+amp+'_'+binn+'.fits']

def load_flat(flat):
    mflat = CCDData.read(flat[0],unit=u.electron)
    return mflat

def create_flat(flat_list,fil,amp,binn,red_path,mbias=None,log=None):
    log.info('Processing files for filter: '+fil)
    log.info(str(len(flat_list))+' files found.')
    flats = []
    flat_scale = []
    for flat in flat_list:
        log.info('Loading file: '+flat)
        raw = CCDData.read(flat,hdu=1,unit=u.adu)
        red = ccdproc.ccd_process(raw, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron)
        log.info('Exposure time of image is '+str(red.header['EXPTIME']))
        log.info('Loading correct master dark')
        mdark = CCDData.read(red_path+'DARK_'+str(red.header['EXPTIME'])+'_'+amp+'_'+binn+'.fits', unit=u.electron)        
        red = ccdproc.subtract_dark(red, mdark, exposure_time='EXPTIME', exposure_unit=u.second)
        red = ccdproc.subtract_overscan(red, overscan=red[:,0:4], overscan_axis=1, model=models.Chebyshev1D(3))
        mask = make_source_mask(red, nsigma=3, npixels=5)
        bkg = Background2D(red, (20, 20), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
        masked = np.array(red)
        masked[mask] = bkg.background[mask]
        log.info('Median flat level: '+str(np.median(masked)))
        norm = 1/np.median(masked[500:1500,500:1500])
        log.info('Flat normalization: '+str(norm))
        flat_scale.append(norm)
        flats.append(CCDData(masked,meta=red.header,unit=u.electron))
    mflat = ccdproc.combine(flats,method='median',scale=flat_scale,sigma_clip=True)
    log.info('Created master flat for filter: '+fil)
    mflat.header['VERSION'] = (__version__, 'Version of setting file.')
    mflat.write(red_path+'mflat_'+fil+'_'+amp+'_'+binn+'.fits',overwrite=True)
    log.info('Master flat written to mflat_'+fil+'_'+amp+'_'+binn+'.fits')
    return

def process_science(sci_list,fil,amp,binn,red_path,mbias=None,mflat=None,proc=None,log=None):
    masks = []
    processed = []
    for sci in sci_list:
        log.info('Loading file: '+sci)
        log.info('Applying gain correction, dark subtraction, overscan correction, flat correction, and trimming image.')
        raw = CCDData.read(sci,hdu=1,unit=u.adu)
        red = ccdproc.ccd_process(raw, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron)
        log.info('Exposure time of science image is '+str(red.header['EXPTIME']))
        log.info('Loading correct master dark')
        mdark = CCDData.read(red_path+'DARK_'+str(red.header['EXPTIME'])+'_'+amp+'_'+binn+'.fits', unit=u.electron)
        red = ccdproc.subtract_dark(red, mdark, exposure_time='EXPTIME', exposure_unit=u.second)#dark_exposure=1*u.second, data_exposure=red.header['EXPTIME']*u.second, scale=True)
        red = ccdproc.subtract_overscan(red, overscan=red[:,0:4], overscan_axis=1, model=models.Chebyshev1D(3))
        red = ccdproc.flat_correct(red, mflat)
        processed_data = ccdproc.ccd_process(red, trim=raw.header['DATASEC'])
        log.info('File proccessed and trimmed.')
        log.info('Cleaning cosmic rays and creating mask.')
        mask = make_source_mask(processed_data, nsigma=3, npixels=5)
        masks.append(mask)
        clean, com_mask = create_mask.create_mask(sci,processed_data,'_mask.fits',static_mask(proc)[0],mask,saturation(red.header),binning(),rdnoise(red.header),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
        processed_data.data = clean
        log.info('Calculating 2D background.')
        bkg = Background2D(processed_data, (510, 510), filter_size=(9, 9),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
        log.info('Median background: '+str(np.median(bkg.background)))
        fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg.fits'),np.array(bkg.background),overwrite=True)
        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found').divide(red.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
        log.info('Background subtracted and image divided by exposure time.')
        final.write(sci.replace('/raw/','/red/'),overwrite=True)
        processed.append(final)
    return processed, masks

def stacked_image(tar,red_path):
    return [red_path+tar+'.fits']

def rdnoise(header):
    return header['RDNOISE']

def binning():
    return [4,4]

def create_dark(cal_list,cal,mbias,red_path,log):
    images = []
    processed = []
    for dark in cal_list[cal]:
        log.info('Creating master dark with exposure time: '+cal.split('_')[1])
        log.info('Loading file: '+dark)
        try:
            raw = CCDData.read(dark, hdu=1, unit='adu')
            images.append(dark)
        except astropy.io.registry.IORegistryError:
            log.error('File not recognized.')
        red = ccdproc.ccd_process(raw, gain=raw.header['GAIN']*u.electron/u.adu, readnoise=raw.header['RDNOISE']*u.electron, master_bias=mbias)
        processed.append(red)
    log.info('Creating master dark.')
    mdark = ccdproc.combine(processed,method='median')
    log.info('Master dark created.')
    mdark_hdu = fits.PrimaryHDU(mdark)
    mdark_hdu.header['VER'] = (__version__, 'Version of telescope parameter file used.')
    for i,im in enumerate(images):
        mdark_hdu.header['FILE'+str(i+1)] = im
    mdark_hdu.header['EXPTIME'] = (np.float(cal.split('_')[1]), 'Exposure time of master dark (sec).')
    mdark_hdu.writeto(red_path+cal+'.fits',overwrite=True)
    log.info('Master dark written to '+cal+'.fits')
    return

def cr_clean_sigclip():
    return 50

def cr_clean_sigcfrac():
    return 0.1

def cr_clean_objlim():
    return 100

def run_phot():
    return True

def catalog_zp():
    return ['2MASS']

def exptime(hdr):
    return hdr['EXPTIME']

def trim(f):
    return False
