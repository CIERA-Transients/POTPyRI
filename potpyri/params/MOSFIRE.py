#parameter file for MOSFIRE/Keck
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

__version__ = 1.6 #last edited 09/11/2021

def static_mask(proc):
    return ['./staticmasks/MF.staticmask.fits']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.1799

def saturation(hdr):
    return hdr['SATURATE']*hdr['SYSGAIN']

def WCS_keywords_old(): #WCS keywords
    return ['WAT0_001', 'WAT1_001', 'WAT2_001']
    
def WCS_keywords(): #WCS keywords
    return ['CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2']

def cal_path():
    return None

def raw_format(proc):
    return 'MF.*.fits.gz'

def dark():
    return False

def bias():
    return False

def flat():
    return True

def raw_header_ext():
    return 0

def science_keyword():
    return ['GRATMODE','MASKNAME']

def science_files():
    return ['imaging','OPEN']

def flat_keyword():
    return ['OBJECT']

def flat_files(): #error
    return ['Flat']

def bias_keyword():
    return []

def bias_files():
    return []

def dark_keyword():
    return ['FILTER']

def dark_files():
    return ['Dark']

def spec_keyword():
    return ['GRATMODE']

def spec_files():
    return ['spectroscopy']

def bad_keyword():
    return ['MASKNAME','PONAME']

def bad_files():
    return ['CLOSED','MIRA']

def target_keyword():
    return 'TARGNAME'

def filter_keyword(hdr):
    filt = hdr['FILTER'].replace(' ','').split('_')[0]
    if filt == 'Ks':
        filt = 'K'
    return filt

def amp_keyword(hdr):
    return '1'

def bin_keyword(hdr):
    return '11'

def time_format(hdr):
    return hdr['MJD-OBS']

def wavelength():
    return 'NIR'

def flat_name(flatpath,fil,amp,binn):
    return [flatpath+'mflat_'+fil+'_'+amp+'_'+binn+'.fits']

def load_flat(flat):
    mflat = CCDData.read(flat[0],unit=u.electron/u.second)
    return mflat

def create_flat(flat_list,fil,amp,binn,red_path,mbias=None,log=None):
    log.info('Processing files for filter: '+fil)
    log.info(str(len(flat_list))+' files found.')
    flats = []
    masks  = []
    flat_scale = []
    for flat in flat_list:
        log.info('Loading file: '+flat)
        raw = CCDData.read(flat,unit=u.adu)
        red = ccdproc.ccd_process(raw, gain=raw.header['SYSGAIN']*u.electron/u.adu, readnoise=rdnoise(raw.header)*u.electron)
        log.info('Exposure time of image is '+str(red.header['TRUITIME']*red.header['COADDONE']))
        with fits.open(static_mask(False)[0]) as hdr:
            mask_bp = -~-hdr[0].data
        mask = make_source_mask(red, nsigma=3, npixels=5)+mask_bp
        bkg = Background2D(red, (20, 20), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
        masked = np.array(red)
        masked[mask] = bkg.background[mask]
        log.info('Median flat level: '+str(np.median(masked)))
        if np.median(masked)==0:
            log.info('Removing file from flat creation.')
        else:
            norm = 1/np.median(masked[224:1824,224:1824])
            log.info('Flat normalization: '+str(norm))
            flat_scale.append(norm)
            flats.append(CCDData(masked,meta=red.header,unit=u.electron))
    mflat = ccdproc.combine(flats,method='median',scale=flat_scale,sigma_clip=True) 
    log.info('Created master flat for filter '+fil+' and '+amp+' amps and '+binn+' binning.')
    mflat.header['VER'] = (__version__, 'Version of telescope parameter file used.')
    mflat.write(red_path+'mflat_'+fil+'_'+amp+'_'+binn+'.fits',overwrite=True)
    log.info('Master flat written to mflat_'+fil+'_'+amp+'_'+binn+'.fits')
    return

def process_science(sci_list,fil,amp,binn,red_path,mbias=None,mflat=None,proc=None,log=None):
    masks = []
    processed = []
    for sci in sci_list:
        log.info('Loading file: '+sci)
        log.info('Applying gain correction and flat correction.')
        with fits.open(sci) as hdr:
            header = hdr[0].header
            data = hdr[0].data
        data[np.isnan(data)] = np.nanmedian(data)
        raw = CCDData(data,meta=header,unit=u.adu)
        red = ccdproc.ccd_process(raw, gain=raw.header['SYSGAIN']*u.electron/u.adu, readnoise=rdnoise(raw.header)*u.electron)
        log.info('Exposure time of science image is '+str(red.header['TRUITIME']*red.header['COADDONE']))
        flat = np.array(ccdproc.flat_correct(red, mflat))
        flat[np.isnan(flat)] = np.nanmedian(flat)
        processed_data = CCDData(flat,unit=u.electron,header=red.header,wcs=red.wcs)
        log.info('File proccessed.')
        log.info('Cleaning cosmic rays and creating mask.')
        with fits.open(static_mask(False)[0]) as hdr:
            mask_bp = -~-hdr[0].data
        mask = make_source_mask(processed_data, nsigma=5, npixels=5)+mask_bp
        masks.append(mask)
        # clean, com_mask = create_mask.create_mask(sci.replace('.gz',''),processed_data,'_mask.fits',static_mask(proc)[0],mask,saturation(red.header),binning(),rdnoise(raw.header),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
        # processed_data.data = clean
        log.info('Calculating 2D background.')
        bkg = Background2D(processed_data, (128, 128), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
        log.info('Median background: '+str(np.median(bkg.background)))
        fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg.fits').replace('.gz',''),np.array(bkg.background),overwrite=True)
        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found').divide(red.header['TRUITIME']*red.header['COADDONE']*u.second,propagate_uncertainties=True,handle_meta='first_found')
        log.info('Background subtracted and image divided by exposure time.')
        final.write(sci.replace('/raw/','/red/').replace('.gz',''),overwrite=True)
        processed.append(final)
    return processed, masks

def stacked_image(tar,red_path):
    return [red_path+tar+'.fits']

def rdnoise(header):
    readnoise = {1:21,4:10.8,8:7.7,16:5.8,32:4.2,64:3.5,128:3.0}
    return readnoise[header['NUMREADS']]

def binning():
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
    return ['2MASS']

def exptime(hdr):
    return hdr['TRUITIME']*hdr['COADDONE']

def trim(f):
    return False