#parameter file for GMOS/Gemini

import os
import astropy
import datetime
import numpy as np
from astropy.coordinates import SkyCoord
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

__version__ = 1.2 #last edited 09/11/2021

def static_mask(proc):
    return ['']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.0803 #0.0807 0.08

def saturation(hdr):
    return 65535

def WCS_keywords_old(): #WCS keywords
    return ['PC1_1','PC1_2','PC2_1','PC2_2',]
    
def WCS_keywords(): #WCS keywords
    return ['CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2']

def cal_path():
    return None

def raw_format(proc):
    return '*.fits.bz2'

def dark():
    return False

def bias():
    return True

def flat():
    return True

def raw_header_ext():
    return 0

def science_keyword():
    return ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']

def science_files():
    return ['OPEN','MIRROR','science','OBJECT']

def flat_keyword():
    return ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']

def flat_files():
    return ['OPEN','MIRROR','dayCal','OBJECT']

def bias_keyword():
    return ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']

def bias_files():
    return ['CLOSED','MIRROR','dayCal','BIAS']

def dark_keyword():
    return []

def dark_files():
    return []

def spec_keyword():
    return ['SHUTTER','FILTER2']

def spec_files():
    return ['OPEN','open2-8']

def bad_keyword():
    return ['RAWGEMQA']

def bad_files():
    return ['FAIL']

def target_keyword():
    return 'OBJECT'

def filter_keyword(hdr):
    return hdr['FILTER2'].replace('_','')

def amp_keyword(hdr):
    nccd = hdr['NCCDS']
    if nccd == '1':
        amp = '4'
    else:
        amp = '12'
    return amp

def bin_keyword(hdr):
    return '22' #hdr['CCDSUM'].replace(' ','') keyword in different ext - check how to calculate

def time_format(hdr):
    return Time(hdr['DATE-OBS']+'T'+hdr['TIME-OBS']).mjd

def wavelength():
    return 'OPT'

def load_bias(red_path,amp,binn):
    bias = red_path+'mbias_'+amp+'_'+binn+'.fits'
    mbias = [CCDData.read(bias,hdu=x+1,unit=u.electron) for x in range(int(amp))]
    return mbias

def create_bias(cal_list,cal,red_path,log):
    amp = cal.split('_')[1]
    binn = cal.split('_')[2]
    log.info('Processing bias files with '+amp+' amps and '+binn+' binning.')
    log.info(str(len(cal_list[cal]))+' files found.')
    for i,bias in enumerate(cal_list[cal]):
        with fits.open(bias) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(bias, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, readnoise=x.header['RDNOISE']*u.electron) for x in raw]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias.replace('/raw/','/red/').replace('.bz2',''),overwrite=True)
        cal_list[cal][i] = bias.replace('/raw/','/red/').replace('.bz2','')    
    mbias = [ccdproc.combine(cal_list[cal],hdu=x+1,unit=u.electron) for x in range(int(amp))]
    mbias_hdu = fits.HDUList([fits.PrimaryHDU()])
    for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
    log.info('Created master bias with '+amp+' amps and '+binn+' binning.')
    mbias_hdu[0].header['VER'] = (__version__, 'Version of telescope parameter file used.')
    mbias_hdu.writeto(red_path+'mbias_'+amp+'_'+binn+'.fits',overwrite=True)
    log.info('Master bias written to mbias_'+amp+'_'+binn+'.fits')
    for bias in cal_list[cal]: os.remove(bias)
    return

def flat_name(flatpath,fil,amp,binn):
    return [flatpath+'mflat_'+fil+'_'+amp+'_'+binn+'.fits']

def load_flat(flat):
    mflat = CCDData.read(flat[0],unit=u.electron)
    return mflat

def create_flat(flat_list,fil,amp,binn,red_path,mbias=None,log=None):
    log.info('Processing files for filter: '+fil)
    log.info(str(len(flat_list))+' files found.')
    scale = []
    flats = []
    for i, flat in enumerate(flat_list):
        log.info('Loading file: '+flat)
        with fits.open(flat) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(flat, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, readnoise=x.header['RDNOISE']*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if header['INSTRUME'] == 'GMOS-S':
            if amp == '12':
                flat_full = CCDData(np.concatenate(red[2:4]+[np.zeros([np.shape(red[3])[0],40])]+red[4:8]+[np.zeros([np.shape(red[7])[0],40])]+red[8:10],axis=1),header=header,unit=u.electron)
            elif amp == '4':
                flat_full = CCDData(np.concatenate(red,axis=1),header=header,unit=u.electron)
        elif header['INSTRUME'] == 'GMOS-N':
            if amp == '12':
                flat_full = CCDData(np.concatenate(red[2:4]+[np.zeros([np.shape(red[3])[0],30])]+red[4:8]+[np.zeros([np.shape(red[7])[0],30])]+red[8:10],axis=1),header=header,unit=u.electron)
            elif amp == '4':
                flat_full = CCDData(np.concatenate(red,axis=1),header=header,unit=u.electron)
        log.info('Exposure time of image is '+str(flat_full.header['EXPTIME']))
        norm = 1/np.median(flat_full[500:1500,700:1400]) #check for binning
        log.info('Flat normalization: '+str(norm))
        scale.append(norm)
        flats.append(flat_full)
    mflat = ccdproc.combine(flats,method='median',scale=scale,sigma_clip=True)
    log.info('Created master flat for filter: '+fil+' and '+amp+' amp '+binn+' biinning.')
    mflat.header['VER'] = (__version__, 'Version of telescope parameter file used.')
    mflat.write(red_path+'mflat_'+fil+'_'+amp+'_'+binn+'.fits',overwrite=True)
    log.info('Master flat written to mflat_'+fil+'_'+amp+'_'+binn+'.fits')
    return

def process_science(sci_list,fil,amp,binn,red_path,mbias=None,mflat=None,proc=None,log=None):
    masks = []
    processed = []
    for sci in sci_list:
        log.info('Loading file: '+sci)
        log.info('Applying bias correction, gain correction and flat correction.')
        with fits.open(sci) as hdr:
            header = hdr[0].header
            if amp == '12':
                header_wcs = hdr[7].header
            elif amp == '4':
                header_wcs = hdr[3].header
        raw = [CCDData.read(sci, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, readnoise=x.header['RDNOISE']*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        header += header_wcs
        if header['INSTRUME'] == 'GMOS-S':
            if amp == '12':
                sci_full = CCDData(np.concatenate(red[2:4]+[np.zeros([np.shape(red[3])[0],40])]+red[4:8]+[np.zeros([np.shape(red[7])[0],40])]+red[8:10],axis=1),header=header,unit=u.electron)
                header['CRPIX1'] = header['CRPIX1']+(np.shape(red[0])[1]*4)+40
            elif amp == '4':
                sci_full = CCDData(np.concatenate(red,axis=1),header=header,unit=u.electron)
        elif header['INSTRUME'] == 'GMOS-N':
            if amp == '12':
                sci_full = CCDData(np.concatenate(red[2:4]+[np.zeros([np.shape(red[3])[0],30])]+red[4:8]+[np.zeros([np.shape(red[7])[0],30])]+red[8:10],axis=1),header=header,unit=u.electron)
                header['CRPIX1'] = header['CRPIX1']+(np.shape(red[0])[1]*4)+40
            elif amp == '4':
                sci_full = CCDData(np.concatenate(red,axis=1),header=header,unit=u.electron)
        log.info('Exposure time of science image is '+str(sci_full.header['EXPTIME']))
        processed_data = ccdproc.flat_correct(sci_full, mflat)
        log.info('File proccessed.')
        log.info('Cleaning cosmic rays and creating mask.')
        mask = make_source_mask(processed_data, nsigma=3, npixels=5)
        masks.append(mask)
        # clean, com_mask = create_mask.create_mask(sci,processed_data,'_mask.fits',static_mask(proc),mask,saturation(red.header),binn(),np.mean(readnoise(amp)),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
        # processed_data.data = clean
        log.info('Calculating 2D background.')
        bkg = Background2D(processed_data, (128, 128), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
        log.info('Median background: '+str(np.median(bkg.background)))
        fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg.fits').replace('.bz2',''),np.array(bkg.background),overwrite=True)
        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found').divide(processed_data.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
        log.info('Background subtracted and image divided by exposure time.')
        final.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(final)[1],np.shape(final)[0]))  
        final.write(sci.replace('/raw/','/red/').replace('.bz2',''),overwrite=True)
        processed.append(final)
    return processed, masks

def stacked_image(tar,red_path):
    return [red_path+tar+'.fits']

def rdnoise(header):
    return header['RDNOISE']

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
    return ['SDSS','PS1']

def exptime(hdr):
    return hdr['EXPTIME']

def fringe_correction(fil):
    if 'z' in fil:
        return True
    else:
        return False
    
def trim(f):
    return False