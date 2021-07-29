#parameter file for LRIS/Keck
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

__version__ = 1.0 #last edited 28/07/2021

def static_mask(proc):
    return ['']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.129 #0.135, 0.123 #new one?

def saturation(hdr):
    return 65535

def WCS_keywords_old(): #WCS keywords
    return []
    
def WCS_keywords(): #WCS keywords
    return ['CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2']

def cal_path():
    return None

def raw_format(proc):
    return 'L*.fits.gz'

def dark():
    return False

def bias():
    return True

def flat():
    return True

def raw_header_ext():
    return 0

def science_keyword():
    return ['KOAIMTYP','SLITNAME','GRANAME','TRAPDOOR']

def science_files():
    return ['object','direct','mirror','open']

def flat_keyword():
    return ['KOAIMTYP']

def flat_files():
    return ['flatlamp']

def bias_keyword():
    return ['KOAIMTYP']

def bias_files():
    return ['bias']

def dark_keyword():
    return []

def dark_files():
    return []

def spec_keyword():
    return ['GRISTRAN']

def spec_files():
    return ['deployed']

def bad_keyword():
    return ['SLITNAME','KOAIMTYP']

def bad_files():
    return ['GOH_LRIS','focus']

def target_keyword():
    return 'TARGNAME'

def filter_keyword(hdr):
    instrument = hdr['INSTRUME']
    if instrument == 'LRISBLUE':
        filt = hdr['BLUFILT']
    if instrument == 'LRIS':
        filt = hdr['REDFILT']
    return filt

def amp_keyword(hdr):
    try:
        amp = str(hdr['NUMAMPS'])
    except:
        amp = str(hdr['NVIDINP']) #check for new - return 1
    instrument = hdr['INSTRUME']
    if instrument == 'LRISBLUE':
        side = 'B'
    if instrument == 'LRIS':
        side = 'R'
    return amp+side

def bin_keyword(hdr):
    return hdr['BINNING'].replace(',','')

def time_format(hdr):
    try:
        mid_time = float(hdr['MJD-OBS'])
    except:
        mid_time = hdr['MJD']
    return mid_time

def wavelength():
    return 'OPT'

def load_bias(red_path,amp,binn):
    bias = red_path+'mbias_'+amp+'_'+binn+'.fits'
    mbias = [CCDData.read(bias,hdu=x+1,unit=u.electron) for x in range(int(amp[0]))]
    return mbias

def create_bias(cal_list,cal,red_path,log):
    amp = cal.split('_')[1]
    binn = cal.split('_')[2]
    gains = gain(amp)
    readnoises = readnoise(amp)
    oscan_reg = overscan_region(amp)
    log.info('Processing bias files with '+amp+' amps and '+binn+' binning.')
    log.info(str(len(cal_list[cal]))+' files found.')
    for i,bias in enumerate(cal_list[cal]):
        with fits.open(bias) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(bias, hdu=x+1, unit='adu') for x in range(int(amp[0]))]
        red = [ccdproc.ccd_process(x, oscan=oscan_reg, oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=gains[j]*u.electron/u.adu, readnoise=readnoises[j]*u.electron) for j,x in enumerate(raw)]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias.replace('/raw/','/red/').replace('.gz',''),overwrite=True)
        cal_list[cal][i] = bias.replace('/raw/','/red/').replace('.gz','')    
    mbias = [ccdproc.combine(cal_list[cal],hdu=x+1,unit=u.electron) for x in range(int(amp[0]))]
    mbias_hdu = fits.HDUList([fits.PrimaryHDU()])
    for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
    log.info('Created master bias with '+amp+' amps and '+binn+' binning.')
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
    gains = gain(amp)
    readnoises = readnoise(amp)
    oscan_reg = overscan_region(amp)
    for i, flat in enumerate(flat_list):
        log.info('Loading file: '+flat)
        with fits.open(flat) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(flat, hdu=x+1, unit='adu') for x in range(int(amp[0]))]
        red = [ccdproc.ccd_process(x, oscan=oscan_reg, oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if amp == '4B':
            flat_full = CCDData(np.concatenate([red[0],np.fliplr(red[1]),np.zeros([np.shape(red[1])[0],111]),red[2],np.fliplr(red[3])],axis=1),header=header,unit=u.electron)
            flat_full = ccdproc.trim_image(flat_full[700:3315,350:3940])
        if amp == '4R':
            flat_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),np.zeros([np.shape(red[0])[0],200]),red[3],np.fliplr(red[2])],axis=1),header=header,unit=u.electron)
        if amp == '1R':
            flat_full = red[0]
        log.info('Exposure time of image is '+str(flat_full.header['ELAPTIME']))
        norm = 1/np.median(flat_full[1200:1600,1200:1600]) #check for binning
        log.info('Flat normalization: '+str(norm))
        scale.append(norm)
        flats.append(flat_full)
    mflat = ccdproc.combine(flats,method='median',scale=scale,sigma_clip=True)
    log.info('Created master flat for filter: '+fil+' and '+amp+' amp '+binn+' biinning.')
    mflat.write(red_path+'mflat_'+fil+'_'+amp+'_'+binn+'.fits',overwrite=True)
    log.info('Master flat written to mflat_'+fil+'_'+amp+'_'+binn+'.fits')
    return

def process_science(sci_list,fil,amp,binn,red_path,mbias=None,mflat=None,proc=None,log=None):
    masks = []
    processed = []
    gains = gain(amp)
    readnoises = readnoise(amp)
    oscan_reg = overscan_region(amp)
    for sci in sci_list:
        log.info('Loading file: '+sci)
        log.info('Applying bias correction, gain correction and flat correction.')
        with fits.open(sci) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(sci, hdu=x+1, unit='adu') for x in range(int(amp[0]))]
        red = [ccdproc.ccd_process(x, oscan=oscan_reg, oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if amp == '4B':
            sci_full = CCDData(np.concatenate([red[0],np.fliplr(red[1]),np.zeros([np.shape(red[1])[0],111]),red[2],np.fliplr(red[3])],axis=1),header=header,unit=u.electron)
            sci_full = ccdproc.trim_image(sci_full[700:3315,350:3940])
        if amp == '4R':
            sci_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),np.zeros([np.shape(red[0])[0],200]),red[3],np.fliplr(red[2])],axis=1),header=header,unit=u.electron)
        if amp == '1R':
            sci_full = red[0]
        log.info('Exposure time of science image is '+str(sci_full.header['ELAPTIME']))
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
        fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg.fits').replace('.gz',''),np.array(bkg.background),overwrite=True)
        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found').divide(processed_data.header['ELAPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
        log.info('Background subtracted and image divided by exposure time.')
        log.info('Writing WCS to file.')
        ra = final.header['RA'].split(':')
        dec = final.header['DEC'].split(':')
        coords = SkyCoord(ra[0]+'h'+ra[1]+'m'+ra[2]+'s',dec[0]+'d'+dec[1]+'m'+dec[2]+'s',frame='icrs')
        final.header['RADECSYS'] = 'ICRS'
        final.header['CUNIT1'] = 'deg'
        final.header['CUNIT2'] = 'deg'
        final.header['CTYPE1'] = 'RA---TAN'
        final.header['CTYPE2'] = 'DEC--TAN'
        final.header['CRVAL1'] = coords.ra.deg
        final.header['CRVAL2'] = coords.dec.deg
        if amp == '4B':
            final.header['CRPIX1'] = 2456
            final.header['CRPIX2'] = 1184
            final.header['CD1_1'] = 0.135/3600*np.sin(np.pi/180.*(final.header['ROTPOSN']+90))
            final.header['CD1_2'] = -0.135/3600*np.cos(np.pi/180.*(final.header['ROTPOSN']+90))
            final.header['CD2_1'] = 0.135/3600*np.cos(np.pi/180.*(final.header['ROTPOSN']+90))
            final.header['CD2_2'] = -0.135/3600*np.sin(np.pi/180.*(final.header['ROTPOSN']+90))
        if amp == '4R':
            final.header['CRPIX1'] = 2400
            final.header['CRPIX2'] = 1200
            final.header['CD1_1'] = 0.123/3600*np.sin(np.pi/180.*(final.header['ROTPOSN']+90))
            final.header['CD1_2'] = -0.123/3600*np.cos(np.pi/180.*(final.header['ROTPOSN']+90))
            final.header['CD2_1'] = 0.123/3600*np.cos(np.pi/180.*(final.header['ROTPOSN']+90))
            final.header['CD2_2'] = -0.123/3600*np.sin(np.pi/180.*(final.header['ROTPOSN']+90))            
        final.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(final)[1],np.shape(final)[0]))  
        final.write(sci.replace('/raw/','/red/').replace('.gz',''),overwrite=True)
        processed.append(final)
    return processed, masks

def stacked_image(tar,red_path):
    return [red_path+tar+'.fits']

def rdnoise(header):
    try:
        amp = str(header['NUMAMPS'])
    except:
        amp = str(header['NVIDINP']) #check for new - return 1
    instrument = header['INSTRUME']
    if instrument == 'LRISBLUE':
        side = 'B'
    if instrument == 'LRIS':
        side = 'R'
    amp = amp+side
    if amp=='4B':
        readnoise = 3.825
    if amp=='4R':
        readnoise = 3.565
    if amp=='1R':
        readnoise = 4
    return readnoise

def binning():
    return [4,4]

def cr_clean_sigclip():
    return 50

def cr_clean_sigcfrac():
    return 0.1

def cr_clean_objlim():
    return 100

def run_phot():
    return False

def catalog_zp():
    return ['SDSS','PS1']

def exptime(hdr):
    return hdr['ELAPTIME']

def gain(amp):
    if  amp=='4B':
        gains = [1.55,1.56,1.63,1.70]
    if amp=='4R':
        gains = [1.71,1.64,1.61,1.67]
    if amp=='1R':
        gains = [1]
    return gains

def readnoise(amp):
    if amp=='4B':
        readnoise = [3.9,4.2,3.6,3.6]
    if amp=='4R':
        readnoise = [3.64,3.45,3.65,3.52]
    if amp=='1R':
        readnoise = [4]
    return readnoise

def overscan_region(amp):
    if amp=='4B':
        oscan_reg = '[1:50,1:4096]'
    if amp=='4R':
        oscan_reg = '[1:7,1:2520]'
    if amp=='1R':
        oscan_reg = '[2065:2170,1:4188]'
    return oscan_reg

def fringe_correction(fil):
    return False