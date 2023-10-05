#parameter file for DEIMOS/Keck
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
import glob

__version__ = 1.5 #last edited 23/11/2021

def static_mask(proc):
    return ['']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.1185

def saturation(hdr):
    return 65000

def WCS_keywords_old(): #WCS keywords
    return []
    
def WCS_keywords(): #WCS keywords
    return ['CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2']

def cal_path():
    return None

def raw_format(proc):
    if proc and str(proc)=='raw':
        return 'd*.fits'
    else:
        return 'DE*.fits.gz'

def dark():
    return False

def bias():
    return True

def flat():
    return True

def raw_header_ext():
    return 0

def science_keyword():
    return ['MOSMODE','OBSMODE','SLMSKNAM','HATCHPOS']

def science_files():
    return ['Direct','imaging','None','open']

def flat_keyword():
    return ['OBSTYPE','OBJECT']

def flat_files():
    return ['DmFlat','lat']

def bias_keyword():
    return ['OBSTYPE']

def bias_files():
    return ['Bias']

def dark_keyword():
    return []

def dark_files():
    return []

def spec_keyword():
    return ['MOSMODE','OBSMODE','SLMSKNAM']

def spec_files():
    return ['Spectral','longslit','LVMslitC']

def bad_keyword():
    return ['SLMSKNAM','PONAME']

def bad_files():
    return ['GOH_X','Mira']

def target_keyword():
    return 'OBJECT'

def filter_keyword(hdr):
    filt = hdr['DWFILNAM']
    return filt

def amp_keyword(hdr):
    return str(hdr['NVIDINP'])

def bin_keyword(hdr):
    return hdr['BINNING'].replace(',','')

def time_format(hdr):
    return float(hdr['MJD-OBS'])

def wavelength():
    return 'OPT'

def load_bias(red_path,amp,binn):
    bias = red_path+'mbias_'+amp+'_'+binn+'.fits'
    mbias = [CCDData.read(bias,hdu=x+1,unit=u.electron) for x in range(int(amp))]
    return mbias

def create_bias(cal_list,cal,red_path,log):
    amp = cal.split('_')[1]
    binn = cal.split('_')[2]
    gains = gain(amp)
    readnoises = readnoise(amp)
    log.info('Processing bias files with '+amp+' amps and '+binn+' binning.')
    log.info(str(len(cal_list[cal]))+' files found.')
    for i,bias in enumerate(cal_list[cal]):
        with fits.open(bias) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(bias, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim='[13:2060,1:2601]', gain=gains[j]*u.electron/u.adu, readnoise=readnoises[j]*u.electron) for j,x in enumerate(raw)]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias.replace('/raw/','/red/').replace('.gz',''),overwrite=True)
        cal_list[cal][i] = bias.replace('/raw/','/red/').replace('.gz','')    
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
    gains = gain(amp)
    readnoises = readnoise(amp)
    for i, flat in enumerate(flat_list):
        log.info('Loading file: '+flat)
        with fits.open(flat) as hdr:
            header = hdr[0].header
        if header['DETSEC01'] == '[1:2048,1:4096]':
            trim_sec = '[13:2060,1321:3921]'
        else:
            trim_sec = '[13:2060,1:2601]'
        raw = [CCDData.read(flat, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim=trim_sec, gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if amp == '4':
            flat_full = CCDData(np.concatenate([red[0],np.zeros([2601,113]),red[1],np.zeros([2601,81]),red[2],np.zeros([2601,113]),red[3]],axis=1),header=header,unit=u.electron/u.second)
        if amp == '8':
            flat_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),red[3],np.fliplr(red[2]),red[5],np.fliplr(red[4]),red[7],np.fliplr(red[6])],axis=1),header=header,unit=u.electron/u.second)            
        log.info('Exposure time of image is '+str(flat_full.header['ELAPTIME']))
        norm = 1/np.median(flat_full[1200:1800,2500:3500])
        log.info('Flat normalization: '+str(norm))
        scale.append(norm)
        flats.append(flat_full)
    mflat = ccdproc.combine(flats,method='median',scale=scale,sigma_clip=True)
    log.info('Created master flat for filter: '+fil)
    mflat.header['VER'] = (__version__, 'Version of telescope parameter file used.')
    mflat.write(red_path+'mflat_'+fil+'_'+amp+'_'+binn+'.fits',overwrite=True)
    log.info('Master flat written to mflat_'+fil+'_'+amp+'_'+binn+'.fits')
    return

def process_science(sci_list,fil,amp,binn,red_path,mbias=None,mflat=None,proc=None,log=None):
    masks = []
    processed = []
    gains = gain(amp)
    readnoises = readnoise(amp)
    for sci in sci_list:
        log.info('Loading file: '+sci)
        log.info('Applying bias correction, gain correction and flat correction.')
        with fits.open(sci) as hdr:
            header = hdr[0].header
        if header['DETSEC01'] == '[1:2048,1:4096]':
            trim_sec = '[13:2060,1321:3921]'
        else:
            trim_sec = '[13:2060,1:2601]'
        raw = [CCDData.read(sci, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], oscan_model=models.Chebyshev1D(3), trim=trim_sec, gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if amp == '4':
            sci_full = CCDData(np.concatenate([red[0],np.zeros([2601,113]),red[1],np.zeros([2601,81]),red[2],np.zeros([2601,113]),red[3]],axis=1),header=header,unit=u.electron)
        if amp == '8':
            sci_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),red[3],np.fliplr(red[2]),red[5],np.fliplr(red[4]),red[7],np.fliplr(red[6])],axis=1),header=header,unit=u.electron)                
        log.info('Exposure time of science image is '+str(sci_full.header['ELAPTIME']))
        processed_data = ccdproc.flat_correct(sci_full, mflat)
        log.info('File proccessed.')
        log.info('Cleaning cosmic rays and creating mask.')
        mask = make_source_mask(processed_data, nsigma=3, npixels=5)
        masks.append(mask)
        # clean, com_mask = create_mask.create_mask(sci,processed_data,'_mask.fits',static_mask(proc),mask,saturation(red.header),binning(),np.mean(readnoise(amp)),cr_clean_sigclip(),cr_clean_sigcfrac(),cr_clean_objlim(),log)
        # processed_data.data = clean
        log.info('Calculating 2D background.')
        bkg = Background2D(processed_data, (128, 128), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
        log.info('Median background: '+str(np.median(bkg.background)))
        fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg.fits').replace('.gz',''),np.array(bkg.background),overwrite=True)
        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found').divide(processed_data.header['ELAPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
        log.info('Background subtracted and image divided by exposure time.')
        log.info('Writing WCS to file.')
        final = ccdproc.trim_image(final[700:2400,300:7607])
        ra = final.header['RA'].split(':')
        dec = final.header['DEC'].split(':')
        coords = SkyCoord(ra[0]+'h'+ra[1]+'m'+ra[2]+'s',dec[0]+'d'+dec[1]+'m'+dec[2]+'s',frame='icrs')
        final.header['RADECSYS'] = 'ICRS'
        final.header['CUNIT1'] = 'deg'
        final.header['CUNIT2'] = 'deg'
        final.header['CTYPE1'] = 'RA---TAN'
        final.header['CTYPE2'] = 'DEC--TAN'
        final.header['CRPIX1'] = 5145
        final.header['CRPIX2'] = 901
        final.header['CRVAL1'] = coords.ra.deg
        final.header['CRVAL2'] = coords.dec.deg
        final.header['CD1_1'] = 0.1185/3600*np.sin(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
        final.header['CD1_2'] = 0.1185/3600*np.cos(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
        final.header['CD2_1'] = 0.1185/3600*np.cos(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
        final.header['CD2_2'] = -0.1185/3600*np.sin(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
        final.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(final)[1],np.shape(final)[0]))        
        final.write(sci.replace('/raw/','/red/').replace('.gz',''),overwrite=True)
        processed.append(final)
    masks = [k[700:2400,300:7607] for k in masks]
    return processed, masks

def stacked_image(tar,red_path):
    return [red_path+tar+'.fits']

def rdnoise(header):
    return 2.60725

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
    return hdr['ELAPTIME']

def gain(amp):
    if amp=='4':
        gains = [1.206, 1.722, 1.211, 1.231]
    if amp=='8':
        gains = [1.232, 1.180, 1.714, 1.173, 1.161, 1.261, 1.246, 1.216]
    return gains

def readnoise(amp):
    if amp=='4':
        readnoise = [2.528, 2.128, 2.5395, 3.2335]
    if amp=='8':
        readnoise = [2.583, 2.473, 1.797, 2.459, 2.434, 2.645, 3.918, 2.549]
    return readnoise

def fringe_correction(fil):
    return False

def trim(f):
    return False

def edit_raw_headers(rawdir):

    for file in glob.glob(os.path.join(rawdir, '*.fits')):

        hdu = fits.open(file)
        h = hdu[0].header

        if ('ELAPTIME' not in h.keys() and 'DATE-BEG' in h.keys() and
            'DATE-END' in h.keys()):

            t1 = Time(h['DATE-BEG'])
            t2 = Time(h['DATE-END'])
            dt = t2 - t1

            hdu[0].header['ELAPTIME']=dt.to_value('sec')

        if (('twi' in h['OBJECT'].lower() and 'flat' in h['OBJECT'].lower()) or
            ('blank' in h['OBJECT'].lower()) or
            ('t_flat' in h['OBJECT'].lower()) or
            ('sky' in h['OBJECT'].lower() and 'flat' in h['OBJECT'].lower())):
            hdu[0].header['OBSTYPE']='DmFlat'

        if 'bias' in h['OBJECT'].lower():
            hdu[0].header['OBSTYPE']='Bias'

        if 'PONAME3' in h.keys() and h['PONAME3'].lower().strip()=='image':
            hdu[0].header['OBSMODE']='imaging'
        else:
            hdu[0].header['OBSMODE']='spectroscopy'


        hdu.writeto(file, overwrite=True, output_verify='silentfix')
