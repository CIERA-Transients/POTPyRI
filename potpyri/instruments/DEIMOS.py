#parameter file for DEIMOS/Keck
import os
import astropy
import datetime
import copy
import ccdproc
import numpy as np

from photutils import Background2D
from photutils import MeanBackground

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import models
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.time import Time

import astropy.units as u
import astropy.wcs as wcs

__version__ = 1.5 #last edited 23/11/2021

def name():
    return 'DEIMOS'

def min_exptime():
    return 1.0

def static_mask(proc):
    return ['']

def wcs_extension():
    return 0

def pixscale(hdr):
    return 0.1185

def saturation(hdr):
    return 65000

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

# Keywords for selecting files from Sort_files object
# This allows a single file type to be used for multiple purposes (e.g., for
# generating a flat-field image from science exposures)
def filetype_keywords():
    return {'SCIENCE':'SCIENCE', 'FLAT':'FLAT', 'DARK':'DARK','BIAS':'BIAS'}

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

# Bias and flat names
def get_mbias_name(red_path, amp, binn):
    return(os.path.join(red_path, f'mbias_{amp}_{binn}.fits'))

def get_mflat_name(red_path, fil, amp, binn):
    return(os.path.join(red_path, f'mflat_{fil}_{amp}_{binn}.fits'))

# Loading bias and flat
def load_bias(red_path, amp, binn):
    bias = get_mbias_name(red_path, amp, binn)
    mbias = [CCDData.read(bias,hdu=x+1,unit=u.electron) 
        for x in range(int(amp))]
    return mbias

def load_flat(flat):
    mflat = CCDData.read(flat)
    return mflat

def create_bias(cal_list, amp, binn, red_path, log):

    if log:
        log.info(f'Processing bias files with {amp} amps and {binn} binning.')
    if log:
        log.info(f'{len(cal_list)} files found.')


    gains = gain(amp)
    readnoises = readnoise(amp)

    outlist = []
    for i,bias in enumerate(sorted(cal_list)):
        if log: log.info(f'Processing bias {bias}')
        basename = os.path.basename(bias)

        with fits.open(bias) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(bias, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], 
            oscan_model=models.Chebyshev1D(3), 
            trim='[13:2060,1:2601]', 
            gain=gains[j]*u.electron/u.adu, 
            readnoise=readnoises[j]*u.electron) for j,x in enumerate(raw)]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        
        outfile = os.path.join(red_path, basename)
        outfile = outfile.replace('.gz', '')
        outfile = outfile.replace('.bz2', '')

        bias_hdu.writeto(outfile, overwrite=True)
        outlist.append(outfile)

    bias_filename = get_mbias_name(red_path, amp, binn)
    
    mbias = [ccdproc.combine(outlist, hdu=x+1, unit=u.electron) 
        for x in range(int(amp))]
    mbias_hdu = fits.HDUList([fits.PrimaryHDU()])
    
    for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
    
    mbias_hdu[0].header['VER'] = (__version__, 
        'Version of telescope parameter file used.')
    mbias_hdu.writeto(bias_filename,overwrite=True)

    if log:
        log.info(f'Created master bias with {amp} amps and {binn} binning.')

    if log:
        log.info(f'Master bias written to {bias_filename}')
    
    for bias in outlist: os.remove(bias)
    
    return

def create_flat(flat_list, fil, amp, binn, red_path, mbias=None, mdark=None, 
    log=None):

    log.info(f'Processing files for filter: {fil}')
    log.info(f'{len(flat_list)} files found.')

    scale = []
    flats = []

    gains = gain(amp)
    readnoises = readnoise(amp)

    for i, flat in enumerate(flat_list):
        flat = os.path.abspath(flat)
        log.info(f'Loading file: {flat}')

        with fits.open(flat) as hdr:
            header = hdr[0].header
        if header['DETSEC01'] == '[1:2048,1:4096]':
            trim_sec = '[13:2060,1321:3921]'
        else:
            trim_sec = '[13:2060,1:2601]'
        raw = [CCDData.read(flat, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], 
            oscan_model=models.Chebyshev1D(3), trim=trim_sec, 
            gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, 
            master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if amp == '4':
            flat_full = CCDData(np.concatenate([red[0],np.zeros([2601,113]),
                red[1],np.zeros([2601,81]),red[2],np.zeros([2601,113]),
                red[3]],axis=1),header=header,unit=u.electron/u.second)
        if amp == '8':
            flat_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),
                red[3],np.fliplr(red[2]),red[5],np.fliplr(red[4]),red[7],
                np.fliplr(red[6])],axis=1),header=header,unit=u.electron/u.second)

        exptime = flat_full.header['ELAPTIME']
        log.info(f'Exposure time of image is {exptime}')
        norm = 1./np.nanmedian(flat_full[1200:1800,2500:3500]) #check for binning
        log.info(f'Flat normalization: {norm}')
        scale.append(norm)
        flats.append(flat_full)
    
    mflat = ccdproc.combine(flats,method='median',scale=scale,sigma_clip=True)
    log.info(f'Created master flat for filter: {fil} and {amp} amp {binn} binning.')
    mflat.header['VER'] = (__version__, 
        'Version of telescope parameter file used.')

    flat_filename = get_mflat_name(red_path, fil, amp, binn)
    mflat.write(flat_filename, overwrite=True)
    log.info(f'Master flat written to {flat_filename}')
    
    return

def process_science(sci_list, fil, amp, binn, red_path, mbias=None,
    mflat=None, mdark=None, proc=None, staticmask=None, skip_skysub=False, 
    log=None):
    
    processed = []
    gains = gain(amp)
    readnoises = readnoise(amp)

    if staticmask is not None and os.path.exists(staticmask[0]):
        hdu = fits.open(staticmask[0])
        imgmask = hdu[1].data.astype(bool)
    else:
        imgmask = None

    for sci in sci_list:
        sci = os.path.abspath(sci)
        if log: log.info(f'Loading file: {sci}')
        if log: log.info('Applying bias, gain, flat correction.')

        with fits.open(sci) as hdr:
            header = hdr[0].header
        if header['DETSEC01'] == '[1:2048,1:4096]':
            trim_sec = '[13:2060,1321:3921]'
        else:
            trim_sec = '[13:2060,1:2601]'
        raw = [CCDData.read(sci, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], 
            oscan_model=models.Chebyshev1D(3), trim=trim_sec, 
            gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, 
            master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
        if amp == '4':
            sci_full = CCDData(np.concatenate([red[0],np.zeros([2601,113]),
                red[1],np.zeros([2601,81]),red[2],np.zeros([2601,113]),
                red[3]],axis=1),header=header,unit=u.electron)
        if amp == '8':
            sci_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),
                red[3],np.fliplr(red[2]),red[5],np.fliplr(red[4]),red[7],
                np.fliplr(red[6])],axis=1),header=header,unit=u.electron)

        exptime = sci_full.header['ELAPTIME']
        if log: log.info(f'Exposure time of science image is {exptime}.')

        # Add SATURATE keyword to sci_full
        sci_full.header['SATURATE'] = saturation(hdr)

        processed_data = ccdproc.flat_correct(sci_full, mflat)
        
        if log: log.info('Calculating 2D background.')
        mean, median, stddev = sigma_clipped_stats(processed_data.data)
        # Add to input mask
        addmask = processed_data.data < median - 3 * stddev

        if imgmask is None:
            input_mask = addmask
        else:
            input_mask = imgmask | addmask
        
        bkg = Background2D(processed_data, (128, 128), filter_size=(3, 3),
            sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(),
            exclude_percentile=80, mask=input_mask, fill_value=median)
        med_background = np.nanmedian(bkg.background)
        if log: log.info(f'Median background: {med_background}')

        bkg_basefile = get_base_science_name(sci).replace('.fits','_bkg.fits')
        bkg_outfile_name = os.path.join(red_path, bkg_basefile)
        fits.writeto(bkg_outfile_name, np.array(bkg.background), overwrite=True)
        if log: log.info(f'Wrote background to: {bkg_outfile_name}')

        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),
            propagate_uncertainties=True, handle_meta='first_found')

        # Apply static mask
        if imgmask is not None:
            final.data[imgmask] = np.nan

        # Mask out excessively negative pixels
        mean, median, stddev  = sigma_clipped_stats(final.data, mask=imgmask)
        mask = final.data < median - 3 * stddev
        final.data[mask] = np.nan

        final.header['SATURATE'] -= med_background.value
        final.header['SKYBKG'] = med_background.value
        if log: log.info('Background subtracted and updated saturation.')

        # Edit additional header values
        final.header['FILTER']=fil
        final.header['AMPS']=amp 
        final.header['BINNING']=binn

        for key in ['CCDSEC','DATASEC','BIASSEC']:
            if key in final.header.keys():
                del final.header[key]

        final_basefile = get_base_science_name(sci)
        final_outfile_name = os.path.join(red_path, final_basefile)

        if log: log.info(f'Writing: {final_outfile_name}')

        final.write(final_outfile_name, overwrite=True)
        
        processed.append(final)
        
    return processed

def number_keyword(hdr):
    elap = Time(hdr['MJD'], format='mjd')-Time('1980-01-01')
    elap = int(np.round(elap.to(u.second).value))

    return elap

def get_base_science_name(image):

    hdu = fits.open(image, mode='readonly')
    hdr = hdu[0].header

    tkwd = target_keyword()
    target = hdr[tkwd].replace(' ','')
    fil = filter_keyword(hdr)[0]
    amp = amp_keyword(hdr)
    binn = bin_keyword(hdr)
    file_time = time_format(hdr)
    number = number_keyword(hdr)

    datestr = Time(file_time, format='mjd').datetime.strftime('ut%y%m%d')

    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.{number}.fits'

    return(filename)

def stacked_image(row, red_path):
    target = row['Target']
    fil = row['Filter'][0]
    amp = row['Amp']
    binn = row['Binning']
    mjd = row['Time']

    datestr = Time(mjd, format='mjd').datetime.strftime('ut%y%m%d')
    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.stk.fits'

    return os.path.join(red_path, filename)

def rdnoise(header):
    return 2.60725

def binning():
    return [4,4]

def catalog_zp(hdr):
    coord = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.deg, u.deg))
    if coord.dec.degree < -30:
        return('SkyMapper')
    else:
        return('PS1')
    
    return('PS1')

def detrend(header):
    return True

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

def edit_raw_headers(files, log=None):

    for file in files:

        basefile = os.path.basename(file)
        if basefile.startswith('DE') and basefile.endswith('.gz'): continue

        if log: log.info(f'Editing headers for: {file}')

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

def edit_stack_headers(stack):

    good_keys = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND',
        'OBSERVAT','TELESCOP','ORIGIN','INSTRUME','EXTNAME','FILTER',
        'DATE-OBS','RA','DEC','POSANG','AZ','EL','ROT','PA',
        'LST','UT','AIRMASS','HA','GSEEING','WSEEING','IMAGETYP',
        'OBJECT','PI','PROPID','GAIN','RADECSYS','SATURATE','SKYBKG',
        'BUNIT','FILENAME','RADISP','DEDISP','MJD-OBS','RADESYS','EXPTOT',
        'RDNOISE','NFILES','OBSTYPE','WCSAXES','CTYPE1','CTYPE2','EQUINOX',
        'LONPOLE','LATPOLE','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CUNIT1',
        'CUNIT2','CD1_1','CD1_2','CD2_1','CD2_2','FWHM','DATE','SKYADU',
        'SKYSIG','NPSFSTAR','NOBJECT','ZPTNSTAR','ZPTMAG','ZPTMUCER',
        'M3SIGMA','M5SIGMA','M10SIGMA','CDELT1','CDELT2']

    # This deletes all keys in the header that are not in good_keys
    for i,h in enumerate(stack):
        while True:
            cont = True
            for key in stack[i].header.keys():
                if key in stack[i].header.keys() and key not in good_keys:
                    del stack[i].header[key]
                    cont = False
            if cont:
                break

    return(stack)

def out_size():
    return 9000

def stack_method(hdr):
    return 'median'
