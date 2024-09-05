#parameter file for GMOS/Gemini
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

__version__ = 1.2 #last edited 09/11/2021

def name():
    return 'GMOS'

def min_exptime():
    return 1.0

def static_mask(paths):
    return [os.path.join(paths['code'], 'data', 'staticmasks', 
        'GMOS-N.12.22.staticmask.fits.fz')]

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.0803*2 #0.0807 0.08

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

# Keywords for selecting files from Sort_files object
# This allows a single file type to be used for multiple purposes (e.g., for
# generating a flat-field image from science exposures)
def filetype_keywords():
    return {'SCIENCE':'SCIENCE', 'FLAT':'FLAT', 'DARK':'DARK','BIAS':'BIAS'}

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
    return '22'

def time_format(hdr):
    return Time(hdr['DATE-OBS']+'T'+hdr['TIME-OBS']).mjd

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

    outlist = []
    for i,bias in enumerate(sorted(cal_list)):
        if log: log.info(f'Processing bias {bias}')
        basename = os.path.basename(bias)

        with fits.open(bias) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(bias, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
            oscan_model=models.Chebyshev1D(3), 
            trim=x.header['DATASEC'], gain=x.header['GAIN']*u.electron/u.adu, 
            readnoise=x.header['RDNOISE']*u.electron) for x in raw]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        
        outfile = os.path.join(red_path, basename)
        outfile = outfile.replace('.gz', '')
        outfile = outfile.replace('.bz2', '')

        bias_hdu.writeto(outfile, overwrite=True)
        outlist.append(outfile)

    bias_filename = get_mbias_name(red_path, amp, binn)
    
    mbias = [ccdproc.combine(outlist, hdu=x+1,unit=u.electron) 
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

def create_flat(flat_list, fil, amp, binn, red_path, mbias=None, log=None):

    log.info(f'Processing files for filter: {fil}')
    log.info(f'{len(flat_list)} files found.')

    scale = []
    flats = []
    for i, flat in enumerate(flat_list):
        flat = os.path.abspath(flat)
        log.info(f'Loading file: {flat}')

        with fits.open(flat) as hdr:
            header = hdr[0].header

        raw = [CCDData.read(flat, hdu=x+1, unit='adu') 
            for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
            oscan_model=models.Chebyshev1D(3), 
            trim=x.header['DATASEC'], 
            gain=x.header['GAIN']*u.electron/u.adu, 
            readnoise=x.header['RDNOISE']*u.electron, 
            master_bias=mbias[k], gain_corrected=True) 
                for k,x in enumerate(raw)]

        if header['INSTRUME'] == 'GMOS-S':
            if amp == '12':
                flat_full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],40])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],40])]+\
                    red[8:10],axis=1),
                    header=header,unit=u.electron)
            elif amp == '4':
                flat_full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)
        elif header['INSTRUME'] == 'GMOS-N':
            if amp == '12':
                flat_full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],30])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],30])]+\
                    red[8:10],axis=1),
                    header=header,unit=u.electron)
            elif amp == '4':
                flat_full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)

        exptime = flat_full.header['EXPTIME']
        log.info(f'Exposure time of image is {exptime}')
        norm = 1./np.nanmedian(flat_full[500:1500,700:1400]) #check for binning
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
    mflat=None, proc=None, staticmask=None, skip_skysub=False, log=None):

    if staticmask is not None and os.path.exists(staticmask[0]):
        hdu = fits.open(staticmask[0])
        imgmask = hdu[1].data.astype(bool)
    else:
        imgmask = None
    
    processed = []
    for sci in sci_list:
        sci = os.path.abspath(sci)
        if log: log.info(f'Loading file: {sci}')
        if log: log.info('Applying bias, gain, flat correction.')

        with fits.open(sci) as hdr:
            header = hdr[0].header
            if amp == '12':
                header_wcs = hdr[7].header
            elif amp == '4':
                header_wcs = hdr[3].header

        raw = [CCDData.read(sci, hdu=x+1, unit='adu') 
            for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
            oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], 
            gain=x.header['GAIN']*u.electron/u.adu, 
            readnoise=x.header['RDNOISE']*u.electron, master_bias=mbias[k], 
            gain_corrected=True) for k,x in enumerate(raw)]

        header += header_wcs
        if header['INSTRUME'] == 'GMOS-S':
            if amp == '12':
                sci_full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],40])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],40])]+\
                    red[8:10],axis=1),header=header,unit=u.electron)
                header['CRPIX1'] = header['CRPIX1']+(np.shape(red[0])[1]*4)+40
            elif amp == '4':
                sci_full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)
        elif header['INSTRUME'] == 'GMOS-N':
            if amp == '12':
                sci_full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],30])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],30])]+\
                    red[8:10],axis=1),header=header,unit=u.electron)
                header['CRPIX1'] = header['CRPIX1']+(np.shape(red[0])[1]*4)+40
            elif amp == '4':
                sci_full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)
        
        exptime = sci_full.header['EXPTIME']
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
            del final.header[key]

        final_basefile = get_base_science_name(sci)
        final_outfile_name = os.path.join(red_path, final_basefile)

        if log: log.info(f'Writing: {final_outfile_name}')

        final.write(final_outfile_name, overwrite=True)
        
        processed.append(final)

    return processed

def stack_method(hdr):
    return('median')

def detrend(header):
    return(True)

def rdnoise(header):
    return header['RDNOISE']

def binning():
    return [4,4]

def catalog_zp(hdr):
    coord = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.deg, u.deg))
    if coord.dec.degree < -30:
        return('SkyMapper')
    else:
        return('PS1')
    
    return('PS1')

def exptime(hdr):
    return hdr['EXPTIME']

def fringe_correction(fil):
    if 'z' in fil:
        return True
    else:
        return False
    
def trim(f):
    return False

def number_keyword(hdr):
    datestr = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
    elap = Time(datestr)-Time('1980-01-01')
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

def edit_stack_headers(stack):

    good_keys = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND',
        'INSTRUME','OBJECT','OBSTYPE','GEMPRGID','OBSERVER','TELESCOP',
        'EPOCH','SSA','RA','DEC','UT','DATE','DETECTOR','DATE-OBS',
        'AIRMASS','MJD-OBS','GAIN','RDNOISE','SATURATE','SKYBKG','BUNIT',
        'CTYPE1','CTYPE2','CUNIT1','CUNIT2','RADECSYS','CD1_1','CD1_2',
        'CD2_1','CD2_2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','NFILES','FILTER',
        'AMPS','BINNING','XTENSION','EXTNAME','TIMESYS','EXPTOT','WCSAXES']

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

def edit_raw_headers(files, log=None):
    pass
