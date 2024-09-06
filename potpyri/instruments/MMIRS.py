#parameter file for MMIRS/MMT
import os
import astropy
import datetime
import numpy as np

from photutils.background import Background2D
from photutils.background import MeanBackground

from astropy.stats import SigmaClip
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
import astropy.units.astrophys as u
import astropy.units as u
import ccdproc
from astropy.modeling import models

__version__ = 1.8 #last edited 09/11/2021

def name():
    return 'MMIRS'

def static_mask(paths):
    return [os.path.join(paths['code'], 'data', 'staticmasks', 
        'MMIRS.staticmask.fits.fz')]

def min_exptime():
    return 0.1

def stack_method(hdr):
    return 'median'

# Detector statistics
def pixscale(hdr):
    return 0.202

def saturation(hdr):
    return 46000.0 * hdr['GAIN']

def gain(hdr):
    return(hdr['GAIN'])

def rdnoise(hdr):
    return(hdr['RDNOISE'])

# Input file format
def raw_format(proc):
    return '*.fits'

# Which calibrations should we perform for this imager
def dark():
    return True

def bias():
    return False

def flat():
    return True

# Parameters for Sort_files
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

# Keywords for selecting files from Sort_files object
# This allows a single file type to be used for multiple purposes (e.g., for
# generating a flat-field image from science exposures)
def filetype_keywords():
    return {'SCIENCE':'SCIENCE', 'FLAT':'SCIENCE', 'DARK':'DARK','BIAS':'BIAS'}

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

def get_mdark_name(red_path, amp, binn):
    return(os.path.join(red_path, f'mdark_{amp}_{binn}.fits'))

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

# Loading bias and flat
def load_dark(red_path, amp, binn):
    dark = get_mdark_name(red_path, amp, binn)
    mdark = CCDData.read(dark)
    return mdark

def load_flat(flat):
    mflat = CCDData.read(flat)
    return mflat

# Takes filename and outputs CCDData object with raw image in units of e-
def import_image(filename, normalize=False):
    hdu = fits.open(filename)

    # Get header names to know how many up the ramp samples there are
    names = [h.name for h in hdu]

    # Get data shape from first image EXT for creating slope frame
    data_shape = hdu['IM1'].data.shape

    all_data = []
    all_times = []
    num_reads = 0
    for i in np.arange(len(hdu)):
        ext=f'IM{i+1}'
        if ext in names:
            all_data.append(hdu[ext].data)
            all_times.append(exptime(hdu[ext].header))
            num_reads+=1

    all_data = np.array(all_data)
    all_times = np.array(all_times)

    # Run a linear regression on 3D array all_data with respect to all_times
    A = np.vstack([all_times, np.ones(len(all_times))]).T
    data = np.linalg.lstsq(A, all_data.reshape(num_reads, -1), rcond=None)

    # Get slope and bias data
    slope = data[0][0].reshape(*data_shape)
    bias = data[0][1].reshape(*data_shape)

    # Rescale the slope image so it's in units of ADU
    # Provide in units of adu instead of adu/second
    if not normalize: slope = slope * np.max(all_times)

    slope = slope.astype(np.float32)

    header = hdu[0].header
    header.update(hdu[f'IM{num_reads}'].header)

    for key in ['BZERO','BSCALE','BITPIX']:
        if key in header.keys():
            del header[key]

    # Multiply by GAIN and apply RDNOISE
    rawdata = CCDData(slope, header=header, unit=u.adu)
    rawdata = ccdproc.ccd_process(rawdata, 
        gain=rawdata.header['GAIN']*u.electron/u.adu,
        readnoise=rawdata.header['RDNOISE']*u.electron)

    # Apply saturation
    rawdata.header['SATURATE'] = saturation(rawdata.header)

    return(rawdata)

def create_flat(flat_list, fil, amp, binn, red_path, mbias=None, mdark=None,
    log=None):

    log.info(f'Processing files for filter: {fil}')
    log.info(f'{len(flat_list)} files found.')

    flats = []
    scale = []

    for flat in flat_list:
        log.info(f'Loading file: {flat}')

        raw = import_image(flat)

        if log: log.info('Subtracting overscan')
        red = ccdproc.subtract_overscan(raw, overscan=raw[0:4,:], 
            overscan_axis=0)

        if log: log.info('Subtracting dark')
        red = ccdproc.subtract_dark(red, mdark, exposure_time='EXPTIME', 
            exposure_unit=u.second)

        exptime = red.header['EXPTIME']

        data_shape = red.data.shape
        norm = 1./np.nanmedian(red[int(data_shape[0]/4):int(3*data_shape[0]/4),
            int(data_shape[1]/4):int(3*data_shape[1]/4)])
        if log: log.info(f'Flat normalization: {norm}, Exptime*norm: {exptime*norm}')
        
        scale.append(norm)
        flats.append(red)

    mflat = ccdproc.combine(flats, method='median', scale=scale,
        sigma_clip=True)

    if log: log.info(f'Created master flat for filter: {fil}')

    mflat.header['VERSION'] = (__version__, 'Version of setting file.')

    flatname = get_mflat_name(red_path, fil, amp, binn)
    if log: log.info(f'Writing master flat to {flatname}')
    mflat.write(flatname, overwrite=True)
    
    return

def process_science(sci_list, fil, amp, binn, red_path, mbias=None, mflat=None,
    mdark=None, proc=None,log=None, staticmask=None, skip_skysub=False):

    processed = []
    for sci in sci_list:
        if log: log.info(f'Loading file: {sci}')
        raw = import_image(sci)

        if log: log.info('Subtracting overscan')
        red = ccdproc.subtract_overscan(raw, overscan=raw[0:4,:], 
            overscan_axis=0)

        if log: log.info('Subtracting dark')
        red = ccdproc.subtract_dark(red, mdark, exposure_time='EXPTIME', 
            exposure_unit=u.second)

        if log: log.info('Applying flat field correction')
        red = ccdproc.flat_correct(red, mflat)

        if log: log.info('Trimming image')
        processed_data = ccdproc.ccd_process(red, trim=raw.header['DATASEC'])

        # Make sure to take out *SEC keywords
        for key in ['DATASEC','CCDSEC','DETSEC','TRIMSEC','DETSIZE']:
            if key in processed_data.header.keys():
                del processed_data.header[key]

        if log: log.info('Calculating 2D background.')
        mean, median, stddev = sigma_clipped_stats(processed_data.data)
        # Add to input mask
        addmask = processed_data.data < median - 3 * stddev
        input_mask = addmask
        
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

        # Mask out excessively negative pixels
        mean, median, stddev  = sigma_clipped_stats(final.data)
        mask = final.data < median - 3 * stddev
        final.data[mask] = np.nan

        final.header['SATURATE'] -= med_background.value
        final.header['SKYBKG'] = med_background.value
        if log: log.info('Background subtracted and updated saturation.')

        final_basefile = get_base_science_name(sci)
        final_outfile_name = os.path.join(red_path, final_basefile)
        if log: log.info(f'Writing: {final_outfile_name}')

        # Edit additional header values
        final.header['FILTER']=fil
        final.header['AMPS']=amp 
        final.header['BINNING']=binn
        final.header['ORGFILE']=sci
        
        final.write(final_outfile_name, overwrite=True)
        
        processed.append(final)

    return processed

def get_base_science_name(image):

    hdu = fits.open(image, mode='readonly')
    hdr = hdu[1].header

    tkwd = target_keyword()
    target = hdr[tkwd].replace(' ','')
    fil = filter_keyword(hdr)[0]
    amp = amp_keyword(hdr)
    binn = bin_keyword(hdr)
    file_time = time_format(hdr)
    number = str(int(Time(hdr['DATE-OBS']).mjd*1e5)-5900000000)

    datestr = Time(file_time, format='mjd').datetime.strftime('ut%y%m%d')

    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.{number}.fits'

    return(filename)


def stacked_image(tar,red_path):
    return [red_path+tar+'.fits']

def rdnoise(header):
    return header['RDNOISE']

def binning():
    return [4,4]

def create_dark(dark_list, amp, binn, red_path, mbias=None, log=None):
    
    processed = []
    exptimes = []
    for dark in dark_list:
        log.info(f'Creating master dark with exposure time: ')
        log.info(f'Loading file: {dark}')
        raw = import_image(dark)

        if log: log.info('Subtracting overscan')
        red = ccdproc.subtract_overscan(raw, overscan=raw[0:4,:], 
            overscan_axis=0)

        processed.append(red)
        exptimes.append(red.header['EXPTIME'])

    exptimes = np.array(exptimes)
    
    if log: log.info('Creating master dark.')
    mdark = ccdproc.combine(processed, method='median', scale=1./exptimes)

    # Rescale to electrons by average of all exposure times
    avg_exptime = np.mean(exptimes)
    mdark.data *= avg_exptime
    mdark.header['EXPTIME']=(avg_exptime, 'Exposure time of master dark (sec).')
    
    # Update other header variables
    mdark.header['VER'] = (__version__, 
        'Version of telescope parameter file used.')    
    for i,im in enumerate(dark_list):
        mdark.header[f'FILE{i+1}']=os.path.basename(im)
    
    darkname = get_mdark_name(red_path, amp, binn)
    if log: log.info(f'Writing master dark to {darkname}')
    mdark.write(darkname, overwrite=True)
    
    return

def cr_clean_sigclip():
    return 50

def cr_clean_sigcfrac():
    return 0.1

def cr_clean_objlim():
    return 100

def run_phot():
    return True

def catalog_zp(hdr):
    return '2MASS'

def min_cat_mag(hdr):
    return 13.5

def exptime(hdr):
    return hdr['EXPTIME']

def trim(f):
    return False

def edit_raw_headers(files, log=None):


    bad_keys = ['PV2_1','PV2_2','PV2_3','PV2_4','PV2_5','CD1_1','CD1_2',
        'CD2_1','CD2_2','DTV1','DTV2','DTM1_1','DTM2_2','LTM1_1','LTM1_2',
        'LTM2_1','LTM2_2','LTV1','LTV2','CTYPE1','CTYPE2','RADECSYS',
        'CRPIX1','CRPIX2','SECPIX1','SECPIX2','DATAMAX','DATAMIN']

    for file in files:
        hdu = fits.open(file)

        for i,h in enumerate(hdu):
            while any([k in bad_keys for k in h.header.keys()]):
                for key in h.header.keys():
                    if key in bad_keys:
                        if key in h.header.keys():
                            del h.header[key]

        hdu.writeto(file, overwrite=True, output_verify='silentfix')

def edit_stack_headers(stack):
    return(stack)

def stacked_image(row, red_path):
    target = row['Target']
    fil = row['Filter'][0]
    amp = row['Amp']
    binn = row['Binning']
    mjd = row['Time']

    datestr = Time(mjd, format='mjd').datetime.strftime('ut%y%m%d')
    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.stk.fits'

    return os.path.join(red_path, filename)

def detrend(header):
    return(True)

def out_size():
    return 2500

