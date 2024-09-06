#parameter file for BINOSPEC/MMT
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

__version__ = 1.4 #last edited 24/08/2021

def name():
    return 'Binospec'

def static_mask(paths):
    return [os.path.join(paths['code'], 'data', 'staticmasks', 
            'bino_proc.trim.staticmask.fits.fz')]

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.24

def saturation(hdr):
    return 65000 #defualt hdr['DATAMAX']*hdr['GAIN']

def min_exptime():
    return 1.0

def cal_path():
    return str(os.getenv("PIPELINE_HOME"))+'/Imaging_pipelines/BINOSPEC_calib/'

def raw_format(proc):
    if proc:
        return 'sci_img_*.fits'
    else:
        return 'sci_img*[!proc].fits'

def stack_method(hdr):
    return 'average'

def dark():
    return False

def bias():
    return False

def flat():
    return True

def raw_header_ext():
    return 1

# Keywords for selecting files from Sort_files object
# This allows a single file type to be used for multiple purposes (e.g., for
# generating a flat-field image from science exposures)
def filetype_keywords():
    return {'SCIENCE':'SCIENCE', 'FLAT':'FLAT', 'DARK':'DARK','BIAS':'BIAS'}

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
    return '2'

def bin_keyword(hdr):
    return hdr['CCDSUM'].replace(' ','')

def time_format(hdr):
    return hdr['MJD']

def wavelength():
    return 'OPT'

def get_mbias_name(red_path, amp, binn):
    return(os.path.join(red_path, f'mbias_{amp}_{binn}.fits'))

def get_mflat_name(red_path, fil, amp, binn):
    return(os.path.join(red_path, f'mflat_{fil}_{amp}_{binn}.fits'))

def load_bias(red_path, amp, binn):
    bias = get_mbias_name(red_path, amp, binn)
    mbias = [CCDData.read(bias,hdu=x+1,unit=u.electron) for x in range(int(amp[0]))]
    return mbias

def load_flat(flat):
    mflat = CCDData.read(flat, hdu=0)
    return mflat

def gain():
    return [1.02, 1.08]

def datasec():
    return ['[1055:3024,217:3911]','[1067:3033,217:3911]']

def biassec():
    return ['[1:1054,1:4112]','[3034:4096,1:4112]']

def create_flat(flat_list, fil, amp, binn, red_path, mbias=None, log=None):

    log.info(f'Processing files for filter: {fil}')
    log.info(f'{len(flat_list)} files found.')

    scale = []
    flats = []
    for i, flat in enumerate(flat_list):
        flat = os.path.abspath(flat)
        log.info(f'Loading file: {flat}')

        with fits.open(flat) as hdr:
            header = hdr[1].header

        raw = [CCDData.read(flat, hdu=x+1, unit='adu') 
            for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=biassec()[k], 
            oscan_model=models.Chebyshev1D(3), 
            trim=datasec()[k], 
            gain=gain()[k]*u.electron/u.adu, 
            readnoise=rdnoise(hdr)*u.electron, 
            gain_corrected=True) for k,x in enumerate(raw)]

        flat_full = CCDData(np.concatenate((red[0], 
            np.empty((red[0].shape[0], 794)), red[1]), axis=1), 
            header=header,unit=u.electron)

        exptime = flat_full.header['EXPTIME']
        log.info(f'Exposure time of image is {exptime}')
        norm = 1./np.nanmedian(flat_full) #check for binning
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

def get_base_science_name(image):

    hdu = fits.open(image, mode='readonly')
    hdr = hdu[1].header

    tkwd = target_keyword()
    target = hdr[tkwd].replace(' ','')
    fil = filter_keyword(hdr)
    amp = amp_keyword(hdr)
    binn = bin_keyword(hdr)
    file_time = time_format(hdr)
    number = str(int(Time(hdr['DATE-OBS']).mjd*1e5)-5900000000)
    datestr = Time(file_time, format='mjd').datetime.strftime('ut%y%m%d')

    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.{number}.fits'

    return(filename)

def process_science(sci_list, fil, amp, binn, red_path, mbias=None,
    mflat=None, proc=None, staticmask=None, skip_skysub=False, log=None):

    processed = []

    if staticmask is not None and os.path.exists(staticmask[0]):
        hdu = fits.open(staticmask[0])
        imgmask = hdu[1].data.astype(bool)
    else:
        imgmask = None

    flat_image = CCDData(np.copy(mflat),unit=u.electron,
        meta=mflat.meta,mask=mflat.mask,uncertainty=mflat.uncertainty)
    
    for sci in sorted(sci_list):
            sci = os.path.abspath(sci)
            if log: log.info(f'Loading file: {sci}')
            if log: log.info('Applying bias correction, gain correction and flat correction.')

            with fits.open(sci) as hdr:
                header = hdr[1].header

            raw = [CCDData.read(sci, hdu=x+1, unit='adu') 
                for x in range(int(amp))]
            red = [ccdproc.ccd_process(x, 
                trim=datasec()[k], 
                gain=gain()[k]*u.electron/u.adu, 
                readnoise=rdnoise(hdr)*u.electron, 
                gain_corrected=True) for k,x in enumerate(raw)]

            sci_full = CCDData(np.concatenate((red[0], 
                np.empty((red[0].shape[0], 794)), red[1]), axis=1), 
                header=header,unit=u.electron)

            sci_full.header['SATURATE'] = saturation(hdr)
            processed_data = ccdproc.flat_correct(sci_full, flat_image)

            if log: log.info('File proccessed.')

            mean, median, stddev = sigma_clipped_stats(processed_data.data)
            # Add to input mask
            addmask = (processed_data.data < median - 3 * stddev) |\
                      np.isnan(processed_data.data) |\
                      np.isinf(processed_data.data) |\
                      (processed_data.data==0.0)

            if imgmask is None:
                input_mask = addmask
            else:
                input_mask = imgmask | addmask
        
            if log: log.info('Calculating 2D background.')
            bkg = Background2D(processed_data, (64, 64), filter_size=(3, 3),
                sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(),
                exclude_percentile=80, mask=input_mask, fill_value=median)
            
            med_background = np.nanmedian(bkg.background)
            if log: log.info(f'Median background: {med_background}')

            bkg_basefile = get_base_science_name(sci).replace('.fits','_bkg.fits')
            bkg_basefile = bkg_basefile.replace('.gz','')
            bkg_outfile_name = os.path.join(red_path, bkg_basefile)
            fits.writeto(bkg_outfile_name,np.array(bkg.background),overwrite=True)
            if log: log.info(f'Wrote background to: {bkg_outfile_name}')

            if not skip_skysub:
                final = processed_data.subtract(CCDData(bkg.background,
                    unit=u.electron), propagate_uncertainties=True, 
                    handle_meta='first_found')
                final.header['SATURATE'] -= med_background.value
                final.header['SKYBKG'] = med_background.value
                if log: log.info('Background subtracted and updated saturation.')
            else:
                final = processed_data
                final.header['SKYBKG'] = 0.0

            # Mask out bad pixels
            mask = np.isnan(final.data) | np.isinf(final.data)
            final.data[mask] = 0.0
            if imgmask is not None:
                final.data[imgmask] = 0.0

            final_basefile = get_base_science_name(sci)
            final_outfile_name = os.path.join(red_path, final_basefile)

            if log: 
                log.info(f'Writing: {final_outfile_name}')
            else:
                print(f'Writing: {final_outfile_name}')
            
            final.write(final_outfile_name,overwrite=True)
            processed.append(final)

    return processed

def stacked_image(row, red_path):
    target = row['Target']
    fil = row['Filter']
    amp = row['Amp']
    binn = row['Binning']
    mjd = row['Time']

    datestr = Time(mjd, format='mjd').datetime.strftime('ut%y%m%d')
    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.stk.fits'

    return os.path.join(red_path, filename)

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

def detrend(header):
    return True

def catalog_zp(hdr):
    return 'PS1'

def exptime(hdr):
    return hdr['EXPTIME']

def fringe_correction(fil):
    if fil == 'z':
        return True
    else:
        return False

def trim(f):
    return False

def edit_raw_headers(files, log=None):

    bad_keys = ['WCSNAMEA','CUNIT1A','CUNIT2A','CTYPE1A','CTYPE2A',
        'CRPIX1A','CRPIX2A','CRVAL1A','CRVAL2A','CD1_1A','CD1_2A',
        'CD2_1A','CD2_2A','SOFTWARE','VER','DATASEC']

    for file in files:
        hdu = fits.open(file)
        for i,h in enumerate(hdu):
            hdr = copy.copy(h.header)
            for key in hdr.keys():
                if key in bad_keys:
                    del hdu[i].header[key]

            # Check values of CRVAL1 and CRVAL2, which can be weird for some 
            # files
            if 'CTYPE1' in hdr.keys() and hdr['CTYPE1']=='RA---TAN' and (hdr['CRVAL1']<0.0 or hdr['CRVAL1']>360.0):
                coord = SkyCoord(hdr['LST'], '31:41:20.04', unit=(u.hour, u.deg))
                hdu[i].header['CRVAL1']=coord.ra.degree

            if 'CTYPE2' in hdr.keys() and hdr['CTYPE2']=='DEC--TAN' and (hdr['CRVAL2']<-90.0 or hdr['CRVAL2']>90.0):
                hdu[i].header['CRVAL2']=31.6889


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
