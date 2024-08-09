#parameter file for LRIS/Keck
import os
import astropy
import datetime
import numpy as np
import ccdproc
import glob

from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.modeling import models
from astropy.stats import SigmaClip
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData
from astropy.wcs import WCS
import astropy.units as u

from photutils import Background2D
from photutils import MeanBackground
from photutils.segmentation import SegmentationImage
from photutils.segmentation import detect_sources
from photutils.segmentation import detect_threshold
from photutils.utils import circular_footprint

__version__ = 1.9 #last edited 18/11/2022


class LRIS():

    def __init__(self, proc=None, ):

        self.name = 'LRIS'
        
        # Detector parameters
        self.pixel_scale = 0.135

        # See LRIS non-linearity: 
        # https://www2.keck.hawaii.edu/inst/lris/lris-red-upgrade-notes_vol2.html
        self.saturation = 55000.0

        self.static_mask = ''

        # Which calibrations do we require
        self.use_dark = False
        self.use_bias = True
        self.use_flat = True

        if str(proc)=='True':
            self.raw_format = '*.fits.gz'
        elif str(proc)=='raw':
            self.raw_format = '*[b,r]*.fits'

def name():
    return 'LRIS'

def static_mask(proc):
    return ''

def pixscale():
    return 0.135 #arcsec/pixel

# Minimum exposure time for science exposures
def min_exptime():
    return 15.0

def saturation(hdr):
    # See LRIS non-linearity: 
    # https://www2.keck.hawaii.edu/inst/lris/lris-red-upgrade-notes_vol2.html
    return 55000

def raw_header_ext():
    return 0

def raw_format(proc):
    if str(proc)=='True':
        return '*.fits.gz'
    elif str(proc)=='raw':
        return '*[b,r]*.fits'

def dark():
    return False

def bias():
    return True

def flat():
    return True

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
        amp = '1' #str(hdr['NVIDINP'])
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

def get_mbias_name(red_path, amp, binn):
    return(os.path.join(red_path, f'mbias_{amp}_{binn}.fits'))

def get_mflat_name(red_path, fil, amp, binn):
    return(os.path.join(red_path, f'mflat_{fil}_{amp}_{binn}.fits'))

def load_bias(red_path, amp, binn):
    bias = get_mbias_name(red_path, amp, binn)
    mbias = [CCDData.read(bias,hdu=x+1,unit=u.electron) for x in range(int(amp[0]))]
    return mbias


def import_file(filename, amp):

    hdu = fits.open(filename)
    raw = []
    if amp=='1R':
        hdu[0].header['CUNIT1'] = 'deg'
        hdu[0].header['CUNIT2'] = 'deg'

        image = CCDData(hdu[0].data, meta=hdu[0].header, 
            wcs=WCS(hdu[0].header), unit='adu')
        raw.append(image)
    elif amp=='4B':
        for i in np.arange(4):
            for key in hdu[0].header.keys():
                if key not in hdu[i+1].header.keys():
                    try:
                        hdu[i+1].header[key]=hdu[0].header[key]
                    except ValueError:
                        pass

            hdu[i+1].header['CUNIT1'] = 'deg'
            hdu[i+1].header['CUNIT2'] = 'deg'
            image = CCDData(hdu[i+1].data, meta=hdu[i+1].header, 
                wcs=WCS(hdu[i+1].header), unit='adu')
            raw.append(image)

    return(raw)

def create_bias(cal_list, amp, binn, red_path, log):

    gains = gain(amp)
    readnoises = readnoise(amp)
    oscan_reg = overscan_region(amp)
    log.info(f'Processing bias files with {amp} amps and {binn} binning.')
    log.info(f'{len(cal_list)} files found.')

    outlist = []
    for i,bias in enumerate(cal_list):
        basename = os.path.basename(bias)

        with fits.open(bias) as hdr:
            header = hdr[0].header
        if amp=='1R':
            raw = import_file(bias, amp)
            red = [ccdproc.ccd_process(x, oscan=oscan_reg, 
                oscan_model=models.Chebyshev1D(3), 
                gain=gains[j]*u.electron/u.adu, 
                readnoise=readnoises[j]*u.electron) for j,x in enumerate(raw)]
        else:
            raw = import_file(bias, amp)
            red = [ccdproc.ccd_process(x, oscan=oscan_reg, 
                oscan_model=models.Chebyshev1D(3), 
                trim=x.header['DATASEC'], gain=gains[j]*u.electron/u.adu, 
                readnoise=readnoises[j]*u.electron) for j,x in enumerate(raw)]
        
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))

        outfile = os.path.join(red_path, basename).replace('.gz','')
        bias_hdu.writeto(outfile, overwrite=True)
        outlist.append(outfile)

    bias_filename = get_mbias_name(red_path, amp, binn)
    mbias = [ccdproc.combine(outlist, hdu=x+1, unit=u.electron, method='median',
        sigma_clip=True, sigma_clip_func=np.ma.median) 
        for x in range(int(amp[0]))]
    mbias_hdu = fits.HDUList([fits.PrimaryHDU()])
    for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
    mbias_hdu[0].header['VER'] = (__version__, 'Version of telescope parameter file used.')
    mbias_hdu.writeto(bias_filename,overwrite=True)

    log.info(f'Created master bias with {amp} amps and {binn} binning.')

    log.info(f'Master bias written to {bias_filename}')
    for bias in outlist: os.remove(bias)
    
    return

def load_flat(flat):
    mflat = CCDData.read(flat,unit=u.electron)
    return mflat

def format_datasec(sec_string, binning=1):
    sec_string = sec_string.replace('[','').replace(']','')
    x,y = sec_string.split(',')
    x1,x2 = x.split(':')
    y1,y2 = y.split(':')

    x1 = float(x1) ; x2 = float(x2) ; y1 = float(y1) ; y2 = float(y2)

    x1 = np.max([int(1.0*x1/binning),1])
    x2 = int(1.0*x2/binning)
    y1 = np.max([int(1.0*y1/binning),1])
    y2 = int(1.0*y2/binning)

    sec_string = f'[{x1}:{x2},{y1}:{y2}]'

    return(sec_string)

def create_flat(flat_list, fil, amp, binn, red_path, mbias=None, log=None):
    log.info(f'Processing files for filter: {fil}')
    log.info(f'{len(flat_list)} files found.')
    scale = []
    flats = []
    gains = gain(amp)
    readnoises = readnoise(amp)
    oscan_reg = overscan_region(amp)
    for i, flat in enumerate(flat_list):
        flat = os.path.abspath(flat)
        log.info(f'Loading file: {flat}')
        with fits.open(flat) as hdr:
            header = hdr[0].header
        
        try:
            n_sat = header['NPIXSAT']
        except KeyError:
            log.error('''No saturation keyword in header, continuing anyway. 
                Please manually check flats and rerun if needed.''')
            n_sat = 0
        if n_sat < 10000:
            if amp == '1R':
                bin1,bin2 = header['BINNING'].split(',')
                bin1 = float(bin1) ; bin2=float(bin2)

                raw = import_file(flat, amp)
                red = [ccdproc.ccd_process(x, oscan=oscan_reg, 
                    oscan_model=models.Chebyshev1D(3), 
                    gain=gains[k]*u.electron/u.adu, 
                    readnoise=readnoises[k]*u.electron, master_bias=mbias[k], 
                    gain_corrected=True) for k,x in enumerate(raw)]
                
                sl1 = str_to_slice(get_1R_datasec(1,binning=bin1))
                sl2 = str_to_slice(get_1R_datasec(2,binning=bin1))
                sl3 = str_to_slice(get_1R_datasec(3,binning=bin2))
                sl4 = str_to_slice(get_1R_datasec(4,binning=bin2))

                flat_full = CCDData(np.concatenate([np.concatenate(
                    [red[0][sl1],red[0][sl2]],axis=1),
                    np.concatenate([red[0][sl3],red[0][sl4]],axis=1)],
                    axis=0),header=header,unit=u.electron)
            else:
                raw = import_file(flat, amp)
                red = [ccdproc.ccd_process(x, oscan=oscan_reg, oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
                if amp == '4B':
                    flat_full = CCDData(np.concatenate([red[0],
                        np.fliplr(red[1]),
                        np.zeros([np.shape(red[1])[0],111]),
                        red[2],np.fliplr(red[3])],axis=1),
                        header=header,unit=u.electron)
                    flat_full = ccdproc.trim_image(flat_full[700:3120,396:3940])
                if amp == '4R':
                    flat_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),np.zeros([np.shape(red[0])[0],200]),red[3],np.fliplr(red[2])],axis=1),header=header,unit=u.electron)
            log.info('Exposure time of image is '+str(flat_full.header['ELAPTIME']))
            norm = 1/np.median(flat_full[1200:1600,1200:1600]) #check for binning
            log.info('Flat normalization: '+str(norm))
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

# Data sections for new red amplifier
# [SS] This needs to be read in from the DSEC keyworkds.
def get_1R_datasec(amp, binning=1):

    if binning==1:
        if amp==1:
            return('[830:2064,284:2064]')
        elif amp==2:
            return('[830:2064,2170:3956]')
        elif amp==3:
            return('[2185:3450,284:2064]')
        elif amp==4:
            return('[2185:3450,2170:3956]')
    elif binning==2:
        if amp==1:
            return('[437:1032,146:1032]')
        elif amp==2:
            return('[437:1032,1177:2069]')
        elif amp==3:
            return('[1152:1784,146:1032]')
        elif amp==4:
            return('[1152:1784,1177:2069]')



def str_to_slice(sec_string):

    sec_string = sec_string.replace('[','').replace(']','')
    x,y = sec_string.split(',')
    x1,x2 = x.split(':')
    y1,y2 = y.split(':')

    x1 = int(x1) ; x2 = int(x2) ; y1 = int(y1) ; y2 = int(y2)

    sl = np.s_[x1:x2,y1:y2]

    return(sl)

def process_science(sci_list, fil, amp, binn, red_path, mbias=None,
    mflat=None, proc=None, log=None):
    
    processed = []
    gains = gain(amp)
    readnoises = readnoise(amp)
    oscan_reg = overscan_region(amp)
    for sci in sorted(sci_list):
        sci = os.path.abspath(sci)
        if log: log.info(f'Loading file: {sci}')
        if log: log.info('Applying bias correction, gain correction and flat correction.')
        with fits.open(sci) as hdr:
            header = hdr[0].header
            if amp == '4B' and hdr[1].header['DATASEC'] != '[52:1075,1:4096]':
                window = True
            elif amp == '4R' and hdr[1].header['DATASEC'] != '[8:1031,1:2500]':
                window = True
            else:
                window = False
            if window:
                if hdr[1].header['DATASEC'] == mbias[0].header['DATASEC']:
                    trim_sec = [mbias[k].header['DATASEC'] for k in range(len(mbias))]
                else:
                    trim_sec = [hdr[k+1].header['DATASEC'] for k in range(len(mbias))]
                    for k,x in enumerate(mbias):
                        x.data = x.data[0:int(trim_sec[k].split(':')[-1].rstrip(']')),0:int(trim_sec[k].split(':')[1].split(',')[0])]
                        mbias[k] = x
            else:
                if amp == '1R':
                    trim_sec = None
                else:
                    trim_sec = [hdr[k+1].header['DATASEC'] for k in range(len(mbias))]
        if amp == '1R':
            bin1,bin2 = header['BINNING'].split(',')
            bin1 = float(bin1) ; bin2=float(bin2)
                
            raw = import_file(sci, amp)
            red = [ccdproc.ccd_process(x, oscan=oscan_reg, oscan_model=models.Chebyshev1D(3), gain=gains[k]*u.electron/u.adu, readnoise=readnoises[k]*u.electron, master_bias=mbias[k], gain_corrected=True) for k,x in enumerate(raw)]
            
            sl1 = str_to_slice(get_1R_datasec(1,binning=bin1))
            sl2 = str_to_slice(get_1R_datasec(2,binning=bin1))
            sl3 = str_to_slice(get_1R_datasec(3,binning=bin2))
            sl4 = str_to_slice(get_1R_datasec(4,binning=bin2))

            sci_full = CCDData(np.concatenate([np.concatenate([red[0][sl1],red[0][sl2]],axis=1),np.concatenate([red[0][sl3],red[0][sl4]],axis=1)],axis=0),header=header,unit=u.electron)
        else:
            raw = import_file(sci, amp)
            red = [ccdproc.ccd_process(x, oscan=oscan_reg, 
                oscan_model=models.Chebyshev1D(3), trim=trim_sec[k], 
                gain=gains[k]*u.electron/u.adu, 
                readnoise=readnoises[k]*u.electron, master_bias=mbias[k], 
                gain_corrected=True) for k,x in enumerate(raw)]
        flat_image = CCDData(np.copy(mflat),unit=u.electron,meta=mflat.meta,mask=mflat.mask,uncertainty=mflat.uncertainty)
        if amp == '4B':
            sci_full = CCDData(np.concatenate([red[0],
                        np.fliplr(red[1]),
                        np.zeros([np.shape(red[1])[0],111]),
                        red[2],np.fliplr(red[3])],axis=1),
                        header=header,unit=u.electron)
            sci_full = ccdproc.trim_image(sci_full[700:3120,396:3940])
        if amp == '4R':
            sci_full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),np.zeros([np.shape(red[0])[0],200]),red[3],np.fliplr(red[2])],axis=1),header=header,unit=u.electron)
            if window and np.shape(sci_full)!=np.shape(flat_image):
                flat_image.data = mflat.data[625:1875,0:3600]
                try:
                    flat_image.mask = mflat.mask[625:1875,0:3600]
                    flat_image.uncertainty = mflat.uncertainty[625:1875,0:3600]
                except:
                    pass

        # Add SATURATE keyword to sci_full
        sci_full.header['SATURATE'] = saturation(hdr)
        
        processed_data = ccdproc.flat_correct(sci_full, flat_image)
        
        if log: log.info('File proccessed.')
        
        if log: log.info('Calculating 2D background.')
        bkg = Background2D(processed_data, (64, 64), filter_size=(3, 3),
            sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(),
            exclude_percentile=80)
        
        med_background = np.median(bkg.background)
        if log: log.info(f'Median background: {med_background}')

        bkg_basefile = os.path.basename(sci).replace('.fits','_bkg.fits')
        bkg_basefile = bkg_basefile.replace('.gz','')
        bkg_outfile_name = os.path.join(red_path, bkg_basefile)
        fits.writeto(bkg_outfile_name,np.array(bkg.background),overwrite=True)
        if log: log.info(f'Wrote background to: {bkg_outfile_name}')
        
        final = processed_data.subtract(CCDData(bkg.background,unit=u.electron),
            propagate_uncertainties=True, handle_meta='first_found')

        final.header['SATURATE'] -= med_background.value
        final.header['SKYBKG'] = med_background.value
        if log: log.info('Background subtracted and updated saturation.')
        
        if window:
            if log: log.info('Image windowed, adding padding.')
            if amp == '4R':
                final = CCDData(np.concatenate([np.zeros([625,
                    np.shape(final)[1]]),final,
                    np.zeros([625,np.shape(final)[1]])],axis=0),
                    header=final.header,unit=u.electron)

        if log: log.info('Writing WCS to file.')
        ra = final.header['RA']
        dec = final.header['DEC']
        coords = SkyCoord(ra, dec, unit=(u.hour, u.deg), frame='icrs')
        final.header['RADECSYS'] = 'ICRS'
        final.header['CUNIT1'] = 'deg'
        final.header['CUNIT2'] = 'deg'
        final.header['CTYPE1'] = 'RA---TAN'
        final.header['CTYPE2'] = 'DEC--TAN'
        final.header['CRVAL1'] = coords.ra.deg
        final.header['CRVAL2'] = coords.dec.deg
        
        pixsc = pixscale()
        if amp == '4B':
            final.header['CRPIX1'] = 2456
            final.header['CRPIX2'] = 1184
            final.header['CD1_1'] = pixsc/3600*np.sin(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
            final.header['CD1_2'] = -pixsc/3600*np.cos(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
            final.header['CD2_1'] = pixsc/3600*np.cos(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
            final.header['CD2_2'] = -pixsc/3600*np.sin(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
        if amp == '4R':
            final.header['CRPIX1'] = 2400
            final.header['CRPIX2'] = 1200
            final.header['CD1_1'] = pixsc/3600*np.sin(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
            final.header['CD1_2'] = -pixsc/3600*np.cos(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
            final.header['CD2_1'] = pixsc/3600*np.cos(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
            final.header['CD2_2'] = -pixsc/3600*np.sin(np.pi/180.*(int(np.round(final.header['ROTPOSN'],0))+90))
        final.header['DATASEC'] = ('[1:%s,1:%s]'%(np.shape(final)[1],np.shape(final)[0]))
        
        final_basefile = os.path.basename(sci).replace('.gz','')
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

def catalog_zp():
    return 'PS1'

def exptime(hdr):
    return hdr['ELAPTIME']

def stack_method(header):
    instrument = header['INSTRUME']
    side=''
    if instrument == 'LRISBLUE':
        side = 'B'
    if instrument == 'LRIS':
        side = 'R'

    if side=='R':
        return('median')
    elif side=='B':
        return('average')
    else:
        return(None)

def detrend(header):
    instrument = header['INSTRUME']
    side=''
    if instrument == 'LRISBLUE':
        side = 'B'
    if instrument == 'LRIS':
        side = 'R'

    if side=='R':
        return(True)
    elif side=='B':
        return(False)
    else:
        return(False)

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

def overscan_region(amp, binning=1):
    if amp=='4B':
        oscan_reg = '[1:50,1:4096]'
    if amp=='4R':
        oscan_reg = '[1:7,1:2520]'
    if amp=='1R':
        b1=int(2065./binning)
        b2=int(2170./binning)
        b3=np.max([int(1./binning),1])
        b4=int(4248./binning)
        oscan_reg = f'[{b1}:{b2},{b3}:{b4}]'
    return oscan_reg

def fringe_correction(fil):
    return False

def trim(f):
    if 'LB' in f:
        return False
    with fits.open(f) as fo:
        hdr = fo[0].header
    amp = amp_keyword(hdr)
    if amp=='4R':
        return True
    else:
        return False

def number_keyword(hdr):
    return hdr['FRAMENO']

def get_base_science_name(image):

    hdu = fits.open(image, mode='readonly')
    hdr = hdu[0].header

    tkwd = target_keyword()
    target = hdr[tkwd].replace(' ','')
    fil = filter_keyword(hdr)
    amp = amp_keyword(hdr)
    binn = bin_keyword(hdr)
    file_time = time_format(hdr)
    number = number_keyword(hdr)
    number = str(number).zfill(5)

    datestr = Time(file_time, format='mjd').datetime.strftime('ut%y%m%d')

    filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.{number}.fits'

    return(filename)

def trim_section(data):
    return data[625:1875,0:3600]

def edit_raw_headers(rawdir, log=None):

    for file in sorted(glob.glob(os.path.join(rawdir, '*.fits'))):

        hdu = fits.open(file)
        h = hdu[0].header

        if ('ELAPTIME' not in h.keys() and 'DATE-BEG' in h.keys() and
            'DATE-END' in h.keys()):

            t1 = Time(h['DATE-BEG'])
            t2 = Time(h['DATE-END'])
            dt = t2 - t1

            hdu[0].header['ELAPTIME']=dt.to_value('sec')

        if ('ELAPTIME' in h.keys() and float(h['ELAPTIME'])==0):
            hdu[0].header['SLITNAME']=''
            hdu[0].header['OBJECT']='bias'
            hdu[0].header['KOAIMTYP']='bias'

        if (('twi' in h['OBJECT'].lower() and 'flat' in h['OBJECT'].lower()) or
            ('blank' in h['OBJECT'].lower()) or
            ('t_flat' in h['OBJECT'].lower()) or
            ('dome' in h['OBJECT'].lower() and 'flat' in h['OBJECT'].lower())):
            hdu[0].header['KOAIMTYP']='flatlamp'

        if ('KOAIMTYP' not in h.keys()):

            if ('twi' in h['OBJECT'] and 'flat' in h['OBJECT']):
                hdu[0].header['KOAIMTYP']='flatlamp'
            elif ('bias' in h['OBJECT'].lower()):
                hdu[0].header['KOAIMTYP']='bias'
            else:
                hdu[0].header['KOAIMTYP']='object'

        if 'NPIXSAT' not in h.keys():
            npixsat = 0
            for nh in hdu:
                if (('VidInp' in nh.name) or
                    (nh.name=='PRIMARY' and nh.data is not None)):
                    npixsat += len(nh.data[np.where(nh.data>60000)])
            hdu[0].header['NPIXSAT']=npixsat

        if 'INSTREV' in h.keys() and h['INSTREV']=='REDMARK4':
            coord = SkyCoord(h['RABASE'],h['DECBASE'],unit=(u.hour, u.deg))

            hdu[0].header['CRPIX1']=2000.0
            hdu[0].header['CRPIX2']=2000.0

            hdu[0].header['CRVAL1']=coord.ra.degree
            hdu[0].header['CRVAL2']=coord.dec.degree

            # TODO: Check that this is the correct dependence on PA
            pa = h['PA']
            bin1,bin2 = h['BINNING'].split(',')

            hdu[0].header['CD1_2']=-np.sin(pa * np.pi/180.0)*3.75E-05*int(bin1)
            hdu[0].header['CD2_1']=np.sin(pa * np.pi/180.0)*3.75E-05*int(bin1)
            hdu[0].header['CD1_1']=np.cos(pa * np.pi/180.0)*3.75E-05*int(bin2)
            hdu[0].header['CD2_2']=-np.cos(pa * np.pi/180.0)*3.75E-05*int(bin2)

        if ('INSTRUME' in h.keys() and h['INSTRUME']=='LRISBLUE' and
            'RA' in h.keys() and 'DEC' in h.keys()):

            coord = SkyCoord(h['RA'], h['DEC'], unit=(u.hour, u.deg))

            hdu[0].header['CRPIX1']=1500.0
            hdu[0].header['CRPIX2']=1300.0

            hdu[0].header['CRVAL1']=coord.ra.degree
            hdu[0].header['CRVAL2']=coord.dec.degree

            pa = h['ROTPOSN']
            hdu[0].header['CD1_2']=np.sin(pa * np.pi/180.0)*3.75E-05
            hdu[0].header['CD2_1']=-np.sin(pa * np.pi/180.0)*3.75E-05
            hdu[0].header['CD1_1']=np.cos(pa * np.pi/180.0)*3.75E-05
            hdu[0].header['CD2_2']=-np.cos(pa * np.pi/180.0)*3.75E-05

            for i in np.arange(len(hdu)):
                if i==0: continue

                hdu[i].header['CUNIT1'] = 'deg'
                hdu[i].header['CUNIT2'] = 'deg'

        hdu.writeto(file, overwrite=True, output_verify='silentfix')

def edit_stack_headers(stack):

    good_keys = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND',
        'DATE-OBS','RA','DEC','TARGRA','TARGDEC','EQUINOX','HA','AIRMASS',
        'PARANG','EL','AZ','LST','TELESCOP','PA','DICHNAME','OBSMODE',
        'CTYPE1','CTYPE2','CUNIT1','CUNIT2','RADECSYS','CD1_1','CD1_2',
        'CD2_1','CD2_2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','SEMESTER',
        'PROGID','PROGPI','PROGINST','OA','SATURATE','SKYBKG','BUNIT',
        'MJD-OBS','EXPTOT','EXPTIME','GAIN','RDNOISE','NFILES','FILTER',
        'OBSTYPE','XTENSION','EXTNAME']

    if 'TARGNAME' in stack[0].header.keys() and 'OBJECT' in stack[0].header.keys():
        stack[0].header['OBJECT'] = stack[0].header['TARGNAME']

    if 'EXPTOT' in stack[0].header.keys() and 'EXPTIME' in stack[0].header.keys():
        stack[0].header['OBJECT'] = stack[0].header['TARGNAME']

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
