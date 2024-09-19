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

__version__ = 1.0 # last edited 09-05-2024

# Generic instrument class for POTPyRI

class Instrument(object):

    def __init__(self):

        # Version
        self.version = __version__

        # Intrument name
        self.name = ''

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0 
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = None
        self.saturation = 65535

        self.min_exptime = 1.0

        # Run dark/bias/flat calibration?
        self.dark = False
        self.bias = False
        self.flat = False

        # Parameters for handling calibration files
        # Run rejection on possible CR pixels in bias
        self.cr_bias = True 

        # How to combine images during stacking
        self.stack_method = 'median'

        self.wavelength = 'OPT'

        self.gain = [1.02, 1.08]
        self.rdnoise = [4.0, 4.0]

        self.datasec = ['[1055:3024,217:3911]','[1067:3033,217:3911]']
        self.biassec = ['[1:1054,1:4112]','[3034:4096,1:4112]']

        # Keywords for selecting files from Sort_files object
        # This allows a single file type to be used for multiple purposes (e.g., for
        # generating a flat-field image from science exposures)
        self.filetype_keywords = {'SCIENCE':'SCIENCE', 'FLAT':'FLAT', 
            'DARK':'DARK','BIAS':'BIAS'}

        # Header keywords
        self.target_keyword = 'OBJECT'
        self.exptime_keyword = 'EXPTIME'
        self.filter_keyword = 'FILTER'
        self.mjd_keyword = 'MJD'
        self.bin_keyword = 'CCDSUM'
        self.amp_keyword = '2'

        # File sorting keywords
        self.science_keywords = ['MASK','SCRN']
        self.science_values = ['imaging','stowed']
        self.flat_keywords = ['MASK','SCRN']
        self.flat_values = ['imaging','deployed']
        self.bias_keywords = []
        self.bias_values = []
        self.dark_keywords = []
        self.dark_values = []
        self.spec_keywords = []
        self.spec_values = []
        self.bad_keywords = []
        self.bad_values = []

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 5000

    # Use these if a single value is needed for gain, rdnoise, etc.
    def get_rdnoise(self, hdr):
        return(self.rdnoise)

    def get_gain(self, hdr):
        return(self.gain)

    # Get specific header keywords from file
    def get_target(self, hdr):
        return(hdr[self.target_keyword].replace(' ',''))

    def get_filter(self, hdr):
        filt = hdr[self.filter_keyword]
        filt = filt.replace(' ','')
        filt = filt.split('_')[0]
        return(filt)

    def get_exptime(self, hdr):
        return(hdr[self.exptime_keyword])

    def get_ampl(self, hdr):
        if self.amp_keyword in hdr.keys():
            return(str(hdr[self.amp_keyword]))
        else:
            return(str(self.amp_keyword))

    def get_binning(self, hdr):
        if self.bin_keyword in hdr.keys():
            binn = str(hdr[self.bin_keyword]).replace(' ','')
            binn = binn.replace(',','')
            return(binn)
        else:
            return(self.bin_keyword)

    def get_time(self, hdr):
        return(float(hdr[self.mjd_keyword]))

    def get_static_mask(self, paths):

        mask_file = os.path.join(paths['code'], 'data', 'staticmasks', 
            f'{self.name.lower()}.staticmask.fits.fz')

        if os.path.exists(mask_file):
            return [mask_file]
        else:
            return [None]

    def get_catalog(self, hdr):
        return(self.catalog_zp)

    def format_datasec(self, sec_string, binning=1):
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

    def raw_format(self, proc):
        if proc:
            return 'sci_img_*.fits'
        else:
            return 'sci_img*[!proc].fits'

    def get_stk_name(self, hdr, red_path):
        
        target = self.get_target(hdr)
        fil = self.get_filter(hdr)
        amp = self.get_ampl(hdr)
        binn = self.get_binning(hdr)
        file_time = self.get_time(hdr)
        datestr = Time(file_time, format='mjd').datetime.strftime('ut%y%m%d')

        filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.stk.fits'
        fullfilename = os.path.join(red_path, filename)

        return(fullfilename)

    def get_sci_name(self, hdr, red_path):

        target = self.get_target(hdr)
        fil = self.get_filter(hdr)
        amp = self.get_ampl(hdr)
        binn = self.get_binning(hdr)
        file_time = self.get_time(hdr)
        number = self.get_number(hdr)
        datestr = Time(file_time, format='mjd').datetime.strftime('ut%y%m%d')

        filename = f'{target}.{fil}.{datestr}.{amp}.{binn}.{number}.fits'
        fullfilename = os.path.join(red_path, filename)

        return(fullfilename)

    def get_bkg_name(self, hdr, red_path):

        sci_name = self.get_sci_name(hdr, red_path)
        bkg_filename = sci_name.replace('.fits','_bkg.fits')

        return(bkg_filename)

    def get_mbias_name(self, red_path, amp, binn):
        return(os.path.join(red_path, f'mbias_{amp}_{binn}.fits'))

    def get_mdark_name(self, red_path, amp, binn):
        return(os.path.join(red_path, f'mdark_{amp}_{binn}.fits'))

    def get_mflat_name(self, red_path, fil, amp, binn):
        return(os.path.join(red_path, f'mflat_{fil}_{amp}_{binn}.fits'))

    def load_bias(self, red_path, amp, binn):
        bias = self.get_mbias_name(red_path, amp, binn)
        mbias = CCDData.read(bias)
        return mbias

    def load_dark(self, red_path, amp, binn):
        dark = self.get_mdark_name(red_path, amp, binn)
        mdark = CCDData.read(dark)
        return mdark

    def load_flat(self, red_path, fil, amp, binn):
        flat = self.get_mflat_name(red_path, fil, amp, binn)
        mflat = CCDData.read(flat)
        return mflat

    def import_image(self, filename, amp, log=None):
        filename = os.path.abspath(filename)
        if log: log.info(f'Loading file: {filename}')

        with fits.open(flat) as hdr:
            header = hdr[1].header

        raw = [CCDData.read(flat, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=self.biassec[k], 
               oscan_model=models.Chebyshev1D(3), 
               trim=self.datasec[k], gain=self.gain[k]*u.electron/u.adu, 
               readnoise=rdnoise(hdr)*u.electron, 
               gain_corrected=True) for k,x in enumerate(raw)]

        frame_full = CCDData(np.concatenate((red[0], 
            np.empty((red[0].shape[0], 794)), red[1]), axis=1), 
            header=header,unit=u.electron)

        frame_full.header['SATURATE'] = saturation(header)

        return(frame_full)

    def create_bias(self, bias_list, amp, binn, red_path, staticmask=None,
        log=None, **kwargs):

        staticmask = self.load_satmask(staticmask)
        
        if log:
            log.info(f'Processing bias files with {amp} amps and {binn} binning.')
            log.info(f'{len(bias_list)} files found.')

        biases = []
        for i, bias in enumerate(bias_list):
            if log: log.info(f'Importing {bias}')
            bias_full = self.import_image(bias, amp, log=log)
            if staticmask is not None:
                if log: log.info('Applying static mask')
                bias_full.data[staticmask]=np.nan
            if self.cr_bias:
                # Add a CR mask to the bias image
                mean, median, stddev = sigma_clipped_stats(bias_full.data)
                mask = bias_full.data > median + 3 * stddev

                if bias_full.mask is not None:
                    bias_full.mask = bias_full.mask | mask
                else:
                    bias_full.mask = mask

                bias_full.data[mask]=median

            biases.append(bias_full)

        mbias = ccdproc.combine(biases, method='median', sigma_clip=True,
            clip_extrema=True)
        mbias.data[np.isnan(mbias.data)]=np.nanmedian(mbias.data)

        if log: log.info(f'Made bias for amp: {amp}, bin: {binn}.')
        mbias.header['VER'] = (__version__, 
            'Version of telescope parameter file used.')

        bias_filename = self.get_mbias_name(red_path, amp, binn)
        mbias.write(bias_filename, overwrite=True)
        log.info(f'Master bias written to {bias_filename}')
        
        return

    def create_dark(self, dark_list, amp, binn, red_path, mbias=None, 
        staticmask=None, log=None):

        staticmask = self.load_satmask(staticmask)

        if log:
            log.info(f'Processing dark files with {amp} amps and {binn} binning.')
            log.info(f'{len(dark_list)} files found.')
    
        darks = []
        exptimes = []
        for dark in dark_list:
            if log: log.info(f'Importing {dark}')
            dark_full = self.import_image(dark, amp, log=log)
            if staticmask is not None:
                if log: log.info('Applying static mask')
                dark_full.data[staticmask]=np.nan
            exptimes.append(self.get_exptime(dark_full.header))

            if mbias is not None:
                if log: log.info('Subtracting bias')
                dark_full = ccdproc.subtract_bias(dark_full, mbias)

            darks.append(dark_full)
        
        if log: log.info('Creating master dark.')
        mdark = ccdproc.combine(darks, method='median', 
            scale=1./np.array(exptimes), weights=np.array(exptimes),
            sigma_clip=True, clip_extrema=True)

        # Rescale to electrons by average of all exposure times
        avg_exptime = np.mean(exptimes)
        mdark.data *= avg_exptime
        mdark.header[self.exptime_keyword]=(avg_exptime, 
            'Exposure time of master dark (sec).')
        
        # Update other header variables
        mdark.header['VER'] = (__version__, 
            'Version of telescope parameter file used.')    
        for i,im in enumerate(dark_list):
            mdark.header[f'FILE{i+1}']=os.path.basename(im)
        
        darkname = self.get_mdark_name(red_path, amp, binn)
        if log: log.info(f'Writing master dark to {darkname}')
        mdark.write(darkname, overwrite=True)
        
        return

    def create_flat(self, flat_list, fil, amp, binn, red_path, mbias=None, 
        mdark=None, is_science=False, staticmask=None, log=None, **kwargs):

        staticmask = self.load_satmask(staticmask)

        if log:
            log.info(f'Processing files for filter: {fil}')
            log.info(f'{len(flat_list)} files found.')

        scale = []
        flats = []
        for i, flat in enumerate(flat_list):
            if log: log.info(f'Importing {flat}')
            flat_full = self.import_image(flat, amp, log=log)
            if staticmask is not None:
                if log: log.info('Applying static mask')
                flat_full.data[staticmask]=np.nan

            if mbias is not None:
                if log: log.info('Subtracting bias')
                flat_full = ccdproc.subtract_bias(flat_full, mbias)

            if mdark is not None:
                if log: log.info('Subtracting dark')
                flat_full = ccdproc.subtract_dark(flat_full, mdark, 
                    exposure_time=self.exptime_keyword, exposure_unit=u.second)

            # Mask flat_full if the image is a science frame
            if is_science:
                mean, median, stddev = sigma_clipped_stats(flat_full.data)
                mask = flat_full.data > 5*stddev + median
                flat_full.mask = mask.astype(np.uint8)

            exptime = self.get_exptime(flat_full.header)
            log.info(f'Exposure time of image is {exptime}')
            
            norm = 1./np.nanmedian(flat_full) #check for binning
            log.info(f'Flat normalization: {norm}')
            scale.append(norm)
            flats.append(flat_full)
        
        mflat = ccdproc.combine(flats, method='median', scale=scale, 
            sigma_clip=True, clip_extrema=True)

        # Mask flat
        mflat.data[np.isinf(mflat.data)]=np.nan
        mflat.data[mflat.data==0.0]=np.nan
        mflat.data[mflat.data<0.0]=np.nan
        mean, median, stddev = sigma_clipped_stats(mflat.data)
        mask = mflat.data > median + 10 * stddev
        mflat.data[mask]=np.nan

        # Remove all nan values from flat
        mask = np.isnan(mflat.data)
        mflat.mask = mask
        mflat.data[mask]=1.0
        
        if log: 
            log.info(f'Made flat for filter: {fil}, amp: {amp}, bin: {binn}.')
        
        mflat.header['VER'] = (__version__, 
            'Version of telescope parameter file used.')

        flat_filename = self.get_mflat_name(red_path, fil, amp, binn)
        mflat.write(flat_filename, overwrite=True)
        log.info(f'Master flat written to {flat_filename}')
        
        return

    def load_satmask(self, satmask_filename):

        if (satmask_filename is not None and 
            len(satmask_filename)>0 and
            satmask_filename[0] is not None and
            os.path.exists(satmask_filename[0])):
            hdu = fits.open(satmask_filename[0])
            satmask = hdu[1].data.astype(bool)
        else:
            satmask = None

        return(satmask)

    def expand_mask(self, input_data, input_mask=None):

        # Get image statistics
        mean, median, stddev = sigma_clipped_stats(input_data.data)
        
        # Add to input mask
        add_mask = np.isnan(input_data.data)
        add_mask = add_mask  | np.isinf(input_data.data)
        add_mask = add_mask  | (input_data.data==0.0) 

        # Expand with input mask if it was provided
        if input_mask is None:
            input_mask = add_mask
        else:
            input_mask = input_mask | add_mask

        return(input_mask)

    def process_science(self, sci_list, fil, amp, binn, red_path, mbias=None,
        mflat=None, mdark=None, proc=None, staticmask=None, skip_skysub=False, 
        log=None):

        staticmask = self.load_satmask(staticmask)

        processed = []
        for sci in sorted(sci_list):
            sci_full = self.import_image(sci, amp, log=log)

            if staticmask is None:
                staticmask = np.zeros(sci_full.data.shape).astype(bool)

            # Subtract bias
            if mbias is not None:
                if log: log.info('Subtracting bias')
                sci_full = ccdproc.subtract_bias(sci_full, mbias)

            # Subtract dark
            if mdark is not None:
                if log: log.info('Subtracting dark')
                sci_full = ccdproc.subtract_dark(sci_full, mdark, 
                    exposure_time=self.exptime_keyword, exposure_unit=u.second)

            # Perform flat-fielding
            if mflat is not None:
                if log: log.info('Flattening image')
                sci_full = ccdproc.flat_correct(sci_full, mflat)
                staticmask = staticmask | (mflat.data < 0.5)

            # Expand input mask
            processed_data = sci_full
            processed_data.mask = self.expand_mask(processed_data, 
                input_mask=staticmask)
            processed_data.data[processed_data.mask]=np.nan

            if not skip_skysub:
                if log: log.info('Calculating 2D background.')
                bkg = Background2D(processed_data, (64, 64), filter_size=(3, 3),
                    sigma_clip=SigmaClip(sigma=3), exclude_percentile=80,
                    bkg_estimator=MeanBackground(), mask=processed_data.mask, 
                    fill_value=np.nanmedian(processed_data.data))
                
                med_background = np.nanmedian(bkg.background)
                if log: log.info(f'Median background: {med_background}')

                bkg_filename = self.get_bkg_name(processed_data.header, red_path)
                if log: log.info(f'Writing background file: {bkg_filename}')
                bkg_hdu = fits.PrimaryHDU(bkg.data)
                bkg_hdu.header = processed_data.header
                bkg_hdu.writeto(bkg_filename, overwrite=True,
                    output_verify='silentfix')

                final = processed_data.subtract(CCDData(bkg.background,
                    unit=u.electron), propagate_uncertainties=True, 
                    handle_meta='first_found')
                final.header['SATURATE'] -= med_background.value
                final.header['SKYBKG'] = med_background.value
                if log: log.info('Updating background and saturation values.')
            else:
                final = processed_data
                final.header['SKYBKG'] = 0.0

            # Convert data format to float32
            final.data = final.data.astype(np.float32)
            final.header['BITPIX']=-32

            # Delete empty header keys
            for key in list(final.header.keys()):
                if not key.strip():
                    if key in final.header.keys():
                        del final.header[key]

            final_filename = self.get_sci_name(final.header, red_path)
            
            # Edit additional header values
            final.header['FILENAME']=final_filename
            final.header['FILTER']=fil
            final.header['AMPS']=amp 
            final.header['BINNING']=binn
            final.header['ORGFILE']=sci
            final.header['EXTNAME']='SCI'

            if log: log.info(f'Writing final file: {final_filename}')
            final.write(final_filename, overwrite=True)
            processed.append(final)

        return processed

    def edit_raw_headers(self, files, log=None):

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
                if ('CTYPE1' in hdr.keys() and 
                    hdr['CTYPE1']=='RA---TAN' and 
                    (hdr['CRVAL1']<0.0 or hdr['CRVAL1']>360.0)):
                    coord = SkyCoord(hdr['LST'], '31:41:20.04', 
                        unit=(u.hour, u.deg))
                    hdu[i].header['CRVAL1']=coord.ra.degree

                if ('CTYPE2' in hdr.keys() and 
                    hdr['CTYPE2']=='DEC--TAN' and 
                    (hdr['CRVAL2']<-90.0 or hdr['CRVAL2']>90.0)):
                    hdu[i].header['CRVAL2']=31.6889


            hdu.writeto(file, overwrite=True, output_verify='silentfix')


    def edit_stack_headers(self, stack):

        good_keys = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND',
            'OBSERVAT','TELESCOP','ORIGIN','INSTRUME','EXTNAME','FILTER',
            'DATE-OBS','RA','DEC','POSANG','AZ','EL','ROT','PA','MJD',
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
