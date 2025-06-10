# Basic instrument class

__version__ = "1.2" # Last edited 09/29/2024

import os
import astropy
import datetime
import copy
import ccdproc
import numpy as np

from photutils.background import Background2D
from photutils.background import MeanBackground

import astropy.units as u
from astropy.io import fits
from astropy.modeling import models
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.time import Time

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

        # Extend header to first file extension for purposes of sort_files
        self.extend_header = False

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

        # Default final image size for binning of 1x1
        self.out_size = 5000

    # Determine whether sky subtraction is needed for the current reduction
    def needs_sky_subtraction(self, filt):
        # Default is to do sky subtraction for all NIR cameras
        if self.wavelength=='NIR':
            return(True)
        else:
            return(False)

    # Get pixel scale for imagers with varying focal ratios (e.g., IMACS)
    def get_pixscale(self):
        return(self.pixscale)

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
            binn = str(hdr[self.bin_keyword]).lower().replace(' ','')
            binn = binn.replace('x','')
            binn = binn.replace(',','')
            return(binn)
        else:
            return(self.bin_keyword)

    def get_out_size(self, hdr):
        binn = int(str(self.get_binning(hdr))[0])
        out_size = int(self.out_size/binn)
        return(out_size)

    def get_time(self, hdr):
        return(float(hdr[self.mjd_keyword]))

    def get_instrument_name(self, hdr):
        return(self.name.lower())

    def get_staticmask_filename(self, hdr, paths):

        instname = self.get_instrument_name(hdr)
        binn = self.get_binning(hdr)

        mask_file = os.path.join(paths['code'], 'data', 'staticmasks', 
            f'{instname}.{binn}.staticmask.fits.fz')

        if os.path.exists(mask_file):
            return([mask_file])
        else:
            return([None])

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
            return('sci_img_*.fits')
        else:
            return('sci_img*[!proc].fits')

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

    def get_mbias_name(self, paths, amp, binn):
        red_path = paths['cal']
        return(os.path.join(red_path, f'mbias_{amp}_{binn}.fits'))

    def get_mdark_name(self, paths, amp, binn):
        red_path = paths['cal']
        return(os.path.join(red_path, f'mdark_{amp}_{binn}.fits'))

    def get_mflat_name(self, paths, fil, amp, binn):
        red_path = paths['cal']
        return(os.path.join(red_path, f'mflat_{fil}_{amp}_{binn}.fits'))

    def get_msky_name(self, paths, fil, amp, binn):
        red_path = paths['cal']
        return(os.path.join(red_path, f'msky_{fil}_{amp}_{binn}.fits'))

    def load_bias(self, paths, amp, binn):
        bias = self.get_mbias_name(paths, amp, binn)
        if os.path.exists(bias):
            mbias = CCDData.read(bias)
        elif os.path.exists(bias+'.fz'):
            hdu = fits.open(bias+'.fz')
            mbias = CCDData(hdu[1].data, header=hdu[1].header, unit=u.electron)
        else:
            raise Exception(f'Could not find bias: {bias}')
        return(mbias)

    def load_dark(self, paths, amp, binn):
        dark = self.get_mdark_name(paths, amp, binn)
        if os.path.exists(dark):
            mdark = CCDData.read(dark)
        elif os.path.exists(dark+'.fz'):
            hdu = fits.open(dark+'.fz')
            mdark = CCDData(hdu[1].data, header=hdu[1].header, unit=u.electron)
        else:
            raise Exception(f'Could not find dark: {dark}')
        return(mdark)

    def load_flat(self, paths, fil, amp, binn):
        flat = self.get_mflat_name(paths, fil, amp, binn)
        if os.path.exists(flat):
            mflat = CCDData.read(flat)
        elif os.path.exists(flat+'.fz'):
            hdu = fits.open(flat+'.fz')
            mflat = CCDData(hdu[1].data, header=hdu[1].header, unit=u.electron)
        else:
            raise Exception(f'Could not find flat: {flat}')
        return(mflat)

    def load_sky(self, paths, fil, amp, binn):
        sky = self.get_msky_name(paths, fil, amp, binn)
        if os.path.exists(sky):
            msky = CCDData.read(sky)
        elif os.path.exists(sky+'.fz'):
            hdu = fits.open(sky+'.fz')
            msky = CCDData(hdu[1].data, header=hdu[1].header, unit=u.electron)
        else:
            raise Exception(f'Could not find sky frame: {sky}')
        return(msky)

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

    def import_sci_image(self, filename, log=None):
        filename = os.path.abspath(filename)
        if log: log.info(f'Loading file: {filename}')

        frame_full = CCDData.read(filename)

        return(frame_full)

    def create_bias(self, bias_list, amp, binn, paths,
        log=None, **kwargs):
        
        if log:
            log.info(f'Processing bias files with {amp} amps and {binn} binning.')
            log.info(f'{len(bias_list)} files found.')

        biases = []
        for i, bias in enumerate(bias_list):
            if log: log.info(f'Importing {bias}')
            bias_full = self.import_image(bias, amp, log=log)

            # Load static mask for this specific file
            staticmask = self.load_staticmask(bias_full.header, paths)

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

        bias_filename = self.get_mbias_name(paths, amp, binn)
        mbias.write(bias_filename, overwrite=True, output_verify='silentfix')
        log.info(f'Master bias written to {bias_filename}')
        
        return

    def create_dark(self, dark_list, amp, binn, paths, mbias=None, log=None):

        if log:
            log.info(f'Processing dark files with {amp} amps and {binn} binning.')
            log.info(f'{len(dark_list)} files found.')
    
        darks = []
        exptimes = []
        for dark in dark_list:
            if log: log.info(f'Importing {dark}')
            dark_full = self.import_image(dark, amp, log=log)

            # Load static mask for this specific file
            staticmask = self.load_staticmask(dark_full.header, paths)

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
            scale=1./np.array(exptimes), sigma_clip=True, clip_extrema=True)

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
        
        darkname = self.get_mdark_name(paths, amp, binn)
        if log: log.info(f'Writing master dark to {darkname}')
        mdark.write(darkname, overwrite=True, output_verify='silentfix')
        
        return

    def create_flat(self, flat_list, fil, amp, binn, paths, mbias=None, 
        mdark=None, is_science=False, log=None, **kwargs):

        if log:
            log.info(f'Processing files for filter: {fil}')
            log.info(f'{len(flat_list)} files found.')

        scale = []
        flats = []
        for i, flat in enumerate(flat_list):
            if log: log.info(f'Importing {flat}')
            flat_full = self.import_image(flat, amp, log=log)
            
            # Load static mask for this specific file
            staticmask = self.load_staticmask(flat_full.header, paths)

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
            if log: log.info(f'Exposure time of image is {exptime}')
            
            norm = 1./np.nanmedian(flat_full.data)
            # Vet the flat normalization - it should not be negative
            if norm > 0.:
                log.info(f'Flat normalization: {norm}')
            else:
                # Skip this file
                log.error(f'Flat normalization: {norm}')
                continue
                
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

        flat_filename = self.get_mflat_name(paths, fil, amp, binn)
        mflat.write(flat_filename, overwrite=True, output_verify='silentfix')
        log.info(f'Master flat written to {flat_filename}')
        
        return

    def create_sky(self, sky_list, fil, amp, binn, paths, log=None, **kwargs):

        if log:
            log.info(f'Processing files for filter: {fil}')
            log.info(f'{len(sky_list)} files found.')

        scale = []
        skys = []
        for i, sky in enumerate(sky_list):
            if log: log.info(f'Importing {sky}')
            sky_full = self.import_sci_image(sky, log=log)

            mean, med, stddev = sigma_clipped_stats(sky_full.data)

            # Mask outliers
            mask = sky_full.data > med + 5 * stddev
            sky_full.data[mask]=np.nan

            # Normalize by median sky background
            mean, med, stddev = sigma_clipped_stats(sky_full.data)
            norm = 1./med
            sky_full = sky_full.multiply(norm)
            
            # Vet the sky normalization - it should not be negative
            if norm > 0.:
                log.info(f'Sky normalization: {norm}')
            else:
                # Skip this file
                log.error(f'Sky normalization: {norm}')
                continue

            sky_full.mask[np.isnan(sky_full.data)]=True
            sky_full.data[np.isnan(sky_full.data)]=1.0
                
            skys.append(sky_full)

        msky = ccdproc.combine(skys, method='median', sigma_clip=True, 
            clip_extrema=True)

        # Mask sky image
        msky.data[np.isinf(msky.data)]=1.0
        msky.data[msky.data==0.0]=1.0
        mean, median, stddev = sigma_clipped_stats(msky.data)
        mask = msky.data > median + 10 * stddev
        msky.data[mask]=1.0
        
        if log: 
            log.info(f'Made sky for filter: {fil}, amp: {amp}, bin: {binn}.')
        
        msky.header['VER'] = (__version__, 
            'Version of telescope parameter file used.')

        sky_filename = self.get_msky_name(paths, fil, amp, binn)
        msky.write(sky_filename, overwrite=True, output_verify='silentfix')
        log.info(f'Master sky written to {sky_filename}')
        
        return

    def load_staticmask(self, hdr, paths):

        satmask_filename = self.get_staticmask_filename(hdr, paths)

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

    def process_science(self, sci_list, fil, amp, binn, paths, mbias=None,
        mflat=None, mdark=None, skip_skysub=False, log=None):

        processed = []
        processed_names = []
        for sci in sorted(sci_list):
            if log: log.info(f'Importing {sci}')
            sci_full = self.import_image(sci, amp, log=log)

            # Load static mask for this specific file
            staticmask = self.load_staticmask(sci_full.header, paths)

            if staticmask is None:
                staticmask = np.zeros(sci_full.data.shape).astype(bool)
            else:
                if log: log.info('Applying static mask')
                sci_full.data[staticmask]=np.nan

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
                # Expand mask wherever sensitivity of the detector is low
                # e.g., due to vignetting
                staticmask = staticmask | (mflat.data < 0.5)

            # Expand input mask
            processed_data = sci_full
            processed_data.mask = self.expand_mask(processed_data, 
                input_mask=staticmask)
            processed_data.data[processed_data.mask]=np.nan

            if log: log.info(f'Wavelength is {self.wavelength}')
            if not skip_skysub and not self.needs_sky_subtraction(fil):
                if log: log.info('Calculating 2D background.')
                bkg = Background2D(processed_data, (64, 64), filter_size=(3, 3),
                    sigma_clip=SigmaClip(sigma=3), exclude_percentile=80,
                    bkg_estimator=MeanBackground(), mask=processed_data.mask, 
                    fill_value=np.nanmedian(processed_data.data))
                
                med_background = np.nanmedian(bkg.background)
                if log: log.info(f'Median background: {med_background}')

                bkg_filename = self.get_bkg_name(processed_data.header, paths['work'])
                if log: log.info(f'Writing background file: {bkg_filename}')
                bkg_hdu = fits.PrimaryHDU(bkg.background.value)
                bkg_hdu.header = processed_data.header
                bkg_hdu.writeto(bkg_filename, overwrite=True,
                    output_verify='silentfix')

                final_img = processed_data.subtract(CCDData(bkg.background,
                    unit=u.electron), propagate_uncertainties=True, 
                    handle_meta='first_found')
                final_img.header['SATURATE'] -= med_background.value
                final_img.header['SKYBKG'] = med_background.value
                if log: log.info('Updating background and saturation values.')
            else:
                final_img = processed_data
                final_img.header['SKYBKG'] = 0.0

            # Apply final masking based on excessively negative values
            mean, median, stddev = sigma_clipped_stats(final_img.data)
            mask = final_img.data < median - 10 * stddev
            final_img.data[mask] = np.nan
            final_img.mask = final_img.mask | mask

            # Convert data format to float32
            final_img.data = final_img.data.astype(np.float32)
            final_img.header['BITPIX']=-32

            # Get read noise for uncertainty image
            rdnoise = np.mean(self.get_rdnoise(final_img.header))

            # Add RMS noise to account for sky background
            data = final_img.data
            rms = 0.5 * (
                np.percentile(data[~np.isnan(data)], 84.13)
                - np.percentile(data[~np.isnan(data)], 15.86)
            )
            rdnoise = rdnoise + rms**2

            if log: log.info('Creating uncertainty image')
            final_img = ccdproc.create_deviation(final_img,
                readnoise=rdnoise*u.electron,
                disregard_nan=True)

            # Delete empty header keys
            for key in list(final_img.header.keys()):
                if not key.strip():
                    if key in final_img.header.keys():
                        del final_img.header[key]

            final_filename = self.get_sci_name(final_img.header, paths['work'])
            
            # Edit additional header values
            final_img.header['FILENAME']=final_filename
            final_img.header['FILTER']=fil
            final_img.header['AMPS']=amp 
            final_img.header['BINNING']=binn
            final_img.header['ORGFILE']=sci
            final_img.header['EXTNAME']='SCI'

            # Get rid of header values if they exist
            for key in ['DATASEC','BIASSEC','CCDSEC']:
                if key in final_img.header.keys():
                    del final_img.header[key]

            if log: log.info(f'Writing final file: {final_filename}')
            final_img.write(final_filename, overwrite=True, output_verify='silentfix')
            
            processed.append(final_img)
            processed_names.append(final_filename)

        # Create sky image and subtract from every science frame
        if self.needs_sky_subtraction(fil):
            sky_frame = self.get_msky_name(paths, fil, amp, binn)
            if not os.path.exists(sky_frame):
                self.create_sky(processed_names, fil, amp, binn, paths, 
                    log=log)
            
            sky_frame = self.load_sky(paths, fil, amp, binn)

            for i,frame in enumerate(processed):

                mean, med, stddev = sigma_clipped_stats(frame.data)
                frame_sky = sky_frame.multiply(med, 
                    propagate_uncertainties=True, handle_meta='first_found')

                processed[i] = frame.subtract(frame_sky, 
                    propagate_uncertainties=True, handle_meta='first_found')
                processed[i].header['SKYBKG']=med
                processed[i].header['SATURATE']-=med

                # Rewrite file
                final_filename = self.get_sci_name(processed[i].header, 
                    paths['work'])

                if log: log.info(f'Writing final file: {final_filename}')
                processed[i].write(final_filename, overwrite=True, output_verify='silentfix')

        return(processed)
