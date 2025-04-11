# Parameter file for Magellan/IMACS

__version__ = "1.0" # Last edited 02/13/2025

import os
import ccdproc
import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData

# Internal dependency
from . import instrument

class FOURSTAR(instrument.Instrument):

    def __init__(self):

        self.version = __version__

        self.name = 'FourStar'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 1
        self.wcs_extension = 1

        # Detector specific characteristics
        self.pixscale = 0.159
        self.saturation = 46000.0

        self.min_exptime = 0.1

        # Run dark/bias/flat calibration?
        self.dark = True
        self.bias = False
        self.flat = True

        # Parameters for handling calibration files
        # Run rejection on possible CR pixels in bias
        self.cr_bias = True 

        # Extend header to first file extension for purposes of sort_files
        self.extend_header = False

        # How to combine images during stacking
        self.stack_method = 'median'

        self.wavelength = 'NIR'

        # https://www.lco.cl/wp-content/uploads/telescopes/magellan/instruments/fourstar/FourStar-Handout.pdf
        self.gain = [2.65, 2.59, 2.51, 2.49]
        self.rdnoise = [25.5, 22.1, 20.5, 18.9]

        self.datasec = ['[1055:3024,217:3911]','[1067:3033,217:3911]']
        self.biassec = ['[1:1054,1:4112]','[3034:4096,1:4112]']

        # Keywords for selecting files from Sort_files object
        # This allows a single file type to be used for multiple purposes (e.g., for
        # generating a flat-field image from science exposures)
        self.filetype_keywords = {'SCIENCE':'SCIENCE', 'FLAT':'[SCIENCE,FLAT]', 
            'DARK':'DARK','BIAS':'BIAS'}

        # Header keywords
        self.target_keyword = 'OBJECT'
        self.exptime_keyword = 'EXPTIME'
        self.filter_keyword = 'FILTER'
        self.mjd_keyword = 'MJD'
        self.bin_keyword = 'CCDSUM'
        self.amp_keyword = 'NAMPS'

        # File sorting keywords
        self.science_keywords = ['OBSMODE','APTYPE']
        self.science_values = ['imaging','open']
        self.flat_keywords = []
        self.flat_values = []
        self.bias_keywords = []
        self.bias_values = []
        self.dark_keywords = ['OBJECT']
        self.dark_values = ['ark']
        self.spec_keywords = ['OBSMODE','APERTURE']
        self.spec_values = ['spectral','open']
        self.bad_keywords = ['MOSID']
        self.bad_values = ['closed']

        self.detrend = True
        self.catalog_zp = '2MASS'

        self.out_size = 2500

    # Get a unique image number that can be derived only from the file header
    def get_number(self, hdr):
        datestr = hdr['DATE-OBS']
        elap = Time(datestr)-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))
        return(elap)

    def raw_format(self, proc):
        return('*.fits')

    def get_rdnoise(self, hdr):
        return(hdr['RDNOISE'])

    def get_gain(self, hdr):
        return(hdr['GAIN'])

    def get_time(self, hdr):
        return Time(hdr['DATE-OBS']).mjd

    # Not currently used - need better implementation of slope measurement to
    # account for jump detection, bias, etc.
    def measure_slope(self, hdu):
        header = hdu[self.raw_header_ext].header

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
                all_times.append(self.get_exptime(hdu[ext].header))
                num_reads+=1

        all_data = np.array(all_data)
        all_times = np.array(all_times)

        # Run a linear regression on 3D array all_data with respect to all_times
        A = np.vstack([all_times, np.ones(len(all_times))]).T
        data = np.linalg.lstsq(A, all_data.reshape(num_reads, -1), rcond=None)

        # Get slope data
        slope = data[0][0].reshape(*data_shape)
        slope = slope.astype(np.float32)
        slope = slope * np.max(all_times)

        return(slope)

    # Takes filename and outputs CCDData object with raw image in units of e-
    def import_image(self, filename, amp, log=None):
        hdu = fits.open(filename)
        header = hdu[self.raw_header_ext].header
        slope = hdu[1].data
        for key in hdu[1].header:
            try:
                header[key] = hdu[1].header[key]
            except ValueError:
                continue

        for key in ['BZERO','BSCALE','BITPIX']:
            if key in header.keys():
                del header[key]

        # Multiply by GAIN and apply RDNOISE
        raw = CCDData(slope, header=header, unit=u.adu)
        red = ccdproc.ccd_process(raw, 
            gain=self.get_gain(header)*u.electron/u.adu,
            readnoise=self.get_rdnoise(header)*u.electron)

        # Overscan correct
        red = ccdproc.subtract_overscan(red, overscan=red[0:4,:], 
            overscan_axis=0)

        # Trim image
        red = ccdproc.ccd_process(red, trim=header['DATASEC'])

        # Make sure to take out *SEC keywords
        for key in ['DATASEC','CCDSEC','DETSEC','TRIMSEC','DETSIZE']:
            if key in red.header.keys():
                del red.header[key]

        # Apply saturation
        saturate = self.saturation * self.get_gain(header)
        red.header['SATURATE'] = saturate

        # Re-apply data to mask
        red.mask = red.data > saturate

        return(red)
