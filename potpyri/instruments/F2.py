# Parameter file for F2/Gemini-S

__version__ = "1.0" # Last edited 11/23/2024

import os
import ccdproc
import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData

# Internal dependency
from . import instrument

class F2(instrument.Instrument):

    def __init__(self):

        self.version = __version__

        self.name = 'F2'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.1792
        self.saturation = 65535.0

        self.min_exptime = 1.0

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

        self.gain = None
        self.rdnoise = None

        self.datasec = None
        self.biassec = None

        # Keywords for selecting files from Sort_files object
        # This allows a single file type to be used for multiple purposes (e.g., for
        # generating a flat-field image from science exposures)
        self.filetype_keywords = {'SCIENCE':'SCIENCE', 'FLAT':'FLAT', 
            'DARK':'DARK','BIAS':'BIAS'}

        # Header keywords
        self.target_keyword = 'OBJECT'
        self.exptime_keyword = 'EXPTIME'
        self.filter_keyword = 'FILTER'
        self.mjd_keyword = 'MJD-OBS'
        self.bin_keyword = '11'
        self.amp_keyword = '1'

        # File sorting keywords
        self.science_keywords = ['DECKER','MASKNAME']
        self.science_values = ['open','none']
        self.flat_keywords = ['OBJECT']
        self.flat_values = ['gcalflat']
        self.bias_keywords = []
        self.bias_values = []
        self.dark_keywords = ['OBSTYPE']
        self.dark_values = ['dark']
        self.spec_keywords = ['GRORDER']
        self.spec_values = ['1']
        self.bad_keywords = ['DECKER']
        self.bad_values = ['closed']

        self.detrend = False
        self.catalog_zp = '2MASS'

        self.out_size = 2500

    def get_saturation(self, hdr):
        # Gives the saturation level in e-
        return(self.saturation*hdr['NREADS']*hdr['GAIN'])

    def raw_format(self, proc):
        return('*.fits.bz2')

    def get_filter(self, hdr):
        filt = hdr['FILTER'].replace(' ','').split('_')[0]
        if filt == 'Ks':
            filt = 'K'
        return(filt)

    def get_rdnoise(self, hdr):
        return(float(hdr['RDNOISE']))

    def get_gain(self, hdr):
        return(float(hdr['GAIN']))

    def get_exptime(self, hdr):
        return(hdr['EXPTIME'])

    # Get a unique image number that can be derived only from the file header
    def get_number(self, hdr):
        datestr = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
        elap = Time(datestr)-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))
        return(elap)

    def import_image(self, filename, amp, log=None):

        hdu = fits.open(filename)

        # Create header
        hdr = hdu[1].header
        for key in hdu[0].header.keys():
            if key not in hdr.keys():
                try:
                    hdr[key]=hdu[0].header[key]
                except ValueError:
                    continue

        # Adjust WCS to be only 2 axes
        for key in ['CRPIX3','CDELT3','CRVAL3','CD3_3']:
            if key in hdr.keys():
                del hdr[key]
        hdr['WCSAXES']=2

        # Create data array from first read
        data = hdu[1].data
        if len(data.shape)==3:
            data = data[0,:,:]

        raw = CCDData(data, header=hdr, unit=u.adu)
        red = ccdproc.ccd_process(raw, 
            gain=self.get_gain(raw.header)*u.electron/u.adu, 
            readnoise=self.get_rdnoise(raw.header)*u.electron)

        red.header['SATURATE'] = self.get_saturation(red.header)

        return(red)
