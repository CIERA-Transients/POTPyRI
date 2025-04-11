# Parameter file for MOSFIRE/Keck

__version__ = "2.0" # Last edited 09/21/2024

import os
import ccdproc
import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData

# Internal dependency
from . import instrument

class MOSFIRE(instrument.Instrument):

    def __init__(self):

        self.version = __version__

        self.name = 'MOSFIRE'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.1799
        self.saturation = 65000.0

        self.min_exptime = 1.0

        # Run dark/bias/flat calibration?
        self.dark = False
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
        self.filetype_keywords = {'SCIENCE':'SCIENCE', 'FLAT':'[SCIENCE,FLAT]', 
            'DARK':'DARK','BIAS':'BIAS'}

        # Header keywords
        self.target_keyword = 'TARGNAME'
        self.exptime_keyword = 'ELAPTIME'
        self.filter_keyword = 'DWFILNAM'
        self.mjd_keyword = 'MJD-OBS'
        self.bin_keyword = '11'
        self.amp_keyword = '1'

        # File sorting keywords
        self.science_keywords = ['GRATMODE','MASKNAME']
        self.science_values = ['imaging','open']
        self.flat_keywords = ['OBJECT']
        self.flat_values = ['flat']
        self.bias_keywords = []
        self.bias_values = []
        self.dark_keywords = ['FILTER']
        self.dark_values = ['dark']
        self.spec_keywords = ['GRATMODE']
        self.spec_values = ['spectroscopy']
        self.bad_keywords = ['MASKNAME','PONAME']
        self.bad_values = ['closed','mira']

        self.detrend = True
        self.catalog_zp = '2MASS'

        self.out_size = 2500

    def get_saturation(self, hdr):
        return(hdr['SATURATE']*hdr['SYSGAIN'])

    def raw_format(self, proc):
        return('MF.*.fits.gz')

    def get_filter(self, hdr):
        filt = hdr['FILTER'].replace(' ','').split('_')[0]
        if filt == 'Ks':
            filt = 'K'
        return(filt)

    def get_rdnoise(self, hdr):
        readnoise = {1:21,4:10.8,8:7.7,16:5.8,32:4.2,64:3.5,128:3.0}
        return readnoise[hdr['NUMREADS']]

    def get_gain(self, hdr):
        return(hdr['SYSGAIN'])

    def get_exptime(self, hdr):
        return(hdr['TRUITIME']*hdr['COADDONE'])

    # Get a unique image number that can be derived only from the file header
    def get_number(self, header):
        number = str(header['FRAMENO']).zfill(5)
        return(number)

    def import_image(self, filename, amp, log=None):

        raw = CCDData.read(filename, unit=u.adu)
        red = ccdproc.ccd_process(raw, 
            gain=self.get_gain(raw.header)*u.electron/u.adu, 
            readnoise=self.get_rdnoise(raw.header)*u.electron)

        red.header['SATURATE'] = self.get_saturation(red.header)

        return(red)
