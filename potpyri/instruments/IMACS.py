# Parameter file for Magellan/IMACS

__version__ = "1.0" # Last edited 02/13/2025

import os
import ccdproc
import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.modeling import models
from astropy.nddata import CCDData
from astropy.time import Time

# Internal dependency
from . import instrument

class IMACS(instrument.Instrument):

    def __init__(self):

        # Version
        self.version = __version__

        # Intrument name
        self.name = 'IMACS'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.111
        self.saturation = 60000

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
        self.bin_keyword = 'BINNING'
        self.amp_keyword = '2'

        # File sorting keywords
        self.science_keywords = ['GISMO','SLITMASK']
        self.science_values = ['none','f/4-imaging']
        self.flat_keywords = ['OBJECT']
        self.flat_values = ['twiflat']
        self.bias_keywords = []
        self.bias_values = []
        self.dark_keywords = []
        self.dark_values = []
        self.spec_keywords = []
        self.spec_values = []
        self.bad_keywords = []
        self.bad_values = []

        self.detrend = False
        self.catalog_zp = 'PS1'

        self.out_size = 4200

    # Get a unique image number that can be derived only from the file header
    def get_number(self, hdr):
        datestr = hdr['DATE-OBS']
        elap = Time(datestr)-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))
        return(elap)

    def get_time(self, hdr):
        return(float(Time(hdr['DATE-OBS']).mjd))

    def get_filter(self, hdr):
        filt = hdr[self.filter_keyword]
        filtmap = {'Bessell_V1':'V',
                   'CTIO-I1': 'I'}
        if filt in filtmap.keys():
            filt = filtmap[filt]
        return(filt)

    def get_ampl(self, hdr):
        try:
            amp = str(hdr['CHIP'])
        except:
            amp = None
        return(amp)

    # Raw image format for ingestion
    def raw_format(self, proc):
        return('iff*.fits*')

    def import_image(self, filename, amp, log=None):
        filename = os.path.abspath(filename)
        if log: log.info(f'Loading file: {filename}')

        hdu = fits.open(filename)

        data = hdu[0].data[:4096,:]
        header = hdu[self.raw_header_ext].header

        raw = CCDData(data, header=header, unit=u.adu)
        red = ccdproc.ccd_process(raw, oscan=header['BIASSEC'], 
               oscan_model=models.Chebyshev1D(3), 
               trim=header['DATASEC'], gain=header['EGAIN']*u.electron/u.adu, 
               readnoise=header['ENOISE']*u.electron, 
               gain_corrected=True)

        frame_full = red
        frame_full.header['SATURATE'] = self.saturation

        return(frame_full)
