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

# Internal dependency
from . import instrument

__version__ = 1.4 #last edited 24/08/2021

class BINOSPEC(instrument.Instrument):

    def __init__(self):

        # Version
        self.version = __version__

        # Intrument name
        self.name = 'BINOSPEC'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 1
        self.wcs_extension = 1

        # Detector specific characteristics
        self.pixscale = 0.24
        self.saturation = 65000

        self.min_exptime = 1.0

        # Run dark/bias/flat calibration?
        self.dark = False
        self.bias = False
        self.flat = True

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
        self.spec_keywords = ['MASK']
        self.spec_values = ['spectroscopy']
        self.bad_keywords = ['MASK']
        self.bad_values = ['mira']

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 5000

    # Get a unique image number that can be derived only from the file header
    def get_number(self, header):
        number = str(int(Time(header['DATE-OBS']).mjd*1e5)-5900000000)
        return(number)

    def import_image(self, filename, amp, log=None):
        filename = os.path.abspath(filename)
        if log: log.info(f'Loading file: {filename}')

        with fits.open(filename) as hdr:
            header = hdr[self.raw_header_ext].header

        raw = [CCDData.read(filename, hdu=x+1, unit='adu') for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=self.biassec[k], 
               oscan_model=models.Chebyshev1D(3), 
               trim=self.datasec[k], gain=self.gain[k]*u.electron/u.adu, 
               readnoise=self.rdnoise[k]*u.electron, 
               gain_corrected=True) for k,x in enumerate(raw)]

        frame_full = CCDData(np.concatenate((red[0], 
            np.empty((red[0].shape[0], 794)), red[1]), axis=1), 
            header=header,unit=u.electron)

        frame_full.header['SATURATE'] = self.saturation

        return(frame_full)
