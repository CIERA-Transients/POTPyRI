# Parameter file for BINOSPEC/MMT

__version__ = "2.0" # Last edited 09/21/2024

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
        self.spec_keywords = ['MASK']
        self.spec_values = ['spectroscopy']
        self.bad_keywords = ['MASK']
        self.bad_values = ['mira']

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 5000

    def raw_format(self, proc):
        if proc:
            return('sci_img_*proc.fits*')
        else:
            return('sci_img*[!proc].fits*')

    # Get a unique image number that can be derived only from the file header
    def get_number(self, hdr):
        datestr = hdr['DATE-OBS']
        elap = Time(datestr)-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))
        return(elap)

    def import_image(self, filename, amp, log=None):
        filename = os.path.abspath(filename)
        if log: log.info(f'Loading file: {filename}')

        hdu = fits.open(filename)
        header = hdu[self.raw_header_ext].header

        raw = []
        for x in range(int(amp)):
            data = hdu[x+1].data
            hdr = hdu[x+1].header

            if 'CRVAL1' in hdr.keys() and hdr['CRVAL1']<0:
                hdr['CRVAL1']=0.0
            if 'CRVAL2' in hdr.keys() and np.abs(hdr['CRVAL2'])>90:
                hdr['CRVAL2']=0.0
            if 'PRESSURE' in hdr.keys():
                del hdr['PRESSURE']
            if 'TEMP' in hdr.keys():
                del hdr['TEMP']

            raw.append(CCDData(data, header=hdr, unit=u.adu))

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
