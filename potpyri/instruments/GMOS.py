#parameter file for GMOS/Gemini
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

__version__ = 1.2 #last edited 09/11/2021

# Internal dependency
from . import instrument

class GMOS(instrument.Instrument):

    def __init__(self, proc=None, ):

        self.version = __version__

        self.name = 'GMOS'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.0803*2
        self.saturation = 65535.0

        self.min_exptime = 1.0

        # Run dark/bias/flat calibration?
        self.dark = False
        self.bias = True
        self.flat = True

        # Parameters for handling calibration files
        # Run rejection on possible CR pixels in bias
        self.cr_bias = True 

        # How to combine images during stacking
        self.stack_method = 'median'

        self.wavelength = 'OPT'

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
        self.filter_keyword = 'FILTER2'
        self.mjd_keyword = 'MJD'
        self.bin_keyword = '22'
        self.amp_keyword = 'NCCDS'

        # File sorting keywords
        self.science_keywords = ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']
        self.science_values = ['OPEN','MIRROR','science','OBJECT']
        self.flat_keywords = ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']
        self.flat_values = ['OPEN','MIRROR','dayCal','OBJECT']
        self.bias_keywords = ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']
        self.bias_values = ['CLOSED','MIRROR','dayCal','BIAS']
        self.dark_keywords = []
        self.dark_values = []
        self.spec_keywords = ['SHUTTER','FILTER2']
        self.spec_values = ['OPEN','open2-8']
        self.bad_keywords = ['RAWGEMQA']
        self.bad_values = ['FAIL']

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 5000

    def raw_format(self, proc):
        return '*.fits.bz2'

    def get_ampl(self, hdr):
        nccd = hdr['NCCDS']
        if nccd == '1':
            amp = '4'
        else:
            amp = '12'
        return amp

    def get_time(self, hdr):
        return Time(hdr['DATE-OBS']+'T'+hdr['TIME-OBS']).mjd

    def get_number(self, hdr):
        datestr = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
        elap = Time(datestr)-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))

        return elap

    def get_catalog(self, hdr):
        coord = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.deg, u.deg))
        if coord.dec.degree < -30:
            return('SkyMapper')
        else:
            return('PS1')
        
        return('PS1')

    def import_image(self, filename, amp, log=None):

        with fits.open(flat) as hdr:
            header = hdr[0].header

        raw = [CCDData.read(filename, hdu=x+1, unit='adu') 
            for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
            oscan_model=models.Chebyshev1D(3), 
            trim=x.header['DATASEC'], 
            gain=x.header['GAIN']*u.electron/u.adu, 
            readnoise=x.header['RDNOISE']*u.electron, 
            gain_corrected=True) 
            for k,x in enumerate(raw)]

        if header['INSTRUME'] == 'GMOS-S':
            if amp == '12':
                full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],40])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],40])]+\
                    red[8:10],axis=1),
                    header=header,unit=u.electron)
            elif amp == '4':
                full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)
        elif header['INSTRUME'] == 'GMOS-N':
            if amp == '12':
                full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],30])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],30])]+\
                    red[8:10],axis=1),
                    header=header,unit=u.electron)
            elif amp == '4':
                full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)

        full.header['SATURATE'] = self.saturation

        return(full)

    def edit_stack_headers(self, stack, log=None):

        good_keys = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND',
            'INSTRUME','OBJECT','OBSTYPE','GEMPRGID','OBSERVER','TELESCOP',
            'EPOCH','SSA','RA','DEC','UT','DATE','DETECTOR','DATE-OBS',
            'AIRMASS','MJD-OBS','GAIN','RDNOISE','SATURATE','SKYBKG','BUNIT',
            'CTYPE1','CTYPE2','CUNIT1','CUNIT2','RADECSYS','CD1_1','CD1_2',
            'CD2_1','CD2_2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','NFILES','FILTER',
            'AMPS','BINNING','XTENSION','EXTNAME','TIMESYS','EXPTOT','WCSAXES']

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

    def edit_raw_headers(self, files, log=None):
        pass
