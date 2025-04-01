# Parameter file for GMOS/Gemini

__version__ = "2.1" # Last edited 09/29/2024

import os
import ccdproc
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import models
from astropy.nddata import CCDData
from astropy.stats import SigmaClip
from astropy.time import Time

# Internal dependency
from . import instrument

class GMOS(instrument.Instrument):

    def __init__(self):

        self.version = __version__

        self.name = 'GMOS'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.0803
        self.saturation = 65535.0

        self.min_exptime = 1.0

        # Run dark/bias/flat calibration?
        self.dark = False
        self.bias = True
        self.flat = True

        # Parameters for handling calibration files
        # Run rejection on possible CR pixels in bias
        self.cr_bias = True 

        # Extend header to first file extension for purposes of sort_files
        self.extend_header = True

        # How to combine images during stacking
        self.stack_method = 'median'

        self.wavelength = 'OPT'

        self.gain = 1.63
        self.rdnoise = 4.14

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
        self.bin_keyword = 'CCDSUM'
        self.amp_keyword = 'NCCDS'

        # File sorting keywords
        self.science_keywords = ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']
        self.science_values = ['open','mirror','science','object']
        self.flat_keywords = ['SHUTTER','GRATING','OBSCLASS','OBSTYPE']
        self.flat_values = ['open','mirror','daycal','object']
        self.bias_keywords = ['SHUTTER','OBSCLASS','OBSTYPE']
        self.bias_values = ['closed','daycal','bias']
        self.dark_keywords = []
        self.dark_values = []
        self.spec_keywords = ['SHUTTER','FILTER2']
        self.spec_values = ['open','open2-8']
        self.bad_keywords = ['RAWGEMQA']
        self.bad_values = ['fail']

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 3200

    # Perform sky subtraction with GMOS z-band data
    def needs_sky_subtraction(self, filt):
        if filt.lower().startswith('z'):
            return(True)
        else:
            return(False)

    def get_rdnoise(self, hdr):
        try:
            return(hdr['RDNOISE'])
        except KeyError:
            return(self.rdnoise)

    def get_gain(self, hdr):
        try:
            return(hdr['GAIN'])
        except KeyError:
            return(self.gain)

    def raw_format(self, proc):
        if str(proc).lower()=='dragons':
            return('*.fits')
        else:
            return('*.fits.bz2')

    def get_ampl(self, hdr):
        nccd = hdr['NCCDS']
        if nccd == '1':
            amp = '4'
        else:
            amp = '12'
        return(amp)

    def get_time(self, hdr):
        return(Time(hdr['DATE-OBS']+'T'+hdr['TIME-OBS']).mjd)

    def get_number(self, hdr):
        datestr = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
        elap = Time(datestr)-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))
        return(elap)

    def get_catalog(self, hdr):
        coord = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.deg, u.deg))
        if coord.dec.degree < -30:
            return('SkyMapper')
        else:
            return('PS1')
        
        return('PS1')

    def get_instrument_name(self, hdr):
        if 'INSTRUME' in hdr.keys():
            return(str(hdr['INSTRUME']).lower())
        else:
            return(self.name.lower())

    def import_image(self, filename, amp, log=None):

        with fits.open(filename) as hdr:
            header = hdr[0].header
            binn = self.get_binning(hdr[1].header)

        raw = [CCDData.read(filename, hdu=x+1, unit='adu') 
            for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
            oscan_model=models.Chebyshev1D(3), 
            trim=x.header['DATASEC'], 
            gain=self.get_gain(x.header)*u.electron/u.adu, 
            readnoise=self.get_rdnoise(x.header)*u.electron, 
            gain_corrected=True) 
            for k,x in enumerate(raw)]

        if header['INSTRUME'] == 'GMOS-S':
            if str(amp) == '12':
                gap = int(80/int(binn[0]))
                full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],gap])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],gap])]+\
                    red[8:10],axis=1),
                    header=header,unit=u.electron)
            elif str(amp) == '4':
                full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)
        elif header['INSTRUME'] == 'GMOS-N':
            if str(amp) == '12':
                gap = int(60/int(binn[0]))
                full = CCDData(np.concatenate(
                    red[2:4]+[np.zeros([np.shape(red[3])[0],gap])]+\
                    red[4:8]+[np.zeros([np.shape(red[7])[0],gap])]+\
                    red[8:10],axis=1),
                    header=header,unit=u.electron)
            elif str(amp) == '4':
                full = CCDData(np.concatenate(red,axis=1),
                    header=header,unit=u.electron)

        full.header['SATURATE'] = self.saturation

        # Add header keyword for binning
        full.header['CCDSUM'] = binn

        return(full)
