#parameter file for DEIMOS/Keck
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

__version__ = 1.5 #last edited 23/11/2021

# Internal dependency
from . import instrument

class DEIMOS(instrument.Instrument):

    def __init__(self, proc=None, ):

        self.version = __version__

        self.name = 'DEIMOS'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.1185
        self.saturation = 65000.0

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
        self.exptime_keyword = 'ELAPTIME'
        self.filter_keyword = 'DWFILNAM'
        self.mjd_keyword = 'MJD-OBS'
        self.bin_keyword = 'BINNING'
        self.amp_keyword = 'NVIDINP'

        # File sorting keywords
        self.science_keywords = ['MOSMODE','OBSMODE','SLMSKNAM','HATCHPOS']
        self.science_values = ['Direct','imaging','None','open']
        self.flat_keywords = ['OBSTYPE','OBJECT']
        self.flat_values = ['DmFlat','lat']
        self.bias_keywords = ['OBSTYPE']
        self.bias_values = ['Bias']
        self.dark_keywords = []
        self.dark_values = []
        self.spec_keywords = ['MOSMODE','OBSMODE','SLMSKNAM']
        self.spec_values = ['Spectral','longslit','LVMslitC']
        self.bad_keywords = ['SLMSKNAM','PONAME']
        self.bad_values = ['GOH_X','Mira']

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 9000

    def raw_format(self, proc):
        if proc and str(proc)=='raw':
            return 'd*.fits'
        else:
            return 'DE*.fits.gz'

    def get_number(self, hdr):
        elap = Time(hdr['MJD'], format='mjd')-Time('1980-01-01')
        elap = int(np.round(elap.to(u.second).value))

        return elap

    def get_ampl(self, hdr):
        return(str(hdr['NVIDINP']))

    def get_gain(self, hdr):
        if amp=='4':
            gains = [1.206, 1.722, 1.211, 1.231]
        if amp=='8':
            gains = [1.232, 1.180, 1.714, 1.173, 1.161, 1.261, 1.246, 1.216]
        return gains

    def get_rdnoise(self, hdr):
        if amp=='4':
            readnoise = [2.528, 2.128, 2.5395, 3.2335]
        if amp=='8':
            readnoise = [2.583, 2.473, 1.797, 2.459, 2.434, 2.645, 3.918, 2.549]
        
        return readnoise

    def import_image(self, filename, amp, log=None):
        with fits.open(filename) as hdr:
            header = hdr[0].header

        gains = self.get_gain(header)
        readnoises = self.get_rdnoise(header)

        if header['DETSEC01'] == '[1:2048,1:4096]':
            trim_sec = '[13:2060,1321:3921]'
        else:
            trim_sec = '[13:2060,1:2601]'
        
        raw = [CCDData.read(flat, hdu=x+1, unit='adu') 
            for x in range(int(amp))]
        red = [ccdproc.ccd_process(x, oscan=x[:,0:13], 
            oscan_model=models.Chebyshev1D(3), trim=trim_sec, 
            gain=gains[k]*u.electron/u.adu, 
            readnoise=readnoises[k]*u.electron, 
            gain_corrected=True) 
            for k,x in enumerate(raw)]

        if amp == '4':
            full = CCDData(np.concatenate([red[0],np.zeros([2601,113]),
                red[1],np.zeros([2601,81]),red[2],np.zeros([2601,113]),
                red[3]],axis=1),
                header=header,
                unit=u.electron/u.second)
        if amp == '8':
            full = CCDData(np.concatenate([red[1],np.fliplr(red[0]),
                red[3],np.fliplr(red[2]),red[5],np.fliplr(red[4]),red[7],
                np.fliplr(red[6])],axis=1),
                header=header,
                unit=u.electron/u.second)

        full.header['SATURATE'] = self.saturation

        return(full)

    def edit_raw_headers(self, files, log=None):

        for file in files:

            basefile = os.path.basename(file)
            if basefile.startswith('DE') and basefile.endswith('.gz'): continue

            if log: log.info(f'Editing headers for: {file}')

            hdu = fits.open(file)
            h = hdu[0].header

            if ('ELAPTIME' not in h.keys() and 'DATE-BEG' in h.keys() and
                'DATE-END' in h.keys()):

                t1 = Time(h['DATE-BEG'])
                t2 = Time(h['DATE-END'])
                dt = t2 - t1

                hdu[0].header['ELAPTIME']=dt.to_value('sec')

            if (('twi' in h['OBJECT'].lower() and 'flat' in h['OBJECT'].lower()) or
                ('blank' in h['OBJECT'].lower()) or
                ('t_flat' in h['OBJECT'].lower()) or
                ('sky' in h['OBJECT'].lower() and 'flat' in h['OBJECT'].lower())):
                hdu[0].header['OBSTYPE']='DmFlat'

            if 'bias' in h['OBJECT'].lower():
                hdu[0].header['OBSTYPE']='Bias'

            if 'PONAME3' in h.keys() and h['PONAME3'].lower().strip()=='image':
                hdu[0].header['OBSMODE']='imaging'
            else:
                hdu[0].header['OBSMODE']='spectroscopy'

            hdu.writeto(file, overwrite=True, output_verify='silentfix')

    def edit_stack_headers(self, stack, log=None):

        good_keys = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND',
            'OBSERVAT','TELESCOP','ORIGIN','INSTRUME','EXTNAME','FILTER',
            'DATE-OBS','RA','DEC','POSANG','AZ','EL','ROT','PA',
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
