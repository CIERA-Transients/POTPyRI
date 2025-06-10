# Parameter file for LRIS/Keck

__version__ = "2.1" # Last edited 05/16/2025

import os
import ccdproc
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData

# Internal dependency
from . import instrument

class LRIS(instrument.Instrument):

    def __init__(self):

        self.version = __version__

        self.name = 'LRIS'

        # Extensions for keeping track of general metadata and WCS
        self.raw_header_ext = 0
        self.wcs_extension = 0

        # Detector specific characteristics
        self.pixscale = 0.135
        # See LRIS non-linearity: 
        # https://www2.keck.hawaii.edu/inst/lris/lris-red-upgrade-notes_vol2.html
        self.saturation = 55000.0

        self.min_exptime = 15.0

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

        self.gain = [1.02, 1.08]
        self.rdnoise = [3.9,4.2,3.6,3.6]

        self.datasec = ['[1055:3024,217:3911]','[1067:3033,217:3911]']
        self.biassec = ['[1:1054,1:4112]','[3034:4096,1:4112]']

        # Keywords for selecting files from Sort_files object
        # This allows a single file type to be used for multiple purposes (e.g., for
        # generating a flat-field image from science exposures)
        self.filetype_keywords = {'SCIENCE':'SCIENCE', 'FLAT':'[SCIENCE,FLAT]', 
            'DARK':'DARK','BIAS':'BIAS'}

        # Header keywords
        self.target_keyword = 'TARGNAME'
        self.exptime_keyword = 'ELAPTIME'
        self.filter_keyword = 'FILTER'
        self.mjd_keyword = 'MJD'
        self.bin_keyword = 'BINNING'
        self.amp_keyword = '2'

        # File sorting keywords
        self.science_keywords = ['KOAIMTYP','SLITNAME','GRANAME','TRAPDOOR']
        self.science_values = ['object','direct','mirror','open']
        self.flat_keywords = ['KOAIMTYP']
        self.flat_values = ['flatlamp']
        self.bias_keywords = ['KOAIMTYP']
        self.bias_values = ['bias']
        self.dark_keywords = []
        self.dark_values = []
        self.spec_keywords = ['GRISTRAN']
        self.spec_values = ['deployed']
        self.bad_keywords = ['SLITNAME','KOAIMTYP']
        self.bad_values = ['goh_lris','focus']

        self.detrend = True
        self.catalog_zp = 'PS1'

        self.out_size = 5000

    # Get a unique image number that can be derived only from the file header
    def get_number(self, header):
        number = str(header['FRAMENO']).zfill(5)
        return(number)

    def raw_format(self, proc):

        if str(proc)=='archive':
            return('*.fits*')
        elif str(proc)=='raw':
            return('*[b,r]*.fits')
        else:
            return('*.fits*')

    def get_instrument_name(self, hdr):
        instrument = hdr['INSTRUME']
        if instrument=='LRISBLUE':
            return('lris.blue')
        elif instrument=='LRIS':
            return('lris.red')
        else:
            raise Exception('Cannot determine instrument name')

    # Specialized procedures for LRIS since we need to deal with red/blue
    def get_filter(self, hdr):
        instrument = hdr['INSTRUME']
        if instrument == 'LRISBLUE':
            filt = hdr['BLUFILT']
        if instrument == 'LRIS':
            filt = hdr['REDFILT']
        return(filt)

    def get_time(self, hdr):
        if 'MJD' in hdr.keys():
            return(float(hdr['MJD']))
        elif 'MJD-OBS' in hdr.keys():
            return(float(hdr['MJD-OBS']))
        else:
            raise Exception('Could not find MJD header keyword.')

    def get_ampl(self, hdr):
        try:
            amp = str(hdr['NUMAMPS'])
        except:
            amp = '1' #str(hdr['NVIDINP'])
        instrument = hdr['INSTRUME']
        if instrument == 'LRISBLUE':
            side = 'B'
        if instrument == 'LRIS':
            side = 'R'

        ampmode = '_'+hdr['AMPMODE'].replace(',','_').replace(':','_')

        fullamp = amp+side+ampmode

        return(fullamp)

    def get_binning(self, header):
        if self.bin_keyword in header.keys():
            binn = str(header[self.bin_keyword]).replace(' ','')
            binn = binn.replace(',','')
            return(binn)
        else:
            return(self.bin_keyword)

    def get_rdnoise(self, hdr):
        amp = self.get_ampl(hdr)
        if amp=='4B_SINGLE_A':
            readnoise = [3.825,3.825,3.825,3.825]
        if amp=='4R_HSPLIT_VSPLIT':
            readnoise = [3.565,3.565,3.565,3.565]
        if amp=='1R_HSPLIT_VUP':
            readnoise = [4]
        if amp=='1R_HSPLIT_VSPLIT':
            readnoise = [4]
        return(readnoise)

    def get_gain(self, hdr):
        amp = self.get_ampl(hdr)
        if  amp=='4B_SINGLE_A':
            gains = [1.55,1.56,1.63,1.70]
        if amp=='4R_HSPLIT_VSPLIT':
            gains = [1.71,1.64,1.61,1.67]
        if amp=='1R_HSPLIT_VUP':
            gains = [1]
        if amp=='1R_HSPLIT_VSPLIT':
            gains = [1]
        return(gains)

    def get_overscan(self, hdr):
        amp = self.get_ampl(hdr)
        binn = int(self.get_binning(hdr)[0])

        if amp=='4B_SINGLE_A':
            oscan_reg = '[1:50,1:4096]'
        if amp=='4R_HSPLIT_VSPLIT':
            oscan_reg = '[1:7,1:2520]'
        if amp=='1R_HSPLIT_VUP':
            b1=int(2065./binn)
            b2=int(2170./binn)
            b3=np.max([int(1./binn),1])
            b4=int(4248./binn)
            oscan_reg = f'[{b1}:{b2},{b3}:{b4}]'
        if amp=='1R_HSPLIT_VSPLIT':
            b1=int(2065./binn)
            b2=int(2170./binn)
            b3=np.max([int(1./binn),1])
            b4=int(4248./binn)
            oscan_reg = f'[{b1}:{b2},{b3}:{b4}]'
        
        return(oscan_reg)

    def get_exptime(self, hdr):
        if self.exptime_keyword in hdr.keys():
            return(float(hdr[self.exptime_keyword]))
        elif 'TTIME' in hdr.keys():
            return(float(hdr['TTIME']))
        else:
            t1 = Time(hdr['DATE-BEG'])
            t2 = Time(hdr['DATE-END'])
            dt = t2 - t1
            return(dt.to_value('sec'))

    def get_1R_datasec(self, amp, binning=1):
        if binning==1:
            if amp==1:
                return(np.s_[830:2064,284:2064])
            elif amp==2:
                return(np.s_[830:2064,2170:3956])
            elif amp==3:
                return(np.s_[2185:3450,284:2064])
            elif amp==4:
                return(np.s_[2185:3450,2170:3956])
        elif binning==2:
            if amp==1:
                return(np.s_[437:1032,146:1032])
            elif amp==2:
                return(np.s_[437:1032,1177:2069])
            elif amp==3:
                return(np.s_[1152:1784,146:1032])
            elif amp==4:
                return(np.s_[1152:1784,1177:2069])

    def get_2R_datasec(self, amp, binning=1):
        if binning==1:
            if amp==1:
                return(np.s_[830:3336,284:2064])
            elif amp==2:
                return(np.s_[830:3336,2170:3956])
        elif binning==2:
            if amp==1:
                return(np.s_[437:1784,146:1032])
            elif amp==2:
                return(np.s_[437:1784,1177:2069])

    def import_image(self, filename, amp, log=None):

        with fits.open(filename) as file_open:
            header = file_open[0].header
            if len(file_open)>1:
                extra_hdr = file_open[1].header
                for key in extra_hdr.keys():
                    if key not in header.keys():
                        header[key] = extra_hdr[key]

        gains = self.get_gain(header)
        readnoises = self.get_rdnoise(header)
        oscan_reg = self.get_overscan(header)

        if amp=='1R_HSPLIT_VUP' or amp=='1R_HSPLIT_VSPLIT':
            with fits.open(filename) as file_open:
                if len(file_open)>1 and file_open[1].name=='COMPRESSED_IMAGE':
                    use_hdu=1
                else:
                    use_hdu=0
                    
            raw = [CCDData.read(filename, hdu=use_hdu, unit='adu')]
            red = [ccdproc.ccd_process(x, oscan=oscan_reg, 
                oscan_model=models.Chebyshev1D(3), 
                gain=gains[j]*u.electron/u.adu, 
                readnoise=readnoises[j]*u.electron) 
                for j,x in enumerate(raw)]
        elif amp=='4B_SINGLE_A' or amp=='4R_HSPLIT_VSPLIT':
            raw = []
            hdu = fits.open(filename)
            for x in range(int(amp[0])):
                hdr = hdu[x+1].header
                hdr['CUNIT1']='deg'
                hdr['CUNIT2']='deg'
                raw.append(CCDData(hdu[x+1].data, header=hdr, unit=u.adu))
            red = [ccdproc.ccd_process(x, oscan=oscan_reg, 
                oscan_model=models.Chebyshev1D(3), trim=x.header['DATASEC'], 
                gain=gains[j]*u.electron/u.adu, 
                readnoise=readnoises[j]*u.electron) 
                for j,x in enumerate(raw)]

        if amp=='1R_HSPLIT_VUP' or amp=='1R_HSPLIT_VSPLIT':
            bin1,bin2 = header['BINNING'].split(',')
            bin1 = float(bin1) ; bin2=float(bin2)
            if amp=='1R_HSPLIT_VUP':
                full = CCDData(
                        np.concatenate(
                        [red[0][self.get_2R_datasec(1, binning=bin1)],
                         red[0][self.get_2R_datasec(2, binning=bin1)]], axis=1),
                            header=header,unit=u.electron)
            elif amp=='1R_HSPLIT_VSPLIT':
                full = CCDData(np.concatenate([np.concatenate(
                        [red[0][self.get_1R_datasec(1, binning=bin1)],
                         red[0][self.get_1R_datasec(2, binning=bin1)]],axis=1),
                        np.concatenate(
                        [red[0][self.get_1R_datasec(3, binning=bin2)],
                         red[0][self.get_1R_datasec(4, binning=bin2)]],axis=1)],
                        axis=0),header=header,unit=u.electron)
            else:
                raise Exception(f'Do not recognize LRISr amp {amp}')
        elif amp=='4B_SINGLE_A':
            full = CCDData(np.concatenate([red[0],
                            np.fliplr(red[1]),
                            np.zeros([np.shape(red[1])[0],111]),
                            red[2],np.fliplr(red[3])],axis=1),
                            header=header,unit=u.electron)
            full = ccdproc.trim_image(full[700:3120,396:3940])
        elif amp=='4R_HSPLIT_VSPLIT':
            full = CCDData(np.concatenate([red[1],
                        np.fliplr(red[0]),
                        np.zeros([np.shape(red[0])[0],200]),
                        red[3],np.fliplr(red[2])],axis=1),
                        header=header, unit=u.electron)

        full.header['SATURATE'] = self.saturation

        return(full)
