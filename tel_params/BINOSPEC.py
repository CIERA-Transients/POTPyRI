#parameter file for BINOSPEC/MMT
import os
import datetime
import numpy as np
from photutils import make_source_mask, Background2D, MeanBackground
from astropy.stats import SigmaClip
from astropy.io import fits
from astropy.time import Time
from astropy.nddata import CCDData
import astropy.units.astrophys as u
import astropy.units as u
import ccdproc
from astropy.modeling import models
import astropy.wcs as wcs

def static_mask(proc):
    if proc:
        return [None,None]#['./staticmasks/bino_proc_left.staticmask.fits','./staticmasks/bino_proc_right.staticmask.fits']
    else:
        return ['./staticmasks/bino_left.staticmask.fits','./staticmasks/bino_right.staticmask.fits']

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.24

def ref_pix(): #place holder
    return 852.492546082, 851.589035034

def WCS_keywords(): #place holder
    WAT0_001 = 'system=image'
    WAT1_001 = 'wtype=tnx axtype=ra lngcor = "3. 4. 4. 2. -0.08078770871202606 0.078'
    WAT1_002 = '07957656086426 -0.07723307309820749 0.08164053570655277 -2.686524069'
    WAT1_003 = '384972E-6 -8.121054526384903E-4 0.01014393940325073 0.37221863445133'
    WAT1_004 = '59 2.360195228674920E-4 5.001808549237395E-4 -0.7035825742463017 -0.'
    WAT1_005 = '006139017761099392 0.3323574770047805 0.03786685758120669 "'
    WAT2_001 = 'wtype=tnx axtype=dec latcor = "3. 4. 4. 2. -0.08078770871202606 0.07'
    WAT2_002 = '807957656086426 -0.07723307309820749 0.08164053570655277 -4.41223056'
    WAT2_003 = '0109198E-6 1.020477159783315E-4 -0.001623611247485439 -0.11558053285'
    WAT2_004 = '78736 -0.001787914383662314 0.00591838599589976 0.4698673160370906 0'
    WAT2_005 = '.002964492811835905 0.3360926419983717 0.7822291841773272 "'
    return WAT0_001, WAT1_001, WAT1_002, WAT1_003, WAT1_004, WAT1_005, WAT2_001, WAT2_002, WAT2_003, WAT2_004, WAT2_005

def cal_path():
    return str(os.getenv("PIPELINE_HOME"))+'/Imaging_pipelines/BINOSPEC_calib/'

def raw_format(proc):
    if proc:
        return 'sci_img_*proc.fits'
    else:
        return 'sci_img*[!proc].fits'

def dark():
    return False

def bias():
    return False

def flat():
    return False

def raw_header_ext():
    return 1

def science_keyword():
    return ['MASK','SCRN']

def science_files():
    return ['imaging','stowed']

def flat_keyword():
    return ['MASK','SCRN']

def flat_files():
    return ['imaging','deployed']

def bias_keyword():
    return ['']

def bias_files():
    return [None]

def dark_keyword():
    return ['']

def dark_files():
    return [None]

def target_keyword():
    return 'OBJECT'

def filter_keyword():
    return 'FILTER'

def time_format(hdr):
    return hdr['MJD']

def wavelength():
    return 'OPT'

def flat_name(cpath,fil):
    return [cpath+'/mflat_'+fil+'_left.fits',cpath+'/mflat_'+fil+'_right.fits']

def load_flat(flat):
    mflat = []
    for f in flat:
        mflat.append(CCDData.read(f,hdu=1,unit=u.electron))
    return mflat

def create_flat(flat_list):
    return None

def gain():
    return [1.085, 1.04649118, 1.04159151, 0.97505369, 1.028, 1.16341855, 1.04742053, 1.0447564]

def process_science(sci_list,fil,cal_path,mdark=None,mbias=None,mflat=None,proc=None):
    masks = []
    processed = []
    flat_left = mflat[0]
    flat_right = mflat[1]
    left_list = []
    right_list = []
    if proc:  
        for j,sci in enumerate(sci_list):
            left = CCDData.read(sci, hdu=1, unit=u.electron)
            left = ccdproc.flat_correct(left,flat_left)
            left = ccdproc.ccd_process(left, trim=left.header['DATASEC'])
            mask = make_source_mask(left, nsigma=3, npixels=5)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_mask_left.fits'),mask.astype(int),overwrite=True)
            bkg = Background2D(left, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_left.fits'),bkg.background,overwrite=True)
            left = left.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            left = left.divide(left.header['EXPTIME'],propagate_uncertainties=True,handle_meta='first_found')
            left.header['DATASEC'] = '[1:'+str(np.shape(left)[1])+',1:'+str(np.shape(left)[0])+']'
            left_list.append(left)
            right = CCDData.read(sci, hdu=2, unit=u.electron)
            right = ccdproc.flat_correct(right,flat_right)
            right = ccdproc.ccd_process(right, trim=right.header['DATASEC'])
            mask = make_source_mask(right, nsigma=3, npixels=5)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_mask_right.fits'),mask.astype(int),overwrite=True)
            bkg = Background2D(right, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_right.fits'),bkg.background,overwrite=True)
            right = right.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            right = right.divide(right.header['EXPTIME'],propagate_uncertainties=True,handle_meta='first_found')
            right.header['DATASEC'] = '[1:'+str(np.shape(right)[1])+',1:'+str(np.shape(right)[0])+']'
            right_list.append(right)
    else:
        for j,sci in enumerate(sci_list):
            with fits.open(sci) as hdr:
                header_left = hdr[1].header
                header_right = hdr[6].header
            data_list = []
            for i in range(8):
                data = ccdproc.CCDData.read(sci,hdu=i+1,unit=u.adu)
                red = ccdproc.ccd_process(data, oscan=data[:,0:50], oscan_model=models.Chebyshev1D(3), trim='[1200:2098,210:2056]', gain=gain()[i]*u.electron/u.adu, readnoise=4*u.electron)
                data_list.append(np.asarray(red).astype(np.float32))
            top_left = np.concatenate([data_list[0],np.fliplr(data_list[1])],axis=1)
            bot_left = np.flipud(np.concatenate([data_list[3],np.fliplr(data_list[2])],axis=1))
            left = CCDData(np.concatenate([top_left,bot_left]),unit=u.electron,header=header_left,wcs=wcs.WCS(header_left))
            left = ccdproc.flat_correct(left,flat_left[209:3903,1149:2947])
            mask = make_source_mask(left, nsigma=3, npixels=5)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_mask_left.fits'),mask.astype(int),overwrite=True)
            bkg = Background2D(left, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_left.fits'),bkg.background,overwrite=True)
            left = left.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            left = left.divide(left.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            left.header['DATASEC'] = '[1:1798,1:3694]'
            left_list.append(left)
            top_right = np.concatenate([data_list[6],np.fliplr(data_list[7])],axis=1)
            bot_right = np.flipud(np.concatenate([data_list[5],np.fliplr(data_list[4])],axis=1))
            right = CCDData(np.concatenate([top_right,bot_right]),unit=u.electron,header=header_right,wcs=wcs.WCS(header_right))
            right = ccdproc.flat_correct(right,flat_right[209:3903,1149:2947])
            mask = make_source_mask(right, nsigma=3, npixels=5)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_mask_right.fits'),mask.astype(int),overwrite=True)
            bkg = Background2D(right, (120, 120), filter_size=(3, 3),sigma_clip=SigmaClip(sigma=3), bkg_estimator=MeanBackground(), mask=mask, exclude_percentile=80)
            fits.writeto(sci.replace('/raw/','/red/').replace('.fits','_bkg_right.fits'),bkg.background,overwrite=True)
            right = right.subtract(CCDData(bkg.background,unit=u.electron),propagate_uncertainties=True,handle_meta='first_found')
            right = right.divide(right.header['EXPTIME']*u.second,propagate_uncertainties=True,handle_meta='first_found')
            right.header['DATASEC'] = '[1:1798,1:3694]'
            right_list.append(right)
    return [left_list,right_list], None

def stacked_image(tar,red_path):
    return [red_path+tar+'_left.fits',red_path+tar+'_right.fits']

def suffix():
    return ['_red_left.fits','_red_right.fits']

def rdnoise(header):
    return 4.0