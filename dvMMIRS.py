#parameter file for MMIRS/MMT
import os
import datetime

def datetime_to_float(d):
    return d.timestamp()

def float_to_datetime(fl):
    return datetime.datetime.fromtimestamp(fl)

def static_mask():
    return './staticmasks/MMIRS.staticmask.fits'

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.202

def ref_pix():
    return 852.492546082, 851.589035034

def WCS_keywords(): #WCS keywords
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
    return os.getenv("HOME")+'/Pipelines/MMIRS_calib/'

def raw_format():
    return '*.fits'

def sky():
    return True

def dark():
    return True

def flat():
    return False

def raw_header_ext():
    return 1

def science_keyword():
    return ['OBSMODE','APTYPE']

def science_files():
    return ['imaging','open']

def flat_keyword():
    return ['']

def flat_files():
    return ['']

def bias_keyword():
    return ['']

def bias_files():
    return ['']

def dark_keyword():
    return ['OBJECT']

def dark_files():
    return ['DARK']

def target_keyword():
    return 'OBJECT'

def filter_keyword():
    return 'FILTER'

def time_format(hdr):
    return (datetime.datetime.strptime(hdr['DATE-OBS'],'%Y-%m-%dT%H:%M:%S')).timestamp()