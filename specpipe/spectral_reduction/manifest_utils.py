from __future__ import print_function
import os,sys,pdb,glob,datetime,shutil
import numpy as np 
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u


### M.S. Stolen from J.B.'s spec browser. Primarily used to determine image types for config files

def make_manifest(dataDir,instrument):
    if instrument == 'LRIS':
        res = make_lris_manifest(dataDir)
    elif instrument == 'KAST':
        res = make_kast_manifest(dataDir)
    else:
        raise ValueError('Instrument {} not supported'.format(instrument))
    return res

def make_lris_manifest(dataDir):
    ''' Make a manifest following the LRIS format '''

    manifest_header = '# Object                  Type Chan '
    manifest_header += 'RA             Dec            '
    manifest_header += 'PA      Exp     UTCDate/Time                  '
    manifest_header += 'SecZ    SlitMask       Dich      Gris/Grat      '
    manifest_header += 'Filter         Filename'
    manifest_text = '{}'.format(manifest_header)

    nullHeaderEntry = 'UNKNOWN'

    # glob all the files in the directory
    allFiles = glob.glob('{}/???????_????.fits'.format(dataDir))
    allFiles.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

    bFiles = []
    rFiles = []
    allSortedFiles = []
    for i,obsFile in enumerate(allFiles):
        if obsFile.split('/')[-1][0].lower() == 'b':
            bFiles.append(obsFile)
        else:
            rFiles.append(obsFile)
    allSortedFiles = bFiles + rFiles

    STANDARD_STAR_LIBRARY = construct_standard_star_library() # dict of std star objects

    #loop over the sorted list
    for i,obsFile in enumerate(allSortedFiles):

        print('Working on {}'.format(obsFile))

        # open the header
        hdu = fits.open(obsFile)
        header = hdu[0].header

        # get the image type
        obsImgType = determine_image_type(header,'LRIS',STANDARD_STAR_LIBRARY)

        # get blue or red channel
        obsInst = header.get('INSTRUME',nullHeaderEntry)
        if obsInst.strip().lower() == 'lrisblue':
            obsChannel = 'BLUE'
        elif obsInst.strip().lower() == 'lris':
            obsChannel = 'RED'
        else:
            obsChannel = nullHeaderEntry

        #get blue or red filter
        if obsChannel == 'BLUE':
            obsFilter = header.get('BLUFILT',nullHeaderEntry)
        elif obsChannel == 'RED':
            obsFilter = header.get('REDFILT',nullHeaderEntry)
        else:
            obsFilter = nullHeaderEntry

        # grating/grism
        if obsChannel == 'BLUE':
            obsGrating = header.get('GRISNAME',nullHeaderEntry).strip().replace(' ','_')
        elif obsChannel == 'RED':
            obsGrating = header.get('GRANAME',nullHeaderEntry).strip().replace(' ','_')
        else:
            obsGrating = nullHeaderEntry

        # get the explicitly named/non-degenerate parameters
        obsDateOnly = header.get('DATE-OBS',nullHeaderEntry)
        obsTimeOnly = header.get('UTC',nullHeaderEntry)
        obsDate = '{}T{}'.format(obsDateOnly,obsTimeOnly)
        obsObject = header.get('TARGNAME',nullHeaderEntry).strip().replace(' ','_')
        obsRA = header.get('RA',nullHeaderEntry)
        obsDec = header.get('DEC',nullHeaderEntry)
        obsPositionAngle = header.get('ROTPOSN',nullHeaderEntry)
        obsAirmass = header.get('AIRMASS',nullHeaderEntry)
        obsExptime = header.get('TTIME',nullHeaderEntry)
        obsSlitmask = header.get('SLITNAME',nullHeaderEntry).strip().replace(' ','_')
        obsDichroic = header.get('DICHNAME',nullHeaderEntry).strip().replace(' ','_')

        #formatting tweaks to values
        if obsDec[0] != '-' and obsDec[0] != '+':
            obsDec = '+{}'.format(obsDec)
        if len(obsDec.split(':')[0]) < 3:
            obsDecDD = obsDec.split(':')[0]
            obsDecMM = obsDec.split(':')[1]
            obsDecSS = obsDec.split(':')[2]
            obsDecDD = '{}0{}'.format(obsDecDD[0],obsDecDD[1])
            obsDec = '{}:{}:{}'.format(obsDecDD,obsDecMM,obsDecSS)

        # prevent empty object names from getting thru
        if obsObject == '':
            obsObject = nullHeaderEntry

        # construct the string
        manifest_text += '\n{:<25} '.format(obsObject)
        manifest_text += '{:<10}'.format(obsImgType)
        manifest_text += '{:<5}'.format(obsChannel)
        manifest_text += '{:<15}{:<15}'.format(obsRA,obsDec)
        if obsPositionAngle != nullHeaderEntry:
            manifest_text += '{:<8.1f}'.format(obsPositionAngle)
        else:
            manifest_text += '{:<8}'.format(nullHeaderEntry)
        if obsExptime != nullHeaderEntry:
            manifest_text += '{:<8.1f}'.format(obsExptime)
        else:
            manifest_text += '{:<8}'.format(nullHeaderEntry)
        manifest_text += '{:<30}'.format(obsDate)
        if obsAirmass != nullHeaderEntry:
            manifest_text += '{:<8.2f}'.format(obsAirmass)
        else:
            manifest_text += '{:<8}'.format(nullHeaderEntry)
        manifest_text += '{:<15}'.format(obsSlitmask)
        manifest_text += '{:<10}'.format(obsDichroic)
        manifest_text += '{:<15}'.format(obsGrating)
        manifest_text += '{:<15}'.format(obsFilter)
        manifest_text += '{:<20}'.format(obsFile.split('/')[-1])



    return manifest_text

def make_kast_manifest(dataDir):
    ''' Make a manifest following the Kast format '''

    manifest_header = '# Object                  Type Chan '
    manifest_header += 'RA             Dec            '
    manifest_header += 'PA      Exp     UTCDate/Time                  '
    manifest_header += 'SecZ    SlitMask       Dich      Gris/Grat      '
    manifest_header += 'Filter         Filename'

    manifest_text = '{}'.format(manifest_header)

    nullHeaderEntry = 'UNKNOWN'

    # glob all the files in the directory
    allFiles = glob.glob('{}/*.fits'.format(dataDir))
    allFiles.sort(key=lambda x: int(x.split('/')[-1].split('.')[0].replace('r','').replace('b','')))

    bFiles = []
    rFiles = []
    allSortedFiles = []
    for i,obsFile in enumerate(allFiles):
        if obsFile.split('/')[-1][0].lower() == 'b':
            bFiles.append(obsFile)
        else:
            rFiles.append(obsFile)
    allSortedFiles = bFiles + rFiles

    STANDARD_STAR_LIBRARY = construct_standard_star_library() # dict of std star objects

    #loop over the sorted list
    for i,obsFile in enumerate(allSortedFiles):

        print('Working on {}'.format(obsFile))

        # open the header
        hdu = fits.open(obsFile)
        header = hdu[0].header

        # get the image type
        try:
            obsImgType = determine_image_type(header,'KAST',STANDARD_STAR_LIBRARY)

        # catch all general exceptions; assign to CAL
        except Exception as e:
            obsImgType = 'CAL'


        # get blue or red channel
        obsVersion = header.get('VERSION',nullHeaderEntry)
        if obsVersion.strip().lower() == 'kastb':
            obsChannel = 'BLUE'
        elif obsVersion.strip().lower() == 'kastr':
            obsChannel = 'RED'
        else:
            obsChannel = nullHeaderEntry

        #get blue or red filter
        if obsChannel == 'BLUE':
            obsFilter = header.get('BLFILT_N',nullHeaderEntry).strip().replace(' ','_')
        elif obsChannel == 'RED':
            obsFilter = header.get('RDFILT_N',nullHeaderEntry).strip().replace(' ','_')
        else:
            obsFilter = nullHeaderEntry

        # get the explicitly named/non-degenerate parameters
        obsDate = header.get('DATE-OBS',nullHeaderEntry)
        obsObject = header.get('OBJECT',nullHeaderEntry).strip().replace(' ','_')
        obsRA = header.get('RA',nullHeaderEntry)
        obsDec = header.get('DEC',nullHeaderEntry)
        obsPositionAngle = header.get('TUB',nullHeaderEntry) # wow
        obsAirmass = header.get('AIRMASS',nullHeaderEntry)
        obsExptime = header.get('EXPTIME',nullHeaderEntry)
        obsGrating = header.get('GRATNG_N',nullHeaderEntry).strip().replace(' ','_')
        obsSlitmask = header.get('SLIT_N',nullHeaderEntry).strip().replace(' ','_')
        obsDichroic = header.get('BSPLIT_N',nullHeaderEntry).strip().replace(' ','_')

        #formatting tweaks to values
        if obsDec[0] != '-':
            obsDec = '+{}'.format(obsDec)
        if len(obsDec.split(':')[0]) < 3:
            obsDecDD = obsDec.split(':')[0]
            obsDecMM = obsDec.split(':')[1]
            obsDecSS = obsDec.split(':')[2]
            obsDecDD = '{}0{}'.format(obsDecDD[0],obsDecDD[1])
            obsDec = '{}:{}:{}'.format(obsDecDD,obsDecMM,obsDecSS)

        # prevent empty object names from getting thru
        if obsObject == '':
            obsObject = nullHeaderEntry

        # construct the string
        manifest_text += '\n{:<25} '.format(obsObject)
        manifest_text += '{:<10}'.format(obsImgType)
        manifest_text += '{:<5}'.format(obsChannel)
        manifest_text += '{:<15}{:<15}'.format(obsRA,obsDec)
        if obsPositionAngle != nullHeaderEntry:
            manifest_text += '{:<8.1f}'.format(obsPositionAngle)
        else:
            manifest_text += '{:<8}'.format(nullHeaderEntry)
        if obsExptime != nullHeaderEntry:
            manifest_text += '{:<8.1f}'.format(obsExptime)
        else:
            manifest_text += '{:<8}'.format(obsExptime)
        manifest_text += '{:<30}'.format(obsDate)
        if obsAirmass != nullHeaderEntry:
            manifest_text += '{:<8.2f}'.format(obsAirmass)
        else:
            manifest_text += '{:<8}'.format(nullHeaderEntry)
        manifest_text += '{:<15}'.format(obsSlitmask)
        manifest_text += '{:<10}'.format(obsDichroic)
        manifest_text += '{:<15}'.format(obsGrating)
        manifest_text += '{:<15}'.format(obsFilter)
        manifest_text += '{:<20}'.format(obsFile.split('/')[-1])

    return manifest_text


def determine_image_type(header,instrument,STANDARD_STAR_LIBRARY):

    # init
    nullHeaderEntry = 'UNKNOWN'

    if instrument == 'KAST':
        ARC_LAMP_KEYS = ['LAMPSTAH','LAMPSTAG','LAMPSTAF','LAMPSTAE',
                         'LAMPSTAD','LAMPSTAC','LAMPSTAB','LAMPSTAA',
                         'LAMPSTAK','LAMPSTAJ','LAMPSTAI']
        FLAT_LAMP_KEYS = ['LAMPSTA5','LAMPSTA4','LAMPSTA3',
                          'LAMPSTA2','LAMPSTA1']
        expTimeKey = 'EXPTIME'

    elif instrument == 'LRIS':
        ARC_LAMP_KEYS = ['MERCURY','NEON','ARGON','CADMIUM',
                         'ZINC','HALOGEN','KRYPTON','XENON',
                         'FEARGON','DEUTERI']
        FLAT_LAMP_KEYS = ['FLAMP1','FLAMP2']
        expTimeKey = 'TTIME'

    else:
        raise ValueError('Instrument {} not supported'.format(instrument))

    pointingTolerance = 60. # arcseconds, for matching against standards
    imageType = ''

    # if lamps are on, it is a calibration image
    # even if they're errantly left on we don't want that labeled SCI
    foc = header.get('OBJECT',nullHeaderEntry).strip().lower()

    for i,key in enumerate(ARC_LAMP_KEYS):
        if header.get(key,nullHeaderEntry).strip().lower() == 'on' and foc != 'focus loop':
            imageType = 'CAL_ARC'
            return imageType
        elif header.get(key,nullHeaderEntry).strip().lower() == 'on' and foc == 'focus loop':
            imageType = 'CAL'
            return imageType
    for i,key in enumerate(FLAT_LAMP_KEYS):
        if header.get(key,nullHeaderEntry).strip().lower() == 'on':
            imageType = 'CAL_FLAT'
            return imageType

    # suppose the lamps were left off or something, but we're not taking
    # on sky exposures, default to general CAL type
    if header.get('ROTMODE',nullHeaderEntry).strip().lower() == 'stationary':
        imageType = 'CAL'
        return imageType

    # if the expTime is <1, its a calibration, but hasn't been caught, default to CAL
    ttime = header.get(expTimeKey,nullHeaderEntry)
    if ttime == nullHeaderEntry or float(ttime) < 1.:
        imageType = 'CAL'
        return imageType

    # if the coordinates are within 10" of a known standard, its probably a standard
    # in the weird case where a true science object is close to a standard, this
    # will require a by-hand fix, but that is rather unlikely.
    pointingCenterStr = '{} {}'.format(header.get('RA',nullHeaderEntry),
                                    header.get('DEC',nullHeaderEntry))
    pointingCenter = SkyCoord(pointingCenterStr, 
                              frame='icrs',
                              unit=(u.hourangle, u.deg))
    for i,standardKey in enumerate(STANDARD_STAR_LIBRARY):
        standard = STANDARD_STAR_LIBRARY[standardKey]
        sep = pointingCenter.separation(standard.coord)
        if sep.arcsecond < pointingTolerance:
            imageType = 'STD'
            return imageType


    # no lamps on, non zero exposure time, not near a standard:
    # it's probably a science frame
    imageType = 'SCI'
    return imageType

def construct_standard_star_library():
    ''' Library of standard stars '''

    ssl = {
        'bd174708': StandardStar(ra='22:11:31.38', dec='+18:05:34.2'),
        'bd262606': StandardStar(ra='14:49:02.36', dec='+25:42:09.1'),
        'bd284211': StandardStar(ra='21:51:11.02', dec='+28:51:50.4'),
        'bd332642': StandardStar(ra='15:51:59.89', dec='+32:56:54.3'),
        'feige15': StandardStar(ra='01:49:09.49', dec='+13:33:11.8'),
        'feige24': StandardStar(ra='02:35:07.59', dec='+03:43:56.8'),
        'feige25': StandardStar(ra='02:38:37.79', dec='+05:28:11.3'),
        'feige34': StandardStar(ra='10:39:36.74', dec='+43:06:09.2'),
        'feige56': StandardStar(ra='12:06:47.24', dec='+11:40:12.7'),
        'feige92': StandardStar(ra='14:11:31.88', dec='+50:07:04.1'),
        'feige98': StandardStar(ra='14:38:15.75', dec='+27:29:32.9'),
        'feige110':StandardStar(ra='23:19:58.4', dec='-05:09:56.2'),
        'g158100': StandardStar(ra='00:33:54', dec='-12:07:57'), ###
        'g191b2b': StandardStar(ra='05:05:30.62', dec='+52:49:51.9'),
        'gd71': StandardStar(ra='05:52:27.62', dec='+15:53:13.2'),
        'gd248': StandardStar(ra='23:26:07', dec='+16:00:21'), ###
        'hd19445': StandardStar(ra='03:08:25.59', dec='+26:19:51.4'),
        'hd84937':StandardStar(ra='09:48:56.1',dec='+13:44:39.3'),
        'hz43':  StandardStar(ra='13:16:21.85', dec='+29:05:55.4'),
        'hz44': StandardStar(ra='13:23:35.26', dec='+36:07:59.5'),
        'ltt1020':StandardStar(ra='01:54:50.27',dec='-27:28:35.7'),
        'ltt1788': StandardStar(ra='03:48:22.61', dec='-39:08:37.0'),
        'ltt2415': StandardStar(ra='05:56:24.74', dec='-27:51:32.4'),
        'ltt3218': StandardStar(ra='08:41:32.43', dec='-32:56:32.9'),
        'ltt3864': StandardStar(ra='10:32:13.62', dec='-35:37:41.7'),
        'ltt4364': StandardStar(ra='11:45:42.92', dec='-64:50:29.5')
        }
    return ssl


class StandardStar():
    def __init__(self,ra='',dec=''):
        self.coord = SkyCoord('{} {}'.format(ra,dec), 
                                frame='icrs',
                                unit=(u.hourangle, u.deg))