from __future__ import print_function
import os,sys,pdb,shutil,glob,argparse,subprocess,shlex
from time import sleep
import numpy as np

from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u

import scipy
from scipy import signal
from scipy import interpolate
from scipy import optimize
from scipy import signal, ndimage

from pyraf import iraf
import pyds9 as pyds9

from flat_utils import combine_flats,inspect_flat

iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.longslit(_doprint=0)
iraf.onedspec(_doprint=0)
iraf.specred(_doprint=0)



class StandardStar():
    ''' A class representing our standard star targets '''
    def __init__(self,ra='',dec=''):
        self.coord = SkyCoord('{} {}'.format(ra,dec), 
                                frame='icrs',
                                unit=(u.hourangle, u.deg))

def show_ds9_list(listFile,instanceName='default'):
    ''' display a list of images in a DS9 instance '''
    
    fileNameArr = np.array([])
    
    # read in the filenames
    with open(listFile,'r') as fin:
        for line in fin:
            if line[0] != '#':
                fileName = line.split()[0]
                fileNameArr = np.append(fileNameArr,fileName)
    
    #Setup the DS9 instance
    disp = pyds9.DS9(instanceName)
    disp.set("frame delete all")
    disp.set("view image no")
    disp.set("view colorbar no")
    disp.set("view panner no")
    disp.set("view info yes")
    disp.set("view magnifier no")
    disp.set("view buttons no")
    disp.set("tile yes")
    disp.set("tile column")
    disp.set("width 1200")
    disp.set("height 275")
      
    #Display the images
    for i in xrange(len(fileNameArr)):
        disp.set("frame {}".format(i))
        ds9cmd = "file fits {}".format(fileNameArr[i])
        disp.set(ds9cmd)
        disp.set("scale mode minmax")
        disp.set("regions delete all")
        disp.set("scale log")
        disp.set("wcs align yes")
        disp.set("cmap invert yes")  
        disp.set(ds9cmd)
    
    return disp


def construct_standard_star_library():
    ''' Construct a library of standard stars '''
    ssl = {
        'BD174708': StandardStar(ra='22:11:31.38', dec='+18:05:34.2'),
        'BD262606': StandardStar(ra='14:49:02.36', dec='+25:42:09.1'),
        'BD284211': StandardStar(ra='21:51:11.02', dec='+28:51:50.4'),
        'BD332642': StandardStar(ra='15:51:59.89', dec='+32:56:54.3'),
        'FEIGE15': StandardStar(ra='01:49:09.49', dec='+13:33:11.8'),
        'FEIGE24': StandardStar(ra='02:35:07.59', dec='+03:43:56.8'),
        'FEIGE25': StandardStar(ra='02:38:37.79', dec='+05:28:11.3'),
        'FEIGE34': StandardStar(ra='10:39:36.74', dec='+43:06:09.2'),
        'FEIGE56': StandardStar(ra='12:06:47.24', dec='+11:40:12.7'),
        'FEIGE66': StandardStar(ra='12:37:23.52',dec='+25:03:59.9'),
        'FEIGE92': StandardStar(ra='14:11:31.88', dec='+50:07:04.1'),
        'FEIGE98': StandardStar(ra='14:38:15.75', dec='+27:29:32.9'),
        'Feige110':StandardStar(ra='23:19:58.4', dec='-05:09:56.2'),
        'G158100': StandardStar(ra='00:33:54', dec='-12:07:57'), ###
        'G191b2b': StandardStar(ra='05:05:30.62', dec='+52:49:51.9'),
        'GD71': StandardStar(ra='05:52:27.62', dec='+15:53:13.2'),
        'GD248': StandardStar(ra='23:26:07', dec='+16:00:21'), ###
        'HD19445': StandardStar(ra='03:08:25.59', dec='+26:19:51.4'),
        'HD84937':StandardStar(ra='09:48:56.1',dec='+13:44:39.3'),
        'HZ43':  StandardStar(ra='13:16:21.85', dec='+29:05:55.4'),
        'HZ44': StandardStar(ra='13:23:35.26', dec='+36:07:59.5'),
        'LTT1020':StandardStar(ra='01:54:50.27',dec='-27:28:35.7'),
        'LTT1788': StandardStar(ra='03:48:22.61', dec='-39:08:37.0'),
        'LTT2415': StandardStar(ra='05:56:24.74', dec='-27:51:32.4'),
        'LTT3218': StandardStar(ra='08:41:32.43', dec='-32:56:32.9'),
        'LTT3864': StandardStar(ra='10:32:13.62', dec='-35:37:41.7'),
        'LTT4364': StandardStar(ra='11:45:42.92', dec='-64:50:29.5')
        }
    return ssl


def get_image_channel(header):
    ''' Determine LRIS camera/channel '''
    chan = header.get('INSTRUME',None)
    if chan == 'LRISBLUE':
        return 'BLUE'
    else:
        return 'RED'


def determine_image_type(header,STANDARD_STAR_LIBRARY):
    ''' Determine image type from LRIS specific header key/values '''

    # init
    nullHeaderEntry = 'UNKNOWN'
    ARC_LAMP_KEYS = ['MERCURY','NEON','ARGON','CADMIUM',
                 'ZINC','HALOGEN','KRYPTON','XENON',
                 'FEARGON','DEUTERI']
    FLAT_LAMP_KEYS = ['FLAMP1','FLAMP2']
    pointingTolerance = 20. # arcseconds

    chan = get_image_channel(header)
    imageType = '{} '.format(chan)

    # if lamps are on, it is a calibration image

    # arcs
    for i,key in enumerate(ARC_LAMP_KEYS):
        if header.get(key,nullHeaderEntry).strip().lower() == 'on':
            imageType += 'ARC'
            return imageType

    # flats
    for i,key in enumerate(FLAT_LAMP_KEYS):
        if header.get(key,nullHeaderEntry).strip().lower() == 'on':
            imageType += 'FLAT'
            return imageType

    # suppose the lamps were left off or something, but we're not taking
    # tracking exposures. Would this capture twilight flats? Probably not.
    if header.get('ROTMODE',nullHeaderEntry).strip().lower() == 'stationary':
        imageType += 'CAL'
        return imageType

    # if the expTime is <1, its a calibration, but of unknown type
    ttime = header.get('TTIME',nullHeaderEntry)
    if ttime == nullHeaderEntry or float(ttime) < 1.:
        imageType += 'CAL'
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
            imageType += 'STD'
            return imageType

    # no lamps on, non zero exposure time, not near a standard:
    # it's probably a science frame
    imageType += 'SCI'
    return imageType


def stack_cal_frames(file_list,outfile='stacked.fits',CLOBBER=False,DISPAXIS=1):
    ''' Stack calibration frames '''
    stackCalFrame = np.array([])
    commentStr = 'keck_prep: stacked files '
    header = {}
    with open(file_list,'r') as fin:
        for line in fin:
            if line[0] != '#' and len(line.split()) > 1:
                infile = line.split()[0]
                hdu = fits.open(infile)
                data = hdu[0].data
                header = hdu[0].header
                expTime = header['ELAPTIME']
                commentStr += '{}'.format(infile)
                if len(stackCalFrame) == 0:
                    stackCalFrame = 1.*data/expTime
                else:
                    stackCalFrame = stackCalFrame + data/expTime

    # adopt the most recent header as a template and add commentStr
    if len(header) > 0:
        header.add_comment(commentStr)
        header['DISPAXIS'] = DISPAXIS

    if len(stackCalFrame) > 0:
        hduOut = fits.PrimaryHDU(stackCalFrame,header)
        if os.path.isfile(outfile):
            if CLOBBER:
                os.remove(outfile)
                hduOut.writeto(outfile)
        else:
            hduOut.writeto(outfile)
    return 0

def reorg_file_list(file_list,obj_list,CLOBBER=False,FULL_CLEAN=False):
    ''' Move files into their appropriate subdirectories '''
    for i,obj in enumerate(obj_list):
    
        if not os.path.isdir(obj):
            os.mkdir(obj)
            
        elif FULL_CLEAN:
            promptStr = 'Do you really want to wipe the dir for {}? [y/n]: '.format(obj)
            usrRespOrig = raw_input(promptStr)
            if usrRespOrig and usrRespOrig[0].strip().upper() == 'Y':
                shutil.rmtree(obj)
                os.mkdir(obj)
        else:
            pass
        
        # now move the appropriate files
        with open(file_list,'r') as fin:
            for line in fin:
                if line[0] != '#' and len(line.split()) > 1:
                    if obj.strip().upper() in line.split()[1].strip().upper():
                        destFile = '{}/{}'.format(obj,line.split()[0])
                        if not os.path.isfile(destFile) or CLOBBER:
                            hdu = fits.open(line.split()[0])
                            header = hdu[0].header
                            imgData = hdu[0].data

                            if CLOBBER:
                                print('Clobbering {}'.format(destFile))

                            # place the file
                            hduOut = fits.PrimaryHDU(imgData,header)
                            hduOut.writeto(destFile,output_verify='ignore')

        # if no files got moved in, remove the dir
        if not os.listdir(obj):
            shutil.rmtree(obj)
    return 0

def add_boolean_arg(parser,name,default=False,help_string=''):
    ''' Add a boolean arg to command line options '''
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--' + name, dest=name, action='store_true',
                        help=help_string)
    group.add_argument('--no_' + name, dest=name, action='store_false',
                        help='DO NOT: {}'.format(help_string))
    parser.set_defaults(**{name:default})


def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Main driver for organizing Keck LRIS data for reduction. '
    #descStr += 'Recommended calling sequence: \n \n'
    #descStr += '$ python keck_basic_2d.py -v -c \n'
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # required args
    #parser.add_argument('requried_arg',type=str,
    #                    help='a required arguement')

    # optional
    parser.add_argument('-v','--verbose',
                        help='print diagnostic info',action='store_true')
    parser.add_argument('-c','--clobber',action='store_true',
                        help='Clobber files already in pre_reduced/ but not subdirs')
    parser.add_argument('-f','--full-clean',action='store_true',
                        help='Do a complete wipe of pre_reduced, including subdirs')
    parser.add_argument('-i', '--inspect', action='store_true',
                        help='Inspect the frames in the file lists (launches DS9 instances)')

    parser.add_argument('--regenerate-all', action='store_true',
                        help='Regenerate allList.txt file list')
    parser.add_argument('--regenerate-arc', action='store_true',
                        help='Regenerate (blue & red)ArcList.txt file lists')
    parser.add_argument('--regenerate-flat', action='store_true',
                        help='Regenerate (blue & red)FlatList.txt file lists')
    parser.add_argument('--regenerate-std', action='store_true',
                        help='Regenerate (blue & red)StdList.txt file lists')
    parser.add_argument('--regenerate-sci', action='store_true',
                        help='Regenerate (blue & red)SciList.txt file lists')


    # mutually exclusives boolean flags
    add_boolean_arg(parser,'reorganize-std', default=True,
                    help_string='Move standard star files into their directories')
    add_boolean_arg(parser,'reorganize-sci', default=True,
                    help_string='Move science files into their directories')
    add_boolean_arg(parser,'stack-cals', default=True,
                    help_string='Stack calibration frames (arcs and flats)')

    # parse
    cmdArgs = parser.parse_args()

    # logic mapping to my args/kwargs
    VERBOSE = cmdArgs.verbose
    CLOBBER = cmdArgs.clobber
    FULL_CLEAN = cmdArgs.full_clean
    REGENERATE_ALL_LIST = cmdArgs.regenerate_all
    REGENERATE_ARC_LIST = cmdArgs.regenerate_arc
    REGENERATE_FLAT_LIST = cmdArgs.regenerate_flat
    REGENERATE_STD_LIST = cmdArgs.regenerate_std
    REGENERATE_SCI_LIST = cmdArgs.regenerate_sci
    REORG_STANDARDS = cmdArgs.reorganize_std
    REORG_SCIENCE = cmdArgs.reorganize_sci
    STACK_CAL_FRAMES = cmdArgs.stack_cals

    # package up
    args = () # no args implemented yet
    kwargs = {}
    kwargs['VERBOSE'] = VERBOSE
    kwargs['CLOBBER'] = CLOBBER
    kwargs['FULL_CLEAN'] = FULL_CLEAN
    kwargs['REGENERATE_ALL_LIST'] = REGENERATE_ALL_LIST
    kwargs['REGENERATE_ARC_LIST'] = REGENERATE_ARC_LIST
    kwargs['REGENERATE_FLAT_LIST'] = REGENERATE_FLAT_LIST
    kwargs['REGENERATE_STD_LIST'] = REGENERATE_STD_LIST
    kwargs['REGENERATE_SCI_LIST'] = REGENERATE_SCI_LIST
    kwargs['REORG_STANDARDS'] = REORG_STANDARDS
    kwargs['REORG_SCIENCE'] = REORG_SCIENCE

    kwargs['STACK_CAL_FRAMES'] = STACK_CAL_FRAMES

    return (args,kwargs)




def main(*args,**kwargs):
    '''
    Run basic 2D CCD reduction on Keck LRIS data

    Parameters
    ----------
    CLOBBER : bool, optional (default=False)
        Overwrite the individual files in pre_reduced, but
        do not wipe subdirectories
    FULL_CLEAN : bool, optional (default=False)
        Completely wipe pre_reduced and all subdirectories
    REGENERATE_ARC_LIST : bool, optional (default=False)
        Regenerate the list of arc files used for wavelength calibration
    REGENERATE_FLAT_LIST : bool, optional (default=False)
        Regenerate the list of flat field files used for flat fielding
    REGENERATE_SCI_LIST : bool, optional (default=False)
        Regenerate the list of science files to be processed
    REGENERATE_STD_LIST : bool, optional (default=False)
        Regenerate the list of standard star files to be processed
    REGENERATE_ALL_LIST : bool, optional (default=False)
        Regenerate the master list of all files in the working directory
    INSPECT_FRAMES : bool, optional (default=False)
        Automatically display (via DS9) the images in the file lists
    STACK_CAL_FRAMES : bool, optional (default=False)
        Stack the calibration frames
    REORG_STANDARDS : bool, optional (default=False)
        Move the standard star files into their respective directories
    REORG_SCIENCE : bool, optional (default=False)
        Move the science files into their respective directories
    TRIM : bool, optional (default=True)
        Trim to some hard coded section of the detector
    MASK_MIDDLE_BLUE : bool (default=False)
        Mask the middle section of the blue images
    MASK_MIDDLE_RED : bool (default=False)
        Mask the middle section of the red images (useful if
        there's wildly disparate values that make the iraf
        windowing tedious)
    
    Returns
    -------
    int : 0, and writes files to disk
    '''


    
    CLOBBER = kwargs.get('CLOBBER',False)
    FULL_CLEAN = kwargs.get('FULL_CLEAN',False)
    REGENERATE_ARC_LIST = kwargs.get('REGENERATE_ARC_LIST',False)
    REGENERATE_FLAT_LIST = kwargs.get('REGENERATE_FLAT_LIST',False)
    REGENERATE_SCI_LIST = kwargs.get('REGENERATE_SCI_LIST',False)
    REGENERATE_STD_LIST = kwargs.get('REGENERATE_STD_LIST',False)
    REGENERATE_ALL_LIST = kwargs.get('REGENERATE_ALL_LIST',False)
    INSPECT_FRAMES = kwargs.get('INSPECT_FRAMES',False)
    STACK_CAL_FRAMES = kwargs.get('STACK_CAL_FRAMES',False)
    REORG_STANDARDS = kwargs.get('REORG_STANDARDS',False)
    REORG_SCIENCE = kwargs.get('REORG_SCIENCE',False)

        
    # This is a dictionary of standard star objects
    STANDARD_STAR_LIBRARY = construct_standard_star_library()

    manifest_file_list = ['blueArcList.txt','redArcList.txt',
                          'blueFlatList.txt','redFlatList.txt',
                          'blueStdList.txt','redStdList.txt',
                          'blueSciList.txt','redSciList.txt',
                          'allList.txt']
    
    cwd = os.getcwd()
    if cwd.split('/')[-1] != 'pre_reduced':
        outStr = 'Looks like you\'re in: \n'
        outStr += '{}\n'.format(cwd) 
        outStr += 'Are you sure you\'re in the right directory?'
        print(outStr)
        pdb.set_trace()

    if FULL_CLEAN:
        promptStr = 'Do you really want to wipe all dirs and txt files from pre_reduced? [y/n]: '
        usrRespOrig = raw_input(promptStr)
        if usrRespOrig and usrRespOrig[0].strip().upper() == 'Y':

            # remove all text files
            for root,dirs,filenames in os.walk('.'):

                if root == '.':
                    for directory in dirs:
                        # remove all subdirectories
                        shutil.rmtree(directory)
                        outStr = 'Removed {}'.format(directory)
                        print(outStr)
                    for filename in filenames:
                        if filename in manifest_file_list:
                            os.remove(filename)
                            outStr = 'Removed {}'.format(filename)
                            print(outStr)


            
    # empties
    SCI_OBJ_LIST = []

    blueArcList = np.array([])
    blueArcAux = np.array([])
    
    redArcList = np.array([])
    redArcAux = np.array([])
    
    blueFlatList = np.array([])
    blueFlatAux = np.array([])
    
    redFlatList = np.array([])
    redFlatAux = np.array([])
    
    blueSciList = np.array([])
    blueSciAux = np.array([])
    
    redSciList = np.array([])
    redSciAux = np.array([])
    
    blueStdList = np.array([])
    blueStdAux = np.array([])
    
    redStdList = np.array([])
    redStdAux = np.array([])
    
    allList = np.array([])
    allAux = np.array([])
            
    # get all the files
    allFiles = sorted(glob.glob('*.fits'))
    
    # parse the headers, adding to lists as you go
    for i in xrange(len(allFiles)):
        
        # read
        inFile = allFiles[i]
        hdu = fits.open(inFile)
        header = hdu[0].header
        
        # parse
        fileType = determine_image_type(header,STANDARD_STAR_LIBRARY)

        # if its a science file, track it
        target_name = header['targname']
        if (('SCI' in fileType.upper()) and (target_name not in SCI_OBJ_LIST)):
            SCI_OBJ_LIST.append(target_name)
                                     
        # get auxilliary header info (object, slitmask)
        auxStr = '{} '.format(header['targname'])
        auxStr += '{} '.format(header['slitname'])
        auxStr += '{} '.format(header['graname'])
        auxStr += '{} '.format(header['ttime'])
        
        # log for the full file list
        allList = np.append(allList,inFile)
        allAux = np.append(allAux,auxStr)
    
        # write to appropriate list
        if fileType.strip().upper() == 'BLUE ARC':
            blueArcList = np.append(blueArcList,inFile)
            blueArcAux = np.append(blueArcAux,auxStr)
        elif fileType.strip().upper() == 'RED ARC':
            redArcList = np.append(redArcList,inFile)
            redArcAux = np.append(redArcAux,auxStr)

        elif fileType.strip().upper() == 'BLUE FLAT':
            blueFlatList = np.append(blueFlatList,inFile)
            blueFlatAux = np.append(blueFlatAux,auxStr)
        elif fileType.strip().upper() == 'RED FLAT':
            redFlatList = np.append(redFlatList,inFile)
            redFlatAux = np.append(redFlatAux,auxStr)

        elif fileType.strip().upper() == 'BLUE SCI':
            blueSciList = np.append(blueSciList,inFile)
            blueSciAux = np.append(blueSciAux,auxStr)
        elif fileType.strip().upper() == 'RED SCI':
            redSciList = np.append(redSciList,inFile)
            redSciAux = np.append(redSciAux,auxStr)

        elif fileType.strip().upper() == 'BLUE STD':
            blueStdList = np.append(blueStdList,inFile)
            blueStdAux = np.append(blueStdAux,auxStr)
        elif fileType.strip().upper() == 'RED STD':
            redStdList = np.append(redStdList,inFile)
            redStdAux = np.append(redStdAux,auxStr)
        else:
            errStr = '{} file type unknown...'.format(inFile)
            print(errStr)
            
            pass
    
    
    
    # write the list files
    # if the rewrite flag was set OR if the files aren't there, then write out
    WRITE_ALL_LIST = REGENERATE_ALL_LIST or not os.path.isfile('allList.txt')
    WRITE_ARC_LISTS = (REGENERATE_ARC_LIST or 
                       not (os.path.isfile('blueArcList.txt') or 
                            os.path.isfile('redArcList.txt')))
    WRITE_FLAT_LISTS = (REGENERATE_FLAT_LIST or
                        not (os.path.isfile('blueFlatList.txt') or 
                             os.path.isfile('redFlatList.txt')))
    WRITE_STD_LISTS = (REGENERATE_STD_LIST or
                        not (os.path.isfile('blueStdList.txt') or 
                             os.path.isfile('redStdList.txt')))
    WRITE_SCI_LISTS = (REGENERATE_SCI_LIST or 
                        not (os.path.isfile('blueSciList.txt') or 
                             os.path.isfile('redSciList.txt')))

    # all files
    if WRITE_ALL_LIST:
        with open('allList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(allList)):
                outStr += '{} {}\n'.format(allList[i],allAux[i])
            fout.write(outStr)

    # arcs
    if WRITE_ARC_LISTS:
        with open('blueArcList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(blueArcList)):
                outStr += '{} {}\n'.format(blueArcList[i],blueArcAux[i])
            fout.write(outStr)
        with open('redArcList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(redArcList)):
                outStr += '{} {}\n'.format(redArcList[i],redArcAux[i])
            fout.write(outStr)
        
    # flats
    if WRITE_FLAT_LISTS:
        with open('blueFlatList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(blueFlatList)):
                outStr += '{} {}\n'.format(blueFlatList[i],blueFlatAux[i])
            fout.write(outStr)
        with open('redFlatList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(redFlatList)):
                outStr += '{} {}\n'.format(redFlatList[i],redFlatAux[i])
            fout.write(outStr)

    # std stars
    if WRITE_STD_LISTS:
        with open('blueStdList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(blueStdList)):
                outStr += '{} {}\n'.format(blueStdList[i],blueStdAux[i])
            fout.write(outStr)
        with open('redStdList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(redStdList)):
                outStr += '{} {}\n'.format(redStdList[i],redStdAux[i])
            fout.write(outStr)
        
    # science
    if WRITE_SCI_LISTS:
        with open('blueSciList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(blueSciList)):
                outStr += '{} {}\n'.format(blueSciList[i],blueSciAux[i])
            fout.write(outStr)
        with open('redSciList.txt','w') as fout:
            outStr = ''
            for i in xrange(len(redSciList)):
                outStr += '{} {}\n'.format(redSciList[i],redSciAux[i])
            fout.write(outStr)
    


    ##### User should do some display/QA/list editing maybe #####
    usrResp = ''
    blueArcDS9 = ''
    redArcDS9 = ''
    blueFlatDS9 = ''
    redFlatDS9 = ''
    blueSciDS9 = ''
    redSciDS9 = ''
    blueStdDS9 = ''
    redStdDS9 = ''
    
    # do visual inspection of frames via ds9 windows
    while usrResp != 'C':
        promptStr = '\nLists are written. First, read the lists, remove '
        promptStr += 'files from the lists you wish to exclude. \n'
        promptStr += 'Then you may: \n'
        promptStr += '  (D)isplay the remaining images, or \n'
        promptStr += '  (C)ontinue with the lists as they are, or \n'
        promptStr += '  (Q)uit the whole thing. \n'
        promptStr += 'I recommend you display, inspect, '
        promptStr += 'remove unwanted frames from lists, then continue.\nCommand: '
        usrRespOrig = raw_input(promptStr)
    
        try:
            usrResp = usrRespOrig.strip().upper()
        except Exception as e:
            usrResp = 'nothing'

        if usrResp == 'Q':
            sys.exit(1)
        
        # display all the images in the lists
        if usrResp == 'D':
            blueArcDS9 = show_ds9_list('blueArcList.txt',instanceName='BlueArcs')
            redArcDS9 = show_ds9_list('redArcList.txt',instanceName='RedArcs')

            blueFlatDS9 = show_ds9_list('blueFlatList.txt',instanceName='BlueFlats')
            redFlatDS9 = show_ds9_list('redFlatList.txt',instanceName='RedFlats')

            blueSciDS9 = show_ds9_list('blueSciList.txt',instanceName='BlueScience')
            redSciDS9 = show_ds9_list('redSciList.txt',instanceName='RedScience')

            blueStdDS9 = show_ds9_list('blueStdList.txt',instanceName='BlueStandards')
            redStdDS9 = show_ds9_list('redStdList.txt',instanceName='RedStandards')        
    
    # the lists are set, now construct the master cal frames
    if STACK_CAL_FRAMES:

        # combine arcs into master arcs
        res = stack_cal_frames('blueArcList.txt',outfile='ARC_blue.fits',CLOBBER=CLOBBER)
        res = stack_cal_frames('redArcList.txt',outfile='ARC_red.fits',CLOBBER=CLOBBER)

        # combine flats into master flats
        res = combine_flats('blueFlatList.txt',READ_FROM_LIST=True,
                            OUTFILE='FLAT_blue.fits',MEDIAN_COMBINE=True)
        res = combine_flats('redFlatList.txt',READ_FROM_LIST=True,
                            OUTFILE='FLAT_red.fits',MEDIAN_COMBINE=True)

        # now do the response curve if the file exists
        if os.path.isfile('FLAT_blue.fits'):
            iraf.specred.response('FLAT_blue.fits', 
                                   normaliz='FLAT_blue.fits', 
                                   response='RESP_blue', 
                                   interac='y', thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=3,high_rej=3, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')
            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['RESP_blue.fits'])

        # now do the response curve if the file exists
        if os.path.isfile('FLAT_red.fits'):
            iraf.specred.response('FLAT_red.fits', 
                                   normaliz='FLAT_red.fits', 
                                   response='RESP_red', 
                                   interac='y', thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=3,high_rej=3, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['RESP_red.fits'])
            
    # calibration frames are all set, now move the std/sci files to their directory
    if REORG_STANDARDS:
        res = reorg_file_list('blueStdList.txt',STANDARD_STAR_LIBRARY,
                                CLOBBER=CLOBBER,FULL_CLEAN=FULL_CLEAN)
        res = reorg_file_list('redStdList.txt',STANDARD_STAR_LIBRARY,
                                CLOBBER=CLOBBER,FULL_CLEAN=FULL_CLEAN)
        
        
    # standards are all set, now move the science files to their directory
    if REORG_SCIENCE:
        res = reorg_file_list('blueSciList.txt',SCI_OBJ_LIST,
                                CLOBBER=CLOBBER,FULL_CLEAN=FULL_CLEAN)
        res = reorg_file_list('redSciList.txt',SCI_OBJ_LIST,
                                CLOBBER=CLOBBER,FULL_CLEAN=FULL_CLEAN)
    
    return 0
    
if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)