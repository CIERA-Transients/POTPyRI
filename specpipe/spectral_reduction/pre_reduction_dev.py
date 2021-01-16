#!/usr/bin/env python
from __future__ import print_function
import sys, os, shutil, pdb, argparse, datetime, json
from optparse import OptionParser
import util
import quick_reduc
import time
import glob
import matplotlib
import instruments
from astropy.io import fits, ascii
from pyraf import iraf
import pyds9 as pyds9
import keck_basic_2d
import manifest_utils as mu
import host_galaxies as host_gals
import numpy as np

matplotlib.use('TkAgg')
from flat_utils import combine_flats,inspect_flat


def show_ds9_listfile(listFile,instanceName='default'):
    ''' reads list of images from a file, displays in a DS9 instance '''
    
    imgList = []
    
    # read in the filenames
    with open(listFile,'r') as fin:
        for line in fin:
            if line[0] != '#':
                fileName = line.split()[0]
                imgList.append(fileName)
    disp = show_ds9_list(imgList,instanceName=instanceName)
    return disp

def show_ds9_list(imgList,instanceName='default'):
    ''' display a list of images in a DS9 instance '''
    
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
    for i,img in enumerate(imgList):
        disp.set("frame {}".format(i))
        ds9cmd = "file fits {}".format(img)
        disp.set(ds9cmd)
        disp.set("scale mode minmax")
        disp.set("regions delete all")
        disp.set("scale log")
        disp.set("cmap invert yes")  
        disp.set(ds9cmd)
    return disp


def config_to_pretty_str(configDict):
    ''' Dummy function for reformatting configDicts '''
    return json.dumps(configDict,indent=4)



def user_adjust_config(configDict,operation=''):
    ''' Allow user to make adjustments to their reduction config interactively '''

    if operation == 'REMOVE':

        usrResp = ''
        while usrResp.upper() not in ['C','Q']:
            configStr = config_to_pretty_str(configDict)

            promptStr = '\n\n\n You\'ve chosen to remove some file(s) from the config.\n'
            promptStr += 'Here\'s the current state:\n'
            promptStr += configStr
            promptStr += '\n\nAt this point you may:\n\n'
            promptStr += '  Enter the filename you wish to remove, or\n'
            promptStr += '  (C)ontinue with these files as is\nCommand: '
            usrRespOrig = raw_input(promptStr)


            try:
                usrResp = usrRespOrig.strip()
            except Exception as e:
                errStr = 'I had some problem with that input: {}'.format(e)
                print(errStr)

            if usrResp.upper() not in ['C','Q']:
                try:
                    removeFilename = usrResp
                    REMOVED = False
                    for imgType,typeDict in configDict.items():
                        for chan,objDict in typeDict.items():
                            for obj,fileList in objDict.items():
                                if removeFilename in fileList:
                                    configDict[imgType][chan][obj].remove(removeFilename)
                                    REMOVED = True
                                    print('i removed {}'.format(removeFilename))
                    if not REMOVED:
                        raise ValueError('Couldn\'t locate {} in config'.format(removeFilename))

                except (Exception,ValueError) as e:
                    errStr = 'I had some problem: {}'.format(str(e))
                    print(errStr)



    if operation == 'ADD':

        usrResp = ''
        while usrResp.upper() not in ['C','Q']:
            configStr = config_to_pretty_str(configDict)

            promptStr = '\n\n\n You\'ve chosen to add some file(s) from the config.\n'
            promptStr += 'Here\'s the current state:\n'
            promptStr += configStr
            promptStr += '\n\nAt this point you may:\n\n'
            promptStr += '  (C)ontinue with these files as is, or\n'
            promptStr += '  Enter the filename you wish to add, \n'
            promptStr += '  according to the format TYPE CHAN OBJECT FILENAME\n\n'
            promptStr += '  (e.g., CAL_ARC BLUE CALIBRATION_ARC r1021.fits)\n'
            promptStr += '  (e.g., CAL_FLAT RED CALIBRATION_FLAT r1022.fits)\n'
            promptStr += '  (e.g., STD BLUE BD28411 b1054.fits)\n'
            promptStr += '  (e.g., SCI BLUE 2019ltw b1055.fits)\n\n'
            promptStr += '  Make sure your input is accurately formatted; this program\n'
            promptStr += '  performs conservative validation and may not add your\n'
            promptStr += '  file if the input is not formatted correctly. In that case\n'
            promptStr += '  quit and just add the file(s) the config file by hand.\nCommand: '
            usrRespOrig = raw_input(promptStr)


            try:
                usrResp = usrRespOrig.strip()
            except Exception as e:
                errStr = 'I had some problem with that input: {}'.format(e)
                print(errStr)

            if usrResp.upper() not in ['C','Q']:
                try:
                    ADDED = False
                    inputList = usrResp.split()
                    imgType = inputList[0]
                    chan = inputList[1]
                    obj = inputList[2]
                    filename = inputList[3]

                    # check if file is available
                    availableFiles = [os.path.basename(file) for file in glob.glob('*.fits')]
                    if filename not in availableFiles:
                        raise ValueError('Couldn\'t find file {}'.format(filename))

                    # try to add
                    if imgType not in ['CAL_ARC','CAL_FLAT','STD','SCI']:
                        raise ValueError('Bad image type supplied!')
                    if chan not in ['BLUE','RED']:
                        raise ValueError('Bad channel supplied!')

                    if obj not in configDict[imgType][chan].keys():
                        configDict[imgType][chan][obj] = [filename]
                        ADDED = True
                    else:
                        if filename not in configDict[imgType][chan][obj]:
                            configDict[imgType][chan][obj].append(filename)
                            ADDED = True

                    if not ADDED:
                        raise ValueError('Couldn\'t add {} to config'.format(filename))

                except Exception as e:
                    errStr = 'I had some problem: {}'.format(str(e))
                    print(errStr)


    return configDict


def basic_2d_proc(rawFile,imgType=None,CLOBBER=False):

    # set up file names based on our convention
    oScanFile = 'pre_reduced/o{}'.format(rawFile)
    toScanFile = 'pre_reduced/to{}'.format(rawFile)
    toScanFlat = 'pre_reduced/to{}'.format(rawFile.split('.')[0]+'_norm.fits')

    # run the basic 2D stuff on this image if necessary
    if (not os.path.isfile(oScanFile)) or CLOBBER:

        # get the instrument configuration
        inst = instruments.blue_or_red(rawFile)[1]
        iraf.specred.dispaxi = inst.get('dispaxis')
        iraf.longslit.dispaxi = inst.get('dispaxis')
        _biassec0 = inst.get('biassec')
        _trimsec0 = inst.get('trimsec')
        _flatsec0 = inst.get('flatsec')

        # remove destination files
        for file in [oScanFile,toScanFile]:
            if os.path.isfile(file) and CLOBBER:
                os.remove(file)

        # check the ultimate desination file, since intermediates get deleted
        if not os.path.isfile(toScanFile):
            if inst.get('observatory')=='lick':

                # do Lick specific bias operations
                util.kastbias(rawFile,oScanFile)

            # do general (IRAF is in such a sorry state I'm not even sure if this works)
            else:
                iraf.ccdproc(rawFile, output=oScanFile, 
                             overscan='yes', trim='yes', 
                             zerocor="no", flatcor="no",readaxi='line',
                             trimsec=str(_trimsec0), biassec=str(_biassec0), 
                             Stdout=1)

            # trim (same trimming operation for all telescopes)
            iraf.ccdproc(oScanFile, output=toScanFile, 
                        overscan='no', trim='yes', zerocor="no", flatcor="no", 
                        readaxi='line',trimsec=str(_trimsec0), Stdout=1)

            #create trimmed flats for norm region
            if imgType == 'CAL_FLAT' and 'kast' in inst.get('name'):
                iraf.ccdproc(oScanFile, output=toScanFlat, 
                        overscan='no', trim='yes', zerocor="no", flatcor="no", 
                        readaxi='line',trimsec=str(_flatsec0), Stdout=1)
            # clean up intermediate files

            # if 'kast_red' in inst.get('name'):
            #     iraf.imtranspose(toScanFile, output=toScanFile)
            #     iraf.wcsreset(toScanFile, wcs='physical')
            #     tfits=fits.open(toScanFile, mode='update')
            #     thead=tfits[0].header
            #     thead.set('DATASEC', '[80:2296, 66:346]')
            #     tfits.flush()
            os.remove(oScanFile)

    return 0



def reorg_files(configDict,CLOBBER=False):
    ''' Move files into their appropriate subdirectories '''
    
    for imgType,typeDict in configDict.items():
        if imgType in ['STD','SCI']:
            for chan,objDict in typeDict.items():
                for obj,fileList in objDict.items():

                    destDir = 'pre_reduced/{}'.format(obj)
                    if not os.path.isdir(destDir):
                        os.mkdir(destDir)

                    for rawFile in fileList:
                        procFile = 'pre_reduced/to{}'.format(rawFile)
                        destFile = '{}/{}'.format(destDir,os.path.basename(procFile))

                        if os.path.isfile(destFile) and CLOBBER:
                            print('Clobbering {}'.format(destFile))
                            os.remove(destFile)

                        if not os.path.isfile(destFile):
                            shutil.copy(procFile,destFile)

    return 0




def pre_reduction_dev(*args,**kwargs):

    # parse kwargs
    VERBOSE = kwargs.get('VERBOSE')
    CLOBBER = kwargs.get('CLOBBER')
    FAKE_BASIC_2D = kwargs.get('FAKE_BASIC_2D')
    FULL_CLEAN = kwargs.get('FULL_CLEAN')
    FAST = kwargs.get('FAST')
    CONFIG_FILE = kwargs.get('CONFIG_FILE')
    MAKE_ARCS = kwargs.get('MAKE_ARCS')
    MAKE_FLATS = kwargs.get('MAKE_FLATS')
    QUICK = kwargs.get('QUICK')
    RED_AMP_BAD = kwargs.get('RED_AMP_BAD')
    HOST = kwargs.get('HOST')

    # init iraf stuff
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)

    iraf.ccdred.verbose = 'no'
    iraf.specred.verbose = 'no'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''

    iraf.longslit.mode = 'h'
    iraf.specred.mode = 'h'
    iraf.noao.mode = 'h'
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"

    # set up config
    if CONFIG_FILE:
        with open(CONFIG_FILE,'r') as fin:
            configDict = json.load(fin)
    else:
        #TODO: Make a better first pass at config file, std have exptime < 250s(?)
        configDict = {
        'SCI': {'BLUE': {}, 
                'RED': {}
                },
        'STD': {'BLUE': {}, 
                'RED': {}
                },
        'CAL_ARC': {'BLUE': {'CALIBRATION_ARC':[]}, 
                    'RED': {'CALIBRATION_ARC':[]}
                },
        'CAL_FLAT': {'BLUE': {'CALIBRATION_FLAT':[]}, 
                     'RED': {'CALIBRATION_FLAT':[]}
                },
        # this last key is a garbage collector; not used
        'CAL': {'BLUE': {}, 
                'RED': {}
                }
        }


        STANDARD_STAR_LIBRARY = mu.construct_standard_star_library()
        observations = sorted(glob.glob('*.fits'))
        arm, inst_dict = instruments.blue_or_red(observations[0])
        if 'lris' in inst_dict.get('name'):
            inst_name = 'LRIS'
        elif 'kast' in inst_dict.get('name'):
            inst_name = 'KAST'
        elif 'goodman' in inst_dict.get('name'):
            inst_name = 'GOODMAN'

        for obsfile in observations:
            header = fits.open(obsfile)[0].header
            imageType = mu.determine_image_type(header, inst_name, STANDARD_STAR_LIBRARY)

            channel, inst_dict = instruments.blue_or_red(obsfile)
            obj = header.get('OBJECT').strip()
            if imageType == 'SCI' or imageType == 'STD':
                if obj in configDict[imageType][channel.upper()].keys():
                    configDict[imageType][channel.upper()][obj].append(obsfile)
                else:
                    configDict[imageType][channel.upper()][obj] = [obsfile]
            if imageType == 'CAL_ARC' and 'foc' not in obsfile:
                configDict[imageType][channel.upper()]['CALIBRATION_ARC'].append(obsfile)
            if imageType == 'CAL_FLAT':
                configDict[imageType][channel.upper()]['CALIBRATION_FLAT'].append(obsfile)


        with open('custom_config.json','w') as fout:
            fout.write(json.dumps(configDict,indent=4))

        outStr = '\n\nOk, not config supplied, so I wrote a first pass custom_config.json\n'
        outStr += 'Use at your own risk! Manually edit if needed and run again with -c custom_config.json\n'
        outStr += 'You can manually add files to the appropriate lists and rerun.\n'
        outStr += 'WARNING: make sure you rename your config file, or it could get overwritten!\n\n'
        print(outStr)

        sys.exit(1)

    if not FAST:
        # do visual inspection of frames via ds9 windows
        usrResp = ''
        while usrResp != 'C':
            promptStr = '\nYou\'ve opted to display images before kicking off the reduction.\n'
            promptStr += 'At this point you may:\n'
            promptStr += '  (D)isplay the current state of the reduction config\n'
            promptStr += '  (C)ontinue with these files as is\n'
            promptStr += '  (R)emove a file from the current config\n'
            promptStr += '  (A)dd a file to the current config\n'
            promptStr += '  (Q)uit the whole thing. \n'
            promptStr += 'I recommend you (D)isplay and remove unwanted frames from your config file,\n'
            promptStr += '(Q)uit, and then rerun with the updated config file.\nCommand: '
            usrRespOrig = raw_input(promptStr)


            try:
                usrResp = usrRespOrig.strip().upper()
            except Exception as e:
                usrResp = 'nothing'

            # (D)isplay all the images in the lists
            if usrResp == 'D':

                blueArcList = configDict['CAL_ARC']['BLUE']['CALIBRATION_ARC']
                redArcList = configDict['CAL_ARC']['RED']['CALIBRATION_ARC']

                blueFlatList = configDict['CAL_FLAT']['BLUE']['CALIBRATION_FLAT']
                redFlatList = configDict['CAL_FLAT']['RED']['CALIBRATION_FLAT']

                blueStdList = []
                redStdList = []

                blueSciList = []
                redSciList = []

                for targ,imgList in configDict['STD']['BLUE'].items():
                    for file in imgList:
                        blueStdList.append(file)
                for targ,imgList in configDict['STD']['RED'].items():
                    for file in imgList:
                        redStdList.append(file)
                for targ,imgList in configDict['SCI']['BLUE'].items():
                    for file in imgList:
                        blueSciList.append(file)
                for targ,imgList in configDict['SCI']['RED'].items():
                    for file in imgList:
                        redSciList.append(file)

                blueArcDS9 = show_ds9_list(blueArcList,instanceName='BlueArcs')
                redArcDS9 = show_ds9_list(redArcList,instanceName='RedArcs')

                blueFlatDS9 = show_ds9_list(blueFlatList,instanceName='BlueFlats')
                redFlatDS9 = show_ds9_list(redFlatList,instanceName='RedFlats')

                blueStdDS9 = show_ds9_list(blueStdList,instanceName='BlueStandards')
                redStdDS9 = show_ds9_list(redStdList,instanceName='RedStandards')

                blueSciDS9 = show_ds9_list(blueSciList,instanceName='BlueScience')
                redSciDS9 = show_ds9_list(redSciList,instanceName='RedScience')

            if usrResp == 'R':
                configDict = user_adjust_config(configDict,operation='REMOVE')
            if usrResp == 'A':
                configDict = user_adjust_config(configDict,operation='ADD')
            if usrResp == 'Q':
                print('Okay, quitting pre_reduction...')
                sys.exit(1)


    # pre_reduced does not exist, needs to be made
    if not os.path.isdir('pre_reduced/'):
        os.mkdir('pre_reduced/')

    if QUICK:
        file =glob.glob('*.fits')[0]
        inst = instruments.blue_or_red(file)[1]
        if 'kast' in inst['name']:
            b_inst = instruments.kast_blue
            r_inst = instruments.kast_red
        if 'lris' in inst['name']:
            b_inst = instruments.lris_blue
            r_inst = instruments.lris_red
        if 'goodman' in inst['name']:
            b_inst = instruments.goodman_blue
            r_inst = instruments.goodman_red
        if not os.path.isdir('pre_reduced/master_files/'):
            os.mkdir('pre_reduced/master_files/')
        b_arcsol = b_inst.get('archive_arc_extracted_id')
        b_resp = b_inst.get('archive_flat_file')
        r_arcsol = r_inst.get('archive_arc_extracted_id')
        r_resp = r_inst.get('archive_flat_file')
        if os.path.isdir('pre_reduced/master_files/'):
            os.system('cp ' + b_arcsol + ' ' + 'pre_reduced/master_files/')
            os.system('cp ' + b_resp + ' ' + 'pre_reduced/')
            os.system('cp ' + r_arcsol + ' ' + 'pre_reduced/master_files/')
            os.system('cp ' + r_resp + ' ' + 'pre_reduced/')


    # pre_reduced exists, but we want to clobber/do a clean reduction
    elif FULL_CLEAN:
        
        promptStr = 'Do you really want to wipe pre_reduced? [y/n]: '
        usrRespOrig = raw_input(promptStr)
        if usrRespOrig and usrRespOrig[0].strip().upper() == 'Y':

            # remove all pre_reduced files
            shutil.rmtree('pre_reduced')
            os.mkdir('pre_reduced/')
      
    # pre_reduced exists, need to document what is there
    else:
        
        # get existing pre_reduced files
        preRedFiles = glob.glob('pre_reduced/*.fits')

    # loop over raw files in configDict, if the destination exists, do nothing
    # # otherwise, do the bias/reorient/trim/output/etc

    MASK_MIDDLE_RED = True #may need to update

    for imgType,typeDict in configDict.items():
        for chan,objDict in typeDict.items():
            for obj,fileList in objDict.items():
                for rawFile in fileList:
                    # try:
                    #     res = basic_2d_proc(rawFile,CLOBBER=CLOBBER)
                    #     if res != 0:
                    #         raise ValueError('Something bad happened in basic_2d_proc on {}'.format(rawFile))
                    # except (Exception,ValueError) as e:
                    #     print('Exception (basic_2d): {}'.format(e))

                    if not FAKE_BASIC_2D:
                        inst = instruments.blue_or_red(rawFile)[1]
                        if inst['name'] == 'lris_blue' or inst['name'] == 'lris_red':
                            # res = keck_basic_2d.main([rawFile])
                            if imgType != 'CAL_FLAT':
                                print (imgType)
                                res = keck_basic_2d.main([rawFile], TRIM=True, ISDFLAT = False, RED_AMP_BAD=RED_AMP_BAD, MASK_MIDDLE_RED=MASK_MIDDLE_RED)
                            else:
                                print (imgType)
                                res = keck_basic_2d.main([rawFile], TRIM=True, ISDFLAT = True, RED_AMP_BAD=RED_AMP_BAD, MASK_MIDDLE_RED=MASK_MIDDLE_RED)
                        else:
                            res = basic_2d_proc(rawFile,imgType=imgType,CLOBBER=CLOBBER)
                    else:
                        # here we're faking the basic 2D reduction because we've done
                        # specialized 2D reduction (e.g., keck_basic_2d)
                        res = 0
                    if res != 0:
                        raise ValueError('Something bad happened in basic_2d_proc on {}'.format(rawFile))
                    
    # move the std and sci files into their appropriate directories
    try:
        res = reorg_files(configDict,CLOBBER=CLOBBER)
        if res != 0:
            raise ValueError('Something bad happened in reorg_files')
    except (Exception,ValueError) as e:
        print('Exception (reorg): {}'.format(e))


    ### some blocks of code from the original pre_reduction ###
    # combine the arcs
    if MAKE_ARCS:
        list_arc_b = configDict['CAL_ARC']['BLUE']['CALIBRATION_ARC']
        list_arc_r = configDict['CAL_ARC']['RED']['CALIBRATION_ARC']
        
        # blue arcs
        if len(list_arc_b) > 0:
            if len(list_arc_b) == 1:
                arc_blue = list_arc_b[0]
                originFile = 'pre_reduced/to{}'.format(arc_blue)
                destFile = 'pre_reduced/ARC_blue.fits'
                shutil.copy(originFile,destFile)
            else:
                arc_str = ''
                for arc in list_arc_b:
                    arc_str += 'pre_reduced/to{},'.format(arc)
                if os.path.isfile('pre_reduced/ARC_blue.fits'):
                    os.remove('pre_reduced/ARC_blue.fits')
                iraf.imcombine(arc_str, output='pre_reduced/ARC_blue.fits')

        # red arcs
        if len(list_arc_r) > 0:
            if len(list_arc_r) == 1:
                arc_red = list_arc_r[0]
                originFile = 'pre_reduced/to{}'.format(arc_red)
                destFile = 'pre_reduced/ARC_red.fits'
                shutil.copy(originFile,destFile)
            else:
                arc_str = ''
                for arc in list_arc_r:
                    arc_str += 'pre_reduced/to{},'.format(arc)
                if os.path.isfile('pre_reduced/ARC_red.fits'):
                    os.remove('pre_reduced/ARC_red.fits')
                iraf.imcombine(arc_str, output='pre_reduced/ARC_red.fits')


    # combine the flats
    if MAKE_FLATS and 'lris' in inst['name']:

        list_flat_b = configDict['CAL_FLAT']['BLUE']['CALIBRATION_FLAT']
        list_flat_r = configDict['CAL_FLAT']['RED']['CALIBRATION_FLAT']
        inter = 'yes'

        b_amp1_list = []
        b_amp2_list = []
        r_amp1_list = []
        r_amp2_list = []
        for flat in list_flat_b:
            b_amp1_list.append(flat.split('.')[0]+'_amp1.'+flat.split('.')[1])
            b_amp2_list.append(flat.split('.')[0]+'_amp2.'+flat.split('.')[1])
        for flat in list_flat_r:
            r_amp1_list.append(flat.split('.')[0]+'_amp1.'+flat.split('.')[1])
            if not RED_AMP_BAD:
                r_amp2_list.append(flat.split('.')[0]+'_amp2.'+flat.split('.')[1])


        # blue flats
        if len(list_flat_b) > 0:
            # br, inst = instruments.blue_or_red(list_flat_b[0])
            br, inst = instruments.blue_or_red('pre_reduced/to{}'.format(b_amp1_list[0]))
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_blue_amp1 = 'pre_reduced/toFlat_blue_amp1.fits'
            Flat_blue_amp2 = 'pre_reduced/toFlat_blue_amp2.fits'

            flat_list_amp1 = []
            for flat in b_amp1_list:
                flat_list_amp1.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_blue_amp1):
                os.remove(Flat_blue_amp1)

            flat_list_amp2 = []
            for flat in b_amp2_list:
                flat_list_amp2.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_blue_amp2):
                os.remove(Flat_blue_amp2)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list_amp1,OUTFILE=Flat_blue_amp1,MEDIAN_COMBINE=False)
            res = combine_flats(flat_list_amp2,OUTFILE=Flat_blue_amp2,MEDIAN_COMBINE=False)
            
            # run iraf response
            iraf.specred.response(Flat_blue_amp1, 
                                   normaliz=Flat_blue_amp1, 
                                   response='pre_reduced/RESP_blue_amp1', 
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=5,high_rej=5, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')
            iraf.specred.response(Flat_blue_amp2, 
                                   normaliz=Flat_blue_amp2, 
                                   response='pre_reduced/RESP_blue_amp2', 
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=5,high_rej=5, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_blue_amp1.fits'], OUTFILE='pre_reduced/RESP_blue_amp1.fits', DISPAXIS=dispaxis)
            res = inspect_flat(['pre_reduced/RESP_blue_amp2.fits'], OUTFILE='pre_reduced/RESP_blue_amp2.fits', DISPAXIS=dispaxis)

            hdu_amp1 = fits.open('pre_reduced/RESP_blue_amp1.fits')
            hdu_amp2 = fits.open('pre_reduced/RESP_blue_amp2.fits')
            amp1_flatten = np.asarray(hdu_amp1[0].data).flatten()
            amp2_flatten = np.asarray(hdu_amp2[0].data).flatten()
            concat_amps = np.concatenate([amp2_flatten, amp1_flatten])
            resp_blue_data = np.reshape(concat_amps, (1000,4096))

            header = hdu_amp1[0].header
            if os.path.isfile('pre_reduced/RESP_blue.fits'):
                os.remove('pre_reduced/RESP_blue.fits')

            hdu = fits.PrimaryHDU(resp_blue_data,header)
            hdu.writeto('pre_reduced/RESP_blue.fits',output_verify='ignore')

            os.remove('pre_reduced/RESP_blue_amp1.fits')
            os.remove('pre_reduced/RESP_blue_amp2.fits')

        # red flats
        if len(list_flat_r) > 0:
            # br, inst = instruments.blue_or_red(list_flat_r[0])
            br, inst = instruments.blue_or_red('pre_reduced/to{}'.format(r_amp1_list[0]))
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_red_amp1 = 'pre_reduced/toFlat_red_amp1.fits'
            Flat_red_amp2 = 'pre_reduced/toFlat_red_amp2.fits'


            flat_list_amp1 = []
            for flat in r_amp1_list:
                flat_list_amp1.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_red_amp1):
                os.remove(Flat_red_amp1)

            flat_list_amp2 = []
            for flat in r_amp2_list:
                flat_list_amp2.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_red_amp2):
                os.remove(Flat_red_amp2)

            # first, combine all the flat files into a master flat
            if not RED_AMP_BAD:
                res = combine_flats(flat_list_amp1,OUTFILE=Flat_red_amp1,MEDIAN_COMBINE=True)
                res = combine_flats(flat_list_amp2,OUTFILE=Flat_red_amp2,MEDIAN_COMBINE=True)
            else:
                res = combine_flats(flat_list_amp1,OUTFILE=Flat_red_amp1,MEDIAN_COMBINE=True)

            #What is the output here? Check for overwrite
            iraf.specred.response(Flat_red_amp1, 
                                  normaliz=Flat_red_amp1, 
                                  response='pre_reduced/RESP_red_amp1', 
                                  interac=inter, thresho='INDEF',
                                  sample='*', naverage=2, function='legendre', 
                                  low_rej=5,high_rej=5, order=80, niterat=20, 
                                  grow=0, graphic='stdgraph')
            if not RED_AMP_BAD:
                iraf.specred.response(Flat_red_amp2, 
                                      normaliz=Flat_red_amp2, 
                                      response='pre_reduced/RESP_red_amp2', 
                                      interac=inter, thresho='INDEF',
                                      sample='*', naverage=2, function='legendre', 
                                      low_rej=5,high_rej=5, order=80, niterat=20, 
                                      grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            if not RED_AMP_BAD:
                res = inspect_flat(['pre_reduced/RESP_red_amp1.fits'], OUTFILE='pre_reduced/RESP_red_amp1.fits', DISPAXIS=dispaxis)
                res = inspect_flat(['pre_reduced/RESP_red_amp2.fits'], OUTFILE='pre_reduced/RESP_red_amp2.fits', DISPAXIS=dispaxis)
            else:
                res = inspect_flat(['pre_reduced/RESP_red_amp1.fits'], OUTFILE='pre_reduced/RESP_red.fits', DISPAXIS=dispaxis)

            if not RED_AMP_BAD:
                hdu_amp1 = fits.open('pre_reduced/RESP_red_amp1.fits')
                hdu_amp2 = fits.open('pre_reduced/RESP_red_amp2.fits')
                amp1_flatten = np.asarray(hdu_amp1[0].data).flatten()
                amp2_flatten = np.asarray(hdu_amp2[0].data).flatten()
                concat_amps = np.concatenate([amp2_flatten, amp1_flatten])
                if not MASK_MIDDLE_RED:
                    resp_red_data = np.reshape(concat_amps, (575,4061))
                    resp_red_data[278:294,:] = 1.
                else:
                    resp_red_data = np.reshape(concat_amps, (575,4061)) #depends on num amps? (500 for 4)

                header = hdu_amp1[0].header
                if os.path.isfile('pre_reduced/RESP_red.fits'):
                    os.remove('pre_reduced/RESP_red.fits')

                hdu = fits.PrimaryHDU(resp_red_data,header)
                hdu.writeto('pre_reduced/RESP_red.fits',output_verify='ignore')

                os.remove('pre_reduced/RESP_red_amp1.fits')
                os.remove('pre_reduced/RESP_red_amp2.fits')
            else:
                os.remove('pre_reduced/RESP_red_amp1.fits')

    elif MAKE_FLATS:
        list_flat_b = configDict['CAL_FLAT']['BLUE']['CALIBRATION_FLAT']
        list_flat_r = configDict['CAL_FLAT']['RED']['CALIBRATION_FLAT']
        inter = 'yes'

        # blue flats
        if len(list_flat_b) > 0:
            # br, inst = instruments.blue_or_red(list_flat_b[0])
            br, inst = instruments.blue_or_red('pre_reduced/to{}'.format(list_flat_b[0]))
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_blue = 'pre_reduced/toFlat_blue.fits'

            flat_list = []
            norm_list = []
            for flat in list_flat_b:
                flat_list.append('pre_reduced/to'+ flat)
                norm_list.append('pre_reduced/to'+ flat.split('.')[0]+'_norm.fits')
            if os.path.isfile(Flat_blue):
                os.remove(Flat_blue)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list,OUTFILE=Flat_blue,MEDIAN_COMBINE=True)
            # combine all the flat files for norm region
            res = combine_flats(norm_list,OUTFILE='pre_reduced/dummy_blue.fits',MEDIAN_COMBINE=True)
        
            #can't get this to work:
            # _flatsec0 = inst.get('flatsec')
            # iraf.ccdproc(Flat_blue, output='pre_reduced/dummy_blue.fits', 
            #         overscan='no', trim='yes', zerocor="no", flatcor="no", 
            #         readaxi='line',trimsec=str(_flatsec0), Stdout=1)

            # run iraf response
            iraf.specred.response(Flat_blue, 
                                   normaliz='pre_reduced/dummy_blue.fits', 
                                   response='pre_reduced/RESP_blue', 
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=5,high_rej=5, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_blue.fits'], OUTFILE='pre_reduced/RESP_blue.fits', DISPAXIS=dispaxis)

            for flat in list_flat_b:
                os.remove('pre_reduced/to'+ flat.split('.')[0]+'_norm.fits')

        # red flats
        if len(list_flat_r) > 0:
            # br, inst = instruments.blue_or_red(list_flat_r[0])
            br, inst = instruments.blue_or_red('pre_reduced/to{}'.format(list_flat_r[0]))
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_red = 'pre_reduced/toFlat_red.fits'

            flat_list= []
            norm_list = []
            for flat in list_flat_r:
                flat_list.append('pre_reduced/to'+ flat)
                norm_list.append('pre_reduced/to'+ flat.split('.')[0]+'_norm.fits')
            if os.path.isfile(Flat_red):
                os.remove(Flat_red)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list,OUTFILE=Flat_red,MEDIAN_COMBINE=True)
            # combine all the flat files for norm region
            res = combine_flats(norm_list,OUTFILE='pre_reduced/dummy_red.fits',MEDIAN_COMBINE=True)

            
            # run iraf response
            iraf.specred.response(Flat_red, 
                                   normaliz='pre_reduced/dummy_red.fits', 
                                   response='pre_reduced/RESP_red', 
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=5,high_rej=5, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_red.fits'], OUTFILE='pre_reduced/RESP_red.fits', DISPAXIS=dispaxis)

            for flat in list_flat_r:
                os.remove('pre_reduced/to'+ flat.split('.')[0]+'_norm.fits')



    if HOST:
        host_gals.make_host_metadata(configDict)

    return 0


def pre_reduction_orig():

    description = "> Performs pre-reduction steps"
    usage = "%prog    \t [option] \n Recommended syntax: %prog -i -c"
  
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option('-n','--nflats',type=int,default=50,
                        help='max number of flats to use (to prevent IRAF crash)')
    option, args = parser.parse_args()
    MAX_N_FLATS = option.nflats
    
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)

    iraf.ccdred.verbose = 'no'
    iraf.specred.verbose = 'no'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''

    iraf.longslit.mode = 'h'
    iraf.specred.mode = 'h'
    iraf.noao.mode = 'h'
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"

    mkarc = raw_input("Make arc? ([y]/n): ")
    mkflat = raw_input("Make flat? ([y]/n): ")

    if len(args) > 1:
        files=[]
        sys.argv.append('--help')
        option, args = parser.parse_args()
        sys.exit()
    elif len(args) == 1:
        files = util.readlist(args[0])
        sys.exit()
    else:
        listfile = glob.glob('*.fits')
        files_science = []
        files_arc = []
        files_dflat = []
        #print 'checking your files ...'
        for img in listfile:
            _type = ''
            hdr0 = util.readhdr(img)
            _type=util.readkey3(hdr0, 'object')
            if 'flat' in _type.lower():
                files_dflat.append(img)
            elif 'arc' not in _type.lower() and 'arc' not in img.lower():
                files_science.append(img)
        if mkarc != 'n':
            mkarc_b = raw_input("List blue arc files to combine (.fits will be added): ").split()
            mkarc_r = raw_input("List red arc files to combine (.fits will be added): ").split()
            for arc in mkarc_b:
                files_arc.append(arc + '.fits')
            for arc in mkarc_r:
                files_arc.append(arc + '.fits')

    if mkarc != 'n':
        list_arc_b = []
        list_arc_r = []
        for arcs in files_arc:
            if instruments.blue_or_red(arcs)[0] == 'blue':
                list_arc_b.append(arcs)
            elif instruments.blue_or_red(arcs)[0] == 'red':
                list_arc_r.append(arcs)
            else:
                sys.exit()

    if mkflat != 'n':
        list_flat_b = []
        list_flat_r = []
        for dflats in files_dflat:
            if instruments.blue_or_red(dflats)[0] == 'blue':
                if len(list_flat_b) < MAX_N_FLATS:
                    list_flat_b.append(dflats)
            elif instruments.blue_or_red(dflats)[0] == 'red':
                if len(list_flat_r) < MAX_N_FLATS:
                    list_flat_r.append(dflats)
            else:
                sys.exit()
                
                
    # make pre_reduced if it doesn't exist
    if not os.path.isdir('pre_reduced/'):
        os.mkdir('pre_reduced/')
        
    # log the existing processed files (need to verify this works if pre_reduced is empty...)
    pfiles = []
    new_files = []
    for root, dirnames, filenames in os.walk('pre_reduced'):
        for file in filenames:
            if file.startswith('to'):
                pfiles.append(file)
    #print(pfiles)

    # loop over each image in pre_reduced
    for img in listfile:
        hdr = util.readhdr(img)
        targ=util.readkey3(hdr, 'object')
        
        # if file is not not a processed file, run the overscan+trim code
        if 'to'+ img not in pfiles:
            
            # if the file is a science file, grab the name for later
            if 'arc' not in targ.lower() and 'flat' not in targ.lower():
                new_files.append(img)
                print ('Adding data for: ' + targ)
                
            inst = instruments.blue_or_red(img)[1]

            iraf.specred.dispaxi = inst.get('dispaxis')
            iraf.longslit.dispaxi = inst.get('dispaxis')

            _biassec0 = inst.get('biassec')
            _trimsec0 = inst.get('trimsec')
            
            ######################################################################
            #
            # JB: this chunk of code needs attention
            # It seems hacky for anything but Kast...
            #
            # overscan
            if not img.startswith('o') and inst.get('observatory')=='lick':
                if os.path.isfile('pre_reduced/o'+img):
                    os.remove('pre_reduced/o'+img)
                util.kastbias(img,'pre_reduced/o'+img)
            elif not img.startswith('o') and inst.get('observatory')!='lick':
                if os.path.isfile('pre_reduced/o'+img):
                    os.remove('pre_reduced/o'+img)
                os.system('cp ' +  img + ' ' + 'pre_reduced/' + img)

                
            # trim
            if not img.startswith('t')and inst.get('observatory')=='lick':
                if os.path.isfile('pre_reduced/to'+img):
                    os.remove('pre_reduced/to'+img)
                iraf.ccdproc('pre_reduced/o'+img, output='pre_reduced/to'+img, 
                             overscan='no', trim='yes', zerocor="no", flatcor="no", 
                             readaxi='line',trimsec=str(_trimsec0), Stdout=1)

            elif not img.startswith('t')and inst.get('observatory')!='lick':
                if os.path.isfile('pre_reduced/to'+img):
                    os.remove('pre_reduced/to'+img)
                # iraf.ccdproc('pre_reduced/'+img, output='pre_reduced/to'+img, 
                #              overscan='yes', trim='yes', zerocor="no", flatcor="no", 
                #              readaxi='line',trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)
                iraf.ccdproc('pre_reduced/'+img, output='pre_reduced/to'+img, 
                             overscan='no', trim='yes', zerocor="no", flatcor="no", 
                             readaxi='line',trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)

        else:
            pass 
            # 'to'+img exists, so do nothing, 




    # combine the arcs
    if mkarc != 'n':
        
        # blue arcs
        if len(list_arc_b) > 0:
            if len(list_arc_b) == 1:
                arc_blue = list_arc_b[0]
                os.system('cp ' + 'pre_reduced/to'+ arc_blue + ' ' + 'pre_reduced/ARC_blue.fits')
            else:
                arc_str = ''
                for arc in list_arc_b:
                    arc_str = arc_str + 'pre_reduced/to'+ arc + ','
                if os.path.isfile('pre_reduced/ARC_blue.fits'):
                    os.remove('pre_reduced/ARC_blue.fits')
                iraf.imcombine(arc_str, output='pre_reduced/ARC_blue.fits')

        # red arcs
        if len(list_arc_r) > 0:
            if len(list_arc_r) == 1:
                arc_red = list_arc_r[0]
                os.system('cp ' + 'pre_reduced/to'+ arc_red + ' ' + 'pre_reduced/ARC_red.fits')
            else:
                arc_str = ''
                for arc in list_arc_r:
                    arc_str = arc_str + 'pre_reduced/to'+ arc + ','
                if os.path.isfile('pre_reduced/ARC_red.fits'):
                    os.remove('pre_reduced/ARC_red.fits')
                iraf.imcombine(arc_str, output='pre_reduced/ARC_red.fits')




    # combine the flats
    if mkflat != 'n':
        inter = 'yes'
        
        # blue flats
        if len(list_flat_b) > 0:
            br, inst = instruments.blue_or_red(list_flat_b[0])
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_blue = 'pre_reduced/toFlat_blue.fits'

            flat_list = []
            for flat in list_flat_b:
                flat_list.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_blue):
                os.remove(Flat_blue)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list,OUTFILE=Flat_blue,MEDIAN_COMBINE=True)
            
            # run iraf response
            iraf.specred.response(Flat_blue, 
                                   normaliz=Flat_blue, 
                                   response='pre_reduced/RESP_blue', 
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=5,high_rej=5, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_blue.fits'],DISPAXIS=dispaxis)

        # red flats
        if len(list_flat_r) > 0:
            br, inst = instruments.blue_or_red(list_flat_r[0])
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_red = 'pre_reduced/toFlat_red.fits'


            flat_list = []
            for flat in list_flat_r:
                flat_list.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_red):
                os.remove(Flat_red)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list,OUTFILE=Flat_red,MEDIAN_COMBINE=True)

            #What is the output here? Check for overwrite
            iraf.specred.response(Flat_red, 
                                  normaliz=Flat_red, 
                                  response='pre_reduced/RESP_red', 
                                  interac=inter, thresho='INDEF',
                                  sample='*', naverage=2, function='legendre', 
                                  low_rej=5,high_rej=5, order=80, niterat=20, 
                                  grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_red.fits'],DISPAXIS=dispaxis)
    









    # science files should have 't' in front now
    # this just gets the base name, to prefix assumed below
    if new_files is not None:
        files_science = new_files

    # get all the science objects for the night
    science_targets = []
    for obj in files_science:
        hdr = util.readhdr(obj)
        _type=util.readkey3(hdr, 'object')
        science_targets.append(_type)

    # make a dir for each sci object
    science_targets = set(science_targets)
    for targ in science_targets:
        if not os.path.isdir('pre_reduced/' + targ + '/'):
            os.mkdir('pre_reduced/'+ targ + '/')

    # copy the files into the obj dir
    for obj in files_science:
        hdr = util.readhdr(obj)
        targ=util.readkey3(hdr, 'object')
        if not obj.startswith('to'):
            os.system('cp ' + 'pre_reduced/to'+ obj + ' ' + 'pre_reduced/' + targ + '/')
        else:
            os.system('cp ' +  'pre_reduced/'+ obj + ' ' + 'pre_reduced/' + targ + '/')

    rawfiles = glob.glob('*.fits')
    ofiles = glob.glob('pre_reduced/o'+ '*.fits')
    tfiles = glob.glob('pre_reduced/to'+ '*.fits')
    
    # delete raw files from the pre_reduced dir
    # there shouldn't be any there though?
    # maybe if the overscan isn't implemented for that detector
    for img in rawfiles:
        util.delete('pre_reduced/' + img)
        
    # delete the ofiles from pre_reduced dir
    for img in ofiles:
        util.delete(img)
    
    
    
def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Pre-reduction for the UCSC Spectroscopic Pipeline'
    parser = argparse.ArgumentParser(description=descStr)
    configCG = parser.add_mutually_exclusive_group()
    basicProcCG = parser.add_mutually_exclusive_group()

    # optional
    parser.add_argument('-q','--quicklook',
                        help='Move archival calibrations to enable quicker reductions', action='store_true')
    parser.add_argument('-v','--verbose',
                        help='Print diagnostic info',action='store_true')
    parser.add_argument('-k','--clobber',
                        help='Clobber existing files/dirs',action='store_true')
    parser.add_argument('-f','--fast',
                        help='Use config as is, don\'t prompt user for anything',action='store_true')


    parser.add_argument('--make-arcs',
                        help='Combine arcs into master arc images', action='store_true')
    parser.add_argument('--make-flats',
                        help='Combine flats into master flat images', action='store_true')
    parser.add_argument('--host',
                        help='Obtain relevant host galaxy metadata', action='store_true')
    parser.add_argument('--red-amp-bad',
                        help='Red side amplifier is bad so trim', action='store_true')

    basicProcCG.add_argument('--fake-basic-2d',
                        help='Fake the basic 2D reductions', action='store_true')
    basicProcCG.add_argument('--full-clean',
                        help='Completely wipe pre_reduction (USE WITH CAUTION)', action='store_true')

    configCG.add_argument('-o','--original',
                    help='Run the original pre_reduction',action='store_true')
    configCG.add_argument('-c','--config-file',
                    help='Config file to use for pre-reduction')




    # parse
    cmdArgs = parser.parse_args()


    args = ()
    kwargs = {}

    kwargs['VERBOSE'] = cmdArgs.verbose
    kwargs['CLOBBER'] = cmdArgs.clobber
    kwargs['FULL_CLEAN'] = cmdArgs.full_clean
    kwargs['FAKE_BASIC_2D'] = cmdArgs.fake_basic_2d
    kwargs['FAST'] = cmdArgs.fast
    kwargs['ORIGINAL'] = cmdArgs.original
    kwargs['CONFIG_FILE'] = cmdArgs.config_file
    kwargs['MAKE_ARCS'] = cmdArgs.make_arcs
    kwargs['MAKE_FLATS'] = cmdArgs.make_flats
    kwargs['QUICK'] = cmdArgs.quicklook
    kwargs['HOST'] = cmdArgs.host
    kwargs['RED_AMP_BAD'] = cmdArgs.red_amp_bad

    return (args,kwargs)

def main(*args,**kwargs):
    '''
    Main driver for running pre-reduction versions

    Parameters
    ----------
    data : list
        List to parse for unique values

    IS_TUPLE : bool, optional
        If True, data is assumed to be a list of tuples

    KEY_INDEX : int, optional
        If IS_TUPLE, KEY_INDEX is used as the tuple index which
        holds the values to construct the unique set from

    Returns
    -------
    set : A set of unique values in data

    '''

    if kwargs.get('ORIGINAL'):
        res = pre_reduction_orig()

    else:
        res = pre_reduction_dev(**kwargs)

    return 0

if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)
