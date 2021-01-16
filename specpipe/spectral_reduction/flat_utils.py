from __future__ import print_function
import os,sys,pdb,argparse,shutil,glob,subprocess,shlex
from time import sleep
import numpy as np

from astropy.io import fits
from astropy.io import ascii

import scipy
from scipy import signal
from scipy import interpolate
from scipy import optimize
from scipy import signal, ndimage

from optparse import OptionParser
import instruments

import matplotlib

matplotlib.use('TkAgg')
matplotlib.rcParams[u'keymap.yscale'].remove(u'l')
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as colors


class fitFlatClass(object):
    '''
    A class to drive flatfield masking operations

    Attributes
    ----------
    rawData : ndarray
        A 2D array representing a stacked color free flatfield image
    flatModelData : ndarray
        A 2D array representing a flatfield model image
    flatCorrData: ndarray
        A 2D array representing a masked/corrected flatfield image
    masterProfile: ndarry
        A 1D array representing average spatial illumination pattern
    fig : matplotlib.figure
        A figure instance displaying the current state of the class
    splineSmooth: float
        Parameter used in smoothing when fitting the telluric features
    regionDict: dictionary
        A dictionary storing masking regions
    dummyRegion: dictionary
        A dictionary representing a prototype masking region

    '''
    def __init__(self,image,fig,inst):
        ''' 
        Initialize the object 

        Parameters
        ----------
        image : ndarray
            A 2D array representing a stacked color free flatfield image
        fig : bool, optional
            A figure instance displaying the current state of the class
        '''
        
        # data
        self.rawData = 1.*image # observed
        self.flatModelData = 1.*image # model
        self.flatCorrData = 1.*image # new flat (obs corrected by model)
        self.masterProfile = 0.*image[:,0] + 1. # single column illumination profile
        
        # figure
        self.fig = fig
        self.inst = inst
        
        # spline smoothing (used in fitting)
        # use larger values to accept worse fits,
        # smaller values to mimic true spline
        self.splineSmooth = 1.
        
        # dict of fit regions
        self.regionDict = {}
        
        
        # dummy region dict
        self.dummyRegion = {'colLo': -1,
                            'colUp': -1,
                            'done': False,
                            'store': False}
        
    def skyRegion_onClick(self,event):
        ''' General function for handling clicks on the canvas '''
        return
        
    def skyRegion_onKeyPress(self,event):
        ''' General function for handling key presses on the canvas '''

        # get the axes for plotting and cursor detection
        ax_list = self.fig.axes
        
        ax1 = ax_list[0]
        ax2 = ax_list[1]
        ax3 = ax_list[2]
        ax4 = ax_list[3]
        ax5 = ax_list[4]
                
        # make sure the cursor was in the correct axis
        if event.inaxes is ax4:
        
            # get the cursor position
            rowPress = event.ydata
            colPress = event.xdata
            
            # match row to nearst existing pixel
            if rowPress < 0:
                row = 0
            elif rowPress >= 0 and rowPress <= self.rawData.shape[0]:
                row = rowPress
            else:
                row = self.rawData.shape[0]
                
            # match col to nearest existing pixel
            if colPress < 0:
                col = 0
            elif colPress >= 0 and colPress <= self.rawData.shape[1]:
                col = colPress
            else:
                col = self.rawData.shape[1]
            
                
            if event.key == 'l':
                outStr = 'Adding lower column'
                print(outStr)
                self.dummyRegion['colLo'] = col
            
            elif event.key == 'u': 
                outStr = 'Adding upper column'
                print(outStr)
                self.dummyRegion['colUp'] = col
            
            elif event.key == 'd':
                outStr = 'Adding region to dict'
                print(outStr)
                self.dummyRegion['done'] = True
                self.dummyRegion['store'] = True
            
            elif event.key == 'q':
                outStr = 'Discarding region'
                print(outStr)
                self.dummyRegion['done'] = True
                self.dummyRegion['store'] = False
                        
            else:
                print('Unrecognized key')
        
        
        return
        
    def add_fit_region(self,name,colLo=-1,colUp=-1):
        '''
        Add a fitting region to the current sky model object

        Parameters
        ----------
        name : string
            Name of the fitting region
        colLo : int
            Lower column of the fitting region
        colUp : int
            Upper column of the fitting region

        Returns
        -------
        int : 0, and modifies self.regionDict
        '''
        
        # reset the dummy region dict
        self.dummyRegion['colLo'] = -1
        self.dummyRegion['colUp'] = -1
        self.dummyRegion['done'] = False
        self.dummyRegion['store'] = False
        
        # init region
        newRegion = flatFitRegion(name)
                
        # get the axes for plotting and cursor detection
        ax_list = self.fig.axes
        ax1 = ax_list[0]
        ax2 = ax_list[1]
        ax3 = ax_list[2]
        ax4 = ax_list[3]
        ax5 = ax_list[4]
        
        # user gave regions
        if colLo > 0 and colUp > 0:
            
            # populate the dummy variable
            self.dummyRegion['colLo'] = colLo
            self.dummyRegion['colUp'] = colUp
            self.dummyRegion['done'] = True
            self.dummyRegion['store'] = True
        
        # enter interactive    
        else:
        
            # print instructions
            outStr = 'Hover over the row you\'d like to add\n'
            outStr += 'Press L to mark a lower column\n'
            outStr += 'Press U to mark an upper column\n'
            outStr += 'Press D to add the region\n'
            outStr += 'Press Q to discard the region\n'
            print(outStr)            
            
            # connect the key press event here
            cid = self.fig.canvas.mpl_connect('button_press_event', self.skyRegion_onClick)
            cid2 = self.fig.canvas.mpl_connect('key_press_event', self.skyRegion_onKeyPress)
            
            lnL = ''
            lnU = ''
            
            while not self.dummyRegion['done']:
                
                if lnL:
                    lnL.remove()
                if lnU:
                    lnU.remove()
                                                
                lnL = ax4.axvline(x=self.dummyRegion['colLo'],c='#0d8202',lw=2.,ls='--')
                lnU = ax4.axvline(x=self.dummyRegion['colUp'],c='#6ef961',lw=2.,ls='--')
                
                # update the plot
                plt.pause(0.05)

            # disconnect
            self.fig.canvas.mpl_disconnect(cid)                
            self.fig.canvas.mpl_disconnect(cid2)                

                
        # flip order if necessary
        if self.dummyRegion['colLo'] > self.dummyRegion['colUp']:
            errStr = 'Had to flip the high/low columns...'
            print(errStr)
            colTemp = self.dummyRegion['colLo']
            self.dummyRegion['colLo'] = self.dummyRegion['colUp']
            self.dummyRegion['colUp'] = colTemp
            
        # assign values
        colLo = int(self.dummyRegion['colLo'])
        colUp = int(self.dummyRegion['colUp'])
        newRegion.colLo = colLo
        newRegion.colUp = colUp
        newRegion.flux = self.rawData[:,colLo:colUp]
        
        # store in the dict
        if self.dummyRegion['store']:
            self.regionDict[name] = newRegion
            
            
        # update the plot
        self.refresh_plot()
        
        return 0

    def remove_fit_region(self,name=''):
        ''' Pops self.regionDict[name] '''
               
        if len(name) > 0:
            trash = self.regionDict.pop(name, None)
                
        elif len(name) == 0:
            
            while True:
                outStr = 'Here\'s the info on the fitting regions:\n'
                for key in self.regionDict.keys():
                    region = self.regionDict[key]
                    outStr += '{} {} {}\n'.format(region.name,region.colLo,region.colUp)
                outStr += 'Enter the name of the region to delete, or Q to quit: '
                usrResp = raw_input(outStr).strip().upper()
            
                # delete the entry
                if usrResp in self.regionDict.keys():
                    trash = self.regionDict.pop(usrResp, None)
                    self.refresh_plot()
                    break
                    
                elif usrResp == 'Q':
                    outStr = 'Quitting region deletion...'
                    print(outStr)
                    break
                    
                else:
                    outStr = 'I do not understand, try again...'
                    print(outStr)
                    
        else:
            outStr = 'Something went wrong...'
            print(outStr)
            
        
        return 0
        
        
    def update_master_profile(self):
        ''' Updates the master profile '''
        
        maskedData = np.array([])
        colArr = np.array([])
        
        # unpack the exclusion regions
        for key in self.regionDict.keys():
            region = self.regionDict[key]
            newCols = np.arange(region.colLo,region.colUp,1)
            colArr = np.append(colArr,newCols)
            
            
        maskedData = self.rawData[:,np.arange(self.rawData.shape[1]) != colArr]
        maskedData = maskedData[:,0,:] # weird
         
        # snippet for sanity checks           
        #fig=plt.figure(figsize=(6,6))    
        #ax1 = plt.subplot2grid((4,4),(0,0),rowspan=2,colspan=4)
        #ax2 = plt.subplot2grid((4,4),(2,0),rowspan=2,colspan=4)
        #ax1.imshow(self.rawData,origin='lower',vmin=0.95,vmax=1.05)
        #ax2.imshow(maskedData,origin='lower',vmin=0.95,vmax=1.05)
        #plt.show(block=True)
        #pdb.set_trace()
        
        updatedProfile = np.median(maskedData,axis=1)
        self.masterProfile = 1.*updatedProfile
        
        return 0
        
    def refresh_plot(self):
        ''' Updates the plot to reflect current state '''
        
        # get the axes for plotting and cursor detection
        ax_list = self.fig.axes
        ax1 = ax_list[0]
        ax2 = ax_list[1]
        ax3 = ax_list[2]
        ax4 = ax_list[3]
        ax5 = ax_list[4]
        
        ax1.cla()
        ax2.cla()
        ax3.cla()
        ax4.cla()
        ax5.cla()
        
        modColArr = np.array([])
        # get cols where sky model is defined
        for key in self.regionDict.keys():
            region = self.regionDict[key]
            newCols = np.arange(region.colLo,region.colUp,1)
            modColArr = np.append(modColArr,newCols)
        modColArr = modColArr.astype(int)

        # ax1 data
        bsd_col_lo = 0
        bsd_col_up = bsd_col_lo + 500
        
        # ax2 data
        # msd_col_lo = self.rawData.shape[1] // 2
        msd_col_lo = 1800
        msd_col_up = msd_col_lo + 300
        
        # ax3 data
        # rsd_col_lo = self.rawData.shape[1] - 500
        rsd_col_lo = 3200
        rsd_col_up = rsd_col_lo + 400
        
        # ax1 data
        blueSkyData = np.median(self.rawData[:,bsd_col_lo:bsd_col_up],axis=0)
        if len(modColArr) > 0:
            blueSkyModel = np.median(self.flatModelData[:,bsd_col_lo:bsd_col_up],axis=0)
        else:
            blueSkyModel = 0.*blueSkyData
        blueSkyX = np.arange(bsd_col_lo,bsd_col_up,1)

    
    
        # ax2 data
        midSkyData = np.median(self.rawData[:,msd_col_lo:msd_col_up],axis=0)
        midSkyX = np.arange(msd_col_lo,msd_col_up,1)
    
        # ax3 data
        redSkyData = np.median(self.rawData[:,rsd_col_lo:rsd_col_up],axis=0)
        if len(modColArr) > 0:
            redSkyModel = np.median(self.flatModelData[:,rsd_col_lo:rsd_col_up],axis=0)
        else:
            redSkyModel = 0.*redSkyData
        redSkyX = np.arange(rsd_col_lo,rsd_col_up,1)

        # plot blue sky
        try:
            pass
            ax1.plot(blueSkyX,blueSkyData,c='k',ls='-',lw=3.)
            ax1.plot(blueSkyX,blueSkyModel,c='r',ls='--',lw=3.)
        except Exception as e:
            print(e)
            pdb.set_trace()

        #plot profile
        # ax2.plot(midSkyX,midSkyData,c='k',ls='-',lw=3.)
        
        # #plot red sky
        # ax3.plot(redSkyX,redSkyData,c='k',ls='-',lw=3.)
        # ax3.plot(redSkyX,redSkyModel,c='r',ls='--',lw=3.)

        # grab the values of the 10 and 90 percentile pixels
        rawData_ravel = np.ravel(self.rawData)
        sortedIndexes = np.argsort(rawData_ravel)
        # vmin = 1.*rawData_ravel[sortedIndexes[int(0.2*len(sortedIndexes))]]
        # vmax = 2.*rawData_ravel[sortedIndexes[int(0.8*len(sortedIndexes))]]

        vmin = 0.9
        vmax = 1.1

        # image
        ax4.imshow(self.rawData,aspect=1.,origin='lower',
                    norm=colors.Normalize(vmin=vmin,vmax=vmax))

        # residuals
        ax5.imshow(self.flatCorrData,aspect=1.,origin='lower',
                    norm=colors.Normalize(vmin=vmin,vmax=vmax))
        
        # over plot the sky regions
        for key in self.regionDict.keys():
            region = self.regionDict[key]
            # lower
            ax4.axvline(x=region.colLo,c='#990000',lw=2,ls='--')
            # upper
            ax4.axvline(x=region.colUp,c='#ff4f4f',lw=2,ls='--')
            
        # sparse axis labels
        ax1.set_xlabel('column')
        ax1.set_ylabel('counts')
        ax2.set_xlabel('column')
        ax3.set_xlabel('column')

        ax1.set_yticklabels([])
        ax2.set_yticklabels([])
        ax3.set_yticklabels([])
        
        # ranges on image plots
        ax4.set_xlim([0,self.rawData.shape[1]])
        ax4.set_ylim([0,self.rawData.shape[0]])
        
        ax5.set_xlim([0,self.rawData.shape[1]])
        ax5.set_ylim([0,self.rawData.shape[0]])
        
        return 0
        

    def fit_sky_background(self,FIT_METHOD='LM'):
        ''' Fit the sky flux and update the sky model '''
        
        # init
        skyImage = np.array([]) # actual sky data
        skyModel = np.array([]) # model of sky data
        skyModelFull = 0.*self.rawData # model across entire chip
        colArr = np.array([])
        
        # update the master profile
        outStr = 'Updating the master profile (takes a bit...)'
        print(outStr)
        self.update_master_profile()
        
        # stack the regions, but preserve the row pixel numbers in rowArr
        # this will be trickier if I implement curved sky regions
        for key in self.regionDict.keys():
            region = self.regionDict[key]
            newCols = np.arange(region.colLo,region.colUp,1)
            colArr = np.append(colArr,newCols)
            
            if len(skyImage) == 0:
                skyImage = 1.*region.flux
            else:
                skyImage = np.hstack((skyImage,1.*region.flux))
                
                
        colArr = colArr.astype(int)
        # for each col, fit flux as a function of pixel number   
        for i in xrange(skyImage.shape[1]): 
            if i % 50 == 0:  
                print('Working on col {} / {}'.format(i,skyImage.shape[1]))
            
            # get the fit data for this column            
            fitX = np.arange(0,skyImage.shape[0]) # abscissa is row pixel number

            # subtract the masterProfile (illumination), and then fit the 
            # result with a smoothed spline. The fit should approximate the sky lines,
            # and the residuals should approximate the pixel-to-pixel sensitivity variations.
            fitY = skyImage[fitX,i] - self.masterProfile
            
            # fit with a smoothed spline
            # note that the results are sensitive to the choice of 's' (the smoothing)
            splineFit = interpolate.UnivariateSpline(fitX, fitY, s=self.splineSmooth)
            skyTheo = splineFit(fitX)
            
            # subtract fit from observed
            residual = self.rawData[:,colArr[i]] - skyTheo
            
            # add residuals (pixel to pixel variations) back in
            skyModelFull[:,colArr[i]] = residual
            
            
            # plot?
            # if colArr[i] == 3348:
            #     PLOT=True
            # else:
            #     PLOT=False
            PLOT = False
            if PLOT:
                fig=plt.figure(figsize=(6,6))    
                axMain = plt.subplot2grid((4,4),(0,0),rowspan=4,colspan=6)
                axMain.plot(fitX,fitY,c='k',ls='',marker='.')
                axMain.plot(fitX,skyTheo,c='r',ls='--',lw=0.5)
                plt.show()
                pdb.set_trace()
                                       
            
            
        # insert into the full object image (no trasposing needed??)
        for key in self.regionDict.keys():
            region = self.regionDict[key]
            newCols = np.arange(region.colLo,region.colUp,1)
            newCols = newCols.astype(int)
            
            self.flatModelData[:,newCols] = skyModelFull[:,newCols] # obviously check this

        self.refresh_plot()
        return 0
        
        return 0
        
        
        
    def subsitute_model_flat(self):
        ''' Substitutes the model in the masked regions '''

        # insert into the full object image (no trasposing needed??)
        for key in self.regionDict.keys():
            region = self.regionDict[key]
            newCols = np.arange(region.colLo,region.colUp,1)
            newCols = newCols.astype(int)
            self.flatCorrData[:,newCols] = 1.*self.flatModelData[:,newCols] 


        if 'blue' in self.inst.get('name'):
            std_tol = 0.03 #subject to change
            count = 0
            found_blue = False
            found_red = False
            for i in range(len(self.flatModelData[0,:])):
                newCols = np.arange(i,i+20,1)
                newCols = newCols.astype(int)
                if i < len(self.flatModelData[0,:]) - 21:
                    std = np.std(np.ravel(self.flatModelData[:,newCols]))
                    print (i, std)
                    if std < std_tol:
                        count+=1
                    if std < std_tol and count > 10 and not found_blue:
                        blue_ind = i - 10
                        found_blue = True
                        count=0
                        if 'blue' in self.inst.get('name'):
                            break

                    #do not substitute in red side of blue detector
                    if 'red' in self.inst.get('name'):
                        if std > std_tol and found_blue:
                            count+=1
                        if std > std_tol and count > 10 and not found_red and found_blue:
                            red_ind = i - 10
                            found_red = True
                            break

            # good_range = [1200,1300]#kast
            good_range = [2500,2700]#lris
            medCols = np.arange(good_range[0],good_range[1],1)
            medCols = medCols.astype(int)

            print ('Substituting median profile at edges up to column ', blue_ind)
            subdata = np.median(self.rawData[:,medCols], axis=1)
            num_blue_cols = blue_ind
            for col in range(num_blue_cols):
                self.flatCorrData[:,col] = subdata


            if 'red' in self.inst.get('name') and found_red:
                num_red_cols = red_ind
                for col in range(num_blue_cols):
                    neg_col = int(-1.*col)
                    self.flatCorrData[:,neg_col] = subdata


        self.refresh_plot()
        return 0
        
        
    #
    # Substitutes the data in a user supplied region with self.masterProfile
    # This is useful if there is a region of the data that is simply trash
    #
    def hard_mask(self,name=''):
        ''' Substitutes the master profile in the specified mask region '''
        
        # sub the master profile
        # inefficient, don't care
        if len(name) > 0:
            colLo = self.regionDict[usrResp].colLo
            colUp = self.regionDict[usrResp].colUp
            for i in xrange(colUp-colLo):
                self.flatCorrData[:,colLo+i] = self.masterProfile
            self.refresh_plot()                
        elif len(name) == 0:
            
            while True:
                outStr = 'Here\'s the info on the fitting regions:\n'
                for key in self.regionDict.keys():
                    region = self.regionDict[key]
                    outStr += '{} {} {}\n'.format(region.name,region.colLo,region.colUp)
                outStr += 'Enter the name of the region to mask, or Q to quit: '
                usrResp = raw_input(outStr).strip().upper()
            
                # sub the master profile
                # inefficient, don't care
                if usrResp in self.regionDict.keys():
                    colLo = self.regionDict[usrResp].colLo
                    colUp = self.regionDict[usrResp].colUp
                    for i in xrange(colUp-colLo):
                        self.flatCorrData[:,colLo+i] = self.masterProfile
                    self.refresh_plot()
                    break
                    
                elif usrResp == 'Q':
                    outStr = 'Quitting region deletion...'
                    print(outStr)
                    break
                    
                else:
                    outStr = 'I do not understand, try again...'
                    print(outStr)
                    
        else:
            outStr = 'Something went wrong...'
            print(outStr)
            
        
        return 0
        
        
        
    def refine(self):
        ''' 
        Substitutes the model flat for raw data and wipes 
        exclusion regions. Use with caution; this is experimental,
        but could be useful for iteratively improving the flat.
        '''
        self.rawData = 1.*self.flatCorrData
        self.regionDict = {}
        self.refresh_plot()
        return 0
        
    def save_flat(self,outFile,header=None):
        ''' Save the corrected flat. Assumes overwriting has been verified '''
        
        # clear space
        if os.path.isfile(outFile):
            os.remove(outFile)

        
        # write correct flat data
        if header is not None:
            hdu = fits.PrimaryHDU(self.flatCorrData,header)
        else:
            hdu = fits.PrimaryHDU(self.flatCorrData)
        hdu.writeto(outFile,output_verify='ignore')  
        
        return 0
        
        

class flatFitRegion(object):
    ''' Class representing an exclusion region '''
    def __init__(self,name):
        self.name = name
        self.colLo = 0
        self.colUp = 1
        self.flux = np.array([])



def combine_flats(flat_list,MEDIAN_COMBINE=False,**kwargs):
    ''' Stacks images in flat_list, returns a master flat and header '''

    # unpack
    outFile = kwargs.get('OUTFILE')
    clobber = kwargs.get('CLOBBER')

    # read data
    flat_comb_image = np.array([])
    median_image_stack = np.array([])
    nImages = 0
    expTime = 0
    nFlatLimit = 15

    # calculate the stack batch size
    if nFlatLimit < len(flat_list):
        batchSize = len(flat_list) // 2 + 1
        while batchSize > len(flat_list):
            batchSize = batchSize // 2 + 1
    else:
        batchSize = len(flat_list)



    # loop over the flat files
    for file in flat_list:
        hdu = fits.open(file)
        br, inst = instruments.blue_or_red(file)
        data = hdu[0].data
        header = hdu[0].header
        nImages += 1
        expTime += header.get('EXPTIME')

        # scale to expTime
        data /= expTime

        if len(flat_comb_image) == 0:
            flat_comb_image = np.copy(data)
        else:

            # if median combining, stack the data
            if MEDIAN_COMBINE:

                # stack images in the z direction
                flat_comb_image = np.dstack((flat_comb_image,data))

                # if flat_comb_image is getting too big, or if we're at the end of flat_list,
                # then squash it along the z axis with a median
                if ((flat_comb_image.shape[2] == batchSize) or 
                    (flat_comb_image.shape[2] == len(flat_list)-1)):

                    # stack the intermediate squashed frames
                    if len(median_image_stack) == 0:
                        median_image_stack = np.median(flat_comb_image,axis=2)
                    else:
                        median_image_stack = np.dstack((median_image_stack,
                                                        np.median(flat_comb_image,axis=2)))

                    # reset the stack
                    flat_comb_image = np.array([])


            # otherwise just sum them
            else:
                flat_comb_image += np.copy(data)

    # if median combining, squash the stack of median
    if MEDIAN_COMBINE:
        # if there are multiple median images in the stack, median them
        if len(median_image_stack.shape) > 2:
            flat_comb_image = np.median(median_image_stack,axis=2)
        # otherwise just return the single frame
        else:
            flat_comb_image = np.copy(median_image_stack)
        header.add_history('combine_flats: median combined {} files'.format(nImages))
    else:
        header.add_history('combine_flats: summed {} files'.format(nImages))

    # counts / sec * expTime
    flat_comb_image *= expTime

    if outFile:
        
        # clear space
        if os.path.isfile(outFile) and clobber:
            os.remove(outFile)

        # write correct flat data
        hdu = fits.PrimaryHDU(flat_comb_image,header)
        hdu.writeto(outFile,output_verify='ignore')  

    return (flat_comb_image,header,inst)





def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Utility for masking bad regions on a flat field '
    parser = argparse.ArgumentParser(description=descStr)

    # required args
    parser.add_argument('flat_list',type=str,nargs='*',default=[],
                       help='The list of flat field files to process')

    # optional
    parser.add_argument('-v','--verbose',
                        help='Print diagnostic info',action='store_true')
    parser.add_argument('-c','--clobber',action='store_true',
                        help='Clobber files already in pre_reduced/ but not subdirs')
    parser.add_argument('-f','--read_from_file',action='store_true',
                        help='flat_list is a file containing a list of flat files to read')
    parser.add_argument('-r','--remove_color',
                        help='Remove color term from flat',action='store_true')

    parser.add_argument('-o','--outfile',
                        help='Output file for the processed flat field')
    parser.add_argument('-d','--dispaxis',type=int,default=1,
                        help='Dispersion axis (default=1)')


    # parse
    cmdArgs = parser.parse_args()

    # logic mapping to my args/kwargs
    FLAT_LIST = cmdArgs.flat_list
    VERBOSE = cmdArgs.verbose
    CLOBBER = cmdArgs.clobber
    READ_FROM_FILE = cmdArgs.read_from_file
    REMOVE_COLOR = cmdArgs.remove_color
    OUTFILE = cmdArgs.outfile
    DISPAXIS = cmdArgs.dispaxis

    # package up
    args = (FLAT_LIST,)
    kwargs = {}
    kwargs['VERBOSE'] = VERBOSE
    kwargs['CLOBBER'] = CLOBBER
    kwargs['READ_FROM_FILE'] = READ_FROM_FILE
    kwargs['REMOVE_COLOR'] = REMOVE_COLOR
    kwargs['OUTFILE'] = OUTFILE
    kwargs['DISPAXIS'] = DISPAXIS


    return (args,kwargs)


def inspect_flat(flat_list,*args,**kwargs):
    '''
    Model a flat field from slitflat data

    Parameters
    ----------
    flat_list : list
        List of flat field image filenames

    OUTFILE : string, optional
        If specified, this will automatically be the file location of 
        the output flat image
    DISPAXIS : int, optional (Default = 1)
        Specifies the dispersion axis of the image. DISPAXIS=1
        corresponds to (spectral, spatial), DISPAXIS=2 corresponds
        to (spatial, spectral)
    REMOVE_COLOR: bool, optional (Default = False)
        If true, remove the color term from the supplied images.

    Returns
    -------
    int : 0, and writes files to disk
    '''
    
    # unpack
    outFile = kwargs.get('OUTFILE',None)
    dispaxis = kwargs.get('DISPAXIS',1)
    remove_color = kwargs.get('REMOVE_COLOR',False)
    read_from_file = kwargs.get('READ_FROM_FILE',False)

    # if user specified a file containing a list, read into flat_list
    if read_from_file:
        flat_list_tmp = []
        with open(flat_list,'r') as fin:
            for line in fin:
                if len(line.split()) > 0 and line.split()[0] != '#':
                    flat_list_tmp.append(line.split()[0].strip())

        # assign to flat_list, then proceed as usual
        flat_list = flat_list_tmp

    # if user passed multiple files, combine them, otherwise, just read in the single file
    if len(flat_list) > 1:
        flat_comb_image, header, inst = combine_flats(flat_list,**kwargs)
    else:
        br, inst = instruments.blue_or_red(flat_list[0])
        hdu = fits.open(flat_list[0])
        flat_comb_image = hdu[0].data
        header = hdu[0].header

    # transpose if we're dealing with cols x rows
    if dispaxis == 2:
        flat_comb_image = flat_comb_image.T

    # this is a hack and inferior to running IRAF response
    if remove_color:
        # for each column, divide by the median
        for i in range(len(flat_comb_image[0,:])):
            flat_comb_image[:,i] /= np.median(flat_comb_image[:,i])
    
    # set up plotting window
    plt.ion()
    
    fig=plt.figure(figsize=(16,8))
    axMain = plt.subplot2grid((36,36), (0,0), rowspan=36, colspan=36)
    ax1 = plt.subplot2grid((36,36), (0,0), rowspan=11, colspan=12)
    ax2 = plt.subplot2grid((36,36), (0,12), rowspan=11, colspan=12)
    ax3 = plt.subplot2grid((36,36), (0,24), rowspan=11, colspan=12)
    ax4 = plt.subplot2grid((36,36), (12,0), rowspan=12, colspan=36)
    ax5 = plt.subplot2grid((36,36), (24,0), rowspan=12, colspan=36)
    
    # ax1 data
    bsd_col_lo = 0
    bsd_col_up = bsd_col_lo + 500
    blueSkyData = np.median(flat_comb_image[:,bsd_col_lo:bsd_col_up],axis=0)
    blueSkyX = np.arange(bsd_col_lo,bsd_col_up,1)
    
    
    # ax2 data
    # msd_col_lo = flat_comb_image.shape[1] // 2
    msd_col_lo = 1800
    msd_col_up = msd_col_lo + 300
    midSkyData = np.median(flat_comb_image[:,msd_col_lo:msd_col_up],axis=0)
    midSkyX = np.arange(msd_col_lo,msd_col_up,1)
    
    # ax3 data
    # rsd_col_lo = flat_comb_image.shape[1] - 500
    rsd_col_lo = 3200
    rsd_col_up = rsd_col_lo + 400
    redSkyData = np.median(flat_comb_image[:,rsd_col_lo:rsd_col_up],axis=0)
    redSkyX = np.arange(rsd_col_lo,rsd_col_up,1)

    # ax1.plot(blueSkyX,blueSkyData,c='k',ls='-',lw=3.)
    # ax2.plot(midSkyX, midSkyData,c='k',ls='-',lw=3.)
    # ax3.plot(redSkyX,redSkyData,c='k',ls='-',lw=3.)

    # grab the values of the 10 and 90 percentile pixels
    flat_comb_image_ravel = np.ravel(flat_comb_image)
    sortedIndexes = np.argsort(flat_comb_image_ravel)
    # vmin = 1.*flat_comb_image_ravel[sortedIndexes[int(0.2*len(sortedIndexes))]]
    # vmax = 2.*flat_comb_image_ravel[sortedIndexes[int(0.8*len(sortedIndexes))]]
    vmin = 0.9
    vmax = 1.1

    # image
    ax4.imshow(flat_comb_image,aspect=1.,origin='lower',
                norm=colors.Normalize(vmin=vmin, vmax=vmax))
    
    # residuals
    ax5.imshow(flat_comb_image,aspect=1.,origin='lower',
                norm=colors.Normalize(vmin=vmin, vmax=vmax))

    # ranges on image plots
    ax4.set_xlim([0,flat_comb_image.shape[1]])
    ax4.set_ylim([0,flat_comb_image.shape[0]])
    
    ax5.set_xlim([0,flat_comb_image.shape[1]])
    ax5.set_ylim([0,flat_comb_image.shape[0]])
    
    # sparse axis labels
    ax1.set_xlabel('column')
    ax1.set_ylabel('counts')
    ax2.set_xlabel('column')
    ax3.set_xlabel('column')
    
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])

    flatFitObj = fitFlatClass(flat_comb_image,fig,inst)

    while True:
        
        # this really should be a dict of key/value pairs
        # and then the prompt is dynamically generated
        validResps = ['A','R','F','S','U','H',  # standard options
                      'AHARD','RHARD','REFINE', # poweruser/hidden options
                      'W','D','Q','Q!']              # stardard ends 
        promptStr = 'Enter (a) to add an exclusion region.\n'
        promptStr += 'Enter (r) to remove a region.\n'
        promptStr += 'Enter (f) to fit the exclusion regions.\n'
        promptStr += 'Enter (s) to substitute model in exclusion regions.\n'
        promptStr += 'Enter (h) to substitute the median profile in a region\n'
        promptStr += 'Enter (u) to undo everything and restart.\n'
        promptStr += 'Enter (w) to write the improved flat to disk.\n'
        promptStr += 'Enter (d) to enter the debugger.\n'
        promptStr += 'Enter (q) to quit and write the current flat to disk.\n'
        promptStr += 'Enter (q!) to quit and do nothing.'
        promptStr += 'Answer: '
        usrResp = raw_input(promptStr).strip().upper()
        
        if usrResp in validResps:
            
            # add region by marking it
            if usrResp == 'A':
                promptStr = 'Enter the name of the sky region (e.g. c1): '
                name = raw_input(promptStr).strip().upper()
                flatFitObj.add_fit_region(name)
                
            # add hardcoded region
            if usrResp == 'AHARD':
                promptStr = 'Enter name colLo colUp (e.g. c1 113 171): '
                usrResp = raw_input(promptStr).upper().strip()
                try:
                    name = usrResp.split()[0]
                    colLo = int(usrResp.split()[1])
                    colUp = int(usrResp.split()[2])
                    flatFitObj.add_fit_region(name,colLo=colLo,colUp=colUp)
                except Exception as e:
                    print(e)
                    
                    
            # remove
            if usrResp == 'R' or usrResp == 'RHARD':
                flatFitObj.remove_fit_region()
                
            # fit sky
            if usrResp == 'F':
                flatFitObj.fit_sky_background()
              
            # substitute model in sketchy regions  
            if usrResp == 'S':
                flatFitObj.subsitute_model_flat()
                
            if usrResp == 'U':
                flatFitObj = fitFlatClass(flat_comb_image,fig,inst)
                flatFitObj.refresh_plot()
                
            if usrResp == 'H':
                flatFitObj.hard_mask()
                
            if usrResp == 'REFINE':
                # this substitutes the corrected
                # flat for the input data...
                # basically an experiment, not intended for use
                flatFitObj.refine()
            
            # write file
            if usrResp == 'W':

                # if cols x rows, transpose before writting
                if dispaxis == 2:
                    flatFitObj.flatCorrData = flatFitObj.flatCorrData.T

                if outFile is None:
                    promptStr = 'Enter name of save file (e.g. RESP_blue): '
                    outFile = raw_input(promptStr).strip()
                    promptStr = 'Write to file {} [y/n]: '.format(outFile)
                    usrResp = raw_input(promptStr).upper().strip()
                    if usrResp == 'Y':
                        flatFitObj.save_flat('pre_reduced/'+outFile,header=header)
                        break
                    else:
                        outFile = None
                        print('Ok, aborting save...')
                else:
                    if os.path.isfile(outFile):
                        os.remove(outFile)
                    flatFitObj.save_flat(outFile,header=header)
            # debug
            if usrResp == 'D':
                pdb.set_trace()
            # quit
            if usrResp == 'Q':

                #mask problematic pixels
                print (flatFitObj.flatCorrData[flatFitObj.flatCorrData <= 0.001])
                flatFitObj.flatCorrData[flatFitObj.flatCorrData <= 0.001] = 1.

                print('Saving and quitting.')
                # if cols x rows, transpose before writting
                if dispaxis == 2:
                    flatFitObj.flatCorrData = flatFitObj.flatCorrData.T
                if outFile is None:
                    promptStr = 'Enter name of save file '
                    promptStr += '(e.g. RESP_blue.fits): '
                    outFile = raw_input(promptStr).strip()
                    promptStr = 'Write to file {} [y/n]: '.format(outFile)
                    usrResp = raw_input(promptStr).upper().strip()
                    if usrResp == 'Y':
                        flatFitObj.save_flat('pre_reduced/'+outFile,header=header)
                        break
                    else:
                        outFile = None
                        print('Ok, aborting save...')
                else:
                    if os.path.isfile(outFile):
                        os.remove(outFile)
                    flatFitObj.save_flat(outFile,header=header)
                break


            # quit
            if usrResp == 'Q!':
                print('Quitting.')
                break
        else:
            errStr = 'I don\'t understand, try again...'
            print(errStr)
    
    return 0
    
if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    inspect_flat(*args,**kwargs)