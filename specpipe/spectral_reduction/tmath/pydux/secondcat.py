def secondcat(wave,obj1,obj2,sig1,sig2,secondtimeflag,wavebeg,waveend,brscale):
    import numpy as np
    import matplotlib.pyplot as plt
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.yesno import yesno
    print('\nCombining two calibrated spectra for second-order correction')
    waveblue=wave
    wavered=wave
    fluxblue=obj1
    fluxred=obj2
    fluxtemp=fluxred
    if (not secondtimeflag):
        print('Plotting blue side as blue, red side as red')
        indexblue=womget_element(waveblue,wavered[0])
        indexred=womget_element(wavered,waveblue[-1])
        meanblue=np.mean(fluxblue[indexblue:])
        meanred=np.mean(fluxred[0:indexred])
        if ((meanblue/meanred) > 1.2) or (meanblue/meanred < 0.8):
            print('Averages very different, scaling red to blue for plot')
            scale=meanblue/meanred
            print('Red multiplied by {}'.format(scale))
            fluxtemp=fluxred*scale
        plt.clf()
        plt.plot(waveblue[indexblue:],fluxblue[indexblue:],drawstyle='steps-mid',color='b')
        plt.plot(wavered[0:indexred],fluxred[:indexred],drawstyle='steps-mid',color='r')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        donescale = False
        while (not donescale):
            print('Change y-scale?')
            answer=yesno('n')
            if (answer == 'y'):
                ymin,ymax=scaleparse()
                plt.ylim([ymin,ymax])
                plt.pause(0.01)
            print('Do again?')
            answer=yesno('n')
            if (answer == 'n'):
                donescale = True
        regiondone = False
        while (not regiondone):
            print('\nOK, mark the two end points of the region to compute average')
            endpoints=plt.ginput(2, timeout=-1)
            element1=endpoints[0][0]
            element2=endpoints[1][0]
            plt.plot(endpoints,'ro')
            if (element1 > element2):
                element1,element2=element2,element1
            binend=womget_element(waveblue,element2)
            binbeg=womget_element(waveblue,element1)
            if (binend > len(waveblue)-1):
                binend=len(waveblue)-1
            if (binbeg < indexblue):
                binbeg=indexblue
            print('\nAre these points OK?')
            answer=yesno('y')
            if (answer == 'y'):
                regiondone=False
            else:
                plt.plot(endpoints,'wo')
        wavebeg=waveblue[binbeg]
        waveend=waveblue[binend]
        binbegred=womget_element(wavered,wavebeg)
        binendred=womget_element(wavered,waveend)
        meanblue=np.mean(fluxblue[binbeg:binend+1])
        meanred=np.mean(fluxred[binbegred:binendred+1])
        print('Average for {}:{}'.format(wavebeg,waveend))
        print('Blue side:  {}'.format(meanblue))
        print('Red side:   {}'.format(meanred))
        brscale=inputter_single('Scale to blue or red? (b/r) ','br')
        if (brscale == 'b'):
            print('Scaling to blue by {}'.format(meanblue/meanred))
            fluxred=fluxred*meanblue/meanred
        else:
            print('Scaling to red by {}'.format(meanred/meanblue))
            fluxblue=fluxblue*meanred/meanblue
        print('\nPlotting blue side as blue, red side as red')
        plt.clf()
        plt.plot(waveblue[binbeg:binend+1],fluxblue[binbeg:binend+1],drawstyle='steps-mid',color='b')
        plt.plot(wavered[binbegred:binendred+1],fluxred[binbegred:binendred+1],drawstyle='steps-mid',color='r')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        donescale = False
        while (not donescale):
            print('Change y-scale?')
            answer=yesno('n')
            if (answer == 'y'):
                ymin,ymax=scaleparse()
                plt.ylim([ymin,ymax])
                plt.pause(0.01)
            print('Do again?')
            answer=yesno('n')
            if (answer == 'n'):
                donescale = True
        regiondone = False
        while (not regiondone):
            print('\nOK, mark the two end points of the region to combine')
            endpoints=plt.ginput(2, timeout=-1)
            element1=endpoints[0][0]
            element2=endpoints[1][0]
            plt.plot(endpoints,'ro')
            if (element1 > element2):
                element1,element2=element2,element1
            binend=womget_element(waveblue,element2)
            binbeg=womget_element(waveblue,element1)
            if (binend > len(waveblue)-1):
                binend=len(waveblue)-1
            if (binbeg < indexblue):
                binbeg=indexblue
            print('\nAre these points OK?')
            answer=yesno('y')
            if (answer == 'y'):
                regiondone=False
            else:
                plt.plot(endpoints,'wo')
        wavebeg=waveblue[binbeg]
        waveend=waveblue[binend]
        binbegred=womget_element(wavered,wavebeg)
        binendred=womget_element(wavered,waveend)
    else:
        binbeg=womget_element(waveblue,wavebeg)
        binend=womget_element(waveblue,waveend)
        binbegred=womget_element(wavered,wavebeg)
        binendred=womget_element(wavered,waveend)
        meanblue=np.mean(fluxblue[binbeg:binend+1])
        meanred=np.mean(fluxred[binbegred:binendred+1])
        if (brscale == 'b'):
            scale=meanblue/meanred
            fluxred=fluxred*scale
            sig2=sig2*scale
        else:
            scale=meanred/meanblue
            fluxblue=fluxblue*scale
            sig1=sig1*scale
    overflux=(fluxblue[binbeg:binend+1] + fluxred[binbegred:binendred+1]) / 2.0
    oversig=np.sqrt(sig1[binbeg:binend+1]**2 + sig2[binbegred:binendred+1]**2)
    newwave=waveblue.copy()
    newflux=fluxblue.copy()
    newsig=sig1.copy()
    newsig[binbeg:binend+1]=oversig
    newsig[binend+1:]=sig2[binendred+1:]
    newflux[binbeg:binend+1]=overflux
    newflux[binend+1:]=fluxred[binendred+1:]
    plt.clf()
    plt.pause(0.01)
    axarr=fig.subplots(2,sharex=True)
    fig.subplots_adjust(hspace=0)
    axarr[0].plot(waveblue[binbeg:binend+1],fluxblue[binbeg:binend+1],drawstyle='steps-mid',color='b')
    axarr[0].plot(wavered[binbegred:binendred+1],fluxred[binbegred:binendred+1],drawstyle='steps-mid',color='r')
    axarr[0].plot(waveblue[binbeg:binend+1],overflux,drawstyle='steps-mid',color='k')
    axarr[0].set_title('Overlap region--with average')
    axarr[0].set_ylabel('Flux')
    plt.pause(0.01)
    axarr[1].plot(newwave,newflux,drawstyle='steps-mid',color='k')
    axarr[1].plot(newwave[binbeg:binend+1],newflux[binbeg:binend+1],drawstyle='steps-mid',color='r')
    axarr[1].set_xlabel('Wavelength')
    axarr[1].set_ylabel('Flux')
    plt.pause(0.01)
    return newflux, newsig, wavebeg, waveend, brscale
    
