def womblo(hop,fig):
    """blotch out bad data"""
    import logging
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import splrep,splev
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.waveparse import waveparse
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.yesno import yesno
    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    print('Select regions to blotch')
    wave=hop[0].wave.copy()
    flux=hop[0].flux.copy()
    done=False
    while (not done):
        plt.cla()
        plt.plot(wave,flux, drawstyle='steps-mid')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title(hop[0].obname)
        wavesub,fluxsub,mode=womwaverange(wave,flux,'none')
        wavebind=womget_element(wave,wavesub[0])
        waverind=womget_element(wave,wavesub[-1])
        plt.cla()
        plt.plot(wave[wavebind:waverind+1],flux[wavebind:waverind+1], \
                 drawstyle='steps-mid')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title(hop[0].obname)
        plt.pause(0.01)
        print('Do you want to enter blotch wavelengths by hand (w),')
        print('mark points (m), fit a spline (s), or quit (q)?')
        choice=inputter_single('(w/m/s/q): ','wmsq')
        if (choice == 'w') or (choice == 'm'):
            blotchgood=False
            while (not blotchgood):
                wavechoicedone=False
                while (not wavechoicedone):
                    if (choice == 'w'):
                        waveselb,waveselr=waveparse()
                    else:
                        print('Mark the two endpoints of the blotch region')
                        endpoints=plt.ginput(2,timeout=-1)
                        waveselb=endpoints[0][0]
                        waveselr=endpoints[1][0]
                    if (waveselb > waveselr):
                        waveselb,waveselr=waveselr,waveselb
                    waveselbind=womget_element(wave,waveselb)
                    waveselrind=womget_element(wave,waveselr)
                    print(waveselb, waveselr,waveselbind,waveselrind)
                    if (waveselbind == 0) or (waveselrind == (len(wave)-1)):
                        print('Wavelengths incorrect--too close to endpoints')
                    else:
                        wavechoicedone=True
                contblue=flux[waveselbind-1]
                contred=flux[waveselrind+1]
                delta=(contred-contblue)/(waveselrind-waveselbind+1)
                fluxcor=flux.copy()
                for i in range(waveselbind,waveselrind+1):
                    fluxcor[i]=contblue+ (i-waveselbind+1)*delta
                plt.plot(wave[wavebind:waverind+1],fluxcor[wavebind:waverind+1], \
                 drawstyle='steps-mid')
                plt.pause(0.01)
                print('Is this acceptable')
                answer=yesno('y')
                if (answer == 'y'):
                    flux=fluxcor.copy()
                    blotchgood=True
                    logging.info('File {} blotched from {} to {}'.format(hop[0].obname, wave[waveselbind], wave[waveselrind]))
        elif (choice == 's'):
            xmin,xmax=plt.xlim()
            ymin,ymax=plt.ylim()
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            nsplinepoints=0
            tmpsplptsx=[]
            tmpsplptsy=[]

            spldone=False
            while (not spldone):
                plt.cla()
                plt.plot(wave[wavebind:waverind+1],flux[wavebind:waverind+1], \
                         drawstyle='steps-mid')
                if (len(tmpsplptsx) > 0):
                    plt.plot(tmpsplptsx,tmpsplptsy,'ro')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.title(hop[0].obname)
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                cid = fig.canvas.mpl_connect('button_press_event', onclick)
                print('\nClick on continuum points for spline fit.')
                print('Spline will replace values between first and last point')
                print('Left button    = add point')
                print('Middle button  = delete point')
                print('Right button   = done\n')
                pflag=''
                while (pflag != 'done'):
                    plt.pause(0.01)
                fig.canvas.mpl_disconnect(cid)

                splptsy=[z for _,z in sorted(zip(tmpsplptsx,tmpsplptsy))]
                splptsx=sorted(tmpsplptsx)
                spline=splrep(splptsx,splptsy,k=3)
                splblueindex=womget_element(wave,splptsx[0])
                splredindex=womget_element(wave,splptsx[-1])
                splwave=wave[splblueindex:splredindex+1].copy()
                splineresult=splev(splwave,spline)
                fluxcor=flux.copy()
                fluxcor[splblueindex:splredindex+1]=splineresult.copy()
                plt.plot(splwave,splineresult,drawstyle='steps-mid')
                print('Is this acceptable')
                answer=yesno('y')
                if (answer == 'y'):
                    flux=fluxcor.copy()
                    spldone=True
                    logging.info('File {} blotched  with spline from {} to {}'.format(hop[0].obname, wave[splblueindex], wave[splredindex]))
        else:
            done=True        
        print('Do another region?')
        another=yesno('n')
        if (another == 'n'):
            done=True
    hop[0].flux=flux.copy()
    return hop

