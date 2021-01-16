def telluric_remove(bstarwave, bstar, bairmass, wave, object, airmass, variance, spectrum, shift=None):
    import numpy as np
    import pdb
    import matplotlib.pyplot as plt
    from tmath.wombat.inputter import inputter
    from tmath.wombat.yesno import yesno
    from tmath.wombat.womscipyrebin import womscipyrebin
    from tmath.wombat.womget_element import womget_element
    from tmath.pydux.xcor import xcor
    from tmath.pydux.finalscaler import finalscaler
    bstartmp=womscipyrebin(bstarwave,bstar,wave)
#    plt.cla()
#    plt.plot(bstarwave,bstartmp)
#    plt.pause(0.01)
#    answer=yesno('y')
    print('\nThe ratio of airmasses (object/B-star) is {}'.format(airmass/bairmass))
    if (airmass/bairmass > 3.0) or (airmass/bairmass < 0.33):
        print('\nWARNING: OBJECT AND B-STAR HAVE WILDLY DIFFERENT')
        print('AIRMASSES: ATMOSPHERIC BAND DIVISION MAY BE LOUSY\n')

    wmin=wave[0]
    wmax=wave[-1]
    npix=len(object)
    wdelt=wave[1]-wave[0]
    print('wdelt',wdelt)
    lag=np.zeros(3)
    lagflag=[False]*3
    xfactor=10
    maxlag=200

    if not shift:
        print('\nCross-correlating object with B-star spectrum\n')
        fig=plt.figure()
        axarr=fig.subplots(2,1)
        if (wmin < 6200) and (wmax > 6400) and (wmax < 6900):
            indblue=womget_element(wave,6200)
            indred=womget_element(wave,6400)
            lag[0]=xcor(object[indblue:indred+1],bstartmp[indblue:indred+1],xfactor,maxlag)
            lagflag[0]=True
            print('The shift at the 6250A band is {} angstroms'.format(lag[0]*wdelt))
        if (wmin < 6800) and (wmax > 6500):
            indblue=womget_element(wave,6800)
            indred=womget_element(wave,6950)
            scale = 1./np.max(object[indblue:indred+1])
            obb=scale*object[indblue:indred+1]
            bb=bstartmp[indblue:indred+1]
            lag[1]=xcor(obb,bb,xfactor,maxlag)
            lagflag[1]=True
            print('The shift at the B band is {} angstroms'.format(lag[1]*wdelt))

            plt.cla()
            # ymin,ymax=finalscaler(object)
            # plt.plot(wave,object,drawstyle='steps-mid',color='r')
            # plt.plot(wave,newobject,drawstyle='steps-mid',color='k')
            ymin,ymax=finalscaler(bstartmp[indblue:indred+1])
            axarr[0].plot(wave[indblue:indred+1], scale*object[indblue:indred+1],drawstyle='steps-mid',color='r')
            axarr[0].plot(wave[indblue:indred+1], bstartmp[indblue:indred+1],drawstyle='steps-mid',color='k')
            axarr[0].plot(wave[indblue:indred+1]+lag[1]*wdelt, bstartmp[indblue:indred+1],drawstyle='steps-mid',color='g')
            plt.pause(0.01)

        if (wmin < 7500) and (wmax > 8000):
            
            indblue=womget_element(wave,7500)
            indred=womget_element(wave,8000)
            scale = 1./np.max(object[indblue:indred+1])
            lag[2]=xcor(scale*object[indblue:indred+1],bstartmp[indblue:indred+1],xfactor,maxlag)
            print('The shift at the A band is {} angstroms'.format(lag[2]*wdelt))
            lagflag[2]=True
            # ymin,ymax=finalscaler(object)
            # plt.plot(wave,object,drawstyle='steps-mid',color='r')
            # plt.plot(wave,newobject,drawstyle='steps-mid',color='k')
            ymin,ymax=finalscaler(bstartmp[indblue:indred+1])
            axarr[1].plot(wave[indblue:indred+1], scale*object[indblue:indred+1],drawstyle='steps-mid',color='r')
            axarr[1].plot(wave[indblue:indred+1], bstartmp[indblue:indred+1],drawstyle='steps-mid',color='k')
            axarr[1].plot(wave[indblue:indred+1]+lag[2]*wdelt, bstartmp[indblue:indred+1],drawstyle='steps-mid',color='g')
            plt.pause(0.01)
            check=inputter('Check plot [enter when done]: ','string',False)

        if (sum(lagflag) > 0):
            avglag=np.sum(lag)/sum(lagflag)
            angshift=avglag*wdelt
            print('The mean shift is {} Angstroms'.format(angshift))
        else:
            angshift=0.0
        plt.close()
    else:
        angshift=shift
    fig = plt.figure(figsize = [9,5])
    telluric_done = False
    bstartmpcopy=bstartmp.copy()
    while (not telluric_done):
        print('Applying a shift of {} Angstroms'.format(angshift))
        bstartmp=bstartmpcopy.copy()
        tmp=womscipyrebin(wave+angshift,bstartmp,wave)
        bstartmp=tmp.copy()
        bstartmp=bstartmp**((airmass/bairmass)**0.55)
        # newobject=object/bstartmp
        newobject=spectrum/bstartmp
        bvar=variance/bstartmp
        print('\nPlotting before and after atmospheric band correction\n')
        plt.cla()
        # ymin,ymax=finalscaler(object)
        # plt.plot(wave,object,drawstyle='steps-mid',color='r')
        # plt.plot(wave,newobject,drawstyle='steps-mid',color='k')
        ymin,ymax=finalscaler(spectrum)
        plt.plot(wave,spectrum,drawstyle='steps-mid',color='r')
        plt.plot(wave,newobject,drawstyle='steps-mid',color='k')
        plt.ylim([ymin,ymax])
        plt.pause(0.01)

        if not shift: 
            print('Is this OK?')
            answer=yesno('y')
            if (answer == 'n'):
                angshift=inputter('Enter B-star shift in Angstroms: ','float',False)
            else:
                telluric_done = True
        else:
            check=inputter('Check plot [enter when done]: ','string',False)
            telluric_done = True
    plt.close()
    return newobject, bvar, angshift
            
