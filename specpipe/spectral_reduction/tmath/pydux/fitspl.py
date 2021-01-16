

def fitspl(wave,flux,airlimit,fig, cal = None):
    """fit spline to spectrum"""
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import splrep,splev
    import tmath.wombat.womconfig as womconfig
    from tmath.wombat.womget_element import womget_element
    from tmath.pydux.wave_telluric import wave_telluric
    from tmath.wombat.yesno import yesno
    from tmath.wombat.onclick import onclick
    from tmath.wombat.onkeypress import onkeypress
    import glob

#    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    
    # starting points for spline
    bandpts = np.array([3000, 3050, 3090, 3200, 3430, 3450, 3500, 3550, 3600, \
               3650, 3700, 3767, 3863, 3945, 4025, 4144, 4200, 4250, \
               4280, 4390, 4450, 4500, 4600, 4655, 4717, 4750, 4908, \
               4950, 5000, 5050, 5100, 5150, 5200, 5250, 5280, 5350, \
               5387, 5439, 5500, 5550, 6100, 6150, 6400, 6430, 6650, \
               6700, 6750, 6800, 7450, 7500, 7550, 8420, 8460, 8520, \
               8570, 8600, 8725, 8770, 9910, 10000, 10200, 10300, \
               10400, 10500, 10600, 10700])
    locsinrange=np.logical_and((bandpts > wave[10]),(bandpts < wave[-10]))
    useband=bandpts[locsinrange]
    for i,_ in enumerate(useband):
        index=womget_element(wave,useband[i])
        useband[i]=index
    # useband now has indices of wavelength positions

    if (min(flux) < 0):
        flux[np.where(flux < 0)]=0.0
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid',color='k')
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    if (airlimit):
        loc=wave_telluric(wave,'high')
    else:
        loc=wave_telluric(wave,'low')
    #invert loc and convert all good sections to np.nan so they won't plot
    loc=np.invert(loc)
    wavetell=wave.copy()
    fluxtell=flux.copy()
    wavetell[loc]=np.nan
    fluxtell[loc]=np.nan
    plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')

    if '../../master_files/' + cal + '_splpts_master.txt' not in glob.glob('../../master_files/*'):
        womconfig.nsplinepoints=len(useband)
        womconfig.tmpsplptsx=wave[useband].copy().tolist()
        womconfig.tmpsplptsy=[]
        for i,_ in enumerate(useband):
            womconfig.tmpsplptsy.append(np.median(flux[useband[i]-2:useband[i]+3]))
        spline=splrep(womconfig.tmpsplptsx,womconfig.tmpsplptsy,k=3)
    else:
        masterx, mastery = np.genfromtxt('../../master_files/' + cal + '_splpts_master.txt')
        womconfig.nsplinepoints=len(masterx)
        womconfig.tmpsplptsx=masterx
        womconfig.tmpsplptsy=mastery
        spline=splrep(womconfig.tmpsplptsx,womconfig.tmpsplptsy,k=3)

    splineresult=splev(wave,spline)
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid',color='k')
    plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
    if (len(womconfig.tmpsplptsx) > 0):
        plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
        plt.plot(wave,splineresult,color='g')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.pause(0.01)
    done=False
    print('Is this OK? ')
    answer=yesno('n')
    if (answer == 'y'):
        done = True
    splptsy=[z for _,z in sorted(zip(womconfig.tmpsplptsx,womconfig.tmpsplptsy))]
    splptsx=sorted(womconfig.tmpsplptsx)
    while (not done):
        plotdone=False
        while (not plotdone):
            plt.cla()
            plt.plot(wave,flux,drawstyle='steps-mid',color='k')
            plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
            if (len(womconfig.tmpsplptsx) > 0):
                plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
                plt.plot(wave,splineresult,color='g')
            plt.xlabel('Wavelength')
            plt.ylabel('Flux')
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.pause(0.01)
            print('Change scale? ')
            answer=yesno('n')
            if (answer == 'y'):
                plt.cla()
                plt.plot(wave,flux,drawstyle='steps-mid',color='k')
                plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
                if (len(womconfig.tmpsplptsx) > 0):
                    plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
                    plt.plot(wave,splineresult,color='g')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                print('Click corners of box to change plot scale')
                newlims=plt.ginput(2, timeout=-1)
                xmin=newlims[0][0]
                ymin=newlims[0][1]
                xmax=newlims[1][0]
                ymax=newlims[1][1]
                plt.cla()
                plt.plot(wave,flux,drawstyle='steps-mid',color='k')
                plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
                if (len(womconfig.tmpsplptsx) > 0):
                    plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
                    plt.plot(wave,splineresult,color='g')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                plotdone=True
            else:
                plotdone = True
            
                
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        cid2 = fig.canvas.mpl_connect('key_press_event', onkeypress)
        print('\nClick on continuum points for spline fit.')
        print('Left button or (a)    = add point')
        print('Middle button or (s) = delete point')
        print('Right button or (d)  = done\n')
        womconfig.pflag=''
        while (womconfig.pflag != 'done'):
            plt.pause(0.01)
        fig.canvas.mpl_disconnect(cid)
        fig.canvas.mpl_disconnect(cid2)

        splptsy=[z for _,z in sorted(zip(womconfig.tmpsplptsx,womconfig.tmpsplptsy))]
        splptsx=sorted(womconfig.tmpsplptsx)
        spline=splrep(splptsx,splptsy,k=3)
        splineresult=splev(wave,spline)
        plt.plot(wave,splineresult,drawstyle='steps-mid',color='g')
        plt.pause(0.01)
        print('Is this fit OK? ')
        answer=yesno('y')
        if (answer == 'y'):
            done=True

    if cal != None:
        np.savetxt('../../master_files/' + cal + '_splpts_master.txt', [splptsx,splptsy])
    return splineresult


def fitspl_dev(wave,flux,airlimit,fig, cal=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import splrep,splev
    import tmath.wombat.womconfig as womconfig
    from tmath.wombat.womget_element import womget_element
    from tmath.pydux.wave_telluric import wave_telluric
    from tmath.wombat.yesno import yesno
    from tmath.wombat.onclick import onclick
    from tmath.wombat.onkeypress import onkeypress
    import glob
    """fit spline to spectrum"""    
    # starting points for spline
    # bandpts = np.array([3000, 3050, 3090, 3200, 3430, 3450, 3500, 3550, 3600, \
    #            3650, 3700, 3767, 3863, 3945, 4025, 4144, 4200, 4250, \
    #            4280, 4390, 4450, 4500, 4600, 4655, 4717, 4750, 4908, \
    #            4950, 5000, 5050, 5100, 5150, 5200, 5250, 5280, 5350, \
    #            5387, 5439, 5500, 5550, 6100, 6150, 6400, 6430, 6650, \
    #            6700, 6750, 6800, 7450, 7500, 7550, 8420, 8460, 8520, \
    #            8570, 8600, 8725, 8770, 9910, 10000, 10200, 10300, \
    #            10400, 10500, 10600, 10700])
    binWidth = 60
    bandpts = np.arange(3000,11000,binWidth)
    locsinrange=np.logical_and((bandpts > wave[10]),(bandpts < wave[-10]))
    #useband=bandpts[locsinrange]

    # now mask based on spectral features
    # you can add mask features here, and the feature name can
    # be anything, they aren't used explicitly
    featureMask = {}

    # empirical masks from Dimitriadis
    featureMask['feature1'] = [3722.56-10.0/2.,3722.56+10.0/2.]
    featureMask['feature2'] = [3736.90-15.0/2.,3736.90+15.0/2.]
    featureMask['feature3'] = [3752.00-15.0/2.,3752.00+15.0/2.]
    featureMask['feature4'] = [3772.00-17.0/2.,3772.00+17.0/2.]
    featureMask['feature5'] = [3800.00-18.0/2.,3800.00+18.0/2.]
    featureMask['feature6'] = [3835.38-20.0/2.,3835.38+20.0/2.]
    featureMask['feature7'] = [3889.05-24.0/2.,3889.05+24.0/2.]
    featureMask['feature8'] = [3933.66-16.0/2.,3933.66+16.0/2.]
    featureMask['feature9'] = [3970.07-18.0/2.,3970.07+18.0/2.]
    featureMask['feature10'] = [4101.74-28.0/2.,4101.74+28.0/2.]
    featureMask['feature11'] = [4340.46-30.0/2.,4340.46+30.0/2.]
    featureMask['feature12'] = [4471.48-20.0/2.,4471.48+20.0/2.]
    featureMask['feature13'] = [4685.70-30.0/2.,4685.70+30.0/2.]
    featureMask['feature14'] = [4861.36-35.0/2.,4861.36+35.0/2.]
    featureMask['feature15'] = [5411.52-35.0/2.,5411.52+35.0/2.]
    featureMask['feature16'] = [5889.95-32.0/2.,5889.95+32.0/2.]
    featureMask['feature17'] = [6562.85-40.0/2.,6562.85+40.0/2.]
    featureMask['feature18'] = [8498.02-25.0/2.,8498.02+25.0/2.]
    featureMask['feature19'] = [8542.09-25.0/2.,8542.09+25.0/2.]
    featureMask['feature20'] = [8662.14-25.0/2.,8662.14+25.0/2.]
    featureMask['feature21'] = [8763.96-20.0/2.,8763.96+20.0/2.]
    featureMask['feature22'] = [8865.75-25.0/2.,8865.75+25.0/2.]
    featureMask['feature23'] = [9010.00-28.0/2.,9010.00+28.0/2.]
    featureMask['feature24'] = [9213.90-30.0/2.,9213.90+30.0/2.]
    featureMask['feature25'] = [9545.97-34.0/2.,9545.97+34.0/2.]
    featureMask['telluric1'] = [3216.0-binWidth/2., 3420.0+binWidth/2.]
    #featureMask['telluric2'] = [5600.0-binWidth/2., 6050.0+binWidth/2.]
    featureMask['telluric3'] = [6250.0-binWidth/2., 6360.0+binWidth/2.]
    featureMask['telluric4'] = [6450.0-binWidth/2., 6530.0+binWidth/2.]
    featureMask['telluric5'] = [6840.0-binWidth/2., 7410.0+binWidth/2.]
    featureMask['telluric6'] = [7560.0-binWidth/2., 8410.0+binWidth/2.]
    featureMask['telluric7'] = [8925.0-binWidth/2., 9900.0+binWidth/2.]

    # inital testing masks from XIDL
    # featureMask['balmer1'] =   [3714.0-binWidth/2., 3723.0+binWidth/2.]
    # featureMask['balmer2'] =   [3650.0-binWidth/2., 3820.0+binWidth/2.]
    # featureMask['balmer3'] =   [3725.0-binWidth/2., 3735.0+binWidth/2.]
    # featureMask['balmer4'] =   [3740.0-binWidth/2., 3755.0+binWidth/2.]
    # featureMask['balmer5'] =   [3760.0-binWidth/2., 3775.0+binWidth/2.]
    # featureMask['balmer6'] =   [3785.0-binWidth/2., 3806.0+binWidth/2.]
    # featureMask['balmer7'] =   [3810.0-binWidth/2., 3820.0+binWidth/2.]
    # featureMask['balmer8'] =   [3824.0-binWidth/2., 3841.0+binWidth/2.]
    # featureMask['balmer9'] =   [3880.0-binWidth/2., 3895.0+binWidth/2.]
    # featureMask['balmer10'] =  [3957.0-binWidth/2., 3979.0+binWidth/2.]
    # featureMask['balmer11'] =  [4000.0-binWidth/2., 4030.0+binWidth/2.]
    # featureMask['balmer12'] =  [4087.0-binWidth/2., 4120.0+binWidth/2.]
    # featureMask['balmer13'] =  [4135.0-binWidth/2., 4145.0+binWidth/2.]
    # featureMask['balmer14'] =  [4328.0-binWidth/2., 4355.0+binWidth/2.]
    # featureMask['balmer15'] =  [4677.0-binWidth/2., 4692.0+binWidth/2.]
    # featureMask['balmer16'] =  [4830.0-binWidth/2., 4931.0+binWidth/2.]
    # featureMask['balmer17'] =  [5402.0-binWidth/2., 5417.0+binWidth/2.]
    # featureMask['balmer18'] =  [6535.0-binWidth/2., 6590.0+binWidth/2.]
    # featureMask['telluric1'] = [3216.0-binWidth/2., 3420.0+binWidth/2.]
    # #featureMask['telluric2'] = [5600.0-binWidth/2., 6050.0+binWidth/2.]
    # featureMask['telluric3'] = [6250.0-binWidth/2., 6360.0+binWidth/2.]
    # featureMask['telluric4'] = [6450.0-binWidth/2., 6530.0+binWidth/2.]
    # featureMask['telluric5'] = [6840.0-binWidth/2., 7410.0+binWidth/2.]
    # featureMask['telluric6'] = [7560.0-binWidth/2., 8410.0+binWidth/2.]
    # featureMask['telluric7'] = [8925.0-binWidth/2., 9900.0+binWidth/2.]

    # generate a mask via a boolean array that excludes regions in stellar features
    maskInds = np.full(len(bandpts),True,dtype=bool)
    for i,key in enumerate(featureMask.keys()):
        featureMaskInds = np.logical_or((bandpts<featureMask[key][0]), 
                                        (bandpts>featureMask[key][1]))
        maskInds *= featureMaskInds
    fullMask = locsinrange * maskInds

    # apply the mask
    useband = bandpts[fullMask]



    for i,_ in enumerate(useband):
        index=womget_element(wave,useband[i])
        useband[i]=index
    # useband now has indices of wavelength positions

    if (min(flux) < 0):
        flux[np.where(flux < 0)]=0.0
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid',color='k')
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    if (airlimit):
        loc=wave_telluric(wave,'high')
    else:
        loc=wave_telluric(wave,'low')
    #invert loc and convert all good sections to np.nan so they won't plot
    loc=np.invert(loc)
    wavetell=wave.copy()
    fluxtell=flux.copy()
    wavetell[loc]=np.nan
    fluxtell[loc]=np.nan
    plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')

    if '../../master_files/' + cal + '_splpts_master.txt' not in glob.glob('../../master_files/*'):
        womconfig.nsplinepoints=len(useband)
        womconfig.tmpsplptsx=wave[useband].copy().tolist()
        womconfig.tmpsplptsy=[]
        for i,_ in enumerate(useband):
            womconfig.tmpsplptsy.append(np.median(flux[useband[i]-2:useband[i]+3]))
        spline=splrep(womconfig.tmpsplptsx,womconfig.tmpsplptsy,k=3)
    else:
        masterx, mastery = np.genfromtxt('../../master_files/' + cal + '_splpts_master.txt')
        womconfig.nsplinepoints=len(masterx)
        womconfig.tmpsplptsx=list(masterx)
        womconfig.tmpsplptsy=list(mastery)
        # print (type(womconfig.tmpsplptsx))
        spline=splrep(womconfig.tmpsplptsx,womconfig.tmpsplptsy,k=3)

    splineresult=splev(wave,spline)
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid',color='k')
    plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
    if (len(womconfig.tmpsplptsx) > 0):
        plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
        plt.plot(wave,splineresult,color='g')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.pause(0.01)
    done=False
    print('Is this OK? ')
    answer=yesno('n')
    if (answer == 'y'):
        done = True
    splptsy=[z for _,z in sorted(zip(womconfig.tmpsplptsx,womconfig.tmpsplptsy))]
    splptsx=sorted(womconfig.tmpsplptsx)
    while (not done):
        plotdone=False
        while (not plotdone):
            plt.cla()
            plt.plot(wave,flux,drawstyle='steps-mid',color='k')
            plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
            if (len(womconfig.tmpsplptsx) > 0):
                plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
                plt.plot(wave,splineresult,color='g')
            plt.xlabel('Wavelength')
            plt.ylabel('Flux')
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.pause(0.01)
            print('Change scale? ')
            answer=yesno('n')
            if (answer == 'y'):
                plt.cla()
                plt.plot(wave,flux,drawstyle='steps-mid',color='k')
                plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
                if (len(womconfig.tmpsplptsx) > 0):
                    plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
                    plt.plot(wave,splineresult,color='g')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                print('Click corners of box to change plot scale')
                newlims=plt.ginput(2, timeout=-1)
                xmin=newlims[0][0]
                ymin=newlims[0][1]
                xmax=newlims[1][0]
                ymax=newlims[1][1]
                plt.cla()
                plt.plot(wave,flux,drawstyle='steps-mid',color='k')
                plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
                if (len(womconfig.tmpsplptsx) > 0):
                    plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
                    plt.plot(wave,splineresult,color='g')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                plotdone=True
            else:
                plotdone = True
            
                
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        cid2 = fig.canvas.mpl_connect('key_press_event', onkeypress)
        print('\nClick on continuum points for spline fit.')
        print('Left button or (a)    = add point')
        print('Middle button or (s) = delete point')
        print('Right button or (d)  = done\n')
        womconfig.pflag=''
        while (womconfig.pflag != 'done'):
            plt.pause(0.01)
        fig.canvas.mpl_disconnect(cid)
        fig.canvas.mpl_disconnect(cid2)

        splptsy=[z for _,z in sorted(zip(womconfig.tmpsplptsx,womconfig.tmpsplptsy))]
        splptsx=sorted(womconfig.tmpsplptsx)
        spline=splrep(splptsx,splptsy,k=3)
        splineresult=splev(wave,spline)
        plt.plot(wave,splineresult,drawstyle='steps-mid',color='g')
        plt.pause(0.01)
        print('Is this fit OK? ')
        answer=yesno('y')
        if (answer == 'y'):
            done=True

    if cal != None:
        np.savetxt('../../master_files/' + cal + '_splpts_master.txt', [splptsx,splptsy])
        
    return splineresult