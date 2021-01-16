def mkbstar(bfile,gratcode):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Cursor
    from astropy.io import fits
    from tmath.wombat.get_screen_size import get_screen_size
    from tmath.wombat.getmswave import getmswave
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.inputter import inputter
    from tmath.wombat.yesno import yesno
    from tmath.pydux.obs_extinction import obs_extinction
    from tmath.pydux.wave_telluric import wave_telluric
    from tmath.pydux.fitspl import fitspl,fitspl_dev
    from tmath.pydux.bblo import bblo
    from tmath.pydux.finalscaler import finalscaler
    from tmath.pydux.pacheck import pacheck

    screen_width, screen_height=get_screen_size()
    plt.ion()
    screenpos='+{}+{}'.format(int(screen_width*0.2),int(screen_height*0.05))
    fig=plt.figure()

    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('B Star')
    fig.set_size_inches(8,5)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    # this should bring the window to the top, but it doesn't
    # wm=plt.get_current_fig_manager()
    #fig.canvas.manager.window.attributes('-topmost', 1)
    # blah=wm.window.attributes('-topmost', 1)

    #plt.pause(0.2)
    #fig.canvas.manager.window.attributes('-topmost', 0)
    # blah=wm.window.attributes('-topmost', 0)

    fitsfile=fits.open(bfile)
    rawdata=fitsfile[0].data
    head=fitsfile[0].header
    num_apertures=rawdata.shape[1]
    wavearr=np.zeros((rawdata.shape[2],rawdata.shape[1]))
    objectname=head['OBJECT']
    airmass=float(head['AIRMASS'])
    exptime=float(head['EXPTIME'])
    
    if (exptime < 1):
        exptime=1.
    head=pacheck(head)
    observat=head['OBSERVAT'].strip().lower()
    sitefactor=obs_extinction(observat)
    for i in range(0,num_apertures):
        wavearr[:,i]=getmswave(head,i)
    if (wavearr[-1,0] < 3000):
        print('************************************************')
        print('Spectrum not wavelength calibrated---bailing out')
        print('************************************************')
        sys.exit(1)
    ap_choice=-1
    if (np.mean(rawdata[0,0,:]) < 1e-7):
        rawdata=rawdata*1e15
    # scale to rational numbers for fit (rational as in sane,
    # not opposed to irrational)
    if (num_apertures != 1):
        axarr=fig.subplots(num_apertures,sharex=True)
        fig.subplots_adjust(hspace=0)
        for i in range(0,num_apertures):
            x=wavearr[:,i]
            y=rawdata[0,i,:]
            axarr[i].clear()

            plt.pause(0.01)
            axarr[i].plot(x,y,drawstyle='steps-mid')
            axarr[i].set_ylabel('Aperture {}'.format(i))
            plt.pause(0.01)
        while (ap_choice < 0) or (ap_choice > num_apertures-1):
            ap_choice=inputter('Which aperture do you want to use for the B star? ','int',False)
        fig.clf()
    else:
        ap_choice=0    
    ymin, ymax=finalscaler(rawdata[0,ap_choice,:])
    fig.clf()
    x1=wavearr[:,ap_choice]
    y1=rawdata[1,ap_choice,:]
    y0=rawdata[0,ap_choice,:]
    plt.plot(x1,y1,drawstyle='steps-mid',color='r')
    plt.plot(x1,y0,drawstyle='steps-mid',color='k')
    plt.pause(0.01)
    try:
        plt.ylim([ymin,ymax])
    except ValueError:
        medValue = np.median(y0)
        plt.ylim([1.e-2*medValue,1.e2*medValue])
    plt.pause(0.01)
    print('\nPlotting optimal as black, normal as red')
    extract=inputter_single('Do you want to use (n)ormal or (o)ptimal extraction? ','no')
    if (extract == 'o'):
        bstar=rawdata[0,ap_choice,:]
    else:
        bstar=rawdata[1,ap_choice,:]
    # fix IRAF weirdness where data can go negative in dispcor interpolation
    bstar[np.where(bstar < 0)] = 0.01
    wave=wavearr[:,ap_choice]
    print('\nAirmass: {}\n'.format(airmass))
    airlimit=0.5
    while (airlimit < 1.0) or (airlimit > 10.):
        airlimit=inputter('Above what airmass is considered high? ','float',False)
        
    print('\nNow fit the continuum manually\n')
    plt.clf()
    ax=fig.add_subplot(111)
    cursor = Cursor(ax, useblit=True, color='k', linewidth=1 )
    airlimit=1.5
    splineresult=fitspl_dev(wave,np.log10(bstar),(airmass>airlimit),fig, cal='bstar{}'.format(gratcode))
    splineresult=10**(splineresult)
    plt.cla()
    plt.plot(wave,splineresult,drawstyle='steps-mid')
    plt.pause(0.01)
    bstar=bstar/splineresult
    if (airmass > airlimit):
        loc=wave_telluric(wave,'high')
    else:
        loc=wave_telluric(wave,'low')
    #invert loc and convert all good sections to np.nan so they won't plot
    loc=np.invert(loc)
    bstar[loc]=1.0
    bstar[np.where(bstar > 1.)]=1.0
    bstar[np.where(bstar <= 0.)]=0.01
    plt.cla()
    plt.plot(wave,bstar,drawstyle='steps-mid')
    plt.pause(0.01)
    print('\nDo you want to blotch the B-star?')
    blotch=yesno('n')
    if (blotch == 'y'):
        bstar=bblo(wave,bstar,(airmass > airlimit),fig)
    delkeylist=['WAT0_001', 'WAT1_001', 'WAT2_001', 'WAT0_002', \
                'WAT1_002', 'WAT2_002', 'WAT3_001', 'WAT2_003', \
                'CTYPE1', 'CTYPE2', 'CTYPE3', 'CD1_1', 'CD2_2', \
                'CD3_3', 'LTM1_1', 'LTM2_2', 'LTM3_3', 'WCSDIM']
    for k in delkeylist:
        try:
            head.remove(k)
        except KeyError:
            pass
    plt.cla()
    plt.plot(wave,bstar,drawstyle='steps-mid')
    plt.pause(0.01)
    head.set('CRPIX1', 1)
    head.set('CRVAL1', wave[0])
    head.set('CDELT1', wave[1]-wave[0])
    head.set('CTYPE1', 'LINEAR')
    outfile='bstar'+gratcode+'.fits'
    print('Writing data to {}'.format(outfile))
    outhdu=fits.PrimaryHDU(bstar)
    hdul=fits.HDUList([outhdu])
    hdul[0].header=head.copy()
    hdul.writeto(outfile,overwrite=True)    
    print('mkbstar')
    print(bfile, gratcode)
    return

