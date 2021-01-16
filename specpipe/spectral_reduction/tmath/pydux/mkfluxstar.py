def mkfluxstar(fluxfile,gratcode):
    import pdb
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Cursor
    from astropy.io import fits
    from tmath.wombat.get_screen_size import get_screen_size
    from tmath.wombat.getmswave import getmswave
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.inputter import inputter
    from tmath.wombat.womscipyrebin import womscipyrebin
    from tmath.pydux.obs_extinction import obs_extinction
    from tmath.pydux.wave_telluric import wave_telluric
    from tmath.pydux.fitspl import fitspl
    from tmath.pydux.fitspl import fitspl_dev
    from tmath.pydux.finalscaler import finalscaler
    from tmath.pydux.abcalc import abcalc
    from tmath.pydux.pacheck import pacheck
    
    screen_width, screen_height=get_screen_size()
    plt.ion()
    screenpos='+{}+{}'.format(int(screen_width*0.2),int(screen_height*0.05))
    fig=plt.figure()
    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('Flux Star')
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

    #extinction terms from Allen, 3rd edition
    extwave= [2400.,2600.,2800.,3000.,3200.,3400.,3600.,3800., \
              4000.,4500.,5000.,5500.,6000.,6500.,7000.,8000., \
              9000.,10000.,12000.,14000.]
    extvals=[68.0,89.0,36.0,4.5,1.30,0.84,0.68,0.55,0.46,0.31, \
             0.23,0.195,0.170,0.126,0.092,0.062,0.048,0.039, \
             0.028,0.021]
                 
    fitsfile=fits.open(fluxfile)
    rawdata=fitsfile[0].data
    head=fitsfile[0].header
    num_apertures=rawdata.shape[1]
        
    wavearr=np.zeros((rawdata.shape[2],rawdata.shape[1]))
    objectname=head['OBJECT']
    airmass=float(head['AIRMASS'])
    exptime=float(head['EXPTIME'])
    head=pacheck(head)
    if (exptime < 1):
        exptime=1.
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
            ap_choice=inputter('Which aperture do you want to use for the flux star? ','int',False)
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
    plt.ylim([ymin,ymax])
    plt.pause(0.01)
    print('\nPlotting optimal as black, normal as red')
    extract=inputter_single('Do you want to use (n)ormal or (o)ptimal extraction? ','no')
    if (extract == 'o'):
        star=rawdata[0,ap_choice,:]
    else:
        star=rawdata[1,ap_choice,:]
    # fix IRAF weirdness where data can go negative in dispcor interpolation
    star[np.where(star < 0)] = 0.01
    wave=wavearr[:,ap_choice]
    extinction=womscipyrebin(extwave,extvals,wave)
    extfactor=np.exp(extinction*sitefactor*airmass)
    star=star*extfactor
    abcurve=abcalc(wave,objectname)
    wdata=10.0**(0.4*abcurve)
    fstar=star*wdata/exptime
    print('Now fit the continuum manually\n')
    plt.clf()
    ax=fig.add_subplot(111)
    cursor = Cursor(ax, useblit=True, color='k', linewidth=1 )
    airlimit=1.5

    splineresult=fitspl_dev(wave,np.log10(fstar),(airmass>airlimit),fig, cal='fluxstar{}'.format(gratcode))
    splineresult=10**(splineresult)
    plt.cla()
    plt.plot(wave,fstar,drawstyle='steps-mid', color='r')
    plt.plot(wave,splineresult,drawstyle='steps-mid')
    plt.pause(0.01)
    delkeylist=['WAT0_001', 'WAT1_001', 'WAT2_001', 'WAT0_002', \
                'WAT1_002', 'WAT2_002', 'WAT3_001', 'WAT2_003', \
                'CTYPE1', 'CTYPE2', 'CTYPE3', 'CD1_1', 'CD2_2', \
                'CD3_3', 'LTM1_1', 'LTM2_2', 'LTM3_3', 'WCSDIM']
    for k in delkeylist:
        try:
            head.remove(k)
        except KeyError:
            pass
    head.set('CRPIX1', 1)
    head.set('CRVAL1', wave[0])
    head.set('CDELT1', wave[1]-wave[0])
    head.set('CTYPE1', 'LINEAR')
    head.set('FLUXFILE',fluxfile)
    outfile='fluxstar'+gratcode+'.fits'
    print('Writing data to {}'.format(outfile))
    outhdu=fits.PrimaryHDU(splineresult)
    hdul=fits.HDUList([outhdu])
    hdul[0].header=head.copy()
    hdul.writeto(outfile,overwrite=True)
    print('mkfluxstar')
    print(fluxfile, gratcode)
    return

