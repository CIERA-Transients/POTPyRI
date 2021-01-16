def jdcnv(year, month, day, hour):
    if (month < 3):
        leap=-1
    else:
        leap=0

    julian=day - 32075 + 1461 * (year +4800 + leap)//4 + \
            367*(month -2 -leap*12)//12 - \
            3*((year + 4900 + leap)//100)//4
    julian=float(julian)+hour/24.0 -0.5
    return julian

def getfitsfile(name,fitstype):
    done=False
    while (not done):
        inputfile=input('Name of fits file for the {}? ({} added if necessary) '.format(name,fitstype))
        inputfile=inputfile.strip()
        if (fitstype not in inputfile):
            inputfile=inputfile+fitstype
        try:
            fitsdata=fits.open(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            fitsdata.close()
            done=True
    return inputfile

def getch():
    import termios
    import sys, tty
    def _getch():
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
    return _getch()

def inputter(statement, datatype, defaultflag,*default):
    done=False
    while (not done):
        answer=input(statement)
        if (defaultflag) and (answer == ''):
            return default[0]
        elif (datatype == 'string'):
            return answer
        elif (datatype == 'float'):
            try:
                answer=float(answer)
            except ValueError:
                print('Please enter a number. \n')
            else:
                return answer
                done=True
        elif (datatype == 'int'):
            try:
                answer=int(answer)
            except ValueError:
                print('Please enter a number. \n')
            else:
                return answer
                done=True


def inputter_single(statement,ans_string):
    reply=' '
    while (reply not in ans_string):
        print(statement)
        reply=getch()
        reply=reply.strip()
        reply=reply.lower()
        if len(reply) > 0:
            reply=reply[0]
        else:
            reply=' '
    return reply

def inputter_single_mix(statement,ans_string):
    reply=' '
    while (reply not in ans_string):
        reply=input(statement)
        reply=reply.strip()
        reply=reply[0]
    return reply

                
def yesno(default):
    answer=''
    while (answer != 'y') and (answer != 'n'):
        print('(y)es or (n)o? (default/return = {}) '.format(default))
        reply=getch()
        reply=reply.strip()
        if len(reply) == 0:
            reply=default
        answer=reply.lower()[0]
    return answer

def waverange(wave,flux,mode):
    if (mode == 'none'):
        print('Enter method to select wavelength range')
        mode=inputter_single('Enter (w)avelengths or mark with the (m)ouse? (w/m) ','wm')
        print(mode)
    print('Spectrum runs from {} to {}.'.format(wave[0],wave[-1]))
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid')
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    wshow()
    done = False
    while (not done):
        if (mode == 'w'):
            waveb,waver=waveparse()
        else:
            pickpoints=plt.ginput(2)
            waveb=pickpoints[0][0]
            waver=pickpoints[1][0]
        if (waveb > waver):
            waveb,waver=waver,waveb
        indexblue=womget_element(wave,waveb)
        indexred=womget_element(wave,waver)
        print('Range selected: {} to {}'.format(wave[indexblue],wave[indexred]))
        plt.plot(wave[indexblue:indexred+1],flux[indexblue:indexred+1], \
                 drawstyle='steps-mid')
        print('Is this range correct?')
        answer=yesno('y')
        if (answer == 'n'):
            plt.cla()
            plt.plot(wave,flux,drawstyle='steps-mid')
            plt.ylabel('Flux')
            plt.xlabel('Wavelength')
            wshow()
        else:
            done=True
    return wave[indexblue:indexred+1],flux[indexblue:indexred+1],mode

def pacheck(head):
    observat=head['OBSERVAT'].strip().lower()
    airmass=float(head['AIRMASS'])
    opt_pa=float(head['OPT_PA'])
    date=head['DATE-OBS'].strip()
    # modern date-obs is YYYY-MM-DD, old version used DD/MM/YY
    if (date[4] == '-'):
        year=int(date[0:4])
        month=int(date[5:7])
    else:
        year=int(date[6:8])
        month=int(date[3:5])
        #last of the Y2K bugs, good until 2050
        if (year > 50):
            year=year+1900
        else:
            year=year+2000
    # all this date stuff is to check for the point when FLWO changed
    # header keys for position angle
    if (year >= 2003) or ((year == 2002) and (month >= 12)):
        flwokey='ROTANGLE'
    else:
        flwokey='POSANGLE'
    posangkey={'keck': ('ROTPOSN',90.),
               'lick': ('TUB', 0.0),
               'palomar': ('CASSPA', 0.0),
               'flwo': (flwokey, 0.0),
               'mmto': ('POSANG', 0.0),
	       'kpno' : ('ROTANGLE', 0.0),
               'sso': ('POSANG', 0.0),
               'vlt': ('POSANG', 0.0),
               'lco': ('ROTANGLE', 60.3),
               'gemini-north': ('PA', 0.0),
               'gemini-south': ('PA', 0.0),
               'gemini-n': ('PA', 0.0),
               'gemini-s': ('PA', 0.0),
               'ctio':  ('PA', 0.0),
               'lapalma': ('FIELD',-90.),
               'soar': ('POSANGLE', 0.0)
               }
    if (observat in posangkey):
        pa=float(head[posangkey[observat][0]])+posangkey[observat][1]
    elif (observat == 'mcdonald'):
        pa=float(head['PARANGLE'])-float(head['RHO_OFFS'])
    else:
        pa=1000.0
    diffpa=abs(optpa-pa)
    if (pa >= 999):
        diffpa=0
    diffpa=diffpa % 180
    if (airmass > 1.1) and (((diffpa > 10) and (diffpa < 170))):
        print('************WARNING***************')
        print('Observed position angle: {}'.format(pa))
        print('Optimal parallactic angle: {}'.format(optpa))
        print('Airmass: {}'.format(airmass))
        print(' ')
        print('Relative flux may be compromised')
        print('Hit any key to indemnify me against any and all')
        print('problems that may arise from this')
        any=getch()
    head.set('OBS_PA',pa,'observed position angle')
    return head

def getmswave(head,aperture):
    # have to add 1 as IRAF is not 0 indexed
    specstr='spec'+str(aperture+1)
    watflag=False
    for i,hcard in enumerate(head.cards):
        if ('WAT2' in hcard[0]):
            if (specstr in hcard[1]):
                watstring=hcard[1]
                specloc=watstring.find(specstr)
                watstringsub=watstring[specloc:]
                quoteloc=watstringsub.find('"')
                specsub=watstringsub[quoteloc+1:]
                specsubsplit=specsub.split()
                crval=float(specsubsplit[3])
                cdelt=float(specsubsplit[4])
                npix=float(specsubsplit[5])
                watflag=True
    if (not watflag):
        crval=float(head['CRVAL1'])
        cdelt=float(head['CD1_1'])
        npix=float(head['NAXIS1'])
    wave=np.arange(npix)*cdelt + crval
    return wave

def finalscaler(flux):
    sortflux=np.sort(flux)
    der=np.gradient(sortflux)
    locs=np.where(np.abs(der) < 1.*np.abs(np.mean(der)))
    sortfluxclip=sortflux[locs]
    ymin=min(sortfluxclip)
    ymax=max(sortfluxclip)
    return ymin, ymax

def scipyrebin(wave,flux,nwave):
    import scipy.interpolate as scint
    inter=scint.interp1d(wave,flux,kind='cubic',bounds_error=False,fill_value='extrapolate')
    nflux=inter(nwave)
    return nflux


def get_element(wavearr,wavelength):
    """get index of given wavelength"""
    index=(np.abs(wavearr-wavelength)).argmin()
    return index

def onclick(event):
#    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
#        event.button, event.x, event.y, event.xdata, event.ydata)
    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    if (event.button == 1):
        plt.plot(event.xdata,event.ydata,'ro')
        tmpsplptsx.append(event.xdata)
        tmpsplptsy.append(event.ydata)
        nsplinepoints=nsplinepoints+1
    if (event.button == 2) and (nsplinepoints != 0):
        delindex=np.sqrt((tmpsplptsx - event.xdata)**2 + (tmpsplptsy - event.ydata)**2).argmin()
        plt.plot(tmpsplptsx[delindex],tmpsplptsy[delindex],'gx',markersize=12)
        del tmpsplptsx[delindex]
        del tmpsplptsy[delindex]
        nsplinepoints=nsplinepoints-1
    if (event.button == 3) and (nsplinepoints <= 3):
        tpl=plt.title('Spline requires FOUR points!',size=24)
        plt.pause(1.0)
        tpl.set_visible(False)
    if (event.button == 3) and (nsplinepoints > 3):
        pflag='done'
        return

def wave_telluric(wave, mode):
    loc1=np.logical_and(wave >= 3216.,wave <= 3420.)
    loc2=np.logical_and(wave >= 5600.,wave <= 6050.)
    loc3=np.logical_and(wave >= 6250.,wave <= 6360.)
    loc4=np.logical_and(wave >= 6450.,wave <= 6530.)
    loc5=np.logical_and(wave >= 6840.,wave <= 7410.)
    loc6=np.logical_and(wave >= 7560.,wave <= 8410.)
    loc7=np.logical_and(wave >= 8925.,wave <= 9900.)
    if (mode == 'high'):
        loc=loc1+loc2+loc3+loc4+loc5+loc6+loc7
    else:
        loc=loc1+loc3+loc5+loc6+loc7
    return loc

def obs_extinction(observat):
    observatory_heights={'keck':4160,
                         'gemini-north': 4213.4,
                         'gemini-south': 2737.0,
                         'gemini-n': 4213.4,
                         'gemini-s': 2737.0,
                         'soar': 2737.0,
                         'kpno': 2064.0,
                         'lick': 1285.0,
                         'palomar': 1706.0 ,
                         'mcdonald': 2075.0,
                         'flwo': 2320.0,
                         'mmto': 2600.0,
                         'sso': 1149.0,
                         'vlt': 2635.0,
                         'lco': 2282.0,
                         'lco-imacs': 2282.0,
                         'ctio': 2399.0
                         }
    if (observat in observatory_heights):
        height=observatory_heights[observat]
    else:
        print('OBSERVATORY UNKNOWN!!!!!!!!')
        height=0.0
    # 8300 meters is scale height of troposphere
    sitefactor=np.exp(-1.0*height/8300.0)
    return sitefactor

def fitspl(wave,flux,airlimit,fig):
    """fit spline to spectrum"""
    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    
    # starting points for spline (ENHANCE!)
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
        index=get_element(wave,useband[i])
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
    nsplinepoints=len(useband)
    tmpsplptsx=wave[useband].copy().tolist()
    tmpsplptsy=[]
    for i,_ in enumerate(useband):
        tmpsplptsy.append(np.mean(flux[useband[i]-2:useband[i]+3]))
    spline=splrep(tmpsplptsx,tmpsplptsy,k=3)
    splineresult=splev(wave,spline)
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid',color='k')
    plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
    if (len(tmpsplptsx) > 0):
        plt.plot(tmpsplptsx,tmpsplptsy,'ro')
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
    while (not done):
        plotdone=False
        while (not plotdone):
            plt.cla()
            plt.plot(wave,flux,drawstyle='steps-mid',color='k')
            plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
            if (len(tmpsplptsx) > 0):
                plt.plot(tmpsplptsx,tmpsplptsy,'ro')
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
                if (len(tmpsplptsx) > 0):
                    plt.plot(tmpsplptsx,tmpsplptsy,'ro')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                print('Click corners of box to change plot scale')
                newlims=plt.ginput(2)
                xmin=newlims[0][0]
                ymin=newlims[0][1]
                xmax=newlims[1][0]
                ymax=newlims[1][1]
                plt.cla()
                plt.plot(wave,flux,drawstyle='steps-mid',color='k')
                plt.plot(wavetell,fluxtell,drawstyle='steps-mid',color='violet')
                if (len(tmpsplptsx) > 0):
                    plt.plot(tmpsplptsx,tmpsplptsy,'ro')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                plotdone=True
            else:
                plotdone = True
            
                
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        print('\nClick on continuum points for spline fit.')
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
        splineresult=splev(wave,spline)
        plt.plot(wave,splineresult,drawstyle='steps-mid',color='g')
        plt.pause(0.01)
        print('Is this fit OK? ')
        answer=yesno('y')
        if (answer == 'y'):
            done=True

    return splineresult


def bblo(wave,bstar,airlimit,fig):
    """blotch out bad data"""
    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    print('Select regions to blotch')
    done=False
    while (not done):
        plt.cla()
        plt.plot(wave,flux, drawstyle='steps-mid')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        wavesub,fluxsub,mode=waverange(wave,flux,'none')
        wavebind=get_element(wave,wavesub[0])
        waverind=get_element(wave,wavesub[-1])
        plt.cla()
        plt.plot(wave[wavebind:waverind+1],flux[wavebind:waverind+1], \
                 drawstyle='steps-mid')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
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
                        endpoints=plt.ginput(2)
                        waveselb=endpoints[0][0]
                        waveselr=endpoints[1][0]
                    if (waveselb > waveselr):
                        waveselb,waveselr=waveselr,waveselb
                    waveselbind=get_element(wave,waveselb)
                    waveselrind=get_element(wave,waveselr)
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
                print('Is this acceptable')
                answer=yesno('y')
                if (answer == 'y'):
                    flux=fluxcor.copy()
                    blotchgood=True
                    logging.info('File {} blotched from {} to {}'.format('bstar', wave[waveselbind], wave[waveselrind]))
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
                splblueindex=get_element(wave,splptsx[0])
                splredindex=get_element(wave,splptsx[-1])
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
                    logging.info('File {} blotched  with spline from {} to {}'.format('bstar', wave[splblueindex], wave[splredindex]))
        else:
            done=True        
        print('Do another region?')
        another=yesno('n')
        if (another == 'n'):
            done=True
    
    return bstar



def mkfluxstar(fluxfile,gratcode,intera):
    root = tk.Tk()
    plt.ion()

    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.withdraw()
    screenpos='+{}+{}'.format(int(screen_width*0.2),int(screen_height*0.05))
    fig=plt.figure()

    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('Flux Star')
    fig.set_size_inches(9,6)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    # this should bring the window to the top, but it doesn't
    wm=plt.get_current_fig_manager()
    #fig.canvas.manager.window.attributes('-topmost', 1)
    blah=wm.window.attributes('-topmost', 1)

    #plt.pause(0.2)
    #fig.canvas.manager.window.attributes('-topmost', 0)
    blah=wm.window.attributes('-topmost', 0)

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
    if (intera == 'y'):        
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
    else:
        star=rawdata[0,ap_choice,:]
    # fix IRAF weirdness where data can go negative in dispcor interpolation
    star[np.where(star < 0)] = 0.01
    wave=wavearr[:,ap_choice]
    extinction=scipyrebin(extwave,extvals,wave)
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
    splineresult=fitspl(wave,np.log10(fstar),(airmass>airlimit),fig)
    splineresult=10**(splineresult)
    plt.cla()
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
    outfile='fluxstar'+gratcode
    print('Writing data to {}'.format(outfile))
    outhdu=fits.PrimaryHDU(splineresult)
    hdul=fits.HDUList([outhdu])
    hdul[0].header=head.copy()
    hdul.writeto(outfile,overwrite=True)
    print('mkfluxstar')
    print(fluxfile, gratcode, intera)
    return

def mkbstar(bfile,gratcode,intera):
    root = tk.Tk()
    plt.ion()

    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.withdraw()
    screenpos='+{}+{}'.format(int(screen_width*0.2),int(screen_height*0.05))
    fig=plt.figure()

    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('B Star')
    fig.set_size_inches(9,6)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    # this should bring the window to the top, but it doesn't
    wm=plt.get_current_fig_manager()
    #fig.canvas.manager.window.attributes('-topmost', 1)
    blah=wm.window.attributes('-topmost', 1)

    #plt.pause(0.2)
    #fig.canvas.manager.window.attributes('-topmost', 0)
    blah=wm.window.attributes('-topmost', 0)

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
    if (intera == 'y'):        
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
            bstar=rawdata[0,ap_choice,:]
        else:
            bstar=rawdata[1,ap_choice,:]
    else:
        bstar=rawdata[0,ap_choice,:]
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
    splineresult=fitspl(wave,np.log10(bstar),(airmass>airlimit),fig)
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
    plt.cla()
    plt.plot(wave,bstar,drawstyle='steps-mid')
    plt.pause(0.01)
    print('\nDo you want to blotch the B-star?')
    blotch=yesno('n')
    if (blotch == 'y'):
        bstar=bblo(wave,bstar,(airmass > airlimit))
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
    outfile='bstar'+gratcode
    print('Writing data to {}'.format(outfile))
    outhdu=fits.PrimaryHDU(splineresult)
    hdul=fits.HDUList([outhdu])
    hdul[0].header=head.copy()
    hdul.writeto(outfile,overwrite=True)    
    print('mkbstar')
    print(bfile, gratcode, intera)
    return

def calibrate(objectlist,gratcode,intera,secondord,gratcode2):
    #extinction terms from Allen, 3rd edition
    extwave= [2400.,2600.,2800.,3000.,3200.,3400.,3600.,3800., \
              4000.,4500.,5000.,5500.,6000.,6500.,7000.,8000., \
              9000.,10000.,12000.,14000.]
    extvals=[68.0,89.0,36.0,4.5,1.30,0.84,0.68,0.55,0.46,0.31, \
             0.23,0.195,0.170,0.126,0.092,0.062,0.048,0.039, \
             0.028,0.021]
    fluxfits=fits.open('fluxstar'+gratcode+'.fits')
    fluxstar=fluxfits[0].data
    fluxhead=fluxfits[0].header
    fluxwavezero=float(fluxhead['CRVAL1'])
    fluxwavedelt=float(fluxhead['CDELT1'])
    fluxwave=np.arange(len(fluxstar))*fluxwavedelt+fluxwavezero
    fluxairmass=float(fluxhead['AIRMASS'])
    fluxname=fluxhead['OBJECT']
    try:
        fluxnum=int(fluxhead['OBSNUM'])
    except KeyError:
        fluxnum=0
    if (secondord):
        fluxfits2=fits.open('fluxstar'+gratcode2+'.fits')
        fluxstar2=fluxfits[0].data
        fluxhead2=fluxfits[0].header
        fluxwavezero2=float(fluxhead2['CRVAL1'])
        fluxwavedelt2=float(fluxhead2['CDELT1'])
        fluxwave2=np.arange(len(fluxstar2))*fluxwavedelt2+fluxwavezero2
        fluxairmass2=float(fluxhead2['AIRMASS'])
        fluxname2=fluxhead2['OBJECT']
        try:
            fluxnum2=int(fluxhead2['OBSNUM'])
        except KeyError:
            fluxnum2=0
    observat=fluxhead['OBSERVAT'].strip().lower()
    sitefactor=obs_extinction(observat)
    infile=open(objectlist,'r')
    for line in infile:
        line=line.strip()
        if ('.ms.fits' not in line):
            inputfile=line+'.ms.fits'
        multifits=fits.open(line)
        multispec=multifits[0].data
        mshead=multifits[0].header
        objectname=mshead['OBJECT']
        print('The object is: {}'.format(objectname))
        airmass=float(mshead['AIRMASS'])
        exptime=float(mshead['EXPTIME'])
        if (exptime < 1):
            exptime=1.0
        num_apertures=multispec.shape[1]
        num_bands=multispec.shape[0]
        wavearr=np.zeros((multispec.shape[2],multispec.shape[1]))
        if (secondord):
            multispec2=multispec.copy()
            mshead2=mshead.copy()
        for i in range(0,num_apertures):
            print('\nAperture {}:'.format(i+1))
            wave=getmswave(mshead,i)
            extinction=scipyrebin(extwave,extvals,wave)
            extfactor=np.exp(extinction*sitefactor*airmass)
            fluxstartmp=scipyrebin(fluxwave,fluxstar,wave)
            wdelt=wave[1]-wave[0]
            for j in range(0,num_bands):
                multispec[j,i,:]=multispec[j,i,:]*extfactor       #extinction
                multispec[j,i,:]=multispec[j,i,:]/fluxstartmp     #flux
                multispec[j,i,:]=multispec[j,i,:]/exptime         #adjust to time
                multispec[j,i,:]=multispec[j,i,:]*10**(-19.44)    #AB->fnu
                multispec[j,i,:]=multispec[j,i,:]*2.99792458e18/wave/wave #fnu->flm
            if (secondord):
                fluxstartmp2=scipyrebin(fluxwave2,fluxstar2,wave)
                for j in range(0,num_bands):
                    multispec2[j,i,:]=multispec2[j,i,:]*extfactor       #extinction
                    multispec2[j,i,:]=multispec2[j,i,:]/fluxstartmp     #flux
                    multispec2[j,i,:]=multispec2[j,i,:]/exptime         #adjust to time
                    multispec2[j,i,:]=multispec2[j,i,:]*10**(-19.44)    #AB->fnu
                    multispec2[j,i,:]=multispec2[j,i,:]*2.99792458e18/wave/wave #fnu->flm
        msfile='c'+gratcode+msfile
        mshead.set('FLUX_Z',fluxairmass,'airmass of flux standard')
        mshead.set('FLUX_NUM',fluxnum,'obsnum of flux standard')
        mshead.set('FLUX_OBJ',fluxname,'id of flux standard')
        outhdu=fits.PrimaryHDU(multispec)
        hdul=fits.HDUList([outhdu])
        hdul[0].header=mshead.copy()
        hdul.writeto(msfile,overwrite=True)
        hdul.close()
        if (secondord):
            msfile='c'+gratcode2+msfile
            mshead2.set('FLX2_Z',fluxairmass,'airmass of flux second ord. standard')
            mshead2.set('FLX2_NUM',fluxnum,'obsnum of flux second ord. standard')
            mshead2.set('FLX2_OBJ',fluxname,'id of flux second ord. standard')
            outhdu=fits.PrimaryHDU(multispec2)
            hdul=fits.HDUList([outhdu])
            hdul[0].header=mshead2.copy()
            hdul.writeto(msfile2,overwrite=True)
            hdul.close()

    infile.close()
    fluxfits.close()
    if (secondord):
        fluxfits2.close()
    print('calibrate')
    print(objectlist,gratcode,intera,secondord,gratcode2)
    return

def final(objectlist,gratcode,intera,secondord,gratcode2,user):
    root = tk.Tk()
    plt.ion()

    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.withdraw()
    screenpos='+{}+{}'.format(int(screen_width*0.2),int(screen_height*0.05))
    fig=plt.figure()

    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('Final Reductions')
    fig.set_size_inches(9,6)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    # this should bring the window to the top, but it doesn't
    wm=plt.get_current_fig_manager()
    #fig.canvas.manager.window.attributes('-topmost', 1)
    blah=wm.window.attributes('-topmost', 1)

    #plt.pause(0.2)
    #fig.canvas.manager.window.attributes('-topmost', 0)
    blah=wm.window.attributes('-topmost', 0)


    bstarfits=fits.open('bstarstar'+gratcode+'.fits')
    bstarstar=bstarfits[0].data
    bstarhead=bstarfits[0].header
    bstarwavezero=float(bstarhead['CRVAL1'])
    bstarwavedelt=float(bstarhead['CDELT1'])
    bstarwave=np.arange(len(bstarstar))*bstarwavedelt+bstarwavezero
    bstarairmass=float(bstarhead['AIRMASS'])
    bstarname=bstarhead['OBJECT']
    try:
        bstarnum=int(bstarhead['OBSNUM'])
    except KeyError:
        bstarnum=0
    if (secondord):
        bstarfits2=fits.open('bstarstar'+gratcode2+'.fits')
        bstarstar2=bstarfits[0].data
        bstarhead2=bstarfits[0].header
        bstarwavezero2=float(bstarhead2['CRVAL1'])
        bstarwavedelt2=float(bstarhead2['CDELT1'])
        bstarwave2=np.arange(len(bstarstar2))*bstarwavedelt2+bstarwavezero2
        bstarairmass2=float(bstarhead2['AIRMASS'])
        bstarname2=bstarhead2['OBJECT']
        try:
            bstarnum2=int(bstarhead2['OBSNUM'])
        except KeyError:
            bstarnum2=0
    observat=bstarhead['OBSERVAT'].strip().lower()
    reduxdir=os.environ(['PY_REDUX'])
    kecksky=['keck','gemini-north','gemini-n','gemini-south','gemini-s', \
             'soar','ctio','vlt','lco','lco-imacs','lapalma']
    licksky=['lick','kpno','flwo','mmto','sso']
    if (observat in kecksky):
        mskyfile=reduxdir+'/kecksky.fits'
    elif (observat in licksky):
        mskyfile=reduxdir+'/licksky.fits'
    else:
        print('\nCannot find mastersky file and observatory unknown\n')
        mskyfile=getfitsfile('master sky','.fits')
    mskyfits=fits.open(mskyfile)
    mskydata=mskyfits[0].data
    mskyhead=mskyfits[0].header
    mskywavezero=float(mskyhead['CRVAL1'])
    mskywavedelt=float(mskyhead['CDELT1'])
    mskywave=np.arange(len(mskystar))*mskywavedelt+mskywavezero
    if (np.abs(np.mean(mskydata)) < 1e-7):
        mskydata=mskydata*1e15
    print('\nHello {}\n'.format(user))
    lat_dict={'keck': 19.8283,
              'gemini-north': 19.8238,
              'gemini-south': -30.228,
              'gemini-n': 19.8238,
              'gemini-s': -30.228,
              'soar': -30.228,
              'kpno': 31.9633,
              'lick': 37.3414,
              'palomar': 33.35611,
              'mcdonald': 30.6717,
              'flwo': 31.681,
              'mmto': 31.688,
              'sso':  -31.2734,
              'vlt':  -24.6254,
              'lco':  -29.01,
              'lco-imacs':  -29.01,
              'lapalma': 28.75833,
              'ctio': -30.16894}
    infile=open(objectlist,'r')
    secondtime = False
    seconddone = False
    for line in infile:
        line=line.strip()
        if ('.ms.fits' not in line):
            inputfile=line+'.ms.fits'
        multifits=fits.open('c'+gratcode+line)
        multispec=multifits[0].data
        mshead=multifits[0].header
        num_apertures=multispec.shape[1]
        num_bands=multispec.shape[0]
        wavearr=np.zeros((multispec.shape[2],multispec.shape[1]))
        if (secondord):
            multifits2=fits.open('c'+gratcode2+line)
            multispec2=multifits2[0].data
            mshead2=multifits2[0].header
        pacheck(mshead)
        objectname=mshead['OBJECT']
        print('The object is: {}'.format(objectname))
        observat=mshead['OBSERVAT'].strip().lower()
        airmass=float(mshead['AIRMASS'])
        if (airmass < 1):
            airmass=1.0
        stringra=mshead['RA']
        stringdec=mshead['DEC']
        geminilist=['gemini-north','gemini-south','gemini-n','gemini-s']
        if (observat) in geminilist:
            ra=float(stringra)
            dec=float(stringdec)
        else:
            coords==SkyCoord(stringra,stringdec,unit=(u.hourangle,u.deg))
            ra=coords.ra.deg
            dec=coords.dec.deg
        ra=ra/RADEG
        dec=dec/RADEG
        stringha=mshead['HA']
        if (stringha[0] not in '+-'):
            stringha='+'+stringha
        halist=stringha.split(':')
        ha=(halist[0]+(halist[1]+halist[2]/60.)/60.)*(360./24.)/RADEG

        if (observat in lat_dict):
            latitude=lat_dict[observat]
        else:
            print('Observatory Unknown!!!!!!!')
            latitude=0.0
        latitude=latitude/RADEG
        # Get Julian date and Earth's velocity
        epoch=float(mshead['EPOCH'])
        date=mshead['DATE-OBS'].strip()
        # get year,month,day from ISO Y2K format YYYY-MM-DDThh:mm:ss.ssss
        # the time following T is optional, so we get T from UTMIDDLE
        if (date[4] == '-'):
            year=int(date.split('-')[0])
            month=int(date.split('-')[1])
            day=int(date.split('-')[2].split('T')[0])
            printdate='{}{:02d}{:02d}'.format(year,month,day)
        else:
            # Old date format DD/MM/YY
            year=int(date[6:8])
            month=int(date[3:5])
            day=int(date[0:2])
            # try to catch old format written after 2000
            # good until 2050, by which time no one should
            # be doing this
            if (year > 50):
                year=year+1900
            else:
                year=year+2000
            printdate='{}{:02d}{:02d}'.format(year,month,day)
        try:
            ut=mshead['UTMIDDLE'].strip()
        except KeyError:
            ut=mshead['UT'].strip()
        hour=int(ut.split(':')[0])
        minute=int(ut.split(':')[1])
        second=float(ut.split(':')[2])
        julian=jdcnv(year,month,day,(hour+(minute/60.)+(second/3600.)))
        vh, vb=baryvel(julian,epoch)
        # Earth's velocity toward object
        v=vb[0]*np.cos(dec)*np.cos(ra) + vb[1]*np.cos(dec)*np.sin(ra) \
           + vb[2]*np.sin(dec)
        # Correct for earth's rotation.
        # Note that this isn't strictly correct because ha is set to
        # the beginning of the observation, while it really should be
        # the middle.  But this is a small difference and doesn't affect
        # the results in any way that we would notice...
        v =  v - ( 0.4651 * cos(latitude) * sin(ha) * cos(dec) )
        print('\nThe velocity of the Earth toward the target is {} km/s'.format(v))
        print('\nThere are {} apertures in the spectrum.\n')
        # clean up header
        delkeylist=['WAT0_001', 'WAT1_001', 'WAT2_001', 'WAT0_002', \
                'WAT1_002', 'WAT2_002', 'WAT3_001', 'WAT2_003', \
                'CTYPE1', 'CTYPE2', 'CTYPE3', 'CD1_1', 'CD2_2', \
                'CD3_3', 'LTM1_1', 'LTM2_2', 'LTM3_3', 'WCSDIM']
        for k in delkeylist:
            try:
                mshead.remove(k)
            except KeyError:
                pass
        for i in range(0,num_apertures):
            print('\nAperture {}:'.format(i+1))
            wave=getmswave(mshead,i)
            wdelt=wave[1]-wave[0]
            mean=np.mean(multispec[0,i,:])
            ymin,ymax=finalscaler(multispec[0,i,:])
            plt.clf()
            gs=gridspec.GridSpec(2,1,height_ratios=[4,1])
            ax0=plt.subplot(gs[0])
            ax1=plt.subplot(gs[1])
            fig.subplots_adjust(hspace=0)
            ax0.plot(wave,multispec[1,i,:],drawstyle='steps-mid',color='r')
            ax0.plot(wave,multispec[0,i,:],drawstyle='steps-mid',color='k')
            plt.pause(0.01)
            ax0.set_ylim((ymin,ymax))
            ax1.semilogy(wave,np.abs(multispec[1,i,:]-multispec[0,i,:])/mean, \
                         drawstyle='steps-mid',color='k')
            ax0.set_ylabel('Flux')
            ax1.set_ylabel('Log of fractional residuals')
            ax1.set_xlabel('Wavelength')
            print('\nPlotting optimal as black, normal as red\n')
            extract=inputter_single('Do you want to use the (n)ormal or the (o)ptimal extraction? ','no')
            if (extract == 'o'):
                object=multispec[0,i,:]
                extractcode='optimal'
                if (secondord):
                    object2=multispec2[0,i,:]
            else:
                object=multispec[1,i,:]
                extractcode='normal'
                if (secondord):
                    object2=multispec2[1,i,:]
            skyband=2
            if (num_bands == 2):
                skyband=1
            if (num_bands > 3):
                sigma=multispec[3,i,:]
            else:
                sigma=np.sqrt(multispec[skyband,i,:])
            if (np.abs(np.mean(object)) < 1e-7):
                print('\nScaling data up by 10^15\n')
                object=object*1.e15
                sigma=sigma*1.e15
            if (secondord):
                if (num_bands > 3):
                    sigma2=multispec2[3,i,:]
                else:
                    sigma2=np.sqrt(multispec2[skyband,i,:])
                    if (np.abs(np.mean(object)) < 1e-7):
                        object2=object2*1.e15
                        sigma2=sigma2*1.e15
            nandslist=['gemini-north','gemini-south','lco-imacs']
            if (observat in nandslist):
                skyfile=msfile.replace('tot','sky')
                skyfits=fits.open(skyfile)
                sky=skyfits[0].data[2,0,:]
            else:
                sky=multispec[skyloc,i,:]
            mx,mn=envelope(sky,envelope_size)
            skycont=congrid(mn,(len(sky),),minusone=True)
            sky=sky-cont
            if (np.abs(np.mean(sky)) < 1.e-7):
                sky=sky*1.e15
            msky=scipyrebin(mskywave,mskydata,wave)
            shift=xcor(msky,sky,xfactor,maxlag)
            angshift=shift*wdelt
            print('The x-cor shift in Angstroms is {}'.format(angshift))
            wave=wave-angshift
            skyshiftdone=False
            npixsky2=len(sky)//2
            while (not skyshiftdone):
                plt.clf()
                axarr=fig.subplots(2,sharex=True)
                fig.subplots_adjust(hspace=0)
                waveplus=wave+angshift
                axarr[0].plot(waveplus[0:npixsky2],msky[0:npixsky2], \
                              drawstyle='steps-mid',color='k')
                axarr[0].plot(wave[0:npixsky2],sky[0:npixsky2], \
                              drawstyle='steps-mid',color='r')
                axarr[0].plot(waveplus[npixsky2-1:],msky[npixsky2-1:], \
                              drawstyle='steps-mid',color='k')
                axarr[0].plot(wave[npixsky2-1:],sky[npixsky2-1:], \
                              drawstyle='steps-mid',color='r')
                plt.pause(0.01)
                print('\nBlack spectrum = master sky')
                print('Red spectrum   = object sky shifted to match master sky')
                print('Is this ok?')
                answer=yesno('y')
                if (answer == 'n'):
                    wave=wave + angshift
                    angshift=inputter('Enter desired shift in Angstroms: ','float',False)
                    wave=wave-angshift
                else:
                    skyshiftdone = True
            # B star removal
            bstarpass=bstar
            bobj, bsig=telluric_remove(bstarwave,bstarpass, bairmass, wave, \
                                       object, airmass, sigma)
            if (secondord):
                bobj2, bsig2=telluric_remove(bstarwave2,bstarpass2, bairmass2, \
                                             wave, object2, airmass, sigma2)

            print('\nRemoving redshift due to motion of Earth...')
            z=-1*v/2.99792458e5
            wave=wave/(1+z)
            # rebin
            plt.clf()
            axarr=fig.subplots(1,2,sharey=True)
            ymin,ymax=finalscaler(bobj[0:100])
            axarr[0].plot(wave[0:100],bobj[0:100],drawstyle='steps-mid')
            axarr[0].set_xlabel('Wavelength')
            axarr[0].set_ylabel('Flux')
            axarr[0].set_ylim((ymin,ymax))
            if (secondtime):
                axarr[0].plot([wavesave0,wavesave0],[ymin,ymax],color='r')
            ymin,ymax=finalscaler(bobj[-100:])
            axarr[1].plot(wave[-100:],bobj[-100:],drawstyle='steps-mid')
            axarr[1].set_xlabel('Wavelength')
            axarr[1].set_ylabel('Flux')
            axarr[1].set_ylim((ymin,ymax))
            if (secondtime):
                axarr[1].plot([wavesaven,wavesaven],[ymin,ymax],color='r')

                newdelt=0.0
            print('\nCurrent A/pix is {}'.format(wave[1]-wave[0]))
            if (secondtime):
                print('\nPrevious resolution choice: {}'.format(deltsave))
            else:
                deltsave = 0.0
            while (newdelt <= 0) or (newdelt > wave[-1]):
                print('Rebin to how many Angstroms per pixel? ')
                newdelt=inputter('         <CR> selects previous choice: ','float',True,deltsave)
                if (newdelt <= 0) or (newdelt > wave[-1]):
                    print('Need positive resoution and smaller than the')
                    print('entire spectrum.  Try again')
            print('\nCurrent range: {} {}'.format(wave[0],wave[-1]))
            if (secondtime):
                print('\nPrevious selection was {} {} (marked in red on plot)'.format(wavesave0,wavesaven))
            else:
                wavesave0 = 0.0
                wavesaven = 0.0
                print('Enter the new wavelength range desired: ')
            print('          <CR> selects previous choice:')
            waveb,waver=waveparse(wave,wavesave0,wavesaven)
            newbin=(waver-waveb)/newdelt +1.0
            frac,whole=np.modf(newbin)
            if (frac > 0.000001):
                print('NON-INTEGER number of bins')
                testwave=newdelt*whole+waveb
                print('Closest match is: {} to {}'.format(waveb,testwave))
                print('Would you like this wavelength range?')
                answer=yesno('y')
                if (answer == 'y'):
                    waver=testwave
            deltsave=newdelt
            wavesave0=waveb
            wavesaven=waver
            waverange=str(waveb) + ' ' + str(waver)
            secondtime = True
            nwave=np.arange(0,whole)*newdelt + waveb
            vartmp=bsig*bsig
            olddeltvec=wave-np.roll(wave,1)
            olddelt=np.mean(olddeltvec)
            finalobj=womashrebin(wave,bobj,nwave)
            finalvar=womashrebin(wave,vartmp,nwave)
            finalsig=np.sqrt(finalvar)*np.sqrt(olddelt/newdelt)
            if (secondord):
                vartmp2=bsig2*bsig2
                finalobj2=womashrebin(wave,bobj2,nwave)
                finalvar2=womashrebin(wave,vartmp2,nwave)
                finalsig2=np.sqrt(finalvar2)*np.sqrt(olddelt/newdelt)
            if (secondord and not seconddone):
                finalobj, finalsig, wavebeg, waveend, brscale=secondcat(nwave,finalobj,finalobj2, finalsig,finalsig2,secondtime, wavebeg, waveend, brscale)
            seconddone = True
            ymin,ymax=finalscaler(finalobj)
            plt.clf()
            plt.plot(nwave,finalobj,drawstyle='steps-mid')
            plt.xlabel('Wavelength')
            plt.ylable('Flux')
            plt.title(objectname)
            outputdone = False
            while (not outputdone):
                print('\nThe file is: {}'.format(msfile))
                print('The object is: {}'.format(objectname))
                print('The DATE-OBS is: {}'.format(date))
                print('The aperture is: {}'.format(i+1))
                print('The previous name was: {}'.format(objname))
                print('\nEnter the object name for the final fits file: ')
                objname=inputter('(UT date and .fits will be added): ','string',False)
                fname=objname+'-'+printdate+'.fits'
                sname=objname+'-'+printdate+'-sigma.fits'
                if (os.path.isfile(fname)):
                    print('{} already exists!!!!')
                    print('Do you wish to overwrite it? ')
                    answer=yesno('y')
                    if (yesno == 'y'):
                        outputdone = True
                else:
                    outputdone = True
            # add to header
            mshead.set('CRPIX1', 1)
            mshead.set('CRVAL1',  nwave[0])
            mshead.set('CDELT1', nwave[1] - nwave[0])
            mshead.set('CTYPE1', 'LINEAR')
            mshead.set('W_RANGE', waverange)
            mshead.set('BSTAR_Z', bairmass)
            mshead.set('BSTARNUM', bstarnum)
            mshead.set('BSTAROBJ', bstarname)
            mshead.set('BARYVEL', v)
            mshead.set('SKYSHIFT', angshift)
            mshead.set('ATMSHIFT', bangshift)
            mshead.set('EXTRACT', extractstring)
            mshead.set('REDUCER', user)
            mshead.set('RED_DATE', rtime, 'EPOCH OF REDUCTION')
            mshead.set('OBJECT', objectname)
            if (secondord):
                mshead.set('SECOND', 'yes',  'Second order correction attempted')
                mshead.set('COMBRANGE',combrange)
                mshead.set('BSTAR_Z2',bairmass2)
                mshead.set('BSTARNU2', bstarnum2)
                mshead.set('BSTAROB2', bstarname2)
                fluxairmass2=mshead2['FLX2_Z']
                fluxnum2=mshead2['FLX2_NUM']
                fluxname2=mshead2['FLX2_OBJ']
                mshead.set('FLX2_Z',fluxairmass2)
                mshead.set('FLX2_NUM', fluxnum2)
                mshead.set('FLX2_OBJ', fluxname2)
            outdata=np.zeros((len(finalobj),2))
            outdata[:,0]=finalobj.copy()
            outdata[:,1]=finalsig.copy()
            outhdu=fits.PrimaryHDU(outdata)
            hdul=fits.HDUList([outhdu])
            mshead.set('NAXIS2',2)
            hdul[0].header=mshead.copy()
            hdul.writeto(fname,overwrite=True)
            hdul.close()
                
    print('final')        
    print(objectlist,gratcode,intera,secondord,gratcode2)
    return

def waveparse(wave,oldwaveb,oldwaver):
    import re
    done=False
    while (not done):
        waveparse=input('Enter wavelength range: ')
        if (waveparse != ''):
            waveparse=waveparse.strip()
            wavesplit=re.split('\W+',waveparse)
            if (len(wavesplit) > 1):
                try:
                    waveb=float(wavesplit[0])
                    waver=float(wavesplit[1])
                except ValueError:
                    print('Please enter numbers, #### ####,\n')
                    print('####,####,or ####-####\n')
                    print('You entered {}\n'.format(waveparse))
                else:
                    done=True
            else:
                print('Enter more than one number\n')
            if (waveb <= 0) or (waver <= 0)  or \
               (waveb < wave[0]) or (waver < wave[0]) or \
               (waveb > wave[-1]) or (waver > wave[-1]) or \
               (waveb > waver):
                print('Wavelengths out of bounds.  Try again')
                done = False
        else:
            waveb=oldwaveb
            waver=oldwaver
    return waveb,waver

def testo():
    x=np.arange(100)
    y=np.sin(x)
    plt.plot(x,y)
    answer=yesno('y')
    return



def abcalc(wave, objectname):

#  creates ab curve to flux standard stars
#  TM 11/21/98


#  EVALUATES AB(NU) MAGNITUDES FOR THE FLUX STAR SPECIFIED BY IDNO.
#
#  IDNO = 1 FOR HD 19445
#         2 FOR HD 84937
#         3 FOR BD +26 2606
#         4 FOR BD +17 4708
#         5 FOR HD 140283
#         6 FOR GD-248
#         7 FOR G158-100
#         8 FOR LTT 1788
#         9 FOR LTT 377
#        10 FOR HZ 4
#        11 FOR HZ 7
#        12 FOR ROSS 627
#        13 FOR LDS 235B
#        14 FOR G99-37
#        15 FOR LB 227
#        16 FOR L745-46A
#        17 FOR FEIGE 24
#        18 FOR G60-54
#        19 FOR G24-9
#        20 FOR BD+28 4211
#        21 FOR G191B2B
#        22 FOR G157-34
#        23 FOR G138-31
#        24 FOR HZ 44
#        25 FOR LTT 9491
#        26 FOR FEIGE 100
#        27 FOR FEIGE 34
#        28 FOR LTT 1020
#        29 FOR LTT 9239
#	 30 FOR HILTNER 600 (APPROXIMATE)
#	 31 FOR BD+25 3941 (APPROXIMATE)
#        32 FOR BD+33 2642 (APPROXIMATE)
#        33 FOR FEIGE 56 (APPROXIMATE)
#        34 FOR G193-74
#        35 FOR eg145
#        36 for FEIGE 25
#        37 for PG0823+546
#        38 for HD 217086
#        39 for HZ 14
#        40 for FEIGE 66
#
#        NUMBERS 25--29 WERE ADDED ON 2 AUGUST 1987# NUMBER 8 WAS
#        MODIFIED SLIGHTLY (UP TO 0.11 MAG) NEAR THE H-K BREAK.
#
#	 NUMBERS 30 AND 31 WERE ADDED ON 24 JANUARY 1988# NUMBERS
#	 ARE ONLY APPROXIMATE NEAR THE BALMER LIMIT, LONGWARD OF
#	 8100 A, AND SHORTWARD OF 3300 A. DITTO FOR NUMBER 32,
#        WHI#H WAS ADDED ON 2 MAY 1988. DITTO FOR NUMBER 33, WHICH
#        WAS ADDED ON 25 AUGUST 1991.
#	THERE ARE 41 WAVELENGTHS.


#  INTERPOLATION IN THIS TABLE SHOULD APPROXIMATE THE MONOCHROMATIC
#  CONTINUUM FLUX DISTRIBUTION OF THE MCSP STANDARD STARS.
#  SYSTEMATIC ERRORS AS LARGE AS 5% MAY OCCUR NEAR THE CONVERGENCE
#  OF THE BALMER SERIES, WHERE THE CONTINUUM IS NOT WELL DEFINED.
#  NO ATTEMPT IS MADE TO FOLLOW THE WINGS OF THE BALMER LINES,
#  A PROCEDURE WHICH SEEMS DANGEROUS AT BEST.
#
#	NOTE THAT ALEX FILIPPENKO HAS MODIFIED ALL AB MAGNITUDES SO
#	THAT ABSORPTION LINES ARE INTERPOLATED OVER SMOOTHLY. THIS
#	MEANS THAT MCSP STANDARDS AND OTHERS CANNOT BE USED IN THE
#	MANNER USED BEFORE (THAT IS, SUMMING ALL COUNTS IN 40-ANG
#	BANDPASSES AND COMPARING WITH TABULATED VALUES).
#
    idstar = 0
    while (idstar < 1) or (idstar > 56):
        print('The object is {}'.format(objectname))
        print(' ')
        print('  WHICH STANDARD STAR IS THIS ?')
        print('  ORIGINAL MCSP FLUX STANDARDS (AB79 SCALE): ')
        print('  (1) HD 19445       (2) HD 84937        (3) BD+262606')
        print('  (4) BD+17 4708       (5) HD 140283 (SOMEWHAT OBSOLETE)')
        print('  OTHERS (AB69, AB79, STONE, OR BALDWIN/STONE SCALE): ')
        print('  (6) GD-248 (AB79)            (7) G158-100 (AB79)')
        print('  (8) LTT 1788 (B/S)           (9) LTT 377 (B/S)')
        print('  (10) HZ 4 (~AB79)            (11) HZ 7 (~AB79)')
        print('  (12) ROSS 627 (~AB79)        (13) LDS 235B (~AB79)')
        print('  (14) G99-37 (~AB79)          (15) LB 227 (~AB79)')
        print('  (16) L745-46A (~AB79)        (17) FEIGE 24 (~AB79)')
        print('  (18) G60-54 (~AB79)          (19) G24-9 (AB79)')
        print('  (20) BD+28 4211 (STONE)      (21) G191B2B (~AB79)')
        print('  (22) G157-34 (AB79)          (23) G138-31 (AB79)')
        print('  (24) HZ 44 (~AB79)           (25) LTT 9491 (B/S)')
        print('  (26) FEIGE 110 (STONE)       (27) FEIGE 34 (STONE)')
        print('  (28) LTT 1020 (B/S)          (29) LTT 9239 (B/S)' )
        print('  (30) HILTNER 600 (STONE)     (31) BD+25 3941 (STONE)')
        print('  (32) BD+33 2642 (STONE)      (33) FEIGE 56 (STONE)'  )
        print('  (34) G193-74 (AB79)          (35) EG145 (OKE,74)  ')
        print('  (36) FEIGE 25 (STONE)        (37) PG0823+546 ')
        print('  (38) HD 217086               (39) HZ 14 ')
        print('  (40) FEIGE 66                (41) Feige 67 (Groot special) ')
        print('  (42) LTT 377                 (43) LTT 2415 ')
        print('  (44) LTT 4364                (45) Feige 15')
        print('  (46) Hiltner 102 (lousy)     (47) LTT 3864 ')
        print('  (48) LTT 3218                (49) CYG OB2  ')
        print('  (50) VMa 2                   (51) GD 71 ')
        print('  (52) HZ 43                   (53) LTT 7379')
        print('  (54) LTT 7987                (55) GD 153 ')
        print('  (56) CD32D9927 ')
        idstar=inputter('Star number? ','int',False)

    #  All numbers from lolita's abcalc.f, we trust that they are ok

    waveab = [3080.0, 3160.0, 3240.0, 3320.0, 3400.0, 3480.0, 3560.0, \
              3640.0, 3680.0, 3720.0, 3760.0, 3780.0, 3820.0, 3860.0, \
              4020.0, 4200.0, 4400.0, 4560.0, 4760.0, 5000.0, 5120.0, \
              5240.0, 5400.0, 5560.0, 5760.0, 6020.0, 6420.0, 6780.0, \
              7100.0, 7460.0, 7780.0, 8100.0, 8380.0, 8780.0, 9300.0, \
              9700.0, 9940.0,10260.0,10820.0,11140.0,12000.0]        
    ab = np.zeros((41, 57))
    ab[:, 1] = [9.476, 9.391, 9.310, 9.232, 9.168, 9.103, 9.039, \
                8.979, 8.945, 8.919, 8.849, 8.814, 8.762, 8.695, \
                8.483, 8.385, 8.293, 8.235, 8.180, 8.118, 8.081, \
                8.048, 8.014, 7.982, 7.948, 7.908, 7.862, 7.837, \
                7.812, 7.790, 7.776, 7.760, 7.759, 7.755, 7.752, \
                7.752, 7.752, 7.752, 7.760, 7.760, 7.763]
    ab[:, 2] = [9.677, 9.610, 9.542, 9.480, 9.422, 9.370, 9.310, \
                9.269, 9.243, 9.205, 9.099, 9.018, 8.926, 8.850, \
                8.650, 8.548, 8.496, 8.448, 8.400, 8.351, 8.332, \
                8.315, 8.290, 8.266, 8.237, 8.202, 8.174, 8.150, \
                8.133, 8.125, 8.123, 8.112, 8.111, 8.115, 8.119, \
                8.123, 8.130, 8.137, 8.150, 8.160, 8.180]
    ab[:, 3] = [11.155, 11.072, 10.998, 10.935, 10.866, 10.806, \
                10.750, 10.698, 10.670, 10.651, 10.540, 10.480, \
                10.400, 10.320, 10.120, 10.021, 9.955, 9.908, 9.843, \
                9.768, 9.741, 9.719, 9.693, 9.660, 9.630, 9.594, \
                9.552, 9.527, 9.502, 9.485, 9.474, 9.463, 9.457, \
                9.457, 9.458, 9.458, 9.461, 9.462, 9.468, 9.473, \
                9.488]
    ab[:, 4] = [10.955, 10.872, 10.798, 10.718, 10.663, 10.601, \
                10.536, 10.487, 10.452, 10.426, 10.307, 10.268, \
                10.190, 10.110, 9.880, 9.765, 9.690, 9.637, 9.573, \
                9.520, 9.491, 9.465, 9.434, 9.405, 9.374, 9.337, \
                9.294, 9.262, 9.233, 9.216, 9.211, 9.200, 9.194, \
                9.194, 9.195, 9.195, 9.198, 9.199, 9.205, 9.210, \
                9.225]
    ab[:, 5] = [9.050, 8.820, 8.630, 8.520, 8.430, 8.360, 8.300, \
                8.230, 8.200, 8.160, 8.050, 7.990, 7.915, 7.870, \
                7.680, 7.560, 7.470, 7.420, 7.340, 7.260, 7.230, \
                7.210, 7.170, 7.120, 7.080, 7.040, 6.975, 6.930, \
                6.890, 6.860, 6.840, 6.830, 6.810, 6.800, 6.810, \
                6.840, 6.860, 6.880, 6.890, 6.900, 6.910]
    ab[:, 6] = [15.295, 15.262, 15.229, 15.200, 15.178, 15.163, \
                15.150, 15.138, 15.133, 15.127, 15.122, 15.119, \
                15.113, 15.108, 15.090, 15.071, 15.058, 15.052, \
                15.050, 15.056, 15.063, 15.070, 15.081, 15.093, \
                15.110, 15.133, 15.175, 15.212, 15.245, 15.283, \
                15.318, 15.353, 15.384, 15.427, 15.484, 15.528, \
                15.554, 15.592, 15.652, 15.688, 15.775]
    ab[:, 7] = [16.725, 16.665, 16.591, 16.502, 16.400, 16.292, \
                16.183, 16.082, 16.034, 15.986, 15.938, 15.914, \
                15.868, 15.824, 15.651, 15.490, 15.346, 15.246, \
                15.139, 15.027, 14.979, 14.932, 14.872, 14.824, \
                14.768, 14.696, 14.609, 14.537, 14.485, 14.430, \
                14.393, 14.361, 14.338, 14.310, 14.281, 14.261, \
                14.248, 14.233, 14.209, 14.199, 14.178]
    ab[:, 8] = [14.785, 14.637, 14.492, 14.348, 14.214, 14.113, \
                14.025, 13.946, 13.911, 13.880, 13.840, 13.818, \
                13.780, 13.737, 13.626, 13.524, 13.440, 13.380, \
                13.310, 13.237, 13.195, 13.166, 13.136, 13.100, \
                13.062, 13.020, 12.976, 12.949, 12.925, 12.898, \
                12.874, 12.850, 12.833, 12.815, 12.809, 12.808, \
                12.804, 12.800, 12.794, 12.790, 12.781]
    ab[:, 9] = [13.180, 13.030, 12.890, 12.660, 12.610, 12.580, \
                12.530, 12.460, 12.375, 12.300, 12.230, 12.200, \
                12.130, 12.020, 11.775, 11.630, 11.570, 11.500, \
                11.390, 11.360, 11.325, 11.285, 11.235, 11.190, \
                11.170, 11.150, 11.085, 11.060, 11.040, 11.035, \
                11.020, 10.985, 10.980, 10.970, 10.965, 10.960, \
                10.960, 10.960, 10.960, 10.960, 10.960]
    ab[:, 10] = [14.660, 14.670, 14.680, 14.690, 14.710, 14.730, \
                 14.750, 14.780, 14.800, 14.820, 14.800, 14.770, \
                 14.770, 14.720, 14.470, 14.360, 14.310, 14.330, \
                 14.380, 14.430, 14.450, 14.480, 14.520, 14.550, \
                 14.590, 14.640, 14.740, 14.840, 14.900, 15.000, \
                 15.080, 15.170, 15.170, 15.170, 15.170, 15.170, \
                 15.170, 15.170, 15.170, 15.170, 15.170]
    ab[:, 11] = [13.910, 13.930, 13.940, 13.960, 13.990, 14.010, \
                 14.035, 14.060, 14.080, 14.100, 14.090, 14.080, \
                 14.080, 14.060, 13.880, 13.850, 13.900, 13.970, \
                 14.040, 14.100, 14.135, 14.170, 14.210, 14.280, \
                 14.340, 14.400, 14.530, 14.630, 14.710, 14.800, \
                 14.880, 14.970, 14.970, 14.970, 14.970, 14.970, \
                 14.970, 14.970, 14.970, 14.970, 14.970]
    ab[:, 12] = [15.080, 14.970, 14.880, 14.825, 14.760, 14.710, \
                 14.652, 14.640, 14.620, 14.600, 14.570, 14.550, \
                 14.540, 14.530, 14.440, 14.350, 14.280, 14.250, \
                 14.220, 14.170, 14.138, 14.120, 14.100, 14.070, \
                 14.045, 14.020, 14.010, 14.005, 13.995, 13.993, \
                 13.991, 13.990, 13.990, 13.990, 13.990, 13.990, \
                 13.990, 13.990, 13.990, 13.990, 13.990]
    ab[:, 13] = [15.410, 15.400, 15.390, 15.390, 15.380, 15.380, \
                 15.390, 15.410, 15.420, 15.430, 15.430, 15.420, \
                 15.410, 15.400, 15.390, 15.415, 15.440, 15.460, \
                 15.500, 15.540, 15.560, 15.570, 15.600, 15.610, \
                 15.630, 15.680, 15.760, 15.850, 15.910, 15.980, \
                 16.060, 16.140, 16.210, 16.310, 16.390, 16.540, \
                 16.620, 16.620, 16.620, 16.620, 16.620]
    ab[:, 14] = [15.795, 15.710, 15.580, 15.490, 15.390, 15.310, \
                 15.220, 15.185, 15.162, 15.140, 15.110, 15.080, \
                 15.050, 15.030, 14.900, 14.810, 14.710, 14.670, \
                 14.610, 14.520, 14.488, 14.460, 14.440, 14.400, \
                 14.360, 14.320, 14.300, 14.290, 14.285, 14.285, \
                 14.290, 14.300, 14.310, 14.320, 14.330, 14.360, \
                 14.380, 14.380, 14.380, 14.380, 14.380]
    ab[:, 15] = [15.270, 15.295, 15.310, 15.320, 15.360, 15.400, \
                 15.415, 15.460, 15.483, 15.490, 15.500, 15.515, \
                 15.530, 15.530, 15.200, 15.050, 15.000, 15.020, \
                 15.090, 15.154, 15.177, 15.215, 15.260, 15.300, \
                 15.355, 15.410, 15.530, 15.630, 15.700, 15.807, \
                 15.910, 16.010, 16.010, 16.010, 16.010, 16.010, \
                 16.010, 16.010, 16.010, 16.010, 16.010]
    ab[:, 16] = [13.820, 13.730, 13.600, 13.520, 13.420, 13.350, \
                 13.300, 13.280, 13.270, 13.265, 13.250, 13.230, \
                 13.210, 13.190, 13.130, 13.070, 13.010, 13.005, \
                 12.985, 12.940, 12.922, 12.920, 12.895, 12.880, \
                 12.860, 12.843, 12.850, 12.870, 12.880, 12.900, \
                 12.920, 12.940, 12.960, 12.990, 13.020, 13.060, \
                 13.090, 13.090, 13.090, 13.090, 13.090]
    ab[:, 17] = [11.540, 11.565, 11.580, 11.640, 11.680, 11.730, \
                 11.760, 11.830, 11.870, 11.890, 11.920, 11.930, \
                 11.950, 11.970, 11.990, 12.010, 12.076, 12.135, \
                 12.210, 12.290, 12.345, 12.338, 12.380, 12.393, \
                 12.415, 12.450, 12.435, 12.414, 12.375, 12.330, \
                 12.310, 12.280, 12.255, 12.210, 12.160, 12.115, \
                 12.100, 12.070, 12.070, 12.070, 12.070]
    ab[:, 18] = [17.749, 17.659, 17.562, 17.460, 17.359, 17.250, \
                 17.146, 17.049, 17.006, 16.945, 16.902, 16.881, \
                 16.828, 16.785, 16.610, 16.436, 16.285, 16.173, \
                 16.068, 15.960, 15.900, 15.860, 15.800, 15.740, \
                 15.680, 15.620, 15.550, 15.508, 15.480, 15.460, \
                 15.450, 15.440, 15.435, 15.430, 15.427, 15.423, \
                 15.425, 15.430, 15.435, 15.440, 15.456]
    ab[:, 19] = [16.780, 16.740, 16.692, 16.634, 16.570, 16.502, \
                 16.431, 16.363, 16.330, 16.300, 16.271, 16.256, \
                 16.229, 16.202, 16.105, 16.020, 15.940, 15.887, \
                 15.829, 15.776, 15.757, 15.740, 15.718, 15.700, \
                 15.681, 15.659, 15.632, 15.615, 15.604, 15.594, \
                 15.586, 15.582, 15.581, 15.586, 15.598, 15.609, \
                 15.617, 15.628, 15.649, 15.662, 15.695]
    ab[:, 20] = [9.500, 9.510, 9.530, 9.560, 9.600, 9.649, 9.695, \
                 9.736, 9.758, 9.772, 9.791, 9.800, 9.820, 9.833, \
                 9.920, 10.011, 10.107, 10.186, 10.261, 10.363, \
                 10.406, 10.449, 10.518, 10.582, 10.660, 10.730, \
                 10.850, 10.937, 11.023, 11.110, 11.165, 11.216, \
                 11.268, 11.310, 11.360, 11.400, 11.427, 11.450, \
                 11.485, 11.500, 11.530]
    ab[:, 21] = [10.750, 10.790, 10.839, 10.883, 10.925, 10.966, \
                 11.011, 11.052, 11.069, 11.089, 11.106, 11.111, \
                 11.129, 11.143, 11.206, 11.283, 11.369, 11.429, \
                 11.508, 11.600, 11.642, 11.686, 11.750, 11.815, \
                 11.888, 11.980, 12.114, 12.225, 12.312, 12.408, \
                 12.488, 12.562, 12.634, 12.726, 12.842, 12.915, \
                 12.955, 12.998, 13.061, 13.090, 13.148]
    ab[:, 22] = [15.889, 15.853, 15.808, 15.766, 15.728, 15.698, \
                 15.670, 15.643, 15.630, 15.618, 15.605, 15.599, \
                 15.586, 15.573, 15.524, 15.478, 15.437, 15.415, \
                 15.393, 15.372, 15.365, 15.358, 15.351, 15.349, \
                 15.350, 15.354, 15.364, 15.373, 15.383, 15.397, \
                 15.413, 15.434, 15.455, 15.487, 15.530, 15.566, \
                 15.587, 15.617, 15.673, 15.705, 15.785]
    ab[:, 23] = [17.252, 17.190, 17.124, 17.060, 17.000, 16.950, \
                 16.891, 16.841, 16.817, 16.793, 16.769, 16.757, \
                 16.732, 16.706, 16.609, 16.515, 16.423, 16.359, \
                 16.288, 16.220, 16.194, 16.170, 16.141, 16.119, \
                 16.099, 16.078, 16.055, 16.041, 16.031, 16.024, \
                 16.023, 16.025, 16.028, 16.033, 16.041, 16.047, \
                 16.051, 16.058, 16.070, 16.078, 16.095]
    ab[:, 24] = [10.737, 10.773, 10.804, 10.830, 10.868, 10.897, \
                 10.930, 10.959, 10.978, 10.989, 11.002, 11.010, \
                 11.029, 11.045, 11.107, 11.163, 11.233, 11.290, \
                 11.359, 11.442, 11.478, 11.516, 11.566, 11.616, \
                 11.680, 11.773, 11.905, 12.023, 12.120, 12.235, \
                 12.321, 12.400, 12.468, 12.540, 12.637, 12.709, \
                 12.760, 12.819, 12.909, 12.954, 13.061]
    ab[:, 25] = [14.156, 14.137, 14.118, 14.099, 14.080, 14.062, \
                 14.045, 14.030, 14.023, 14.017, 14.012, 14.009, \
                 14.005, 14.002, 13.996, 13.997, 14.003, 14.011, \
                 14.026, 14.053, 14.066, 14.079, 14.096, 14.113, \
                 14.139, 14.173, 14.212, 14.261, 14.321, 14.380, \
                 14.412, 14.436, 14.450, 14.462, 14.467, 14.469, \
                 14.471, 14.475, 14.481, 14.485, 14.496]
    ab[:, 26] = [10.883, 10.922, 10.960, 10.998, 11.031, 11.065, \
                 11.102, 11.138, 11.154, 11.171, 11.187, 11.195, \
                 11.210, 11.226, 11.283, 11.350, 11.438, 11.505, \
                 11.588, 11.681, 11.723, 11.764, 11.822, 11.879, \
                 11.937, 12.014, 12.146, 12.236, 12.321, 12.415, \
                 12.487, 12.539, 12.585, 12.650, 12.735, 12.800, \
                 12.839, 12.892, 12.983, 13.035, 13.176]
    ab[:, 27] = [10.126, 10.165, 10.204, 10.242, 10.283, 10.324, \
                 10.367, 10.409, 10.428, 10.446, 10.463, 10.471, \
                 10.487, 10.502, 10.568, 10.653, 10.750, 10.824, \
                 10.920, 11.024, 11.070, 11.116, 11.179, 11.241, \
                 11.311, 11.406, 11.551, 11.653, 11.721, 11.796, \
                 11.852, 11.892, 11.920, 11.960, 12.011, 12.051, \
                 12.075, 12.106, 12.162, 12.194, 12.279]
    ab[:, 28] = [13.481, 13.272, 13.063, 12.872, 12.730, 12.617, \
                 12.528, 12.461, 12.430, 12.399, 12.365, 12.348, \
                 12.313, 12.276, 12.131, 11.991, 11.858, 11.775, \
                 11.700, 11.626, 11.594, 11.564, 11.522, 11.477, \
                 11.428, 11.379, 11.304, 11.260, 11.233, 11.192, \
                 11.158, 11.133, 11.113, 11.082, 11.061, 11.056, \
                 11.056, 11.057, 11.058, 11.058, 11.059]
    ab[:, 29] = [14.338, 14.094, 13.850, 13.625, 13.455, 13.330, \
                 13.225, 13.134, 13.093, 13.053, 13.016, 12.997, \
                 12.962, 12.926, 12.772, 12.605, 12.460, 12.364, \
                 12.264, 12.179, 12.140, 12.101, 12.051, 12.008, \
                 11.964, 11.913, 11.839, 11.784, 11.740, 11.696, \
                 11.675, 11.646, 11.631, 11.599, 11.566, 11.553, \
                 11.549, 11.548, 11.547, 11.546, 11.544]
    ab[:, 30] = [11.010, 10.970, 10.930, 10.850, 10.830, 10.810, \
                 10.800, 10.790, 10.780, 10.750, 10.720, 10.680, \
                 10.580, 10.500, 10.470, 10.460, 10.460, 10.450, \
                 10.450, 10.440, 10.430, 10.430, 10.420, 10.420, \
                 10.430, 10.440, 10.480, 10.500, 10.520, 10.550, \
                 10.570, 10.590, 10.610, 10.640, 10.670, 10.700, \
                 10.720, 10.740, 10.780, 10.800, 10.860]
    ab[:, 31] = [12.190, 12.110, 12.060, 11.950, 11.850, 11.790, \
                 11.740, 11.660, 11.630, 11.580, 11.500, 11.450, \
                 11.350, 11.260, 11.160, 11.070, 10.990, 10.870, \
                 10.790, 10.690, 10.590, 10.500, 10.410, 10.360, \
                 10.320, 10.280, 10.210, 10.130, 10.070, 10.020, \
                 9.980, 9.960, 9.970, 9.980, 9.990, 10.000, 10.010, \
                 10.020, 10.030, 10.040, 10.050]
    ab[:, 32] = [10.780, 10.730, 10.670, 10.628, 10.641, 10.663, \
                 10.688, 10.700, 10.687, 10.656, 10.585, 10.546, \
                 10.483, 10.434, 10.456, 10.516, 10.573, 10.617, \
                 10.672, 10.746, 10.777, 10.802, 10.840, 10.882, \
                 10.930, 10.991, 11.080, 11.132, 11.169, 11.219, \
                 11.249, 11.312, 11.415, 11.511, 11.623, 11.730, \
                 11.799, 11.890, 12.028, 12.098, 12.283]
    ab[:, 33] = [11.510, 11.450, 11.332, 11.290, 11.270, 11.270, \
                 11.288, 11.290, 11.290, 11.250, 11.125, 11.020, \
                 10.900, 10.830, 10.700, 10.764, 10.820, 10.845, \
                 10.902, 10.980, 11.010, 11.035, 11.070, 11.110, \
                 11.144, 11.205, 11.307, 11.375, 11.430, 11.518, \
                 11.593, 11.648, 11.690, 11.750, 11.890, 11.960, \
                 12.015, 12.090, 12.210, 12.280, 12.420]
    ab[:, 34] = [16.350, 16.300, 16.263, 16.195, 16.160, 16.100, \
                 16.070, 16.025, 16.010, 15.980, 15.960, 15.958, \
                 15.930, 15.890, 15.858, 15.807, 15.737, 15.732, \
                 15.702, 15.670, 15.648, 15.635, 15.624, 15.613, \
                 15.600, 15.587, 15.582, 15.583, 15.589, 15.602, \
                 15.619, 15.633, 15.649, 15.682, 15.752, 15.825, \
                 15.860, 15.940, 16.020, 16.080, 16.250]
    ab[:, 35] = [14.50, 14.50, 14.41, 14.45, 14.42, 14.45, 14.44, \
                 14.49, 14.46, 14.49, 14.47, 14.53, 14.60, 14.58, \
                 14.55, 14.49, 14.53, 14.53, 14.56, 14.66, 14.67, \
                 14.63, 14.68, 14.70, 14.77, 14.83, 14.86, 14.97, \
                 14.99, 15.07, 15.12, 15.16, 15.20, 15.27, 15.35, \
                 15.44, 15.38, 15.36, 15.35, 15.34, 15.30]
    ab[:, 36] = [12.590, 12.590, 12.598, 12.556, 12.540, 12.494, \
                 12.515, 12.509, 12.497, 12.427, 12.270, 12.192, \
                 12.035, 11.878, 11.752, 11.768, 11.815, 11.830, \
                 11.883, 11.950, 11.959, 11.968, 11.989, 12.010, \
                 12.024, 12.080, 12.189, 12.209, 12.240, 12.288, \
                 12.340, 12.361, 12.400, 12.400, 12.400, 12.400, \
                 12.400, 12.400, 12.400, 12.400, 12.400]
    ab[:, 37] = [13.5600, 13.5600, 13.4400, 13.4600, 13.5800, 13.5760, \
                 13.5480, 13.6440, 13.6780, 13.6900, 13.7000, 13.7200, \
                 13.7480, 13.7680, 13.8420, 13.9500, 14.0000, 14.1100, \
                 14.1540, 14.2500, 14.3040, 14.3260, 14.3900, 14.3660, \
                 14.4340, 14.5520, 14.6627, 14.7020, 14.7900, 14.8200, \
                 14.8780, 14.9218, 14.9320, 14.8800, 14.9777, 15.0173, \
                 15.0647, 15.1137, 15.1137, 15.1137, 15.1137]
    ab[:, 38] = [8.95000, 8.95000, 8.87000, 8.73600, 8.65000, 8.62000, \
                 8.52800, 8.49600, 8.47800, 8.46200, 8.43400, 8.40200, \
                 8.35400, 8.32600, 8.22000, 8.15000, 8.08000, 7.96400, \
                 7.88800, 7.79000, 7.70800, 7.61600, 7.58000, 7.48000, \
                 7.45400, 7.39600, 7.34200, 7.23422, 7.20000, 7.15864, \
                 7.11746, 7.09000, 7.09000, 7.09000, 7.09000, 7.09000, \
                 7.09000, 7.09000, 7.09000, 7.09000, 7.09000]
    ab[:, 39] = [13.3800, 13.3800, 13.3912, 13.3576, 13.3700, 13.3724, \
                 13.3819, 13.3912, 13.3959, 13.4007, 13.4054, 13.4077, \
                 13.4123, 13.4170, 13.4360, 13.4600, 13.5064, 13.5620, \
                 13.6630, 13.7901, 13.8228, 13.8848, 13.9700, 14.0092, \
                 14.1128, 14.1708, 14.3481, 14.4030, 14.5000, 14.5856, \
                 14.5978, 14.6716, 14.6716, 14.6716, 14.6716, 14.6716, \
                 14.6716, 14.6716, 14.6716, 14.6716, 14.6716]
    ab[:,40]= [9.71000, 9.71000, 9.69960, 9.71940, 9.75000, 9.78720, \
               9.80216, 9.81480, 9.82124, 9.82851, 9.83550, 9.83918, \
               9.84723, 9.85676, 9.89451, 9.96000, 10.0500, 10.1132, \
               10.2140, 10.3200, 10.3714, 10.4180, 10.4700, 10.5340, \
               10.6116, 10.7332, 10.9026, 10.9693, 11.0300, 11.1424, \
               11.2164, 11.2752, 11.3143, 11.2954, 11.4150, 11.5300, \
               11.6312, 11.7700, 11.7700, 11.7700, 11.7700]
    ab[:, 41] = [10.7060, 10.7060, 10.7070, 10.8412, 10.8751, 10.9346, \
                 10.9988, 11.0508, 11.0876, 11.1117, 11.1198, 11.1309, \
                 11.1536, 11.1253, 11.1920, 11.2856, 11.3744, 11.4628, \
                 11.5659, 11.6423, 11.6924, 11.7352, 11.8133, 11.8565, \
                 11.9284, 12.0111, 12.1394, 12.2520, 12.3412, 12.4391, \
                 12.5198, 12.5946, 12.6563, 12.7585, 12.8850, 12.8850, \
                 12.8850, 12.8850, 12.8850, 12.8850, 12.8850]
    ab[:, 42] = [12.6650, 12.6650, 12.6650, 12.6689, 12.6060, 12.5745, \
                 12.4804, 12.3900, 12.3252, 12.2909, 12.2230, 12.1593, \
                 12.0765, 12.0277, 11.7515, 11.6150, 11.5000, 11.4408, \
                 11.3709, 11.3090, 11.2841, 11.2483, 11.2300, 11.1848, \
                 11.1531, 11.1270, 11.0906, 11.0527, 11.0500, 11.0405, \
                 11.0333, 11.0334, 11.0415, 11.0607, 11.0660, 11.0623, \
                 11.0590, 11.0628, 11.1040, 11.1040, 11.1040]
    ab[:, 43] = [13.4480, 13.4480, 13.4480, 13.4304, 13.3760, 13.3201, \
                 13.2367, 13.1485, 13.1038, 13.0542, 12.9668, 12.9094, \
                 12.8061, 12.7308, 12.5630, 12.4730, 12.4300, 12.3768, \
                 12.3296, 12.2690, 12.2535, 12.2343, 12.2000, 12.1752, \
                 12.1656, 12.1231, 12.0979, 12.0637, 12.0610, 12.0447, \
                 12.0358, 12.0350, 12.0426, 12.0589, 12.0848, 12.0996, \
                 12.1069, 12.1220, 12.1220, 12.1220, 12.1220]
    ab[:, 44] = [11.9290, 11.9290, 11.9290, 11.9046, 11.8710, 11.8479, \
                 11.8000, 11.7477, 11.7263, 11.7090, 11.6945, 11.6854, \
                 11.6673, 11.6534, 11.6223, 11.5960, 11.5740, 11.5590, \
                 11.5371, 11.5154, 11.5068, 11.5009, 11.5000, 11.5009, \
                 11.5019, 11.5053, 11.5110, 11.5168, 11.5520, 11.5775, \
                 11.5920, 11.6202, 11.6513, 11.6907, 11.7590, 11.8070, \
                 11.8338, 11.8690, 11.8690, 11.8690, 11.8690]
    ab[:, 45] = [11.7000, 11.7000, 11.6640, 11.5452, 11.5300, 11.5080, \
                 11.4916, 11.3878, 10.6773, 10.2642, 10.2610, 10.2566, \
                 10.2443, 10.2311, 10.1995, 10.2156, 10.2294, 10.2487, \
                 10.2753, 10.3380, 10.3431, 10.3561, 10.3827, 10.4104, \
                 10.4257, 10.4725, 10.5476, 10.5986, 10.6300, 10.6840, \
                 10.7300, 10.7407, 10.7600, 10.7600, 10.7600, 10.7600, \
                 10.7600, 10.7600, 10.7600, 10.7600, 10.7600]
    ab[:, 46] = [11.9800, 11.9800, 11.9368, 11.7878, 11.7300, 11.6804, \
                 11.6335, 11.5881, 11.5662, 11.5266, 11.4671, 11.4369, \
                 11.3756, 11.3132, 11.2430, 11.1481, 11.0673, 10.9546, \
                 10.8776, 10.7300, 10.6385, 10.5474, 10.4362, 10.3487, \
                 10.2950, 10.2300, 10.1340, 10.0134, 9.93000, 9.87329, \
                 9.81000, 9.73858, 9.72000, 9.72000, 9.72000, 9.72000, \
                 9.72000, 9.72000, 9.72000, 9.72000, 9.72000]
    ab[:, 47] = [13.7300, 13.7300, 13.6836, 13.5524, 13.3700, 13.3659, \
                 13.2409, 13.0914, 13.0160, 12.9800, 12.9800, 12.9800, \
                 12.9312, 12.8542, 12.6721, 12.5692, 12.4777, 12.4015, \
                 12.3132, 12.2700, 12.2428, 12.1974, 12.1997, 12.1187, \
                 12.0558, 12.0370, 11.9921, 11.9808, 11.9500, 11.9355, \
                 11.9100, 11.8901, 11.8895, 11.8679, 11.8612, 11.8664, \
                 11.8605, 11.8813, 11.9400, 11.9400, 11.9400]
    ab[:, 48] = [12.4400, 12.4400, 12.3899, 12.3234, 12.3017, 12.2882, \
                 12.2779, 12.2677, 12.2625, 12.2733, 12.2800, 12.2691, \
                 12.1489, 11.9884, 11.9632, 11.8777, 11.8386, 11.8197, \
                 11.8198, 11.8300, 11.8369, 11.8437, 11.8516, 11.8599, \
                 11.8679, 11.8826, 11.9407, 11.9594, 11.9600, 12.0066, \
                 12.0400, 12.0421, 12.1004, 12.1041, 12.1510, 12.2075, \
                 12.2522, 12.3002, 12.3000, 12.3000, 12.3000]
    ab[:, 49] = [15.1000, 15.1000, 15.0400, 14.9527, 14.8300, 14.5632, \
                 14.3972, 14.2112, 14.1330, 14.0266, 13.9064, 13.8616, \
                 13.7702, 13.6812, 13.4022, 13.0700, 12.7400, 12.3942, \
                 12.1244, 11.7000, 11.4582, 11.2224, 11.0100, 10.6924, \
                 10.4476, 10.1932, 9.86280, 9.46970, 9.21000, 9.01880, \
                 8.80220, 8.55000, 8.55000, 8.55000, 8.55000, 8.55000, \
                 8.55000, 8.55000, 8.55000, 8.55000, 8.55000]
    ab[:, 50] = [14.1000, 14.1000, 14.1000, 13.9600, 13.7600, 13.5200, \
                 13.4500, 13.3600, 13.2671, 13.2056, 13.1505, 13.1267, \
                 13.0874, 13.0555, 12.9277, 12.8200, 12.6900, 12.5800, \
                 12.5000, 12.4500, 12.4100, 12.4200, 12.4100, 12.3600, \
                 12.3556, 12.3556, 12.3100, 12.3000, 12.2756, 12.2987, \
                 12.3100, 12.2544, 12.2514, 12.3034, 12.2867, 12.2867, \
                 12.3830, 12.3622, 12.3500, 12.3500, 12.3500]
    ab[:,51] = [12.194600, 12.228300, 12.270300, 12.321600, 12.354300, \
                12.393000, 12.426700, 12.465400, 12.479000, 12.496700, \
                12.515500, 12.523100, 12.534000, 12.538500, 12.523600, \
                12.551700, 12.626000, 12.682800, 12.769200, 12.860300, \
                12.903400, 12.945100, 12.996200, 13.052800, 13.126000, \
                13.209800, 13.338300, 13.440400, 13.530800, 13.625200, \
                13.706700, 13.789400, 13.856400, 13.955600, 14.046900, \
                14.128700, 14.173900, 14.244800, 14.352100, 14.404500, \
                14.556200]
    ab[:,52]= [11.890830, 11.953902, 12.001304, 12.038769, 12.068400, \
               12.117474, 12.158472, 12.193123, 12.204728, 12.247052, \
               12.256253, 12.264238, 12.288051, 12.299415, 12.346685, \
               12.400056, 12.476015, 12.550985, 12.628268, 12.718414, \
               12.767287, 12.800667, 12.874192, 12.927598, 13.002908, \
               13.089129, 13.226055, 13.322492, 13.418968, 13.513610, \
               13.603105, 13.685889, 13.728714, 13.841738, 13.964003, \
               14.030111, 14.105518, 14.132800, 14.132800, 14.132800, \
               14.132800]
    ab[:,53]=[ 11.850000, 11.850000, 11.850000, 11.862280, 11.753000, \
               11.710168, 11.583240, 11.455560, 11.391720, 11.327880, \
               11.264040, 11.232120, 11.168280, 11.104440, 10.878012, \
               10.764100, 10.590400, 10.528520, 10.427080, 10.354000, \
               10.322270, 10.293700, 10.238000, 10.182120, 10.128600, \
               10.090520, 10.032000, 9.9762198, 9.9617100, 9.9360398, \
               9.9185802, 9.9072905, 9.9031467, 9.9140937, 9.9183197, \
               9.9221601, 9.9244636, 9.9169998, 9.9169998, 9.9169998, \
               9.9169998]                                           
    ab[:,54]=[12.373000, 12.373000, 12.373000, 12.356440, 12.367000, \
              12.397760, 12.395280, 12.386800, 12.403720, 12.419600, \
              12.422360, 12.421640, 12.412320, 12.403216, 12.223200, \
              11.962500, 11.955000, 11.989480, 12.050820, 12.113800, \
              12.140126, 12.165320, 12.207000, 12.247120, 12.309800, \
              12.370320, 12.459800, 12.549230, 12.609700, 12.681920, \
              12.752000, 12.815400, 12.880920, 12.958500, 13.082100, \
              13.161000, 13.220040, 13.175000, 13.175000, 13.175000, \
              13.175000]
    ab[:,55] = [12.431100, 12.467100, 12.510300, 12.547300, 12.601100, \
                12.643000, 12.672900, 12.715400, 12.725200, 12.746000, \
                12.767500, 12.777500, 12.793400, 12.797200, 12.814100, \
                12.861100, 12.933300, 12.994100, 13.078300, 13.172600, \
                13.212700, 13.258900, 13.314200, 13.377200, 13.450300, \
                13.528800, 13.662500, 13.764600, 13.859700, 13.955800, \
                14.044400, 14.125700, 14.184500, 14.283600, 14.377200, \
                14.483600, 14.542900, 14.592500, 14.699500, 14.761600, \
                14.911200]
    ab[:,56]=[12.024000, 12.024000, 12.024000, 12.010360, 11.966999, \
              11.948019, 11.845119, 11.752359, 11.703140, 11.617559, \
              11.453240, 11.356260, 11.181759, 11.091120, 10.767380, \
              10.643000, 10.570000, 10.488140, 10.471140, 10.456500, \
              10.454224, 10.448364, 10.446000, 10.416080, 10.418639, \
              10.420500, 10.441330, 10.444580, 10.465000, 10.494888, \
              10.510241, 10.525400, 10.539100, 10.558599, 10.585599, \
              10.609200, 10.623360, 10.631999, 10.631999, 10.631999, \
              10.631999]



    # use scipy interp here as the very non-linear wavelength scale
    # precludes the use of womashrebin --- this is going finer anyway
    abcurve=scipyrebin(waveab, ab[:, idstar], wave)
    return abcurve


def baryvel(dje):
#       Based on Fortran routine of P. STUMPFF (IBM-VERSION 1979):
#                 ASTRON. ASTROPHYS.  SUPPL. SER. 41, 1 (1980)
#                 M. H. SLOVAK (VAX 11/780 IMPLEMENTATION 1986)
#                 G. TORRES (1989)


    au=1.4959787e8
    dct0=2415020.0
    dcjul=36525.0
    dcbes=0.313
    dctrop=365.24219572
    dc1900=1900.0
    dc2pi=2.0*np.pi
    forbel=[None]*7
    sorbel=[None]*17
    sn=[None]*4
    dvelh=[None]*3
    dvelb=[None]*3

#  constants dcfel(i,k) of fast-changing elements

#                  i = 1           i = 2         i = 3
    dcfel=[[1.7400353e+00, 6.2833195099091e+02, 5.2796e-06], \
                 [6.2565836e+00, 6.2830194572674e+02,-2.6180e-06], \
                 [4.7199666e+00, 8.3997091449254e+03,-1.9780e-05], \
                 [1.9636505e-01, 8.4334662911720e+03,-5.6044e-05], \
                 [4.1547339e+00, 5.2993466764997e+01, 5.8845e-06], \
                 [4.6524223e+00, 2.1354275911213e+01, 5.6797e-06], \
                 [4.2620486e+00, 7.5025342197656e+00, 5.5317e-06], \
                 [1.4740694e+00, 3.8377331909193e+00, 5.6093e-06]]

#  constants dceps and ccsel(i,k) of slowly changing elements

#                  i = 1      i = 2       i = 3

    dceps=[4.093198e-01,-2.271110e-04,-2.860401e-08]

    ccsel=[[1.675104e-02,-4.179579e-05,-1.260516e-07], \
                 [2.220221e-01, 2.809917e-02, 1.852532e-05], \
                 [1.589963e+00, 3.418075e-02, 1.430200e-05], \
                 [2.994089e+00, 2.590824e-02, 4.155840e-06], \
                 [8.155457e-01, 2.486352e-02, 6.836840e-06], \
                 [1.735614e+00, 1.763719e-02, 6.370440e-06], \
                 [1.968564e+00, 1.524020e-02,-2.517152e-06], \
                 [1.282417e+00, 8.703393e-03, 2.289292e-05], \
                 [2.280820e+00, 1.918010e-02, 4.484520e-06], \
                 [4.833473e-02, 1.641773e-04,-4.654200e-07], \
                 [5.589232e-02,-3.455092e-04,-7.388560e-07], \
                 [4.634443e-02,-2.658234e-05, 7.757000e-08], \
                 [8.997041e-03, 6.329728e-06,-1.939256e-09], \
                 [2.284178e-02,-9.941590e-05, 6.787400e-08], \
                 [4.350267e-02,-6.839749e-05,-2.714956e-07], \
                 [1.348204e-02, 1.091504e-05, 6.903760e-07], \
                 [3.106570e-02,-1.665665e-04,-1.590188e-07]]

#  constants of the arguments of the short-period perturbations by
#  the planets:  dcargs(i,k)

#                 i = 1           i = 2
    
    dcargs=[[5.0974222e+00,-7.8604195454652e+02], \
                  [3.9584962e+00,-5.7533848094674e+02], \
                  [1.6338070e+00,-1.1506769618935e+03], \
                  [2.5487111e+00,-3.9302097727326e+02], \
                  [4.9255514e+00,-5.8849265665348e+02], \
                  [1.3363463e+00,-5.5076098609303e+02], \
                  [1.6072053e+00,-5.2237501616674e+02], \
                  [1.3629480e+00,-1.1790629318198e+03], \
                  [5.5657014e+00,-1.0977134971135e+03], \
                  [5.0708205e+00,-1.5774000881978e+02], \
                  [3.9318944e+00, 5.2963464780000e+01], \
                  [4.8989497e+00, 3.9809289073258e+01], \
                  [1.3097446e+00, 7.7540959633708e+01], \
                  [3.5147141e+00, 7.9618578146517e+01], \
                  [3.5413158e+00,-5.4868336758022e+02]]
    
#  amplitudes ccamps(n,k) of the short-period perturbations
    
#       n = 1      n = 2      n = 3      n = 4      n = 5
    
    ccamps= \
        [[-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5,-2.490817e-7], \
               [-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5,-1.823138e-7], \
               [ 6.593466e-7, 1.322572e-5, 9.258695e-6,-4.674248e-7,-3.646275e-7], \
               [ 1.140767e-5,-2.049792e-5,-4.747930e-6,-2.638763e-6,-1.245408e-7], \
               [ 9.516893e-6,-2.748894e-6,-1.319381e-6,-4.549908e-6,-1.864821e-7], \
               [ 7.310990e-6,-1.924710e-6,-8.772849e-7,-3.334143e-6,-1.745256e-7], \
               [-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6,-1.655307e-7], \
               [-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6,-3.736225e-7], \
               [ 3.442177e-7, 2.671323e-6, 1.832858e-6,-2.394688e-7,-3.478444e-7], \
               [ 8.702406e-6,-8.421214e-6,-1.372341e-6,-1.455234e-6,-4.998479e-8], \
               [-1.488378e-6,-1.251789e-5, 5.226868e-7,-2.049301e-7, 0.0], \
               [-8.043059e-6,-2.991300e-6, 1.473654e-7,-3.154542e-7, 0.0], \
               [ 3.699128e-6,-3.316126e-6, 2.901257e-7, 3.407826e-7, 0.0], \
               [ 2.550120e-6,-1.241123e-6, 9.901116e-8, 2.210482e-7, 0.0], \
               [-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.0]]

#  constants of the secular perturbations in longitude ccsec3 and
#  ccsec(n,k)
    
#                 n = 1      n = 2        n = 3

    ccsec3=-7.757020e-08

    ccsec=[[1.289600e-06, 5.550147e-01, 2.076942e+00], \
                 [3.102810e-05, 4.035027e+00, 3.525565e-01], \
                 [9.124190e-06, 9.990265e-01, 2.622706e+00], \
                 [9.793240e-07, 5.508259e+00, 1.559103e+01]]
    
#  sidereal rate dcsld in longitude, rate ccsgd in mean anomaly

    dcsld=1.990987e-07
    ccsgd=1.990969e-07

#  some constants used in the calculation of the lunar contribution

    cckm=3.122140e-05
    ccmld=2.661699e-06
    ccfdi=2.399485e-07

#  constants dcargm(i,k) of the arguments of the perturbations of the
#  motion of the moon

#                  i = 1           i = 2

    dcargm=[[5.1679830e+00, 8.3286911095275e+03], \
       [5.4913150e+00,-7.2140632838100e+03], \
       [5.9598530e+00, 1.5542754389685e+04]]
    
#  amplitudes ccampm(n,k) of the perturbations of the moon
    
#        n = 1       n = 2       n = 3       n = 4
    
    ccampm = \
        [[1.097594e-01, 2.896773e-07, 5.450474e-02, 1.438491e-07], \
               [-2.223581e-02, 5.083103e-08, 1.002548e-02,-2.291823e-08], \
               [1.148966e-02, 5.658888e-08, 8.249439e-03, 4.063015e-08]]

#  ccpamv = a*m*dl/dt (planets); dc1mme = 1 - mass(earth+moon)
    
    ccpamv=[8.326827e-11,1.843484e-11,1.988712e-12,1.881276e-12]
    dc1mme=0.99999696
    
#  program execution begins
    
#  time-arguments
    
    dt = (dje - dct0)/dcjul
    t = dt
    dtsq = dt * dt
    tsq = dtsq
#  values of all elements for the instant dje

    for k in range(0,8):
        dlocal=np.fmod(dcfel[k][0]+dt*dcfel[k][1]+dtsq*dcfel[k][2], dc2pi)
        if (k == 0):
            dml = dlocal
        if (k != 0):
            forbel[k-1]=dlocal

    g=forbel[0]
    deps=np.fmod(dceps[0]+dt*dceps[1]+dtsq*dceps[2], dc2pi)

    for k in range(0,17):
        sorbel[k] = np.fmod(ccsel[k][0]+t*ccsel[k][1]+tsq*ccsel[k][2], dc2pi)
    e=sorbel[0]
#  secular perturbations in longitude
    for k in range(0,4):
        a =np.fmod(ccsec[k][1]+t*ccsec[k][2], dc2pi)
        sn[k] = np.sin(a)

#  periodic perturbations of the emb (earth-moon barycenter)

    pertl = ccsec[0][0]*sn[0] +ccsec[1][0]*sn[1] \
        +(ccsec[2][0]+t*ccsec3)*sn[2] +ccsec[3][0]*sn[3]

    pertld = 0.0
    pertr =  0.0
    pertrd = 0.0
    
    for k in range(0,15):
        a = np.fmod(dcargs[k][0]+dt*dcargs[k][1], dc2pi)
        cosa = np.cos(a)
        sina = np.sin(a)
        pertl = pertl + ccamps[k][0]*cosa + ccamps[k][1]*sina
        pertr = pertr + ccamps[k][2]*cosa + ccamps[k][3]*sina
        if (k < 11):
            pertld=pertld+(ccamps[k][1]*cosa-ccamps[k][0]*sina)*ccamps[k][4]
            pertrd=pertrd+(ccamps[k][3]*cosa-ccamps[k][2]*sina)*ccamps[k][4]


#  elliptic part of the motion of the emb
    esq = e * e
    dparam = 1.0 - esq
    param = dparam
    twoe = e + e
    twog = g + g
    phi = twoe*((1.0 - esq*(1.0/8.0))*np.sin(g) + e*(5.0/8.0)*np.sin(twog) \
                      + esq*0.5416667*np.sin(g+twog) )
    f = g + phi
    sinf = np.sin(f)
    cosf = np.cos(f)
    dpsi = dparam/(1.0 + e*cosf)
    phid = twoe*ccsgd*((1.0+esq*1.5)*cosf+e*(1.25-sinf*sinf*0.5))
    psid = ccsgd*e*sinf/np.sqrt(param)

#  perturbed heliocentric motion of the emb

    d1pdro = (1.0 + pertr)
    drd = d1pdro*(psid+dpsi*pertrd)
    drld = d1pdro*dpsi*(dcsld+phid+pertld)
    dtl = np.fmod(dml+phi+pertl, dc2pi)
    dsinls = np.sin(dtl)
    dcosls = np.cos(dtl)
    dxhd = drd*dcosls - drld*dsinls
    dyhd = drd*dsinls + drld*dcosls

#  influence of eccentricity, evection and variation of the geocentric
#  motion of the moon

    pertl =  0.0
    pertld = 0.0
    pertp =  0.0
    pertpd = 0.0
    
    for k in range(0,3):
        a = np.fmod(dcargm[k][0] + dt*dcargm[k][1], dc2pi)
        sina = np.sin(a)
        cosa = np.cos(a)
        pertl   = pertl  + ccampm[k][0]*sina
        pertld  = pertld + ccampm[k][1]*cosa
        pertp   = pertp  + ccampm[k][2]*cosa
        pertpd  = pertpd - ccampm[k][3]*sina


#  heliocentric motion of the earth

    tl =  forbel[1] + pertl
    sinlm = np.sin(tl)
    coslm = np.cos(tl)
    sigma = cckm/(1.0 + pertp)
    a = sigma*(ccmld+pertld)
    b = sigma*pertpd
    dxhd = dxhd + a*sinlm + b*coslm
    dyhd = dyhd - a*coslm + b*sinlm
    dzhd =      - sigma*ccfdi*np.cos(forbel[2])

#  barycentric motion of the earth

    dxbd = dxhd*dc1mme
    dybd = dyhd*dc1mme
    dzbd = dzhd*dc1mme

    for k in range(0,4):
        plon = forbel[k+3]
        pomg = sorbel[k+1]
        pecc = sorbel[k+9]
        tl = np.fmod(plon+2.0*pecc*np.sin(plon-pomg), dc2pi)
#        sinlp[k] = sin(tl)
#        coslp[k] = np.cos(tl)
        dxbd = dxbd + ccpamv[k]*(np.sin(tl) + pecc*np.sin(pomg))
        dybd = dybd - ccpamv[k]*(np.cos(tl) + pecc*np.cos(pomg))
        dzbd = dzbd - ccpamv[k]*sorbel[k+13]*np.cos(plon-sorbel[k+5])

#  transition to mean equator of date

    dcosep = np.cos(deps)
    dsinep = np.sin(deps)
    dyahd = dcosep*dyhd - dsinep*dzhd
    dzahd = dsinep*dyhd + dcosep*dzhd
    dyabd = dcosep*dybd - dsinep*dzbd
    dzabd = dsinep*dybd + dcosep*dzbd
#
    dvelh[0] = dxhd*au
    dvelh[1] = dyahd*au
    dvelh[2] = dzahd*au
    dvelb[0] = dxbd*au
    dvelb[1] = dyabd*au
    dvelb[2] = dzabd*au

    return dvelh,dvelb

def envelope(array, window_size):
    win_num=len(array)//window_size
    maxarr=np.zeros(win_num)
    minarr=np.zeros(win_num)
    for i in range(0,win_num):
        maxarr[i]=max(array[i*window_size:(i+1)*window_size])
        minarr[i]=min(array[i*window_size:(i+1)*window_size])
    return maxarr, minarr

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    import scipy.interpolate
    import scipy.ndimage
    """Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    """
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print('[congrid] dimensions error.')
        print('This routine currently only support')
        print('rebinning to the same number of dimensions.')
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + list(range( ndims - 1 ))
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = np.mgrid[nslices]

        newcoords_dims = list(range(np.rank(newcoords)))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print('Congrid error: Unrecognized interpolation type.')
        print("Currently only 'neighbour', 'nearest','linear'") 
        print("and 'spline' are supported.")
        return None

def crosscorrelation(x, y, maxlag):
    """
    Cross correlation with a maximum number of lags.

    `x` and `y` must be one-dimensional numpy arrays with the same length.

    This computes the same result as
        numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag]

    The return vaue has length 2*maxlag + 1.
    np.argmax of the return correlation vector gives the value to add to y
    np.argmax-maxlag

    """
    from numpy.lib.stride_tricks import as_strided
    py = np.pad(y.conj(), 2*maxlag, mode='constant')
    T = as_strided(py[2*maxlag:], shape=(2*maxlag+1, len(y) + 2*maxlag),
                   strides=(-py.strides[0], py.strides[0]))
    px = np.pad(x, maxlag, mode='constant')
    return T.dot(px)

def xcor(x, y, xfactor, maxlag):
    npix=len(x)
    xbig=congrid(x,(npix*xfactor,),minusone=True)
    ybig=congrid(y,(npix*xfactor,),minusone=True)
    cc=crosscorrelation(xbig,ybig,maxlag)
    shift=np.argmax(cc)-maxlag
    shift=shift/xfactor
    if (np.argmax(cc) == 0) or (np.argmax(cc) == npix - 1):
        print('Cross-correlation failed!  Try larger maxlag?')
        print('Setting shift to zero...\n')
        shift=0.0
    return shift
    
def secondcat(wave,obj1,obj2,sig1,sig2,secondtimeflag,wavebeg,waveend,brscale):
    print('\nCombining two calibrated spectra for second-order correction')
    waveblue=wave
    wavered=wave
    fluxblue=obj1
    fluxred=obj2
    fluxtemp=fluxred
    if (not secondtimeflag):
        print('Plotting blue side as blue, red side as red')
        indexblue=get_element(waveblue,wavered[0])
        indexred=get_element(wavered,waveblue[-1])
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
            endpoints=plt.ginput(2)
            element1=endpoints[0][0]
            element2=endpoints[1][0]
            plt.plot(endpoints,'ro')
            if (element1 > element2):
                element1,element2=element2,element1
            binend=get_element(waveblue,element2)
            binbeg=get_element(waveblue,element1)
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
        binbegred=get_element(wavered,wavebeg)
        binendred=get_element(wavered,waveend)
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
            endpoints=plt.ginput(2)
            element1=endpoints[0][0]
            element2=endpoints[1][0]
            plt.plot(endpoints,'ro')
            if (element1 > element2):
                element1,element2=element2,element1
            binend=get_element(waveblue,element2)
            binbeg=get_element(waveblue,element1)
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
        binbegred=get_element(wavered,wavebeg)
        binendred=get_element(wavered,waveend)
    else:
        binbeg=get_element(waveblue,wavebeg)
        binend=get_element(waveblue,waveend)
        binbegred=get_element(wavered,wavebeg)
        binendred=get_element(wavered,waveend)
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
    
def scaleparse():
    import re
    done=False
    while (not done):
        scaleparse=input('Enter new ymin and ymax values: ')
        if (scaleparse != ''):
            scaleparse=scaleparse.strip()
            scalesplit=re.split('\W+',scaleparse)
            if (len(scalesplit) > 1):
                try:
                    ymin=float(scalesplit[0])
                    ymax=float(scalesplit[1])
                except ValueError:
                    print('Please enter numbers, ### ###,\n')
                    print('###,###,or ###-###\n')
                    print('You entered {}\n'.format(scaleparse))
                else:
                    done=True
            else:
                print('Enter more than one number\n')
    return  ymin, ymax
