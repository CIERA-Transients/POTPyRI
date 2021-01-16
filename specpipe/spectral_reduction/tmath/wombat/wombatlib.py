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
        question='(y)es or (n)o? (default/return = {}) '.format(default)
        reply=input(question)
        reply=reply.strip()
        if len(reply) == 0:
            reply=default
        answer=reply.lower()[0]
    return answer




def airtovac(wave):
    """
    AIR to VAC from Goddard IDL library

    PURPOSE:
       Convert air wavelengths to vacuum wavelengths 
    EXPLANATION:
       Wavelengths are corrected for the index of refraction of air under 
       standard conditions.  Wavelength values below 2000 A will not be 
       altered.  Uses the IAU standard for conversion given in Morton 
       (1991 Ap.J. Suppl. 77, 119)

    INPUT/OUTPUT:
       WAVE - Wavelength in Angstroms, numpy vector
               WAVE should be input as air wavelength(s)
               returns vacuum

    EXAMPLE:
       If the air wavelength is  W = 6056.125 (a Krypton line), then 
       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019

    METHOD:
       See Morton (Ap. J. Suppl. 77, 119) for the formula used
    """
    import numpy as np
    wavenumber_squared=(10000./wave)**2

    factor=1. + 6.4328e-5 + (2.94981e-2)/(146. - wavenumber_squared) \
            + (2.5540e-4)/(41.-wavenumber_squared)
    type(wave)
    vaclocs=np.where(wave < 2000)
    factor[vaclocs]=1.
    wave=wave*factor
    return wave


def wommkatmdisp(hop):
    light_speed=2.99792458e10
    h_planck=6.6260755e-27
    k_boltzmann=1.380658e-16
    print('This routine will create a curve of the dispersion of light')
    print('by the atmosphere in arc seconds per 1 Angstrom bin')
    print('\n')
    print('Enter airmass for the calculation: ')
    airmass=inputter('(default = 1.5): ','float',True,1.5)
    if (airmass < 1.0) or (airmass > 7.0):
        airmass=1.5
    print('Enter temperature at telescope in degrees Celsius: ')
    temp=inputter('(default = 7C [45F]): ','float',True,7.0)
    if (temp <= -100) or (temp >= 100):
        temp = 7.0
    print('Enter barometric pressure at telescope in mm of Hg: ')
    press=inputter('(default = 600 mm Hg): ','float',True,600.0)
    if (press <= 0):
        press=600.0
    print('Enter water vapor pressure at telescope in mm of Hg: ')
    water=inputter('(default = 8 mm Hg): ','float',True,8.0)
    if (water <= 0):
        water=8.0

    waveb,waver=waveparse()
    if (waver < waveb):
        waveb,waver=waver,waveb
    if (waver == waveb):
        waver=waver+1.
    print('\nOK, calculating dispersion of light in arc seconds over the')
    print('range {}A to {}A at temperature {}C, pressure {} mm Hg,'.format(waveb,waver,temp,press))
    print('and water vapor pressure {} mm Hg at airmass {}.'.format(water,airmass)) 
    print('Zero is set at 5000A.')

    wave=np.arange(waveb,waver+1.)
    vacuum=airtovac(wave)
    lfactor=(1.e4/vacuum)**2
    waterfactor = water*((0.0624-0.000680*lfactor)/(1.0+0.003661*temp))
    nstp = 1E-6*(64.328+(29498.1/(146.0-lfactor))+(255.4/(41-lfactor)))
    n = (nstp)*((press*(1.0+(1.049-0.0157*temp)*1E-6*press))/  \
              (720.883*(1.0+0.003661*temp)))-waterfactor*1.0E-6
    n = n+1.0
    five=airtovac(np.array([5000.0]))
    fivefact=(1.e4/five)**2
    wfive = water*((0.0624-0.000680*fivefact)/(1.0+0.003661*temp))
    nstpfive = 1E-6*(64.328+(29498.1/(146.0-fivefact))+(255.4/(41-fivefact)))
    nfive = (nstpfive)*((press*(1.0+(1.049-0.0157*temp)*1E-6*press))/ \
                        (720.883*(1.0+0.003661*temp)))-wfive*1.0E-6
    nfive = nfive+1.0
    cosz=1./airmass
    tanz=(np.sqrt(1-cosz**2))/cosz
    flux=206265.*(n-nfive)*tanz
    spectxt='Atmospheric dispersion curve at z = {}'.format(airmass)
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid')
    plt.title(spectxt)
    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].obname=spectxt
    hop[0].var=np.ones(wave.shape)
    hop[0].header=''
    return hop

def womhop(hop):
    hopnum=0
    while (hopnum < 1) or (hopnum > HOPSIZE):
        hopnum=inputter('Store in which hopper: ','int',False)
    hop[hopnum]=copy.deepcopy(hop[0])
    return hop

def womrdhop(hop):
    hopnum=0
    while (hopnum < 1) or (hopnum > HOPSIZE):
        hopnum=inputter('Read from which hopper: ','int',False)
    if (len(hop[hopnum].flux) == 0):
        print('Nothing in hopper {}'.format(hopnum))
        print('Active spectrum still {}'.format(hop[0].obname))
    else:
        hop[0]=copy.deepcopy(hop[hopnum])
    print('Object is {}'.format(hop[0].obname))
    return hop
    
def wombuf(hop):
    for i,_ in enumerate(hop):
        print('Hopper {}: {}'.format(i,hop[i].obname))
    return hop

def womrpl(hop):
    done=False
    while (not done):
        inputfile=input('Name of file to be read? ')
        inputfile=inputfile.strip()
        if (inputfile == ''):
            return
        try:
            data=np.loadtxt(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            done=True
    wave=data[:,0]
    flux=data[:,1]
    if (data.shape[1] > 2):
        var=data[:,2]
    else:
        var=np.ones(wave.shape)

    var=var**2
    print('Found {} bins.'.format(len(wave)))

    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].var=var.copy()
    hop[0].obname=inputfile
    hop[0].header=''
    return hop

def womwpl(hop):
    print('\nObject is {}\n'.format(hop[0].obname))
    spectxt=input('Enter the name for the output file: ')
    if (spectxt == ''):
        return hop
    spectxt=spectxt.strip()
    np.savetxt(spectxt,np.transpose([hop[0].wave,hop[0].flux,np.sqrt(hop[0].var)]))
    return hop

def womplot(hop):
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    if (hop[0].wave[0] < 0):
        plt.xlabel('Velocity')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    wshow()
    print('Change scale?')
    answer=yesno('n')
    if (answer == 'y'):
        xmin,xmax=plt.xlim()
        ymin,ymax=plt.ylim()
        done=False
        while (not done):
            print('Click corners of box to change plot scale')
            newlims=plt.ginput(2)
            xmin=newlims[0][0]
            ymin=newlims[0][1]
            xmax=newlims[1][0]
            ymax=newlims[1][1]
            plt.cla()
            plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
            plt.xlabel('Wavelength')
            if (hop[0].wave[0] < 0):
                plt.xlabel('Velocity')
            plt.ylabel('Flux')
            plt.title(hop[0].obname)
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            print('Is this OK?')
            loopanswer=yesno('y')
            if (loopanswer == 'y'):
                done=True
            
    return hop

def womapl(hop):
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    return hop


def waveparse():
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
    return waveb,waver

def wommkbb(hop):
    light_speed=2.99792458e10
    h_planck=6.6260755e-27
    k_boltzmann=1.380658e-16
    print('This routine will create a blackbody curve for a given temperature')
    print('with 1 Angstrom bins')
    print('\n')
    temp=inputter('Enter temperature in degrees Kelvin: ','float',False)
    if (temp <= 0):
        temp=1.0
    waveb,waver=waveparse()
    if (waver < waveb):
        waveb,waver=waver,waveb
    if (waver == waveb):
        waver=waver+1.
    wave=np.arange(waveb,waver+1)
    wavecm=wave/1.e8
    answer=inputter_single('Calculate B_nu(n) or B_lambda(l) ','nl')
    if (answer.lower()[0] == 'l'):
        flux=(2.0*h_planck*light_speed*light_speed/wavecm**5)/ \
              (np.exp((h_planck*light_speed)/(wavecm*k_boltzmann*temp))-1.0)
        flux=np.pi*flux/(1.e8)
        ext='flm'
    else:
        nu=light_speed/wavecm
        flux=(2.0*h_planck*nu**3/light_speed/light_speed)/ \
              (np.exp((h_planck*nu)/(k_boltzmann*temp))-1.0)
        flux=np.pi*flux*1.e11
        ext='fnu'
    spectxt='bb{}.{}'.format(temp,ext)
    plt.cla()
    plt.plot(wave,flux,drawstyle='steps-mid')
    plt.title(spectxt)
    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].obname=spectxt
    hop[0].var=np.ones(wave.shape)
    hop[0].header=''
    return hop

def womrdfits(hop):
    done=False
    while (not done):
        inputfile=input('Name of fits file to be read? (.fits added if necessary) ')
        inputfile=inputfile.strip()
        if (inputfile == ''):
            return
        if ('.fits' not in inputfile):
            inputfile=inputfile+'.fits'
        try:
            fitsdata=fits.open(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            done=True
    if (len(fitsdata) > 1):
        print('Multi-Extension File!  Bailing out!')
        return hop
    flux=fitsdata[0].data
    flux=flux.astype(float)
    crval1=fitsdata[0].header['CRVAL1']
    wdelt=fitsdata[0].header['CDELT1']
    if (wdelt < 0.000001):
        wdelt=fitsdata[0].header['CD1_1']
    objectname=fitsdata[0].header['OBJECT']
    if (not isinstance(objectname, str)):
        objectname=inputfile
    print('Object is {}'.format(objectname))
    wave=np.arange(1,len(flux)+1)*wdelt+crval1

    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].obname=objectname
    hop[0].var=np.ones(wave.shape)
    hop[0].header=fitsdata[0].header
    fitsdata.close()
    return hop

def wshow():
    wm=plt.get_current_fig_manager()
    blah=wm.window.attributes('-topmost',1)
    blah=wm.window.attributes('-topmost',0)
    return

def womwaverange(wave,flux,mode):
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

def womget_element(wavearr,wavelength):
    """get index of given wavelength"""
    index=(np.abs(wavearr-wavelength)).argmin()
    return index

def womcho(hop):
    """choose region of spectrum with wavelengths or mouse"""
    wave=hop[0].wave.copy()
    flux=hop[0].flux.copy()
    var=hop[0].var.copy()
    print('Current A/pix is {}'.format(wave[1]-wave[0]))
    wave,flux,mode=womwaverange(wave,flux,'none')
    indexblue=womget_element(hop[0].wave,wave[0])
    indexred=womget_element(hop[0].wave,wave[-1])
    var=var[indexblue:indexred+1].copy()
    print('Overplotting chosen section')
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.plot(wave,flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].var=var.copy()
    ##FIX header for fits
    return hop

def womashrebinselector(hop):
    hop=womrebindriver(hop,'ash')
    return hop

def womscipyrebinselector(hop):
    hop=womrebindriver(hop,'scipy')
    return hop
    

def womrebindriver(hop,mode):
    """input info to call womrebin functions"""
    print('Current A/pix is {}'.format(hop[0].wave[1]-hop[0].wave[0]))
    print('Wavelength rage: {} to {}'.format(hop[0].wave[0],hop[0].wave[-1]))
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    done = False
    while (not done):
        newdelt=inputter('Rebin to how many Angstroms per pixel? ', \
                         'float',False)
        waveb,waver=waveparse()
        if (waveb > waver):
            waveb,waver=waver,waveb
        if (newdelt > 0):
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
                    done=True
            else:
                print('Rebinning to {} A/pix, {} to {}'.format(newdelt,waveb,\
                                                       waver))
                print('Is this OK?')
                answer=yesno('y')
                if (answer == 'y'):
                    done=True
    
    nwave=np.arange(0,whole)*newdelt+waveb
    if (mode == 'ash'):
        nflux=womashrebin(hop[0].wave, hop[0].flux, nwave)
        nvar=womashrebin(hop[0].wave, hop[0].var, nwave)
    else:
        nflux=womscipyrebin(hop[0].wave, hop[0].flux, nwave)
        nvar=womscipyrebin(hop[0].wave, hop[0].var, nwave)
    print('Overplotting rebinned spectrum')
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.plot(nwave,nflux,drawstyle='steps-mid')
    
    #FIX header
    hop[0].wave=nwave.copy()
    hop[0].flux=nflux.copy()
    hop[0].var=nvar.copy()
    return hop

def womply(a,i,x):
    """poly evaluator for womashrebin
;
;             I+1
;  EVALUATES  SUM A(J) * X**(I+1-J)
;             J=1
;
;  I IS THUS THE ORDER OF THE POLY, WHICH HAS I+1 COEFFICIENTS,
;  AND A(I+1) IS THE CONSTANT TERM.
;

    x is not array-like

    """
    if (i < 0):
        polyval=0.
        print('Value of I out of bounds in womply')
        return polyval
    polyval=a[0]
    if (i == 0):
        return polyval
    for j in range(2,i+2):
        polyval=x*polyval+a[j-1]
    return polyval

def womrebin_tform(rx,xin,arctwo,darcone,in21,acc,rslope,in11,in12,arcone):
    """needed for ashrebin"""
    n=1
    x=xin
    rl=womply(arctwo,in21,rx)
    while True:
        l=womply(arcone,in11,x)
        dx=(rl-l)*rslope
        x=x+dx
        if (abs(dx) < acc):
            break
        if (n > 100):
            break
        n=n+1
        rslope=1./(womply(darcone,in12,x))
    return x

def womashrebin(wave,flux,nwave):
    """moo"""
    ACCBIN=0.0001
    arcone=np.zeros(2)
    arctwo=np.zeros(2)
    darcone=np.zeros(2)
    npix=len(wave)
    nrbin=len(nwave)
    rbin=np.zeros(nrbin)
    arc1=np.array([wave[0]-(wave[-1]-wave[0])/(npix-1), \
                   ((wave[-1]-wave[0])*npix)/(npix-1)])
    arc2=np.array([nwave[0]-(nwave[-1]-nwave[0])/(nrbin-1), \
                   ((nwave[-1]-nwave[0])*nrbin)/(nrbin-1)])
    inarc1=2
    inarc2=2
    in11=inarc1-1
    in12=inarc1-2
    in21=inarc2-1
    lstart=0
    acc=ACCBIN/npix
    rnpix=1./npix
    rnrbin=1./nrbin
    rx2=0.5*rnrbin
    for i in range(0,inarc1):
        arcone[inarc1-1-i]=arc1[i]
        darcone[inarc1-1-i]=i*arc1[i]
    for i in range(0, inarc2):
        arctwo[inarc2-1-i]=arc2[i]
    rslope=1./womply(darcone,in12,0.2)
    x1in=0.2
    x1=womrebin_tform(rx2,x1in,arctwo,darcone,in21,acc,rslope,in11,in12,arcone)
    x1=x1*npix
    dx=0
    nsgn=1
    if((womply(arctwo,in21,1)-arc2[1])*rslope < 0):
        nsgn=-1
    nstop=nrbin
    j1=np.round(x1)-1
    for k in range(0, nstop):
        rx2=rx2+rnrbin
        x2in=(x1+dx)*rnpix
        x2=womrebin_tform(rx2,x2in,arctwo,darcone,in21,acc,rslope,in11,in12,arcone)
        x2=x2*npix
        dx=x2-x1
        j2=np.round(x2)-1
        d=0
        if (lstart == 0):
            lstart=1
            m1=int(max([min([j1-1,npix-1]),1]))
            m2=int(max([min([j1,npix-1]),1]))
            m3=int(max([min([j1+1,npix-1]),1]))
            a=(flux[m1]+flux[m3])*0.5
            b=(a-flux[m1])*0.5
            c=(13./12.)*flux[m2]-a/12.0
            a=(a-flux[m2])/3.0
            y=x1-j1-1
            dd=nsgn*((((a*y)+b)*y+c)*y-b*0.25)+a*0.125+c*0.5

        m1=int(max([min([j2-1,npix-1]),1]))
        m2=int(max([min([j2,npix-1]),1]))
        m3=int(max([min([j2+1,npix-1]),1]))
        a=(flux[m1]+flux[m3])*0.5
        b=(a-flux[m1])*0.5
        c=1.0833333333333333*flux[m2]-a*0.0833333333333333
        a=(a-flux[m2])*0.333333333333333
        y=x2-j2-1
        d=d-dd
        dd=nsgn*((((a*y)+b)*y+c)*y-b*0.25)
        ddd=a*0.125+c*0.5
        d=d+dd-ddd
        dd=dd+ddd
      #  print(d,dd,ddd)
        for kk in range(int(j1), int(j2)+1*int(nsgn), int(nsgn)):
            d=d+flux[int(max([min([kk,npix-1]),1]))]
       # print(d)

        rbin[k]=d/np.abs(dx)
        x1=x2
        j1=j2
    nflux=rbin.copy()
    return nflux

def womscipyrebin(wave,flux,nwave):
    import scipy.interpolate as scint
    inter=scint.interp1d(wave,flux,kind='cubic',bounds_error=False,fill_value='extrapolate')
    nflux=inter(nwave)
    return nflux

def womstat(hop):
    """print statistics of region of spectrum"""
    import scipy.stats as st
    print('\nObject is {}'.format(hop[0].obname))
    print('\nEnter range for statistics: ')
    wave,flux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
    print('\nMean: {}'.format(np.mean(flux)))
    print('Variance: {}'.format(np.var(flux,ddof=1)))
    print('Std. Dev.: {}'.format(np.std(flux,ddof=1)))
    print('Mean Dev.: {}'.format(np.mean(np.abs(flux-np.mean(flux)))))
    print('S/N: {}'.format(np.mean(flux)/np.std(flux,ddof=1)))
    #this is IRAF definition for S/N
    print('Skewness: {}'.formt(st.skew(flux)))
    print('Kurtosis: {}'.format(st.kurtosis(flux)))
    print('Median: {}'.format(np.median(flux)))
    print('No. points: {}'.format(len(flux)))
    return hop
    
def womexam(hop):
    """print spectrum bin values"""
    done = False
    while (not done):
        answer='y'
        wave,flux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
        indexblue=womget_element(hop[0].wave,wave[0])
        indexred=womget_element(hop[0].wave,wave[-1])
        var=hop[0].var[indexblue:indexred+1]
        if (len(wave) > 30):
            print('This will print out {} values'.format(len(wave)))
            print('Do you really want to do that?')
            answer=yesno('n')
        if (answer == 'y'):
            done=True
    print('\nWave    Flux')
    for i,_ in enumerate(wave):
        print('{} {} {}'.format(wave[i],flux[i], var[i]))
    return hop
            
def womsmo(hop):
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    mode=inputter_single('(b)oxcar or (s)avitzky-Golay smoothing? (b/s) ', \
                         'bs')
    if (mode == 'b'):
        print('\nSmoothing with a boxcar\n')
        from astropy.convolution import convolve, Box1DKernel
        width=inputter('Enter the width of the boxcar: ','int',False)
        if (width > len(hop[0].wave)):
            width = 1
        box_k=Box1DKernel(width)
        sflux=convolve(hop[0].flux,box_k,boundary='extend')
    else:
        print('\nThe Savitzky-Golay filter smooths the data while conserving')
        print('flux and retaining the dynamic range of variations.  The ')
        print('routine suggests a width of 1-2 times the FWHM of the ')
        print('desired features, assuming you know what features you ')
        print('do desire.  This currently uses a polynomial of degree')
        print('2 to create the filter, so the width must be at least 3.')
        print('Good luck.\n')
        from scipy.signal import savgol_filter
        width=inputter('Enter the width for the filter: ','int',False)
        if (width < 3) or (width > len(hop[0].wave)):
            width=3
        sflux=savgol_filter(hop[0].flux,width,2)
    print('Overplotting smoothed spectrum')
    plt.plot(hop[0].wave,sflux,drawstyle='steps-mid')
    hop[0].flux=sflux.copy()
    return hop

def womplotlog(hop):
    """plot with log scale for x, y, or both"""
    
    print('Object is {}'.format(hop[0].obname))
    print('Spectrum runs from {} to {}.\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    mode=inputter_single('Do you want log scale on x, y, or (b)oth axes? (x/y/b) ','xyb')
    plt.cla()
    if (mode == 'x'):
        plt.semilogx(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    elif (mode == 'y'):
        plt.semilogy(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    else:
        plt.loglog(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    return hop
    
def womscale(hop):
    """scale a spectrum by multiplicative value"""
    factor=inputter('Enter multiplicative scale factor: ','float',False)
    hop[0].flux=hop[0].flux*factor
    return hop

def womscale2match(hop):
    """scale one spectrum to match another"""
    print('This will scale one hopper to match another')
    done = False
    while (not done):
        hopnum1=inputter('Enter first hopper: ','int',False)
        hopnum2=inputter('Enter second hopper: ','int',False)
        if (hopnum1 < 1) or (hopnum1 > HOPSIZE) or (hopnum2 < 1) or \
           (hopnum2 > HOPSIZE):
            print('Hopper numbers must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    if (hop[hopnum1].wave[0] != hop[hopnum2].wave[0]) or \
       (hop[hopnum1].wave[1] != hop[hopnum2].wave[1]) or \
       (hop[hopnum1].wave[-1] != hop[hopnum2].wave[-1]):
        print('Hoppers to not have the same wavelength scale!')
        return hop
    print('\nSpectra run from {} to {}'.format(hop[hopnum1].wave[0], \
                                             hop[hopnum1].wave[-1]))

    plt.cla()
    plt.plot(hop[hopnum1].wave,hop[hopnum1].flux,drawstyle='steps-mid', \
             color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[hopnum1].obname)
    plt.plot(hop[hopnum2].wave,hop[hopnum2].flux,drawstyle='steps-mid', \
             color='r')
    print('\nHopper A in black, hopper B in red')
    print('\nChoose region to compute average')
    wave,flux,mode=womwaverange(hop[hopnum1].wave,hop[hopnum1].flux,'none')
    indexblue=womget_element(hop[hopnum1].wave,wave[0])
    indexred=womget_element(hop[hopnum1].wave,wave[-1])
    avg1=np.mean(hop[hopnum1].flux[indexblue:indexred+1])
    avg2=np.mean(hop[hopnum2].flux[indexblue:indexred+1])
    print('\n Hopper A: {}'.format(avg1))
    print('\n Hopper B: {}'.format(avg2))
    choice=inputter_single('Scale to (A) or (B) ','ab')
    print(' ')
    done = False
    while (not done):
        hopnum3=inputter('Store scaled spec in which hopper? ','int',False)
        if (hopnum3 < 1) or (hopnum3 > HOPSIZE):
            print('Hopper numbers must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    if (choice == 'a'):
        hop[hopnum3].flux=hop[hopnum2].flux*(avg1/avg2).copy()
        print('Scale: {}'.format(avg1/avg2))
        hop[hopnum3].obname=hop[hopnum2].obname
        hop[hopnum3].header=hop[hopnum2].header
        hop[hopnum3].var=hop[hopnum2].var*(avg1/avg2).copy()
        lstring='File: {} scaled by {} to match {}'.format(hop[hopnum2], \
                                                           avg1/avg2, \
                                                           hop[hopnum1])
    else:
        hop[hopnum3].flux=hop[hopnum1].flux*(avg2/avg1).copy()
        print('Scale: {}'.format(avg1/avg2))
        hop[hopnum3].obname=hop[hopnum1].obname
        hop[hopnum3].header=hop[hopnum1].header
        hop[hopnum3].var=hop[hopnum1].var*(avg2/avg1).copy()
        lstring='File: {} scaled by {} to match {}'.format(hop[hopnum1].obname, \
                                                           avg2/avg1, \
                                                           hop[hopnum2].obname)
    hop[hopnum3].wave=hop[hopnum1].wave.copy()
    logging.debug(lstring)
    return hop


def womscalemany(hop):
    """scale many hoppers to one"""
    print('This will scale many hoppers to one')
    done = False
    while (not done):
        hopnum1=inputter('Enter fiducial hopper: ','int',False)
        if (hopnum1 < 1) or (hopnum1 > HOPSIZE):
            print('Hopper number must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    hoplist=[]
    hopnum2=0
    while (hopnum2 != 99):
        done = False
        while (not done):
            hopnum2=inputter('Enter another hopper: (99 to end) ','int',False)
            if (hopnum2 != 99):
                if (hopnum2 < 1) or (hopnum2 > HOPSIZE):
                    print('Hopper numbers must be in range 1-{}'.format(HOPSIZE))
                elif (len(hop[hopnum2].wave) == 0):
                    print('Nothing in that hopper')
                else:
                    hoplist.append(hopnum2)
                    done=True
                if (hopnum2 > 0) and (hopnum2 < 21) and \
                   (len(hop[hopnum2].wave) != 0):
                    if (hop[hopnum1].wave[0] != hop[hopnum2].wave[0])  \
                       or (hop[hopnum1].wave[1] != hop[hopnum2].wave[1]) \
                       or (hop[hopnum1].wave[-1] != hop[hopnum2].wave[-1]):
                        print('Hoppers to not have the same wavelength scale!')
                        return hop
            else:
                done=True
    print('\nSpectra run from {} to {}'.format(hop[hopnum1].wave[0], \
                                             hop[hopnum1].wave[-1]))

    plt.cla()
    plt.plot(hop[hopnum1].wave,hop[hopnum1].flux,drawstyle='steps-mid', \
             color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[hopnum1].obname)
    print('\nPlotting fiducial spectrum')
    print('\nChoose region to compute average')
    wave,flux,mode=womwaverange(hop[hopnum1].wave,hop[hopnum1].flux,'none')
    indexblue=womget_element(hop[hopnum1].wave,wave[0])
    indexred=womget_element(hop[hopnum1].wave,wave[-1])
    avg1=np.mean(hop[hopnum1].flux[indexblue:indexred+1])
    print('Fiducial spectrum mean in range: {}'.format(avg1))
    for i,_ in enumerate(hoplist):
        print('{} {}'.format(i,hoplist[i]))
        avgloop=np.mean(hop[hoplist[i]].flux[indexblue:indexred+1])
        print('Hopper {} mean in range: {}'.format(hoplist[i], avgloop))
        hop[hoplist[i]].flux=hop[hoplist[i]].flux*avg1/avgloop
        hop[hoplist[i]].var=hop[hoplist[i]].var*avg1/avgloop
        logging.debug('File: {} scaled by {} to match {}'.format(hop[hoplist[i]].obname, avg1/avgloop, hop[hopnum1].obname))
    return hop
                     
def womredshift(hop):
    """remove redshift from spectrum"""
    light_speed=2.99792458e5
    print('Current A/pix is {}'.format(hop[0].wave[1]-hop[0].wave[0]))
    print('\nWavelength range: {} to {}\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    answer=inputter_single('Remove redshift in (z) or (k)m/s? (z/k) ','zk')
    z=inputter('Enter the redshift (positive is away from us): ','float',False)
    if (answer == 'k'):
        z=np.sqrt((1.0+z/c)/(1.0-z/c)) - 1.0
    hop[0].wave=hop[0].wave/(1.0+z)
    print('\nNew wavelength range: {} to {}\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    print('Rebin spectrum?')
    rebin=yesno('n')
    if (rebin == 'y'):
        hop=womrebindriver(hop,'scipy')
    else:
        pass
        #FIX header
    return hop

def womblueshift(hop):
    """add redshift to spectrum"""
    light_speed=2.99792458e5
    print('Current A/pix is {}'.format(hop[0].wave[1]-hop[0].wave[0]))
    print('\nWavelength range: {} to {}\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    answer=inputter_single('Remove redshift in (z) or (k)m/s? (z/k) ','zk')
    z=inputter('Enter the desired blueshift (positive is away from us): ','float',False)
    if (answer == 'k'):
        z=np.sqrt((1.0+z/c)/(1.0-z/c)) - 1.0
    hop[0].wave=hop[0].wave*(1.0+z)
    print('\nNew wavelength range: {} to {}\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    print('Rebin spectrum?')
    rebin=yesno('n')
    if (rebin == 'y'):
        hop=womrebindriver(hop,'scipy')
    else:
        pass
        #FIX header
    return hop

def womms(hop):
    """arithmetic on spectra"""
    print('\nThis routine will perfrom arithmatic operations with two hoppers\n')
    print('Do you want to (a)dd, (s)ubtract, (m)ultiply, or (d)ivide spectra?\n')
    choice=inputter_single('(a/s/m/d) ','asmd')
    if (choice == 's'):
        print('Second spectrum will be subtracted from the first.')
    if (choice == 'd'):
        print('First spectrum will be divided by the first.')
    done = False
    while (not done):
        hopnum1=inputter('\nEnter the first hopper: ','int',False)
        print(' ')
        hopnum2=inputter('Enter the second hopper: ','int',False)
        if (hopnum1 < 1) or (hopnum1 > HOPSIZE) or (hopnum2 < 1) or \
           (hopnum2 > HOPSIZE):
            print('\nHopper must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    if (hop[hopnum1].wave[0] != hop[hopnum2].wave[0])  \
       or (hop[hopnum1].wave[1] != hop[hopnum2].wave[1]) \
       or (hop[hopnum1].wave[-1] != hop[hopnum2].wave[-1]):
        print('Hoppers to not have the same wavelength scale!')
        return hop
    if (choice == 'a'):
        hop[0].flux=hop[hopnum1].flux + hop[hopnum2].flux
        hop[0].var=hop[hopnum1].var + hop[hopnum2].var
    elif (choice == 's'):
        hop[0].flux=hop[hopnum1].flux - hop[hopnum2].flux
        hop[0].var=hop[hopnum1].var + hop[hopnum2].var       
    elif (choice == 'm'):
        hop[0].flux=hop[hopnum1].flux * hop[hopnum2].flux
        hop[0].var=(hop[hopnum2].flux)**2*hop[hopnum1].var + \
                    (hop[hopnum1].flux)**2*hop[hopnum2].var
    elif (choice == 'd'):
        hop[0].flux=hop[hopnum1].flux / hop[hopnum2].flux
        hop[0].var=(hop[hopnum2].flux)**2*hop[hopnum1].var + \
                    (hop[hopnum1].flux)**2*hop[hopnum2].var
    hop[0].wave=hop[hopnum1].wave
    hop[0].obname=hop[hopnum1].obname
    hop[0].header=hop[hopnum1].header
    return hop

def womfnuflam(hop):
    """convert flux densities"""
    light_speed=2.99792458e18
    print('Converting flux units.  Is the input spectrum in f-(n)u,')
    print('f-(l)ambda, (a)b magnitudes, or (s)t magnitudes?')
    flun=inputter_single('(n/l/a/s)','nlas')
    fluxmean=np.mean(hop[0].flux)
    scale=False
    #convert input to f-lambda, remove any arbitrary scaling
    if (flun == 'n'):
        if (fluxmean > 1.e-6):
            hop[0].flux=hop[0].flux/1.e26
            scale=True
        hop[0].flux=hop[0].flux*light_speed/hop[0].wave/hop[0].wave
    elif (flun == 'a'):
        hop[0].flux=10**(-0.4*hop[0].flux - 19.44)
        hop[0].flux=hop[0].flux*light_speed/hop[0].wave/hop[0].wave
    elif (flun == 's'):
        hop[0].flux=10**(-0.4*hop[0].flux - 8.44)
    else:
        if (fluxmean > 1.e-6):
            hop[0].flux=hop[0].flux/1.e15
            scale=True

    print('\nConvert to f-(n)u, f-(l)ambda, (a)b magnitudes,')
    print('or (s)t magnitudes?')
    flout=inputter_single('(n/l/a/s)','nlas')
    if (flout == 'n'):
        hop[0].flux=hop[0].flux*hop[0].wave*hop[0].wave/light_speed
        if (scale):
            hop[0].flux=hop[0].flux*1.e26
    elif (flout == 'a'):
        hop[0].flux=hop[0].flux*hop[0].wave*hop[0].wave/light_speed
        hop[0].flux=hop[0].flux*10**19.44
        negloc=np.where(hop[0].flux <= 0)
        hop[0].flux[negloc]=1.0
        hop[0].flux=-2.5*np.log10(hop[0].flux)
    elif (flout == 's'):
        negloc=np.where(hop[0].flux <= 0)
        hop[0].flux[negloc]=3.63e-19
        hop[0].flux=-2.5*np.log10(hop[0].flux) - 21.10
    else:
        if (scale):
            hop[0].flux=hop[0].flux*1.e15
    return hop

def womlup(hop):
    """scales to asinh luptitudes"""
    print('\nTaking asinh of flux data\n')

    hop[0].flux=np.arcsinh(hop[0].flux)

    return hop

def wombluen(hop):
    """bluen spectrum with scattering law"""
    print('\nThis routine will bluen a spectrum with a power law')
    print('of lambda^{-a}.  You will supply the value of a.\n')

    factor=inputter('Enter exponential factor: ','float',False)

    wavefac=hop[0].wave**(-1.0*factor)
    wavefac=wavefac/wavefac[-1]
    hop[0].flux=hop[0].flux*wavefac

    return hop

def womrelvel(hop):
    """relativistic velocity calculation"""
    light_speed=2.99792458e5
    print('Relativistic velocity calculation\n')
    
    lambda1=inputter('Observed wavelength ','float',False)
    print(' ')
    lambda0=inputter('Rest wavelength ','float',False)
    if (lambda0 <= 0):
        print('Invalid rest wavelength')
        return hop
    z=(lambda1-lambda0)/lambda0
    sq=(z+1.0)**2
    vel=z*light_speed
    relvel=((sq-1.0)/(sq+1.0))*light_speed

    print('\nz: {}'.format(z))
    print('Velocity: {}'.format(vel))
    print('Relativistic velocity {}'.format(relvel))
    
    return hop

def womxshift(hop):
    """linearly shift wavelength by arbitrary amount"""
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    wshow()

    print('Routine to linearly shift wavelength scale\n')

    shift=inputter('Enter wavelength shift in Angstroms: ','float',False)

    hop[0].wave=hop[0].wave+shift

    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')

    logging.debug('File {} wavelength scale shifted by {} A'.format\
                  (hop[0].obname,shift))

    #FIX header
    return hop

def womyshift(hop):
    """linearly shift flux by arbitrary amount"""
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    wshow()

    print('Routine to linearly shift flux scale\n')

    shift=inputter('Enter flux shift: ','float',False)

    hop[0].flux=hop[0].flux+shift

    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')

    logging.debug('File {} flux scale shifted by {} A'.format\
                  (hop[0].obname,shift))

    #FIX header
    return hop

def womzap(hop):
    """routine to zap outliers (CRs)"""
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    wshow()
    print('\nRoutine to zap outliers\n')
    done = False
    while (not done):
        nsig=inputter('Enter zapping threshold in sigmas (0=replace all with median): ','float',False)
        print(' ')
        boxsize=inputter('Enter box size for computing statistics: (odd integer < 45) ','int',False)
        if (boxsize < 3) or (boxsize > 45) or (nsig < 0):
            print('Invalid boxsize or sigma')
        else:
            done = True
    if (boxsize % 2 == 0):
        boxsize=boxsize+1
    half=int(boxsize/2.)
    newflux=hop[0].flux.copy()
    if (nsig > 0):
        mode=inputter_single('Use inter(q)uartile or (m)edian variance (q/m)? ','qm')
        if (mode == 'm'):
            for i in range(half,len(newflux)-half):
                sample=newflux[i-half:i+half+1]
                medval=np.median(sample)
                varpix=(sample[half]-medval)**2
                sum=-1.0*varpix
                for k in range(0,boxsize):
                    diff=sample[k]-medval
                    sum=sum+diff*diff
                sum=sum/(boxsize-1)
                if (varpix > nsig*nsig*sum):
                    newflux[i]=medval
        if (mode == 'q'):
            for i in range(half,len(newflux)-half):
                sample=newflux[i-half:i+half+1]
                medval=np.median(sample)
                q25=np.percentile(sample,25.)
                q75=np.percentile(sample,75.)
                qsig=0.7414272*(q75-q25)
                diff=abs(sample[half]-medval)
                if (diff > nsig*qsig):
                    newflux[i]=medval
    if (nsig == 0):
        from scipy.ndimage.filters import median_filter
        newflux=median_filter(newflux,boxsize)
    plt.plot(hop[0].wave,newflux,drawstyle='steps-mid')
    print('Does this zapping look good?')
    good=yesno('y')
    if (good == 'y'):
        hop[0].flux=newflux.copy()
        logging.debug('File {} zapped with sigma {} and boxsize {}'.format\
                      (hop[0].obname,nsig,boxsize))
    else:
        print('OK, active spectrum unchanged')
    #FIX var
    return hop
                                                                           
def womcomselector(hop):
    """select combine 2 equally"""
    hop=womcom(hop,2,False)
    return hop

def womcmwselector(hop):
    """select combine 2 with weights"""
    hop=womcom(hop,2,True)
    return hop

def womcom3selector(hop):
    """select combine 3 equally"""
    hop=womcom(hop,3,False)
    return hop

def womcommanyselector(hop):
    """select combine many equally"""
    hop=womcom(hop,HOPSIZE+1,False)
    return hop

def womcommanywselector(hop):
    """select combine many with weights"""
    hop=womcom(hop,HOPSIZE+1,True)
    return hop

def womcom(hop,num,weights):
    """combine spectra, with or without weights"""
    hopnum=[0]
    weight=[0]
    while (hopnum[0] < 1) or (hopnum[0] > HOPSIZE):
        hopnum[0]=inputter('Enter first hopper: ','int',False)
    if (weights):
        weight[0]=inputter('Enter weight for first hopper: ','float',False)
    for i in range(1,num):
        hoploop=0
        weightloop=0
        if (num > 3):
            print('Use hopper 99 to end')
        while (hoploop < 1) or (hoploop > HOPSIZE):
            hoploop=inputter('Enter next hopper: ','int',False)
            if (hoploop == 99):
                break
        if (hoploop == 99):
            break
        if (hop[hopnum[0]].wave[0] != hop[hoploop].wave[0]) or \
       (hop[hopnum[0]].wave[1] != hop[hoploop].wave[1]) or \
       (hop[hopnum[0]].wave[-1] != hop[hoploop].wave[-1]):
            print('Hoppers to not have the same wavelength scale!')
            return hop
        hopnum.append(hoploop)
        if (weights):
            weightloop=inputter('Enter next weight: ','float',False)
        weight.append(weightloop)
    if (weights) and (abs(sum(weight)-1.0) > 0.00001):
        print('Weights do not add to 1.0')
        return hop
    if (not weights):
        weight=[1./len(hopnum)]*len(hopnum)
    newflux=np.zeros(len(hop[hopnum[0]].flux))
    logging.debug('Combining spectra:')
    
    for i in range(0,len(hopnum)):
        newflux=newflux+hop[hopnum[i]].flux*weight[i]
        logging.debug('Combining {} with weight {}'.format(hop[hopnum[i]].obname,\
                                                           weight[i]))
    hopout=0
    while (hopout < 1) or (hopout > HOPSIZE):
        hopout=inputter('Enter hopper to store combined spectrum: ','int',False)
    hop[hopout].wave=hop[hopnum[0]].wave.copy()
    hop[hopout].flux=newflux.copy()
    hop[hopout].obname=hop[hopnum[0]].obname
    hop[hopout].header=hop[hopnum[0]].header
    hop[hopout].var=hop[hopnum[0]].var.copy()
    plt.cla()
    plt.plot(hop[hopout].wave,hop[hopout].flux,drawstyle='steps-mid',color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[hopout].obname)
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    print('\nPlotting combined spectrum in black, components in color\n')
    for i in range(0,len(hopnum)):
        plt.plot(hop[hopnum[i]].wave,hop[hopnum[i]].flux,drawstyle='steps-mid')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    
    #FIX var
    return hop

def womjoin(hop):
    """join two spectra that abut in wavelength"""
    print('\nThis will join two hoppers (with no overlap)\n')
    hopnum1=0
    hopnum2=0
    while (hopnum1 < 1) or (hopnum1 > HOPSIZE):
        hopnum1=inputter('Enter first hopper: ','int',False)
    while (hopnum2 < 1) or (hopnum2 > HOPSIZE):
        hopnum2=inputter('Enter first hopper: ','int',False)
    if (hop[hopnum1].wave[0] > hop[hopnum2].wave[0]):
        hopnum1,hopnum2=hopnum2,hopnum1
    wdelt1=hop[hopnum1].wave[1]-hop[hopnum1].wave[0]
    wdelt2=hop[hopnum2].wave[1]-hop[hopnum2].wave[0]
    # check if wavelength dispersion same
    if (abs(wdelt1 -wdelt2) > 0.00001):
        print('Spectra do not have same Angstrom/pixel')
        print('Blue side: {}'.format(wdelt1))
        print('Red side: {}'.format(wdelt2))
        return hop
    #check if spectra abut
    if (abs(hop[hopnum2].wave[0] - (hop[hopnum1].wave[-1] + wdelt1)) \
        > 0.00001):
        print('\nSpectra do not abut\n')
        print('Red end of blue: {}'.format(hop[hopnum1].wave[-1]))
        print('Blue end of red: {}\n'.format(hop[hopnum2].wave[0]))
        return hop
    print('Joining from {} to {}'.format(hop[hopnum1].wave[-1], \
                                         hop[hopnum2].wave[0]))
    hopout=0
    while (hopout < 1) or (hopout > HOPSIZE):
        hopout=inputter('Enter hopper to store combined spectrum: ','int',False)
    hop[hopout].wave=np.concatenate([hop[hopnum1].wave,hop[hopnum2].wave])
    hop[hopout].flux=np.concatenate([hop[hopnum1].flux,hop[hopnum2].flux])
    hop[hopout].var=np.concatenate([hop[hopnum1].var,hop[hopnum2].var])
    hop[hopout].obname=hop[hopnum1].obname
    hop[hopout].header=hop[hopnum1].header
    logging.debug('Files: {} and {} joined from {} to {}'.format(hop[hopnum1].obname, hop[hopnum2].obname, hop[hopnum1].wave[-1],hop[hopnum2].wave[0]))
    #FIX header
    return hop

def womreddening(hop):
    """redden or deredden with various reddening laws"""
    import extinction
    r_v=3.1
    print('Redden or deredden a spectrum')
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid',color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    flux=hop[0].flux.copy()
    action=inputter_single('(r)edden or (d)eredden the spectrum? (r/d) ', 'rd')
    print(' ')
    type=inputter_single('Do you want to enter the (c)olor excess, or (v)isual extinction? ','cv')
    print(' ')
    if (type == 'v'):
        av=inputter('Enter A_V in magnitudes: ','float',False)
    else:
        ebv=inputter('Enter E(B-V) in magnitudes: ','float',False)
        av=r_v*ebv
    print(' ')
    print('Do you want to use: ')
    print('(c)ardelli, Clayton, Mathis 1989')
    print("(o)'donnell 1994")
    print('(f)itzpatrick 1999\n')
    
    method=inputter_single('(c/o/f) ','cof')
    if (action == 'r'):
        if (method == 'c'):
            newflux=extinction.apply(extinction.ccm89(hop[0].wave,av,r_v),flux)
        elif (method == 'o'):
            newflux=extinction.apply(extinction.odonnell94(hop[0].wave,av,r_v),flux)
        else:
            newflux=extinction.apply(extinction.fitzpatrick99(hop[0].wave,av,r_v),flux)
    else:
        if (method == 'c'):
            ext=extinction.ccm89(hop[0].wave,av,r_v)
        elif (method == 'o'):
            ext=extinction.odonnell94(hop[0].wave,av,r_v)
        else:
            ext=extinction.fitzpatrick99(hop[0].wave,av,r_v)
        newflux=flux*10**(0.4*ext)
    plt.plot(hop[0].wave,newflux,drawstyle='steps-mid',color='r')
    print('\nOriginal spectrum in black, red/dered in red\n')
    print('Is this OK?\n')
    answer=yesno('y')
    if (answer == 'y'):
        hop[0].flux=newflux.copy()
        print('\nActive spectrum now changed')
    else:
        print('\nSorry to disappoint you, active spectrum unchanged')
    return hop

def womheader(hop):
    """print header to screen"""
    if (len(hop[0].header)) < 2:
        print('No header!\n')
        return hop
    for i,_ in enumerate(hop[0].header):
        print(str(hop[0].header.cards[i]))
        if (i > 0) and (i % 19 == 0):
            print('Hit any key to continue')
            any=getch()
    return hop

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

def womrmsfits(hop):
    """read multispec fits files"""
    done=False
    while (not done):
        inputfile=input('Name of fits file to be read? (.ms.fits added if necessary) ')
        inputfile=inputfile.strip()
        if (inputfile == ''):
            return
        if ('.ms.fits' not in inputfile):
            inputfile=inputfile+'.ms.fits'
        try:
            fitsdata=fits.open(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            done=True
    if (len(fitsdata) > 1):
        print('Multi-Extension File!  Bailing out!')
        return hop
    flux=fitsdata[0].data
    flux=flux.astype(float)
    fhdr=fitsdata[0].header
    objectname=fhdr['OBJECT']
    if (not isinstance(objectname, str)):
        objectname=inputfile
    print('Object is {}'.format(objectname))
    #FIX for aps with different wavelength solutions, you need to look
    # at other keywords, but I don't have an example file right now
    # this is a long-term fix
    num_apertures=flux.shape[1]
    print('\nThere are {} apertures in this file.'.format(num_apertures))
    ap_choice=0
    if (num_apertures > 1):
        while (ap_choice < 1) or (ap_choice > num_apertures):
            ap_choice=inputter('Which aperture do you want? ','int',False)
    else:
        ap_choice=1
    ap_choice=ap_choice - 1
    wave=getmswave(fhdr,ap_choice)
    
    num_bands=flux.shape[0]
    print('\nThere are {} bands in this aperture.'.format(num_bands))
    # clean up keywords
    delkeylist=['WAT0_001', 'WAT1_001', 'WAT2_001', 'WAT0_002', \
                'WAT1_002', 'WAT2_002', 'WAT3_001', 'WAT2_003', \
                'CTYPE1', 'CTYPE2', 'CTYPE3', 'CD1_1', 'CD2_2', \
                'CD3_3', 'LTM1_1', 'LTM2_2', 'LTM3_3', 'WCSDIM']
    for k in delkeylist:
        try:
            fhdr.remove(k)
        except KeyError:
            pass
    fhdr.set('CRPIX1', 1)
    fhdr.set('CRVAL1', wave[0])
    fhdr.set('CDELT1', wdelt)
    fhdr.set('CTYPE1', 'LINEAR')
    if (num_bands > 2):
        var=flux[3,ap_choice,:]**2
    for i in range(0,num_bands):
        hopnum=0
        while (hopnum < 1) or (hopnum > HOPSIZE):
            hopnum=inputter('Store band {} in which hopper? '.format(i+1),'int',False)
        hop[hopnum].wave=wave.copy()
        hop[hopnum].flux=flux[i,ap_choice,:].copy()
        hop[hopnum].obname=objectname
        # OK, strictly speaking this var only applies to the optimal extraction
        # in band 1 (0 in python index).  We still have access to 'var' in
        # the hopper we put it in.
        hop[hopnum].var=var.copy()
        hop[hopnum].header=fhdr
    fitsdata.close()
    return hop

def womwfits(hop):
    """write fits file"""
    print('Object is {}'.format(hop[0].obname))
    fileout=inputter('Enter the name for the output file: ' + \
                     '(.fits will be added,if necessary)','string',False)
    outhdu=fits.PrimaryHDU(hop[0].flux)
    hdul=fits.HDUList([outhdu])
#    hdul[0].header.set('EXTEND','F')
    hdul[0].header.set('OBJECT',hop[0].obname)
    if ('.fits' not in fileout):
            fileout=fileout+'.fits'
    if (len(hop[0].header) > 1):
        hdul[0].header=hop[0].header.copy()
  #  hdul[0].header['SIMPLE']=('T','Written by Wombat')
    hdul[0].header.set('CRPIX1',1)
    hdul[0].header.set('CRVAL1',hop[0].wave[0])
    hdul[0].header.set('CDELT1',hop[0].wave[1]-hop[0].wave[0])
    hdul[0].header.set('CTYPE1','LINEAR')
    hdul[0].header.set('COMMENT','Written by Wombat')
    hdul.writeto(fileout,overwrite=True)
    hdul.close()
    return hop

def filtermag(wave,flux,filtw,filtt,zp):
    """calculate magnitude in given filter"""
    if (wave[0] > filtw[-1]) or (wave[-1] < filtw[0]):
        print('Filter and spectrum do not overlap')
        print('Bailing out and returning mag of 999')
        mag = 999
        flux = 999
        coverage = 0
        efflambda = 0
        totflux=0
        return mag, flux, coverage, efflambda, totflux
    blueover=np.where(wave < filtw[0])
    if (len(blueover[0]) < 1):
        blueindex=0
    else:
        blueindex=blueover[0][-1] + 1
    redover=np.where(wave > filtw[-1])
    if (len(redover[0]) < 1):
        redindex=len(wave)-1
    else:
        redindex=redover[0][0]-1
    wavematch=wave[blueindex:redindex+1]
    fluxmatch=flux[blueindex:redindex+1]
    spline=splrep(filtw, filtt, k=3)
    splineresult=splev(wavematch,spline)
    flux,totflux,efflambda,coverage=filterconv(wavematch,fluxmatch,splineresult,filtw,filtt)
    mag=-2.5*np.log10(flux/zp)
    return mag, flux, coverage, efflambda, totflux

def filterconv(wavematch,fluxmatch,splineresult,filtw,filtt):
    """convolve filter, find eff lambda"""
    wavematchplusone=np.roll(wavematch,-1)
    wavematchminusone=np.roll(wavematch,1)
    wavematchplusone[-1]=wavematchplusone[-2] + wavematchplusone[-2] - \
                          wavematchplusone[-3]
    wavematchminusone[0]=wavematchminusone[1]-(wavematchminusone[2] - \
                          wavematchminusone[1])
    totflux=np.sum(fluxmatch*splineresult*(wavematchplusone-wavematchminusone)/2.0)
    normalization=np.sum(splineresult*(wavematchplusone-wavematchminusone)/2.0)
    flux=totflux/normalization
    totwave=np.sum(wavematch*fluxmatch*splineresult*(wavematchplusone-wavematchminusone)/2.0)
    efflambda=totwave/totflux
    filtwmatchplusone=np.roll(filtw,-1)
    filtwmatchminusone=np.roll(filtw,1)
    filtwmatchplusone[-1]=filtwmatchplusone[-2] + filtwmatchplusone[-2] - \
                          filtwmatchplusone[-3]
    filtwmatchminusone[0]=filtwmatchminusone[1]-(filtwmatchminusone[2] - \
                          filtwmatchminusone[1])
    com=np.sum(filtt*(filtwmatchplusone-filtwmatchminusone)/2.0)
    coverage=normalization/com
    return flux,totflux,efflambda,coverage


    

def womfilters(hop):
    """calculate photometric values from spectra"""
    print('NOTE:  The routine expects an f_lambda spectrum')
    print('       I will try to guess if the spectrum')
    print('       has been scaled by 1E15')
    print(' ')
    print('       Check this before believing fluxes')
    print(' ')
    flux=hop[0].flux.copy()
    if (np.mean(flux) > 0.00001):
        flux = flux *1.e-15

    filtwave=np.zeros((141,10))
    filttran=np.zeros((141,10))
    filtwave[0:21, 1] = [3600.00, 3700.00, 3800.00, 3900.00, \
                         4000.00, 4100.00, 4200.00, 4300.00, \
                         4400.00, 4500.00, 4600.00, 4700.00, \
                         4800.00, 4900.00, 5000.00, 5100.00, \
                         5200.00, 5300.00, 5400.00, 5500.00, \
                         5600.00] 
    filttran[0:21, 1] = [0.00000, 0.03000, 0.13400, 0.56700, \
                         0.92000, 0.97800, 1.00000, 0.97800, \
                         0.93500, 0.85300, 0.74000, 0.64000, \
                         0.53600, 0.42400, 0.32500, 0.23500, \
                         0.15000, 0.09500, 0.04300, 0.00900, \
                         0.00000]
    filtwave[0:24, 2] = [4700.00, 4800.00, 4900.00, 5000.00, \
                      5100.00, 5200.00, 5300.00, 5400.00, \
                      5500.00, 5600.00, 5700.00, 5800.00, \
                      5900.00, 6000.00, 6100.00, 6200.00, \
                      6300.00, 6400.00, 6500.00, 6600.00, \
                      6700.00, 6800.00, 6900.00, 7000.00] 
    filttran[0:24:, 2] = [0.00000, 0.03000, 0.16300, 0.45800, \
                      0.78000, 0.96700, 1.00000, 0.97300, \
                      0.89800, 0.79200, 0.68400, 0.57400, \
                      0.46100, 0.35900, 0.27000, 0.19700, \
                      0.13500, 0.08100, 0.04500, 0.02500, \
                      0.01700, 0.01300, 0.00900, 0.00000]
    filtwave[0:24, 3] = [5500.00, 5600.00, 5700.00, 5800.00, \
                      5900.00, 6000.00, 6100.00, 6200.00, \
                      6300.00, 6400.00, 6500.00, 6600.00, \
                      6700.00, 6800.00, 6900.00, 7000.00, \
                      7100.00, 7200.00, 7300.00, 7400.00, \
                      7500.00, 8000.00, 8500.00, 9000.00]
    filttran[0:24:, 3] = [0.00000, 0.23000, 0.74000, 0.91000, \
                      0.98000, 1.00000, 0.98000, 0.96000, \
                      0.93000, 0.90000, 0.86000, 0.81000, \
                      0.78000, 0.72000, 0.67000, 0.61000, \
                      0.56000, 0.51000, 0.46000, 0.40000, \
                      0.35000, 0.14000, 0.03000, 0.00000]
    filtwave[0:23, 4] = [7000.00, 7100.00, 7200.00, 7300.00, \
                         7400.00, 7500.00, 7600.00, 7700.00, \
                         7800.00, 7900.00, 8000.00, 8100.00, \
                         8200.00, 8300.00, 8400.00, 8500.00, \
                         8600.00, 8700.00, 8800.00, 8900.00, \
                         9000.00, 9100.00, 9200.00]
    filttran[0:23, 4] = [0.00000, 0.02400, 0.23200, 0.55500, \
                         0.78500, 0.91000, 0.96500, 0.98500, \
                         0.99000, 0.99500, 1.00000, 1.00000, \
                         0.99000, 0.98000, 0.95000, 0.91000, \
                         0.86000, 0.75000, 0.56000, 0.33000, \
                         0.15000, 0.03000, 0.00000]

    filtwave[0:24, 0] = [3050.00, 3100.00, 3150.00, 3200.00, \
                      3250.00, 3300.00, 3350.00, 3400.00, \
                      3450.00, 3500.00, 3550.00, 3600.00, \
                      3650.00, 3700.00, 3750.00, 3800.00, \
                      3850.00, 3900.00, 3950.00, 4000.00, \
                      4050.00, 4100.00, 4150.00, 4200.00]
    filttran[0:24, 0] = [0.00000, 0.02000, 0.07700, 0.13500, \
                      0.20400, 0.28200, 0.38500, 0.49300, \
                      0.60000, 0.70500, 0.82000, 0.90000, \
                      0.95900, 0.99300, 1.00000, 0.97500, \
                      0.85000, 0.64500, 0.40000, 0.22300, \
                      0.12500, 0.05700, 0.00500, 0.00000]
    
    filtwave[0:47,5]=[2980., 3005., 3030., 3055., 3080., 3105., 3130.,  \
                      3155., 3180.,3205., 3230., 3255., 3280., 3305., \
                      3330., 3355., 3380., 3405., 3430., 3455., 3480., \
                      3505., 3530., 3555., 3580., 3605., 3630., 3655., \
                      3680., 3705., 3730., 3755., 3780., 3805., 3830., \
                      3855., 3880., 3905., 3930., 3955., 3980., 4005., \
                      4030., 4055., 4080., 4105., 4130.]
    filttran[0:47,5]=[0.    , 0.0014, 0.0071, 0.0127, 0.0198, 0.0314, \
                      0.0464, 0.0629, 0.0794, 0.0949, 0.1093, 0.1229, \
                      0.1352, 0.1458, 0.1545, 0.1617, 0.1679, 0.1737, \
                      0.1786, 0.1819, 0.1842, 0.186 , 0.187 , 0.1868, \
                      0.1862, 0.1858, 0.1853, 0.1841, 0.1812, 0.1754, \
                      0.1669, 0.1558, 0.1419, 0.1247, 0.1054, 0.0851, \
                      0.0634, 0.0405, 0.0216, 0.011 , 0.0062, 0.0032, \
                      0.0015, 0.0008, 0.0006, 0.0003, 0.    ]
    filtwave[0:89,6]=[3630., 3655., 3680., 3705., 3730., 3755., 3780., \
                      3805., 3830., 3855., 3880., 3905., 3930., 3955., \
                      3980., 4005., 4030., 4055., 4080., 4105., 4130., \
                      4155., 4180., 4205., 4230., 4255., 4280., 4305., \
                      4330., 4355., 4380., 4405., 4430., 4455., 4480., \
                      4505., 4530., 4555., 4580., 4605., 4630., 4655., \
                      4680., 4705., 4730., 4755., 4780., 4805., 4830., \
                      4855., 4880., 4905., 4930., 4955., 4980., 5005., \
                      5030., 5055., 5080., 5105., 5130., 5155., 5180., \
                      5205., 5230., 5255., 5280., 5305., 5330., 5355., \
                      5380., 5405., 5430., 5455., 5480., 5505., 5530., \
                      5555., 5580., 5605., 5630., 5655., 5680., 5705., \
                      5730., 5755., 5780., 5805., 5830.]
    filttran[0:89,6]=[0.000e+00, 5.000e-04, 1.300e-03, 2.200e-03, \
                      3.000e-03, 3.900e-03, 5.500e-03, 8.700e-03, \
                      1.620e-02, 3.010e-02, 5.000e-02, 7.450e-02, \
                      1.024e-01, 1.324e-01, 1.629e-01, 1.924e-01, \
                      2.191e-01, 2.419e-01,2.609e-01, 2.767e-01, \
                      2.899e-01, 3.010e-01, 3.105e-01, 3.186e-01, \
                      3.258e-01, 3.324e-01, 3.385e-01, 3.442e-01, \
                      3.496e-01, 3.548e-01, 3.596e-01, 3.640e-01, \
                      3.678e-01, 3.709e-01, 3.736e-01, 3.763e-01, \
                      3.792e-01, 3.827e-01, 3.863e-01, 3.899e-01, \
                      3.931e-01, 3.955e-01, 3.973e-01, 3.986e-01, \
                      3.997e-01, 4.008e-01, 4.019e-01, 4.030e-01, \
                      4.043e-01, 4.057e-01, 4.073e-01, 4.091e-01, \
                      4.110e-01, 4.129e-01, 4.147e-01, 4.165e-01, \
                      4.181e-01, 4.194e-01, 4.201e-01, 4.201e-01, \
                      4.191e-01, 4.169e-01, 4.147e-01, 4.115e-01, \
                      3.988e-01, 3.684e-01, 3.233e-01, 2.690e-01, \
                      2.112e-01, 1.550e-01, 1.043e-01, 6.270e-02, \
                      3.370e-02, 1.900e-02, 1.280e-02, 8.700e-03, \
                      5.700e-03, 3.700e-03, 2.400e-03, 1.700e-03, \
                      1.400e-03, 1.200e-03, 1.000e-03, 9.000e-04, \
                      7.000e-04, 5.000e-04, 3.000e-04, 1.000e-04, 0.000e+00]
    filtwave[0:75,7]=[5380., 5405., 5430., 5455., 5480., 5505., 5530., \
                      5555., 5580., 5605., 5630., 5655., 5680., 5705., \
                      5730., 5755., 5780., 5805., 5830., 5855., 5880., \
                      5905., 5930., 5955., 5980., 6005., 6030., 6055., \
                      6080., 6105., 6130., 6155., 6180., 6205., 6230., \
                      6255., 6280., 6305., 6330., 6355., 6380., 6405., \
                      6430., 6455., 6480., 6505., 6530., 6555., 6580., \
                      6605., 6630., 6655., 6680., 6705., 6730., 6755., \
                      6780., 6805., 6830., 6855., 6880., 6905., 6930., \
                      6955., 6980., 7005., 7030., 7055., 7080., 7105., \
                      7130., 7155., 7180., 7205., 7230.]
    filttran[0:75,7]=[0.000e+00, 1.600e-03, 1.130e-02, 2.970e-02, \
                      5.680e-02, 9.230e-02, 1.356e-01, 1.856e-01, \
                      2.390e-01, 2.917e-01, 3.395e-01, 3.794e-01, \
                      4.116e-01, 4.371e-01, 4.570e-01, 4.723e-01, \
                      4.839e-01, 4.925e-01, 4.990e-01, 5.040e-01, \
                      5.080e-01, 5.112e-01, 5.141e-01, 5.169e-01, \
                      5.194e-01, 5.213e-01, 5.222e-01, 5.220e-01, \
                      5.212e-01, 5.202e-01, 5.197e-01, 5.202e-01, \
                      5.215e-01, 5.233e-01, 5.254e-01, 5.275e-01, \
                      5.294e-01, 5.310e-01, 5.319e-01, 5.320e-01, \
                      5.316e-01, 5.310e-01, 5.305e-01, 5.302e-01, \
                      5.299e-01, 5.290e-01, 5.271e-01, 5.241e-01, \
                      5.211e-01, 5.176e-01, 5.057e-01, 4.775e-01, \
                      4.341e-01, 3.792e-01, 3.162e-01, 2.488e-01, \
                      1.824e-01, 1.225e-01, 7.470e-02, 4.300e-02, \
                      2.470e-02, 1.550e-02, 1.120e-02, 8.300e-03, \
                      5.900e-03, 4.100e-03, 2.900e-03, 2.100e-03, \
                      1.600e-03, 1.300e-03, 1.000e-03, 8.000e-04, \
                      5.000e-04, 2.000e-04, 0.000e+00]
    filtwave[0:89,8]=[6430., 6455., 6480., 6505., 6530., 6555., 6580., \
                      6605., 6630., 6655., 6680., 6705., 6730., 6755., \
                      6780., 6805., 6830., 6855., 6880., 6905., 6930., \
                      6955., 6980., 7005., 7030., 7055., 7080., 7105., \
                      7130., 7155., 7180., 7205., 7230., 7255., 7280., \
                      7305., 7330., 7355., 7380., 7405., 7430., 7455., \
                      7480., 7505., 7530., 7555., 7580., 7605., 7630., \
                      7655., 7680., 7705., 7730., 7755., 7780., 7805., \
                      7830., 7855., 7880., 7905., 7930., 7955., 7980., \
                      8005., 8030., 8055., 8080., 8105., 8130., 8155., \
                      8180., 8205., 8230., 8255., 8280., 8305., 8330., \
                      8355., 8380., 8405., 8430., 8455., 8480., 8505., \
                      8530., 8555., 8580., 8605., 8630.]
    filttran[0:89,8]=[0.000e+00, 1.000e-04, 3.000e-04, 4.000e-04, 5.000e-04, \
                      4.000e-04, 3.000e-04, 5.000e-04, 1.000e-03, 2.100e-03, \
                      3.600e-03, 6.000e-03, 1.110e-02, 2.080e-02, 3.660e-02, \
                      5.970e-02, 9.130e-02, 1.317e-01, 1.779e-01, 2.260e-01, \
                      2.719e-01, 3.125e-01, 3.470e-01, 3.755e-01, 3.978e-01, \
                      4.142e-01, 4.256e-01, 4.331e-01, 4.377e-01, 4.405e-01, \
                      4.416e-01, 4.411e-01, 4.392e-01, 4.358e-01, 4.315e-01, \
                      4.265e-01, 4.214e-01, 4.165e-01, 4.119e-01, 4.077e-01, \
                      4.039e-01, 4.006e-01, 3.975e-01, 3.943e-01, 3.906e-01, \
                      3.862e-01, 3.812e-01, 3.757e-01, 3.700e-01, 3.641e-01, \
                      3.583e-01, 3.526e-01, 3.473e-01, 3.424e-01, 3.379e-01, \
                      3.337e-01, 3.297e-01, 3.259e-01, 3.224e-01, 3.194e-01, \
                      3.169e-01, 3.150e-01, 3.132e-01, 3.111e-01, 3.081e-01, \
                      3.039e-01, 2.996e-01, 2.945e-01, 2.803e-01, 2.493e-01, \
                      2.060e-01, 1.578e-01, 1.118e-01, 7.430e-02, 4.580e-02, \
                      2.570e-02, 1.340e-02, 7.700e-03, 5.500e-03, 3.700e-03, \
                      2.300e-03, 1.500e-03, 1.100e-03, 1.100e-03, 1.100e-03, \
                      9.000e-04, 6.000e-04, 3.000e-04, 0.000e+00]
    filtwave[0:141,9]=[ 7730.,  7755.,  7780.,  7805.,  7830.,  7855.,  7880., \
                       7905., 7930.,  7955.,  7980.,  8005.,  8030.,  8055., \
                       8080.,  8105., 8130.,  8155.,  8180.,  8205.,  8230., \
                       8255.,  8280.,  8305., 8330.,  8355.,  8380.,  8405., \
                       8430.,  8455.,  8480.,  8505., 8530.,  8555.,  8580., \
                       8605.,  8630.,  8655.,  8680.,  8705., 8730.,  8755., \
                       8780.,  8805.,  8830.,  8855.,  8880.,  8905., 8930., \
                       8955.,  8980.,  9005.,  9030.,  9055.,  9080.,  9105.,\
                       9130.,  9155.,  9180.,  9205.,  9230.,  9255.,  9280.,\
                       9305., 9330.,  9355.,  9380.,  9405.,  9430.,  9455., \
                       9480.,  9505., 9530.,  9555.,  9580.,  9605.,  9630., \
                       9655.,  9680.,  9705., 9730.,  9755.,  9780.,  9805.,  \
                       9830.,  9855.,  9880.,  9905., 9930.,  9955.,  9980., \
                       10005., 10030., 10055., 10080., 10105., 10130., 10155.,\
                       10180., 10205., 10230., 10255., 10280., 10305., 10330.,\
                       10355., 10380., 10405., 10430., 10455., 10480., 10505.,\
                       10530., 10555., 10580., 10605., 10630., 10655., 10680.,\
                       10705., 10730., 10755., 10780., 10805., 10830., 10855.,\
                       10880., 10905., 10930., 10955., 10980., 11005., 11030.,\
                       11055., 11080., 11105., 11130., 11155., 11180., 11205.,\
                       11230.]
    filttran[0:141,9]=[0.00e+00, 0.00e+00, 1.00e-04, 1.00e-04, 1.00e-04, \
                       2.00e-04, 2.00e-04, 3.00e-04, 5.00e-04, 7.00e-04, \
                       1.10e-03, 1.70e-03, 2.70e-03, 4.00e-03, 5.80e-03, \
                       8.20e-03, 1.14e-02, 1.55e-02, 2.02e-02, 2.55e-02, \
                       3.11e-02, 3.69e-02, 4.28e-02, 4.84e-02, 5.36e-02, \
                       5.83e-02, 6.25e-02, 6.61e-02, 6.93e-02, 7.20e-02, \
                       7.44e-02, 7.63e-02, 7.79e-02, 7.92e-02, 8.01e-02, \
                       8.08e-02, 8.12e-02, 8.13e-02, 8.12e-02, 8.07e-02, \
                       8.01e-02, 7.91e-02, 7.79e-02, 7.66e-02, 7.50e-02, \
                       7.34e-02, 7.16e-02, 6.98e-02, 6.79e-02, 6.61e-02, \
                       6.42e-02, 6.24e-02, 6.07e-02, 5.90e-02, 5.74e-02, \
                       5.59e-02, 5.46e-02, 5.35e-02, 5.24e-02, 5.15e-02, \
                       5.05e-02, 4.96e-02, 4.85e-02, 4.74e-02, 4.62e-02, \
                       4.50e-02, 4.38e-02, 4.26e-02, 4.15e-02, 4.04e-02, \
                       3.93e-02, 3.83e-02, 3.73e-02, 3.63e-02, 3.53e-02, \
                       3.42e-02, 3.31e-02, 3.19e-02, 3.07e-02, 2.94e-02, \
                       2.80e-02, 2.67e-02, 2.53e-02, 2.40e-02, 2.27e-02, \
                       2.13e-02, 2.01e-02, 1.88e-02, 1.76e-02, 1.65e-02, \
                       1.53e-02, 1.43e-02, 1.32e-02, 1.22e-02, 1.12e-02, \
                       1.03e-02, 9.40e-03, 8.60e-03, 7.80e-03, 7.10e-03, \
                       6.40e-03, 5.80e-03, 5.20e-03, 4.70e-03, 4.20e-03, \
                       3.80e-03, 3.50e-03, 3.10e-03, 2.80e-03, 2.60e-03, \
                       2.40e-03, 2.20e-03, 2.00e-03, 1.90e-03, 1.80e-03, \
                       1.60e-03, 1.50e-03, 1.40e-03, 1.30e-03, 1.20e-03, \
                       1.10e-03, 1.00e-03, 9.00e-04, 8.00e-04, 8.00e-04, \
                       7.00e-04, 7.00e-04, 6.00e-04, 6.00e-04, 5.00e-04, \
                       5.00e-04, 4.00e-04, 4.00e-04, 3.00e-04, 3.00e-04, \
                       2.00e-04, 2.00e-04, 1.00e-04, 1.00e-04, 0.00e+00, \
                       0.00e+00]

    filtsize = [24, 21, 24, 24, 23, 47, 89, 75, 89, 141]
#		Holds the filter zero-points as determined from
#		Vega model by Dreiling & Bell (ApJ, 241,736, 1980)
#
#		B	6.268e-9   erg cm-2 s-1 A-1
#		V	3.604e-9
#		R	2.161e-9
#		I	1.126e-9
#
#		The following zero-points are from Lamla
#		(Landolt-Boernstein Vol. 2b, eds. K. Schaifer & 
#		H.H. Voigt, Berlin: Springer, p. 73, 1982 QC61.L332)
#
#		U	4.22e-9   erg cm-2 s-1 A-1
#
#		J	3.1e-10
#		H	1.2e-10
#		K	3.9e-11
#
#               U        B          V        R         I
    zp_johnson = np.array([4.22e-9, 6.268e-9, 3.604e-9, 2.161e-9, 1.126e-9])
    leff_sloan=np.array([3560.,  4830.,  6260.,  7670.,  9100])
    zp_sloan=3.631e-20*2.99792458e18/(leff_sloan*leff_sloan)
    zeropoint=np.concatenate([zp_johnson,zp_sloan])
    mag=np.zeros(10)
    filtflux=mag.copy()
    coverage=mag.copy()
    efflambda=mag.copy()
    totflux=mag.copy()
    filtername = ['U', 'B', 'V', 'R', 'I','u','g','r','i','z']
    for i,_ in enumerate(filtername):
        filtw=filtwave[0:filtsize[i],i]
        filtt=filttran[0:filtsize[i],i]
        mag[i], filtflux[i], coverage[i], efflambda[i], totflux[i]= \
              filtermag(hop[0].wave,flux, filtw, filtt, \
              zeropoint[i])                                                            
    logging.info('For object {}'.format(hop[0].obname))
    logging.info('Filter magnitude  Flux(erg/s/cm^2/A) Flux(erg/s/cm^2)  Coverage(%)  Eff. Lambda')
    for i in range(0,5):
        if (mag[i] > 99):
            logging.info('  {:1s}        FILTER AND SPECTRUM DO NOT OVERLAP'.format(filtername[i]))
        else:
            logging.info('  {:1s}     {:6.3f}      {:10.4e}        {:10.4e}         {:5.1f}         {:6.1f}'.format(filtername[i],mag[i],filtflux[i],totflux[i],coverage[i]*100.,efflambda[i]))

    for i in range(5,10):
        if (mag[i] > 99):
            logging.info('  {:1s}        FILTER AND SPECTRUM DO NOT OVERLAP'.format(filtername[i]))
        else:
            logging.info('  {:1s}     {:6.3f}      {:10.4e}        {:10.4e}         {:5.1f}         {:6.1f}'.format(filtername[i],mag[i],filtflux[i],totflux[i],coverage[i]*100.,efflambda[i]))

            
    print(' ')
    logging.info('Colors')
    colortab=[[0,1],[1,2],[2,3],[2,4],[5,6],[6,7],[7,8],[8,9]]
    for i in range(0,4):
        if (mag[colortab[i][0]] > 99) or (mag[colortab[i][1]] > 99):
            logging.info('{}-{}    ONE OR BOTH FILTERS DO NOT OVERLAP SPECTRUM'.format(filtername[colortab[i][0]],filtername[colortab[i][1]]))
        else:
            logging.info('{:1s}-{:1s}    {:12.4f}'.format(filtername[colortab[i][0]],filtername[colortab[i][1]],mag[colortab[i][0]]-mag[colortab[i][1]]))
    for i in range(4,8):
        if (mag[colortab[i][0]] > 99) or (mag[colortab[i][1]] > 99):
            logging.info('{}-{}    ONE OR BOTH FILTERS DO NOT OVERLAP SPECTRUM'.format(filtername[colortab[i][0]],filtername[colortab[i][1]]))
        else:
            logging.info('{:1s}-{:1s}    {:12.4f}'.format(filtername[colortab[i][0]],filtername[colortab[i][1]],mag[colortab[i][0]]-mag[colortab[i][1]]))

    print('\nWould you like to scale the spectrum to match photometry?\n')
    answer=yesno('n')
    if (answer == 'y'):
        print('\nWhich filter do you have?')
        scalefilt=inputter_single_mix('U/B/V/R/I/u/g/r/i/z: ','UBVRIugriz')
        filtindex=filtername.index(scalefilt)
        scalemag=inputter('Enter your value for filter {}: '.format(filtername[filtindex]),'float',False)
        print(' ')
        logging.info('Scaling {} from {}={:6.3f} to {}={}'.format(hop[0].obname,filtername[filtindex],mag[filtindex],filtername[filtindex],scalemag))
        logging.info('Multiplying by {:.3f}'.format(10**(0.4*(mag[filtindex]-scalemag))))
        hop[0].flux=hop[0].flux*10**(0.4*(mag[filtindex]-scalemag))
    

    return hop


def womirfilters(hop):
    """calculate photometric values from spectra"""
    print('NOTE:  The routine expects an f_lambda spectrum')
    print('       I will try to guess if the spectrum')
    print('       has been scaled by 1E15')
    print(' ')
    print('       Check this before believing fluxes')
    print(' ')
    print('NOTE Also:  These are the 2MASS filter curves')
    print(' ')
    flux=hop[0].flux.copy()
    if (np.mean(flux) > 0.00001):
        flux = flux *1.e-15

    filtwave=np.zeros((109,3))
    filttran=np.zeros((109,3))

    filtwave[:,0]=[1.050, 1.051, 1.062, 1.066, 1.070, 1.075, 1.078, 1.082, \
        1.084, 1.087, 1.089, 1.093, 1.096, 1.102, 1.105, 1.107, 1.109, 1.112, \
        1.116, 1.117, 1.120, 1.123, 1.128, 1.129, 1.132, 1.134, 1.138, 1.140, \
        1.143, 1.147, 1.154, 1.159, 1.164, 1.167, 1.170, 1.173, 1.175, 1.179, \
        1.182, 1.186, 1.188, 1.192, 1.195, 1.199, 1.202, 1.209, 1.216, 1.221, \
        1.227, 1.231, 1.236, 1.240, 1.244, 1.247, 1.253, 1.255, 1.258, 1.260, \
        1.265, 1.270, 1.275, 1.279, 1.286, 1.292, 1.297, 1.302, 1.305, 1.307, \
        1.310, 1.313, 1.316, 1.319, 1.323, 1.326, 1.330, 1.333, 1.334, 1.336, \
        1.339, 1.343, 1.346, 1.349, 1.353, 1.355, 1.360, 1.363, 1.370, 1.373, \
        1.377, 1.383, 1.388, 1.392, 1.395, 1.396, 1.397, 1.398, 1.400, 1.401, \
        1.402, 1.404, 1.406, 1.407, 1.410, 1.412, 1.416, 1.421, 1.426, 1.442, \
        1.450]


    filttran[:,0]=[0.0000, 0.0000, 0.0000, 0.0023, 0.0087, 0.0150, 0.0309, 0.0690, \
        0.1136, 0.1709, 0.2282, 0.2886, 0.3491, 0.4255, 0.4668, 0.5209, \
        0.5687, 0.6228, 0.6546, 0.6864, 0.7150, 0.7437, 0.7595, 0.7595, \
        0.7435, 0.7276, 0.6861, 0.6575, 0.6224, 0.5873, 0.5649, 0.5840, \
        0.6157, 0.6571, 0.6857, 0.7271, 0.7685, 0.8162, 0.8416, 0.8511, \
        0.8447, 0.8256, 0.7937, 0.7554, 0.7172, 0.6757, 0.6629, 0.6883, \
        0.7391, 0.7869, 0.8505, 0.8823, 0.8950, 0.8854, 0.8471, 0.8184, \
        0.7802, 0.7324, 0.6845, 0.6239, 0.5889, 0.5729, 0.5728, 0.5918, \
        0.6172, 0.6681, 0.6968, 0.7286, 0.7667, 0.7954, 0.8431, 0.8813, \
        0.9194, 0.9353, 0.9257, 0.9225, 0.9129, 0.8906, 0.8524, 0.8141, \
        0.7854, 0.7599, 0.7439, 0.7375, 0.7247, 0.7183, 0.7087, 0.7023, \
        0.7022, 0.7181, 0.7339, 0.7147, 0.6829, 0.6446, 0.6160, 0.5873, \
        0.5172, 0.4662, 0.3770, 0.2305, 0.1350, 0.1126, 0.0712, 0.0362, \
        0.0170, 0.0042, 0.0009, 0.0007, 0.0000]


    filtwave[0:57,1]=[1.315, 1.341, 1.368, 1.397, 1.418, 1.440, 1.462, 1.478, \
        1.486, 1.493, 1.504, 1.515, 1.528, 1.539, 1.546, 1.551, 1.556, 1.565, \
        1.572, 1.577, 1.583, 1.592, 1.597, 1.602, 1.613, 1.619, 1.628, 1.633, \
        1.642, 1.648, 1.657, 1.659, 1.671, 1.684, 1.701, 1.715, 1.727, 1.739, \
        1.746, 1.751, 1.753, 1.756, 1.764, 1.775, 1.785, 1.790, 1.796, 1.803, \
        1.810, 1.813, 1.818, 1.828, 1.835, 1.850, 1.871, 1.893, 1.914]


    filttran[0:57,1]=[0.0014, 0.0014, 0.0000, 0.0000, 0.0014, 0.0028, 0.0070, \
        0.0252, 0.0700, 0.1807, 0.3529, 0.4972, 0.6527, 0.7591, 0.8109, \
        0.8319, 0.8403, 0.8389, 0.8305, 0.8235, 0.8193, 0.8277, 0.8347, \
        0.8375, 0.8319, 0.8193, 0.8081, 0.8053, 0.8095, 0.8165, 0.8263, \
        0.8305, 0.8375, 0.8431, 0.8501, 0.8529, 0.8543, 0.8529, 0.8445, \
        0.8305, 0.8151, 0.7927, 0.7255, 0.6275, 0.5084, 0.4258, 0.3291, \
        0.2101, 0.1275, 0.0882, 0.0560, 0.0294, 0.0154, 0.0070, 0.0028, \
        0.0014, 0.0000]


    filtwave[0:76,2]=[1.900, 1.915, 1.927, 1.934, 1.939, 1.948, 1.957, 1.962, \
        1.969, 1.976, 1.981, 1.989, 1.990, 1.998, 2.008, 2.014, 2.019, 2.028, \
        2.037, 2.045, 2.061, 2.072, 2.075, 2.082, 2.089, 2.099, 2.106, 2.113, \
        2.120, 2.124, 2.138, 2.145, 2.155, 2.169, 2.176, 2.185, 2.197, 2.208, \
        2.213, 2.218, 2.232, 2.237, 2.248, 2.256, 2.260, 2.263, 2.265, 2.270, \
        2.272, 2.276, 2.277, 2.281, 2.284, 2.286, 2.291, 2.293, 2.295, 2.297, \
        2.299, 2.306, 2.311, 2.316, 2.320, 2.325, 2.328, 2.335, 2.339, 2.344, \
        2.346, 2.352, 2.361, 2.363, 2.370, 2.375, 2.384, 2.399]

    filttran[0:76,2]=[0.0000, 0.0013, 0.0027, 0.0040, 0.0082, 0.0153, 0.0293, \
        0.0462, 0.0743, 0.1222, 0.1714, 0.2672, 0.3517, 0.4263, 0.6262, \
        0.6797, 0.7487, 0.7853, 0.8120, 0.8303, 0.8485, 0.8513, 0.8583, \
        0.8597, 0.8667, 0.8751, 0.8765, 0.8835, 0.8891, 0.8863, 0.8848, \
        0.8819, 0.8805, 0.8748, 0.8804, 0.8818, 0.8902, 0.8986, 0.9014, \
        0.8999, 0.8999, 0.8956, 0.8913, 0.8969, 0.8997, 0.8997, 0.9053, \
        0.9109, 0.9166, 0.9109, 0.9025, 0.8870, 0.8686, 0.8433, 0.7714, \
        0.7292, 0.6650, 0.5950, 0.5333, 0.4094, 0.3108, 0.2234, 0.1544, \
        0.1234, 0.0896, 0.0599, 0.0416, 0.0320, 0.0300, 0.0162, 0.0063, \
        0.0007, 0.0034, 0.0020, 0.0006, 0.0000]

    filtwave=filtwave*10000.0
    
    filtsize = [109, 57, 76]
    #		Holds the filter zero-points as determined from
#		Vega model by Dreiling & Bell (ApJ, 241,736, 1980)
#
#		B	6.268e-9   erg cm-2 s-1 A-1
#		V	3.604e-9
#		R	2.161e-9
#		I	1.126e-9
#
#		The following zero-points are from Lamla
#		(Landolt-Boernstein Vol. 2b, eds. K. Schaifer & 
#		H.H. Voigt, Berlin: Springer, p. 73, 1982 QC61.L332)
#
#		U	4.22e-9   erg cm-2 s-1 A-1
#
#		J	3.1e-10
#		H	1.2e-10
#		K	3.9e-11
#
#               U        B          V        R         I

    zeropoint = [3.1e-10, 1.2e-10,3.9e-11]

    mag=np.zeros(3)
    filtflux=mag.copy()
    coverage=mag.copy()
    efflambda=mag.copy()
    totflux=mag.copy()
    filtername = ['J', 'H', 'K']
    for i,_ in enumerate(filtername):
        filtw=filtwave[0:filtsize[i],i]
        filtt=filttran[0:filtsize[i],i]
        mag[i], filtflux[i], coverage[i], efflambda[i], totflux[i]= \
              filtermag(hop[0].wave,flux, filtw, filtt, \
              zeropoint[i])                                                            
    logging.info('For object {}'.format(hop[0].obname))
    logging.info('Filter magnitude  Flux(erg/s/cm^2/A) Flux(erg/s/cm^2)  Coverage(%)  Eff. Lambda')
    for i in range(0,3):
        if (mag[i] > 99):
            logging.info('  {:1s}        FILTER AND SPECTRUM DO NOT OVERLAP'.format(filtername[i]))
        else:
            logging.info('  {:1s}     {:6.3f}      {:10.4e}        {:10.4e}         {:5.1f}         {:7.1f}'.format(filtername[i],mag[i],filtflux[i],totflux[i],coverage[i]*100.,efflambda[i]))


            
    print(' ')
    logging.info('Colors')
    colortab=[[0,1],[1,2]]
    for i in range(0,2):
        if (mag[colortab[i][0]] > 99) or (mag[colortab[i][1]] > 99):
            logging.info('{}-{}    ONE OR BOTH FILTERS DO NOT OVERLAP SPECTRUM'.format(filtername[colortab[i][0]],filtername[colortab[i][1]]))
        else:
            logging.info('{:1s}-{:1s}    {:12.4f}'.format(filtername[colortab[i][0]],filtername[colortab[i][1]],mag[colortab[i][0]]-mag[colortab[i][1]]))


    print('\nWould you like to scale the spectrum to match photometry?\n')
    answer=yesno('n')
    if (answer == 'y'):
        print('\nWhich filter do you have?')
        scalefilt=inputter_single_mix('J/H/K: ','JHK')
        filtindex=filtername.index(scalefilt)
        scalemag=inputter('Enter your value for filter {}: '.format(filtername[filtindex]),'float',False)
        print(' ')
        logging.info('Scaling {} from {}={:6.3f} to {}={}'.format(hop[0].obname,filtername[filtindex],mag[filtindex],filtername[filtindex],scalemag))
        logging.info('Multiplying by {:.3f}'.format(10**(0.4*(mag[filtindex]-scalemag))))
        hop[0].flux=hop[0].flux*10**(0.4*(mag[filtindex]-scalemag))
    

    return hop

def womhertz(hop):
    """converts wavelength to hertz for f_nu spectrum"""
    print('NOTE:  The routine expects an f_nu spectrum')
    print('       I will try to guess if the spectrum')
    print("       has been scaled by 1E26 (it shouldn't be)")
    print(' ')
    print('       Check this before believing any result')
    print(' ')


    if (np.mean(hop[0].flux) > 0.00001):
        hop[0].flux=hop[0].flux*1.e-26

    hop[0].wave = hop[0].wave*1.e-10
    hop[0].wave = 2.99792458e8/hop[0].wave
    #FIX head

    print('Active spectrum now in hertz vs. f_nu')
    return hop

def womsqr(hop):
    """square flux"""
    hop[0].flux=hop[0].flux * hop[0].flux
    print('\nSpectrum has been squared\n')
    return hop

def womsqrt(hop):
    """square root of flux"""
    hop[0].flux=np.sqrt(hop[0].flux)
    print('\nSquare root of spectrum now active\n')
    return hop

def womblo(hop,fig):
    """blotch out bad data"""
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


def womspl(hop,fig):
    """fit spline to spectrum"""
    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    print('\nObject is {}\n'.format(hop[0].obname))
    womplot(hop)
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    nsplinepoints=0
    tmpsplptsx=[]
    tmpsplptsy=[]

    done=False
    while (not done):
        plt.cla()
        plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
        if (len(tmpsplptsx) > 0):
            plt.plot(tmpsplptsx,tmpsplptsy,'ro')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title(hop[0].obname)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
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
        splineresult=splev(hop[0].wave,spline)
        plt.plot(hop[0].wave,splineresult,drawstyle='steps-mid')
        print('Is this fit OK? ')
        answer=yesno('y')
        if (answer == 'y'):
            done=True
    print('\nSubtract spline fit from flux?\n')
    sub=yesno('n')
    if (sub == 'y'):
        hop[0].flux=hop[0].flux - splineresult
    print('\nStore spline in hopper?\n')
    store=yesno('y')
    if (store == 'y'):
        hopnum=0
        while (hopnum < 1) or (hopnum > HOPSIZE):
            hopnum=inputter('Store in which hopper: ','int',False)
        hop[hopnum]=copy.deepcopy(hop[0])
        hop[hopnum].flux=splineresult.copy()
        hop[hopnum].obname=hop[hopnum].obname+'spline'
        hop[hopnum].var=np.zeros(len(hop[0].wave))
    return hop

def womvelcho(hop):
    """select region of spectrum, put onto velocity scale"""
    wave=hop[0].wave.copy()
    flux=hop[0].flux.copy()
    var=hop[0].var.copy()
    nwave,nflux,mode=womwaverange(wave,flux,'none')
    wavebindex=womget_element(wave,nwave[0])
    waverindex=womget_element(wave,nwave[-1])
    nvar=var[wavebindex:waverindex+1].copy()
    binvec=np.arange(len(nwave))
    wavelog=np.log(nwave[0])+((np.log(nwave[-1])-np.log(nwave[0]))/len(nwave))*binvec
    wavelog=np.exp(wavelog)
    fluxlog=womscipyrebin(nwave,nflux,wavelog)
    varlog=womscipyrebin(nwave,nvar,wavelog)
    zp=wavelog[0]-1.0
    while (zp < wavelog[0]) or (zp > wavelog[-1]):
        zp=inputter('Enter the zero point for velocity (in angstroms): ','float',False)
    indexzp=womget_element(wavelog, zp)
    print(' ')
    logging.info('Zero at bin {}'.format(indexzp))
    logging.info('with lambda {}'.format(wavelog[indexzp]))
    print(' ')
    z=(wavelog - zp)/zp
    square=(z+1)*(z+1)
    wavekmsrel=((square-1.)/(square+1.))*299792.458
    kmsperbin=np.zeros(len(nwave))
    for i in range(1,len(nwave)):
        kmsperbin[i]=2.0*2.99792458e5*(wavelog[i]-wavelog[i-1])/\
                      (wavelog[i]+wavelog[i-1])
    kmsmean=np.mean(kmsperbin[1:])
    logging.info('Average km/s per bin: {}'.format(kmsmean))
    logging.info('km/s at bin 1: {}'.format(kmsperbin[1]))
    logging.info('km/s at bin n: {}'.format(kmsperbin[-1]))
    print(' ')
    wavekms=kmsmean*binvec+kmsmean/2.0
    indexzp=womget_element(wavekms,zp)
    wavekms=wavekms-wavekms[indexzp]
    hop[0].wave=wavekmsrel.copy()
    hop[0].flux=fluxlog.copy()
    hop[0].var=varlog.copy()
    return hop

def womint(hop):
    """calculate intensity in given wavelength range"""
    print(' ')
    logging.info('Object is {}'.format(hop[0].obname))
    print(' ')
    print('Spectrum runs from {} to {}'.format(hop[0].wave[0],hop[0].wave[-1]))
    print(' ')
    print('This routine expects the spectrum to be in flambda units.')
    print('It also expects a linear wavelength scale.')
    print(' ')
    print('Choose general region of spectrum\n')
    nwave,nflux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
    print('\nNow pick the exact range for the intensity calculation')
    wavint,fluxint,mode=womwaverange(nwave,nflux,mode)
    indexblue=womget_element(nwave, wavint[0])
    indexred=womget_element(nwave,wavint[-1])
    wdelt=nwave[1]-nwave[0]
    lineflux=np.sum(nflux[indexblue:indexred+1])*wdelt
    linefluxin=np.sum(nflux[indexblue+1:indexred])*wdelt
    linefluxout=np.sum(nflux[indexblue-1:indexred+2])*wdelt
    print(' ')
    logging.info('FWZI (approximate): {}'.format(nwave[indexred]-nwave[indexblue]))
    logging.info('Line flux (ergs/sec/cm^2): {}'.format(lineflux))
    logging.info('Line flux one pixel in: {}'.format(linefluxin))
    logging.info('Line flux one pixel out: {}'.format(linefluxout))
    logging.info('Note that flux may need to be scaled by 1e-15')
    logging.info('Average difference (between line flux and one pixel')
    avgdiff=(np.abs(linefluxin - lineflux) + np.abs(linefluxout-lineflux))/2.0
    logging.info('in or out): {}'.format(avgdiff))
    logging.info('As a percentage of line flux: {}'.format(100.0*avgdiff/lineflux))
    return hop


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gauss_cont(x, *p):
    A, mu, sigma, m, b = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + m*x + b


def womgau(hop):
    """fit gaussian to line"""
    from scipy.optimize import curve_fit
    print(' ')
    logging.info('Object is {}'.format(hop[0].obname))
    print(' ')
    print('Spectrum runs from {} to {}'.format(hop[0].wave[0],hop[0].wave[-1]))
    print(' ')
    print('This routine expects the spectrum to be in flambda units.')
    print('It also expects a linear wavelength scale.')
    print(' ')
    print('Choose general region of spectrum\n')
    nwave,nflux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
    print('\nNow pick the exact range for the fit')
    waveint,fluxint,mode=womwaverange(nwave,nflux,mode)
    indexblue=womget_element(nwave, waveint[0])
    indexred=womget_element(nwave,waveint[-1])
    if (mode == 'w'):
        done = False
        while (not done):
            print(' ')
            wavecenter=inputter('Enter approximate center of Gaussian : ','float',False)
            indexcenter=womget_element(waveint,wavecenter)
            if (indexcenter <= 0) or (wavecenter > waveint[-1]):
                print('Bad central wavelength, try again')
            else:
                done = True
    else:
        done=False
        while (not done):
            print('Mark the approximate center of the Gaussian')
            pickcent=plt.ginput(1)
            indexcenter=womget_element(waveint,pickcent[0][0])
            print('\nApproximate center at {}'.format(waveint[indexcenter]))
            print('\nIs this OK?')
            answer=yesno('y')
            if (answer == 'y'):
                done=True
    weights=np.sqrt(hop[0].var[indexblue:indexred+1])
    print(' ')
    continuum=inputter_single('Do you want to fit gaussian with (c)ontinuum, or (n)o continuum? ','cn')
    if (continuum == 'c'):
        p=[fluxint[indexcenter], waveint[indexcenter],3.0,1.0,waveint[0]]
        result=curve_fit(gauss_cont,waveint,fluxint,sigma=weights,p0=p,absolute_sigma=True,full_output=True)
    else:
        p=[fluxint[indexcenter], waveint[indexcenter],3.0]
        result=curve_fit(gauss,waveint,fluxint,sigma=weights,p0=p,absolute_sigma=True,full_output=True)
    coefferr=np.sqrt(np.diag(result[1]))
    coeff=result[0]
        # make 'finer-grained' version of fit, 0.2A/pix for calculations
    wavecalc=np.arange(2*5*50*abs(coeff[2]))*0.2+coeff[1]-0.2*5*50*abs(coeff[2])
    calccenter=womget_element(wavecalc,coeff[1])
    if (continuum == 'c'):
        fluxcalc=gauss_cont(wavecalc,*coeff)
        fluxcont=wavecalc*coeff[3]+coeff[4]
        fluxgaussian=fluxcalc-fluxcont
        linecont=fluxcont[calccenter]
    else:
        fluxcalc=gauss(wavecalc,*coeff)
    
    
    deltafit=wavecalc[1]-wavecalc[0]
    calcindexblue=womget_element(wavecalc,waveint[0])
    calcindexred=womget_element(wavecalc,waveint[-1])
    sumfluxcalc=np.sum(fluxcalc[calcindexblue:calcindexred+1]*deltafit)
    sumallfluxcalc=np.sum(fluxcalc*deltafit)
    chi=(result[2]['fvec']**2).sum()
    redchi=chi/(len(waveint)-len(coeff))
    if (continuum == 'c'):
        sumfluxgaussian=np.sum(fluxgaussian[calcindexblue:calcindexred+1]*deltafit)
        sumallfluxgaussian=np.sum(fluxgaussian*deltafit)
        sumfluxcont=np.sum(fluxcont[calcindexblue:calcindexred+1]*deltafit)
        sumallfluxcont=np.sum(fluxcont*deltafit)
        # propagate uncertainty (from old version) not sure this is correct
        height_pct=coefferr[0]/coeff[0]
        sigma_pct=coefferr[2]/coeff[2]
        flux_pct=np.sqrt(height_pct**2+sigma_pct**2)
        sumfluxgaussiansig=sumfluxgaussian*flux_pct
        sumallfluxgaussiansig=sumallfluxgaussian*flux_pct
    plt.cla()
    plt.plot(nwave,nflux,drawstyle='steps-mid',color='k')
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.plot(wavecalc,fluxcalc,drawstyle='steps-mid',color='b')
    if (continuum == 'c'):
        plt.plot(wavecalc,fluxgaussian,drawstyle='steps-mid',color='r')
        plt.plot(wavecalc,fluxcont,drawstyle='steps-mid',color='g')
    plt.plot([waveint[0],waveint[0]],[ymin,ymax],color='k',linestyle='--')
    plt.plot([waveint[-1],waveint[-1]],[ymin,ymax],color='k',linestyle='--')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    logging.info('For object {} Gaussian fit'.format(hop[0].obname))
    if (continuum == 'c'):
        print('\nData = Black, Fit = Blue, Continuum = Green, Fit-Continuum = Red\n')
    else:
        print('\nData = Black, Fit = Blue\n')
    logging.info('Height      {:16.8f}+/-{:16.8f}'.format(coeff[0],coefferr[0]))
    logging.info('Center      {:16.8f}+/-{:16.8f}'.format(coeff[1],coefferr[1]))
    logging.info('Sigma       {:16.8f}+/-{:16.8f}'.format(coeff[2],coefferr[2]))
    if (continuum == 'c'):
        logging.info('Slope       {:16.8f}+/-{:16.8f}'.format(coeff[3],coefferr[3]))
        logging.info('Y-intercept {:16.8f}+/-{:16.8f}'.format(coeff[4],coefferr[4]))
        logging.info('FWHM        {:16.8f}+/-{:16.8f}'.format(2.35482*np.abs(coeff[2]),2.35482*coefferr[2]))
        logging.info('Flux between dotted lines (Gaussian): {:16.8f}+/-{:16.8f}'.format(sumfluxgaussian, sumfluxgaussiansig))
        logging.info('EW between dotted lines (Gaussian): {:16.8f}'.format(sumfluxgaussian/linecont))
        logging.info('Flux for full (Gaussian): {:16.8f}+/-{:16.8f}'.format(sumallfluxgaussian, sumallfluxgaussiansig))
        logging.info('EW for full (Gaussian): {:16.8f}'.format(sumallfluxgaussian/linecont))
        logging.info('Continuum flux at line center: {:16.8f}'.format(linecont))
    logging.info('Chi^2: {}'.format(chi))
    logging.info('Reduced chi^2: {}'.format(redchi))
    logging.info('All fluxes might need to be scaled by 1e-15')
    print(' ')
    return hop

def womwavescale(hop):
    """multiplicative change to wavelength scale, convert units"""
    print('Routine does a multiplicative shift to the wavelength scale')
    shift=inputter('Enter wavelength factor: ','float',False)
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    wshow()
    hop[0].wave=hop[0].wave*shift
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    logging.info('File {} wavelength scale multiplied by {}'.format(hop[0].obname, shift))
    #FIX header
    return hop


def womhelp(hop):
    """list function of each command"""
    print(' ')
    print('Wombat.py is the python replacement for Sarah.  She does many of the ')
    print('same things, but usually in a faster and simpler way.  The ')
    print('organization is quite simple.  Once a spectrum is read in, ')
    print('it is the active object for all subsequent operations, until ')
    print('another spectrum replaces it as the active one.  There is a ')
    print('twenty (20) slot hopper that can store the spectra for later ')
    print('use.  Wombat will maintain header information from a fits file, ')
    print('and adjust it as operations are performed.  She is unfortunately ')
    print('very finicky about input, so don''t try to use letters when you ')
    print('should be using numbers.  Questions should be directed to ')
    print('Tom Matheson.  Good luck.  You''ll need it.')
    print(' ')
    print('Hit any key')
    a = getch()
    print(' ')
    print('rpl:  read two-column ascii spectrum, (wavelength, flux)')
    print('wpl:  write two-column ascii spectrum, (wavelength, flux)')
    print('rfits:  read fits format spectrum (1-D only)')
    print('wfits:  write fits format spectrum (1-D only)')
    print('rmsfits: read multi-spec fits file')
    print('hop:  store spectrum in hopper (20 available)')
    print('rh, rdhop: read spectrum from hopper to active status')
    print('buf, buffer: list contents of hopper')
    print('e, examine: print flux values for given wavelength range')
    print('p, plot: plot spectrum (asks for wavelength range, y-scale adjustable')
    print('ph, hard: OBSOLETE, use matplotlib save button')
    print('pl, plotl, plotlog:  plot with x, y, or both on log-scale')
    print('apl, oplot: add plot to current plot window')
    print('cho, choose: select subset of spectrum by wavelength range')
    print('stat: compute statistics on given wavelength range')
    print('bin: change wavelength binning scale (linear binning, ')
    print('      uses scipy interpolate with cubic')
    print('ashrebin: Ashley quadratic interpolation (classic mode)')
    print('b, blotch: blotch out bad pixels interactively')
    print('com, combine: add two spectra with equal weight')
    print('cmw: add two spectra with unequal weights')
    print('cat: interactively splice together blue and red halves of spectra')
    print('join: combine red and blue pieces of abutting spectra')
    print(' ')
    print('Hit any key')
    a = getch()
    print(' ')
    print('spl, spline: fit spline to spectrum, option to subtract, save spline')
    print('w, redshift: remove redshift in z or km/s, rebins automatically')
    print('ms, arith: spectrum arithmetic, add, subtract, multiply, divide')
    print('fnu, flux: convert flux scale')
    print('smo, smooth: smooth spectrum with boxcar or Savitzky-Golay')
    print('sca, scale: multiply spectrum by constant')
    print('avt, scale2: scale one spectrum to match another')
#    print('win, window: eliminate flux values beyond given ranges')
    print('velcho: choose subset, put on velocity scale, zero at given value')
    print('gau: fit gaussian+line to emission lines, returns flux, EW')
    print('int: calculate intensity in emission lines')
    print('linesub: use one line as profile, subtract elsewhere')
    print('red: (de)redden spectra')
    print('bluen: bluen spectrum with scattering law')
    print('com3: combine three spectra equally')
    print('head, fitshead: page through the header')
    print('zap: median filter spectrum')
    print('xshift: add constant to abcissa (wavelength)')
    print('yshift: add constant to ordinate (flux)')
    print('bb, planck: calculate blackbody curve for given temp.')
    print('atmdisp: calculate atmospheric dispersion effects')
#    print('zcalc: determine redshift between two spectra (BETA)')
    print('$[command]: executes unix (shell) [command]')
    print('help, ?: print these help pages')
    print('q, quit: quit the program')
    return hop
