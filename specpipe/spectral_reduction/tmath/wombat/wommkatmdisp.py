def wommkatmdisp(hop):
    import numpy as np
    import matplotlib.pyplot as plt
    from tmath.wombat.waveparse import waveparse
    from tmath.wombat.inputter import inputter
    from tmath.wombat.airtovac import airtovac
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

