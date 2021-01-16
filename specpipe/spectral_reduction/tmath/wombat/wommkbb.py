def wommkbb(hop):
    import numpy as np
    import matplotlib.pyplot as plt
    from tmath.wombat.waveparse import waveparse
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
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

