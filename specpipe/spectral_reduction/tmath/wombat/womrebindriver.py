def womrebindriver(hop,mode):
    """input info to call womrebin functions"""
    import matplotlib.pyplot as plt
    import numpy as np
    from tmath.wombat.inputter import inputter
    from tmath.wombat.waveparse import waveparse
    from tmath.wombat.yesno import yesno
    from tmath.wombat.womashrebin import womashrebin
    from tmath.wombat.womscipyrebin import womscipyrebin
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

