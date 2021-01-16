def womwaverange(wave,flux,mode):
    import matplotlib.pyplot as plt
    from tmath.wombat.waveparse import waveparse
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.yesno import yesno
    from tmath.wombat.wshow import wshow
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
            pickpoints=plt.ginput(2,timeout=-1)
            waveb=pickpoints[0][0]
            waver=pickpoints[1][0]
        if (waveb > waver):
            waveb,waver=waver,waveb
        indexblue=womget_element(wave,waveb)
        indexred=womget_element(wave,waver)
        print('Range selected: {} to {}'.format(wave[indexblue],wave[indexred]))
        plt.plot(wave[indexblue:indexred+1],flux[indexblue:indexred+1], \
                 drawstyle='steps-mid',color='r')
        plt.pause(0.01)
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

