def womwaverange2(waveb,fluxb,waver,fluxr,mode):
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
    plt.cla()
    plt.plot(waveb,fluxb,drawstyle='steps-mid',color='b')
    plt.plot(waver,fluxr,drawstyle='steps-mid',color='r')
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    wshow()
    done = False
    while (not done):
        if (mode == 'w'):
            waveblue,wavered=waveparse()
        else:
            pickpoints=plt.ginput(2,timeout=-1)
            waveblue=pickpoints[0][0]
            wavered=pickpoints[1][0]
        if (waveblue > wavered):
            waveblue,wavered=wavered,waveblue
        indexblue=womget_element(waveb,waveblue)
        indexred=womget_element(waveb,wavered)
        print('Range selected: {} to {}'.format(waveb[indexblue],waveb[indexred]))
        plt.plot(waveb[indexblue:indexred+1],fluxb[indexblue:indexred+1], \
                 drawstyle='steps-mid',color='k')
        plt.pause(0.01)
        print('Is this range correct?')
        answer=yesno('y')
        if (answer == 'n'):
            plt.cla()
            plt.plot(waveb,fluxb,drawstyle='steps-mid',color='b')
            plt.plot(waver,fluxr,drawstyle='steps-mid',color='r')
            plt.ylabel('Flux')
            plt.xlabel('Wavelength')
            wshow()
        else:
            done=True
    return waveblue,wavered,mode

