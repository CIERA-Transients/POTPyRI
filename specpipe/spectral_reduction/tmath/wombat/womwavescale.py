def womwavescale(hop):
    """multiplicative change to wavelength scale, convert units"""
    import matplotlib.pyplot as plt
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat.wshow import wshow
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


