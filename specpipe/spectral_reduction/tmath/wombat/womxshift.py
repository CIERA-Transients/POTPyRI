def womxshift(hop):
    """linearly shift wavelength by arbitrary amount"""
    import matplotlib.pyplot as plt
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat.wshow import wshow
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

