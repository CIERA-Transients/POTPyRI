def womsmo(hop):
    import matplotlib.pyplot as plt
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.inputter import inputter
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

