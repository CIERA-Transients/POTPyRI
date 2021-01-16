def womplotlog(hop):
    """plot with log scale for x, y, or both"""
    import matplotlib.pyplot as plt
    from tmath.wombat.inputter_single import inputter_single
    
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
    
