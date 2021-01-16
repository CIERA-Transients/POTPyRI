def womcho(hop):
    """choose region of spectrum with wavelengths or mouse"""
    import matplotlib.pyplot as plt
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
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

