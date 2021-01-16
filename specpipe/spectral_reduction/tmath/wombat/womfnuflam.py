def womfnuflam(hop):
    """convert flux densities"""
    import numpy as np
    from tmath.wombat.inputter_single import inputter_single
    light_speed=2.99792458e18
    print('Converting flux units.  Is the input spectrum in f-(n)u,')
    print('f-(l)ambda, (a)b magnitudes, or (s)t magnitudes?')
    flun=inputter_single('(n/l/a/s)','nlas')
    fluxmean=np.mean(hop[0].flux)
    scale=False
    #convert input to f-lambda, remove any arbitrary scaling
    if (flun == 'n'):
        if (fluxmean > 1.e-6):
            hop[0].flux=hop[0].flux/1.e26
            scale=True
        hop[0].flux=hop[0].flux*light_speed/hop[0].wave/hop[0].wave
    elif (flun == 'a'):
        hop[0].flux=10**(-0.4*hop[0].flux - 19.44)
        hop[0].flux=hop[0].flux*light_speed/hop[0].wave/hop[0].wave
    elif (flun == 's'):
        hop[0].flux=10**(-0.4*hop[0].flux - 8.44)
    else:
        if (fluxmean > 1.e-6):
            hop[0].flux=hop[0].flux/1.e15
            scale=True

    print('\nConvert to f-(n)u, f-(l)ambda, (a)b magnitudes,')
    print('or (s)t magnitudes?')
    flout=inputter_single('(n/l/a/s)','nlas')
    if (flout == 'n'):
        hop[0].flux=hop[0].flux*hop[0].wave*hop[0].wave/light_speed
        if (scale):
            hop[0].flux=hop[0].flux*1.e26
    elif (flout == 'a'):
        hop[0].flux=hop[0].flux*hop[0].wave*hop[0].wave/light_speed
        hop[0].flux=hop[0].flux*10**19.44
        negloc=np.where(hop[0].flux <= 0)
        hop[0].flux[negloc]=1.0
        hop[0].flux=-2.5*np.log10(hop[0].flux)
    elif (flout == 's'):
        negloc=np.where(hop[0].flux <= 0)
        hop[0].flux[negloc]=3.63e-19
        hop[0].flux=-2.5*np.log10(hop[0].flux) - 21.10
    else:
        if (scale):
            hop[0].flux=hop[0].flux*1.e15
    return hop

