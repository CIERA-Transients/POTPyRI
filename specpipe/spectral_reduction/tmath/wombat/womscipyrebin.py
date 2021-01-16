def womscipyrebin(wave,flux,nwave):
    import scipy.interpolate as scint
    inter=scint.interp1d(wave,flux,kind='cubic',bounds_error=False,fill_value='extrapolate')
    nflux=inter(nwave)
    return nflux

