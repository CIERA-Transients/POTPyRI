def filtermag(wave,flux,filtw,filtt,zp):
    """calculate magnitude in given filter"""
    import numpy as np
    from scipy.interpolate import splrep,splev
    from tmath.wombat.filterconv import filterconv
    if (wave[0] > filtw[-1]) or (wave[-1] < filtw[0]):
        print('Filter and spectrum do not overlap')
        print('Bailing out and returning mag of 999')
        mag = 999
        flux = 999
        coverage = 0
        efflambda = 0
        totflux=0
        return mag, flux, coverage, efflambda, totflux
    blueover=np.where(wave < filtw[0])
    if (len(blueover[0]) < 1):
        blueindex=0
    else:
        blueindex=blueover[0][-1] + 1
    redover=np.where(wave > filtw[-1])
    if (len(redover[0]) < 1):
        redindex=len(wave)-1
    else:
        redindex=redover[0][0]-1
    wavematch=wave[blueindex:redindex+1]
    fluxmatch=flux[blueindex:redindex+1]
    spline=splrep(filtw, filtt, k=3)
    splineresult=splev(wavematch,spline)
    flux,totflux,efflambda,coverage=filterconv(wavematch,fluxmatch,splineresult,filtw,filtt)
    mag=-2.5*np.log10(flux/zp)
    return mag, flux, coverage, efflambda, totflux

