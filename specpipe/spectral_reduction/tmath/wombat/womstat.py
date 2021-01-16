def womstat(hop):
    """print statistics of region of spectrum"""
    import scipy.stats as st
    import numpy as np
    from tmath.wombat.womwaverange import womwaverange
    print('\nObject is {}'.format(hop[0].obname))
    print('\nEnter range for statistics: ')
    wave,flux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
    print('\nMean: {}'.format(np.mean(flux)))
    print('Variance: {}'.format(np.var(flux,ddof=1)))
    print('Std. Dev.: {}'.format(np.std(flux,ddof=1)))
    print('Mean Dev.: {}'.format(np.mean(np.abs(flux-np.mean(flux)))))
    print('S/N: {}'.format(np.mean(flux)/np.std(flux,ddof=1)))
    #this is IRAF definition for S/N
    print('Skewness: {}'.format(st.skew(flux)))
    print('Kurtosis: {}'.format(st.kurtosis(flux)))
    print('Median: {}'.format(np.median(flux)))
    print('No. points: {}'.format(len(flux)))
    return hop
    
