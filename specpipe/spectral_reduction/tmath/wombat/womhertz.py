def womhertz(hop):
    """converts wavelength to hertz for f_nu spectrum"""
    import numpy as np
    print('NOTE:  The routine expects an f_nu spectrum')
    print('       I will try to guess if the spectrum')
    print("       has been scaled by 1E26 (it shouldn't be)")
    print(' ')
    print('       Check this before believing any result')
    print(' ')


    if (np.mean(hop[0].flux) > 0.00001):
        hop[0].flux=hop[0].flux*1.e-26

    hop[0].wave = hop[0].wave*1.e-10
    hop[0].wave = 2.99792458e8/hop[0].wave
    #FIX head

    print('Active spectrum now in hertz vs. f_nu')
    return hop

