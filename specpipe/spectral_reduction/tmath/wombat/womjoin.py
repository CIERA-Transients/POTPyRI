def womjoin(hop):
    """join two spectra that abut in wavelength"""
    import numpy as np
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat import HOPSIZE
    print('\nThis will join two hoppers (with no overlap)\n')
    hopnum1=0
    hopnum2=0
    while (hopnum1 < 1) or (hopnum1 > HOPSIZE):
        hopnum1=inputter('Enter first hopper: ','int',False)
    while (hopnum2 < 1) or (hopnum2 > HOPSIZE):
        hopnum2=inputter('Enter first hopper: ','int',False)
    if (hop[hopnum1].wave[0] > hop[hopnum2].wave[0]):
        hopnum1,hopnum2=hopnum2,hopnum1
    wdelt1=hop[hopnum1].wave[1]-hop[hopnum1].wave[0]
    wdelt2=hop[hopnum2].wave[1]-hop[hopnum2].wave[0]
    # check if wavelength dispersion same
    if (abs(wdelt1 -wdelt2) > 0.00001):
        print('Spectra do not have same Angstrom/pixel')
        print('Blue side: {}'.format(wdelt1))
        print('Red side: {}'.format(wdelt2))
        return hop
    #check if spectra abut
    if (abs(hop[hopnum2].wave[0] - (hop[hopnum1].wave[-1] + wdelt1)) \
        > 0.00001):
        print('\nSpectra do not abut\n')
        print('Red end of blue: {}'.format(hop[hopnum1].wave[-1]))
        print('Blue end of red: {}\n'.format(hop[hopnum2].wave[0]))
        return hop
    print('Joining from {} to {}'.format(hop[hopnum1].wave[-1], \
                                         hop[hopnum2].wave[0]))
    hopout=0
    while (hopout < 1) or (hopout > HOPSIZE):
        hopout=inputter('Enter hopper to store combined spectrum: ','int',False)
    hop[hopout].wave=np.concatenate([hop[hopnum1].wave,hop[hopnum2].wave])
    hop[hopout].flux=np.concatenate([hop[hopnum1].flux,hop[hopnum2].flux])
    hop[hopout].var=np.concatenate([hop[hopnum1].var,hop[hopnum2].var])
    hop[hopout].obname=hop[hopnum1].obname
    hop[hopout].header=hop[hopnum1].header
    logging.debug('Files: {} and {} joined from {} to {}'.format(hop[hopnum1].obname, hop[hopnum2].obname, hop[hopnum1].wave[-1],hop[hopnum2].wave[0]))
    #FIX header
    return hop

