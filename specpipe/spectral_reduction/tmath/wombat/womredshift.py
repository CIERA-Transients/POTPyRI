def womredshift(hop):
    """remove redshift from spectrum"""
    import numpy as np
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.inputter import inputter
    from tmath.wombat.yesno import yesno
    from tmath.wombat.womrebindriver import womrebindriver
    light_speed=2.99792458e5
    print('Current A/pix is {}'.format(hop[0].wave[1]-hop[0].wave[0]))
    print('\nWavelength range: {} to {}\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    answer=inputter_single('Remove redshift in (z) or (k)m/s? (z/k) ','zk')
    z=inputter('Enter the redshift (positive is away from us): ','float',False)
    if (answer == 'k'):
        z=np.sqrt((1.0+z/c)/(1.0-z/c)) - 1.0
    hop[0].wave=hop[0].wave/(1.0+z)
    print('\nNew wavelength range: {} to {}\n'.format(hop[0].wave[0],hop[0].wave[-1]))
    print('Rebin spectrum?')
    rebin=yesno('n')
    if (rebin == 'y'):
        hop=womrebindriver(hop,'scipy')
    else:
        pass
        #FIX header
    return hop

