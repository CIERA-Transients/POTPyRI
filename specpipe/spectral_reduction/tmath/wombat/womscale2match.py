def womscale2match(hop):
    """scale one spectrum to match another"""
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat import HOPSIZE
    print('This will scale one hopper to match another')
    done = False
    while (not done):
        hopnum1=inputter('Enter first hopper: ','int',False)
        hopnum2=inputter('Enter second hopper: ','int',False)
        if (hopnum1 < 1) or (hopnum1 > HOPSIZE) or (hopnum2 < 1) or \
           (hopnum2 > HOPSIZE):
            print('Hopper numbers must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    if (hop[hopnum1].wave[0] != hop[hopnum2].wave[0]) or \
       (hop[hopnum1].wave[1] != hop[hopnum2].wave[1]) or \
       (hop[hopnum1].wave[-1] != hop[hopnum2].wave[-1]):
        print('Hoppers to not have the same wavelength scale!')
        return hop
    print('\nSpectra run from {} to {}'.format(hop[hopnum1].wave[0], \
                                             hop[hopnum1].wave[-1]))

    plt.cla()
    plt.plot(hop[hopnum1].wave,hop[hopnum1].flux,drawstyle='steps-mid', \
             color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[hopnum1].obname)
    plt.plot(hop[hopnum2].wave,hop[hopnum2].flux,drawstyle='steps-mid', \
             color='r')
    print('\nHopper A in black, hopper B in red')
    print('\nChoose region to compute average')
    wave,flux,mode=womwaverange(hop[hopnum1].wave,hop[hopnum1].flux,'none')
    indexblue=womget_element(hop[hopnum1].wave,wave[0])
    indexred=womget_element(hop[hopnum1].wave,wave[-1])
    avg1=np.mean(hop[hopnum1].flux[indexblue:indexred+1])
    avg2=np.mean(hop[hopnum2].flux[indexblue:indexred+1])
    print('\n Hopper A: {}'.format(avg1))
    print('\n Hopper B: {}'.format(avg2))
    choice=inputter_single('Scale to (A) or (B) ','ab')
    print(' ')
    done = False
    while (not done):
        hopnum3=inputter('Store scaled spec in which hopper? ','int',False)
        if (hopnum3 < 1) or (hopnum3 > HOPSIZE):
            print('Hopper numbers must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    if (choice == 'a'):
        hop[hopnum3].flux=hop[hopnum2].flux*(avg1/avg2).copy()
        print('Scale: {}'.format(avg1/avg2))
        hop[hopnum3].obname=hop[hopnum2].obname
        hop[hopnum3].header=hop[hopnum2].header
        hop[hopnum3].var=hop[hopnum2].var*(avg1/avg2).copy()
        lstring='File: {} scaled by {} to match {}'.format(hop[hopnum2], \
                                                           avg1/avg2, \
                                                           hop[hopnum1])
    else:
        hop[hopnum3].flux=hop[hopnum1].flux*(avg2/avg1).copy()
        print('Scale: {}'.format(avg1/avg2))
        hop[hopnum3].obname=hop[hopnum1].obname
        hop[hopnum3].header=hop[hopnum1].header
        hop[hopnum3].var=hop[hopnum1].var*(avg2/avg1).copy()
        lstring='File: {} scaled by {} to match {}'.format(hop[hopnum1].obname, \
                                                           avg2/avg1, \
                                                           hop[hopnum2].obname)
    hop[hopnum3].wave=hop[hopnum1].wave.copy()
    logging.debug(lstring)
    return hop


