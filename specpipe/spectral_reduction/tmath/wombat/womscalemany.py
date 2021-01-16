def womscalemany(hop):
    """scale many hoppers to one"""
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat import HOPSIZE
    print('This will scale many hoppers to one')
    done = False
    while (not done):
        hopnum1=inputter('Enter fiducial hopper: ','int',False)
        if (hopnum1 < 1) or (hopnum1 > HOPSIZE):
            print('Hopper number must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    hoplist=[]
    hopnum2=0
    while (hopnum2 != 99):
        done = False
        while (not done):
            hopnum2=inputter('Enter another hopper: (99 to end) ','int',False)
            if (hopnum2 != 99):
                if (hopnum2 < 1) or (hopnum2 > HOPSIZE):
                    print('Hopper numbers must be in range 1-{}'.format(HOPSIZE))
                elif (len(hop[hopnum2].wave) == 0):
                    print('Nothing in that hopper')
                else:
                    hoplist.append(hopnum2)
                    done=True
                if (hopnum2 > 0) and (hopnum2 < 21) and \
                   (len(hop[hopnum2].wave) != 0):
                    if (hop[hopnum1].wave[0] != hop[hopnum2].wave[0])  \
                       or (hop[hopnum1].wave[1] != hop[hopnum2].wave[1]) \
                       or (hop[hopnum1].wave[-1] != hop[hopnum2].wave[-1]):
                        print('Hoppers to not have the same wavelength scale!')
                        return hop
            else:
                done=True
    print('\nSpectra run from {} to {}'.format(hop[hopnum1].wave[0], \
                                             hop[hopnum1].wave[-1]))

    plt.cla()
    plt.plot(hop[hopnum1].wave,hop[hopnum1].flux,drawstyle='steps-mid', \
             color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[hopnum1].obname)
    print('\nPlotting fiducial spectrum')
    print('\nChoose region to compute average')
    wave,flux,mode=womwaverange(hop[hopnum1].wave,hop[hopnum1].flux,'none')
    indexblue=womget_element(hop[hopnum1].wave,wave[0])
    indexred=womget_element(hop[hopnum1].wave,wave[-1])
    avg1=np.mean(hop[hopnum1].flux[indexblue:indexred+1])
    print('Fiducial spectrum mean in range: {}'.format(avg1))
    for i,_ in enumerate(hoplist):
        print('{} {}'.format(i,hoplist[i]))
        avgloop=np.mean(hop[hoplist[i]].flux[indexblue:indexred+1])
        print('Hopper {} mean in range: {}'.format(hoplist[i], avgloop))
        hop[hoplist[i]].flux=hop[hoplist[i]].flux*avg1/avgloop
        hop[hoplist[i]].var=hop[hoplist[i]].var*avg1/avgloop
        logging.debug('File: {} scaled by {} to match {}'.format(hop[hoplist[i]].obname, avg1/avgloop, hop[hopnum1].obname))
    return hop
                     
