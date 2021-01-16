def womcom(hop,num,weights):
    """combine spectra, with or without weights"""
    import logging
    import matplotlib.pyplot as plt
    import numpy as np
    from tmath.wombat.inputter import inputter
    from tmath.wombat import HOPSIZE
    hopnum=[0]
    weight=[0]
    while (hopnum[0] < 1) or (hopnum[0] > HOPSIZE):
        hopnum[0]=inputter('Enter first hopper: ','int',False)
    if (weights):
        weight[0]=inputter('Enter weight for first hopper: ','float',False)
    for i in range(1,num):
        hoploop=0
        weightloop=0
        if (num > 3):
            print('Use hopper 99 to end')
        while (hoploop < 1) or (hoploop > HOPSIZE):
            hoploop=inputter('Enter next hopper: ','int',False)
            if (hoploop == 99):
                break
        if (hoploop == 99):
            break
        if (hop[hopnum[0]].wave[0] != hop[hoploop].wave[0]) or \
       (hop[hopnum[0]].wave[1] != hop[hoploop].wave[1]) or \
       (hop[hopnum[0]].wave[-1] != hop[hoploop].wave[-1]):
            print('Hoppers to not have the same wavelength scale!')
            return hop
        hopnum.append(hoploop)
        if (weights):
            weightloop=inputter('Enter next weight: ','float',False)
        weight.append(weightloop)
    if (weights) and (abs(sum(weight)-1.0) > 0.00001):
        print('Weights do not add to 1.0')
        return hop
    if (not weights):
        weight=[1./len(hopnum)]*len(hopnum)
    newflux=np.zeros(len(hop[hopnum[0]].flux))
    logging.debug('Combining spectra:')
    
    for i in range(0,len(hopnum)):
        newflux=newflux+hop[hopnum[i]].flux*weight[i]
        logging.debug('Combining {} with weight {}'.format(hop[hopnum[i]].obname,\
                                                           weight[i]))
    hopout=0
    while (hopout < 1) or (hopout > HOPSIZE):
        hopout=inputter('Enter hopper to store combined spectrum: ','int',False)
    hop[hopout].wave=hop[hopnum[0]].wave.copy()
    hop[hopout].flux=newflux.copy()
    hop[hopout].obname=hop[hopnum[0]].obname
    hop[hopout].header=hop[hopnum[0]].header
    hop[hopout].var=hop[hopnum[0]].var.copy()
    plt.cla()
    plt.plot(hop[hopout].wave,hop[hopout].flux,drawstyle='steps-mid',color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[hopout].obname)
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    print('\nPlotting combined spectrum in black, components in color\n')
    for i in range(0,len(hopnum)):
        plt.plot(hop[hopnum[i]].wave,hop[hopnum[i]].flux,drawstyle='steps-mid')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    
    #FIX var
    return hop

