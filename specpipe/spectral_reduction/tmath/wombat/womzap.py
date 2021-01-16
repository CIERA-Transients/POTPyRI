def womzap(hop):
    """routine to zap outliers (CRs)"""
    import numpy as np
    import matplotlib.pyplot as plt
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.yesno import yesno
    from tmath.wombat.wshow import wshow
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    wshow()
    print('\nRoutine to zap outliers\n')
    done = False
    while (not done):
        nsig=inputter('Enter zapping threshold in sigmas (0=replace all with median): ','float',False)
        print(' ')
        boxsize=inputter('Enter box size for computing statistics: (odd integer < 45) ','int',False)
        if (boxsize < 3) or (boxsize > 45) or (nsig < 0):
            print('Invalid boxsize or sigma')
        else:
            done = True
    if (boxsize % 2 == 0):
        boxsize=boxsize+1
    half=int(boxsize/2.)
    newflux=hop[0].flux.copy()
    if (nsig > 0):
        mode=inputter_single('Use inter(q)uartile or (m)edian variance (q/m)? ','qm')
        if (mode == 'm'):
            for i in range(half,len(newflux)-half):
                sample=newflux[i-half:i+half+1]
                medval=np.median(sample)
                varpix=(sample[half]-medval)**2
                sum=-1.0*varpix
                for k in range(0,boxsize):
                    diff=sample[k]-medval
                    sum=sum+diff*diff
                sum=sum/(boxsize-1)
                if (varpix > nsig*nsig*sum):
                    newflux[i]=medval
        if (mode == 'q'):
            for i in range(half,len(newflux)-half):
                sample=newflux[i-half:i+half+1]
                medval=np.median(sample)
                q25=np.percentile(sample,25.)
                q75=np.percentile(sample,75.)
                qsig=0.7414272*(q75-q25)
                diff=abs(sample[half]-medval)
                if (diff > nsig*qsig):
                    newflux[i]=medval
    if (nsig == 0):
        from scipy.ndimage.filters import median_filter
        newflux=median_filter(newflux,boxsize)
    plt.plot(hop[0].wave,newflux,drawstyle='steps-mid')
    print('Does this zapping look good?')
    good=yesno('y')
    if (good == 'y'):
        hop[0].flux=newflux.copy()
        logging.debug('File {} zapped with sigma {} and boxsize {}'.format\
                      (hop[0].obname,nsig,boxsize))
    else:
        print('OK, active spectrum unchanged')
    #FIX var
    return hop
                                                                           
