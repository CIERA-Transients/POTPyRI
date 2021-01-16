def womvelcho(hop):
    """select region of spectrum, put onto velocity scale"""
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    from tmath.wombat.inputter import inputter
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.womscipyrebin import womscipyrebin
    wave=hop[0].wave.copy()
    flux=hop[0].flux.copy()
    var=hop[0].var.copy()
    nwave,nflux,mode=womwaverange(wave,flux,'none')
    wavebindex=womget_element(wave,nwave[0])
    waverindex=womget_element(wave,nwave[-1])
    nvar=var[wavebindex:waverindex+1].copy()
    binvec=np.arange(len(nwave))
    wavelog=np.log(nwave[0])+((np.log(nwave[-1])-np.log(nwave[0]))/len(nwave))*binvec
    wavelog=np.exp(wavelog)
    fluxlog=womscipyrebin(nwave,nflux,wavelog)
    varlog=womscipyrebin(nwave,nvar,wavelog)
    zp=wavelog[0]-1.0
    while (zp < wavelog[0]) or (zp > wavelog[-1]):
        zp=inputter('Enter the zero point for velocity (in angstroms): ','float',False)
    indexzp=womget_element(wavelog, zp)
    print(' ')
    logging.info('Zero at bin {}'.format(indexzp))
    logging.info('with lambda {}'.format(wavelog[indexzp]))
    print(' ')
    z=(wavelog - zp)/zp
    square=(z+1)*(z+1)
    wavekmsrel=((square-1.)/(square+1.))*299792.458
    kmsperbin=np.zeros(len(nwave))
    for i in range(1,len(nwave)):
        kmsperbin[i]=2.0*2.99792458e5*(wavelog[i]-wavelog[i-1])/\
                      (wavelog[i]+wavelog[i-1])
    kmsmean=np.mean(kmsperbin[1:])
    logging.info('Average km/s per bin: {}'.format(kmsmean))
    logging.info('km/s at bin 1: {}'.format(kmsperbin[1]))
    logging.info('km/s at bin n: {}'.format(kmsperbin[-1]))
    print(' ')
    wavekms=kmsmean*binvec+kmsmean/2.0
    indexzp=womget_element(wavekms,zp)
    wavekms=wavekms-wavekms[indexzp]
    hop[0].wave=wavekmsrel.copy()
    hop[0].flux=fluxlog.copy()
    hop[0].var=varlog.copy()
    return hop

