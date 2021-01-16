def filterconv(wavematch,fluxmatch,splineresult,filtw,filtt):
    """convolve filter, find eff lambda"""
    import numpy as np
    wavematchplusone=np.roll(wavematch,-1)
    wavematchminusone=np.roll(wavematch,1)
    wavematchplusone[-1]=wavematchplusone[-2] + wavematchplusone[-2] - \
                          wavematchplusone[-3]
    wavematchminusone[0]=wavematchminusone[1]-(wavematchminusone[2] - \
                          wavematchminusone[1])
    totflux=np.sum(fluxmatch*splineresult*(wavematchplusone-wavematchminusone)/2.0)
    normalization=np.sum(splineresult*(wavematchplusone-wavematchminusone)/2.0)
    flux=totflux/normalization
    totwave=np.sum(wavematch*fluxmatch*splineresult*(wavematchplusone-wavematchminusone)/2.0)
    efflambda=totwave/totflux
    filtwmatchplusone=np.roll(filtw,-1)
    filtwmatchminusone=np.roll(filtw,1)
    filtwmatchplusone[-1]=filtwmatchplusone[-2] + filtwmatchplusone[-2] - \
                          filtwmatchplusone[-3]
    filtwmatchminusone[0]=filtwmatchminusone[1]-(filtwmatchminusone[2] - \
                          filtwmatchminusone[1])
    com=np.sum(filtt*(filtwmatchplusone-filtwmatchminusone)/2.0)
    coverage=normalization/com
    return flux,totflux,efflambda,coverage


    

