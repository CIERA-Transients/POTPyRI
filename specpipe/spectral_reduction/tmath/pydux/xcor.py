def xcor(x, y, xfactor, maxlag):
    import numpy as np
    import matplotlib.pyplot as plt
    from tmath.pydux.congrid import congrid
    from tmath.pydux.crosscorrelation import crosscorrelation
    npix=len(x)
    xbig=congrid(x,(npix*xfactor,),minusone=True)
    ybig=congrid(y,(npix*xfactor,),minusone=True)
    cc=crosscorrelation(xbig,ybig,maxlag)
    shift=np.argmax(cc)-maxlag
    shift=shift/xfactor
    if (np.argmax(cc) == 0) or (np.argmax(cc) == npix - 1):
        print('Cross-correlation failed!  Try larger maxlag?')
        print('Setting shift to zero...\n')
        shift=0.0
    return shift
    
