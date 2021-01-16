def crosscorrelation(x, y, maxlag):
    """
    https://stackoverflow.com/questions/30677241/how-to-limit-cross-correlation-window-width-in-numpy?rq=1
    User Warren Weckesser
    https://stackoverflow.com/users/1217358/warren-weckesser

    Cross correlation with a maximum number of lags.

    `x` and `y` must be one-dimensional numpy arrays with the same length.

    This computes the same result as
        numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag]

    The return vaue has length 2*maxlag + 1.
    np.argmax of the return correlation vector gives the value to add to y
    np.argmax-maxlag
    """
    import numpy as np
    from numpy.lib.stride_tricks import as_strided
    # need to normalize here so that spectra with wildly different fluxes will
    # correlate properly (zero-normalized cross correlation)
    xx = (x-np.mean(x))/np.std(x)
    yy = (y-np.mean(y))/np.std(y)/len(y)
    py = np.pad(yy.conj(), 2*maxlag, mode='constant')
    T = as_strided(py[2*maxlag:], shape=(2*maxlag+1, len(yy) + 2*maxlag),
                   strides=(-py.strides[0], py.strides[0]))
    px = np.pad(xx, maxlag, mode='constant')
    return T.dot(px)

