def gauss(x, *p):
    import numpy as np
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

