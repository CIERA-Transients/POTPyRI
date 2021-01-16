def gauss_cont(x, *p):
    import numpy as np
    A, mu, sigma, m, b = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + m*x + b


