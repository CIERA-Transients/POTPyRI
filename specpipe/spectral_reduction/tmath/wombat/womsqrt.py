def womsqrt(hop):
    """square root of flux"""
    import numpy as np
    hop[0].flux=np.sqrt(hop[0].flux)
    print('\nSquare root of spectrum now active\n')
    return hop

