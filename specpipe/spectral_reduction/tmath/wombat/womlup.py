def womlup(hop):
    """scales to asinh luptitudes"""
    import numpy as np
    print('\nTaking asinh of flux data\n')

    hop[0].flux=np.arcsinh(hop[0].flux)

    return hop

