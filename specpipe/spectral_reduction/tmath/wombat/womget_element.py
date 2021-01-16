def womget_element(wavearr,wavelength):
    """get index of given wavelength"""
    import numpy as np
    index=(np.abs(wavearr-wavelength)).argmin()
    return index

