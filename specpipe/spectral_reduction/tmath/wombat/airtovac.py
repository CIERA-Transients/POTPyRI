def airtovac(wave):
    """
    AIR to VAC from Goddard IDL library

    PURPOSE:
       Convert air wavelengths to vacuum wavelengths 
    EXPLANATION:
       Wavelengths are corrected for the index of refraction of air under 
       standard conditions.  Wavelength values below 2000 A will not be 
       altered.  Uses the IAU standard for conversion given in Morton 
       (1991 Ap.J. Suppl. 77, 119)

    INPUT/OUTPUT:
       WAVE - Wavelength in Angstroms, numpy vector
               WAVE should be input as air wavelength(s)
               returns vacuum

    EXAMPLE:
       If the air wavelength is  W = 6056.125 (a Krypton line), then 
       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019

    METHOD:
       See Morton (Ap. J. Suppl. 77, 119) for the formula used
    """
    import numpy as np
    wavenumber_squared=(10000./wave)**2

    factor=1. + 6.4328e-5 + (2.94981e-2)/(146. - wavenumber_squared) \
            + (2.5540e-4)/(41.-wavenumber_squared)
    type(wave)
    vaclocs=np.where(wave < 2000)
    factor[vaclocs]=1.
    wave=wave*factor
    return wave


