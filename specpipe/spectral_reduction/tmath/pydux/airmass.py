def airmass(el):
    """
    Calculate airmass from elevation

    :param el: float
    :returns: airmass
    :rtype: float
   
    Notes
    -----
    Airmass calculation:  the old log program used
   airmass=(1/cos(zdrad))-0.0018167*((1/cos(zdrad))-1)-
   0.002875*((1/cos(zdrad))-1)^2-0.0008083*((1/cos(dzrad))-1)^3, which
   is the empirical fit of Hardie, R. H. 1962, in Astronomical
   Techniques. Hiltner, W. A., ed. Chicago: University of Chicago
   Press.  This blows up at zenith distances of 85 degrees and up.
   It isn't really critical for this program, but I'm going to
   substitute the formulation of Rozenberg (Rozenberg,
   G. V. 1966. Twilight: A Study in Atmospheric Optics. New York:
   Plenum Press, 160. Translated from the Russian by R. B. Rodman).
   This does much better at higher airmasses and is entirely
   consistent at lower zenith distances (better than 0.1% agreement).
   there are plenty of other good fits at high airmass, but they
   don't all give 1.0 at zenith!
    """
    import numpy as np
    degrad=180.0/np.pi
    zd=90.0-el
    zdrad=zd/degrad
    z=1.0/(np.cos(zdrad)+0.025*np.exp(-11.0*np.cos(zdrad)))
    return z
