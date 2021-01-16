def getmswave(head,aperture):
    import numpy as np
    # have to add 1 as IRAF is not 0 indexed
    specstr='spec'+str(aperture+1)
    watflag=False
    for i,hcard in enumerate(head.cards):
        if ('WAT2' in hcard[0]):
            if (specstr in hcard[1]):
                watstring=hcard[1]
                specloc=watstring.find(specstr)
                watstringsub=watstring[specloc:]
                quoteloc=watstringsub.find('"')
                specsub=watstringsub[quoteloc+1:]
                specsubsplit=specsub.split()
                crval=float(specsubsplit[3])
                cdelt=float(specsubsplit[4])
                npix=float(specsubsplit[5])
                watflag=True
    if (not watflag):
        crval=float(head['CRVAL1'])
        cdelt=float(head['CD1_1'])
        npix=float(head['NAXIS1'])
    wave=np.arange(npix)*cdelt + crval
    return wave

