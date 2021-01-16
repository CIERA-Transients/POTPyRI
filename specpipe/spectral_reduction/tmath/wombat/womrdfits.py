def womrdfits(hop):
    from astropy.io import fits
    import numpy as np
    done=False
    while (not done):
        inputfile=input('Name of fits file to be read? (.fits added if necessary) ')
        inputfile=inputfile.strip()
        if (inputfile == ''):
            return hop
        if ('.fits' not in inputfile):
            inputfile=inputfile+'.fits'
        try:
            fitsdata=fits.open(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            done=True
    if (len(fitsdata) > 1):
        print('Multi-Extension File!  Bailing out!')
        return hop

    crval1=fitsdata[0].header['CRVAL1']
    try:
        wdelt=fitsdata[0].header['CDELT1']
    except KeyError:
        wdelt = 0.0000001
    if (wdelt < 0.000001):
        wdelt=fitsdata[0].header['CD1_1']
    try:
        objectname=fitsdata[0].header['OBJECT']
        if (not isinstance(objectname, str)):
            objectname=inputfile
    except KeyError:
        objectname = inputfile
    print('Object is {}'.format(objectname))
    flux=fitsdata[0].data
    # if two vectors, then this is flux/sigma from pydux (let's hope)
    if (len(flux.shape) == 2):
        var=flux[:,1]
        var=var.astype(float)
        var=var*var
        flux=flux[:,0]
        flux=flux.astype(float)
    else:
        flux=flux.astype(float)
        var=np.ones(flux.shape)
    wave=np.arange(1,len(flux)+1)*wdelt+crval1
    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].obname=objectname
    hop[0].var=var.copy()
    hop[0].header=fitsdata[0].header
    fitsdata.close()
    return hop

