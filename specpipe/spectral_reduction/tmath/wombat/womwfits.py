def womwfits(hop):
    """write fits file"""
    from astropy.io import fits
    from tmath.wombat.inputter import inputter
    print('Object is {}'.format(hop[0].obname))
    fileout=inputter('Enter the name for the output file: ' + \
                     '(.fits will be added,if necessary)','string',False)
    outhdu=fits.PrimaryHDU(hop[0].flux)
    hdul=fits.HDUList([outhdu])
#    hdul[0].header.set('EXTEND','F')
    hdul[0].header.set('OBJECT',hop[0].obname)
    if ('.fits' not in fileout):
            fileout=fileout+'.fits'
    if (len(hop[0].header) > 1):
        hdul[0].header=hop[0].header.copy()
  #  hdul[0].header['SIMPLE']=('T','Written by Wombat')
    hdul[0].header.set('CRPIX1',1)
    hdul[0].header.set('CRVAL1',hop[0].wave[0])
    hdul[0].header.set('CDELT1',hop[0].wave[1]-hop[0].wave[0])
    hdul[0].header.set('CTYPE1','LINEAR')
    hdul[0].header.set('COMMENT','Written by Wombat')
    hdul.writeto(fileout,overwrite=True)
    hdul.close()
    return hop

