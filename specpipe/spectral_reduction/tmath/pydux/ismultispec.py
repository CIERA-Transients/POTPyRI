def ismultispec(inputfile):
    """tests if fits file is multispec using presence of NAXIS3 
    header keyword and value of 3 or 4"""
    from astropy.io import fits

    infits=fits.open(inputfile)
    try:
        naxis3 = infits[0].header['NAXIS3']
    except KeyError:
        return False
    if (naxis3 == 3 or naxis3 == 4):
        return True
    else:
        return False
