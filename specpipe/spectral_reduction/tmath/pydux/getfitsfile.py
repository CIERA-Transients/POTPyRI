def getfitsfile(name,fitstype, gratcode=None):
    from astropy.io import fits
    import glob
    done=False
    while (not done):
        if gratcode:
            listfile = glob.glob('*.fits')
            for f in listfile:
                if 'ex' in f and gratcode in f:
                    inputfile = f
                    fitsdata=fits.open(inputfile)
                    break
            done=True
        else:
            inputfile=input('Name of fits file for the {}? ({} added if necessary) '.format(name,fitstype))
            inputfile=inputfile.strip()
            if (fitstype not in inputfile):
                inputfile=inputfile+fitstype
        try:
            fitsdata=fits.open(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            fitsdata.close()
            done=True
    return inputfile

