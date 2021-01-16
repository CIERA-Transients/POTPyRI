def womrmsfits(hop):
    """read multispec fits files"""
    from astropy.io import fits
    from tmath.wombat.inputter import inputter
    from tmath.wombat import HOPSIZE
    from tmath.wombat.getmswave import getmswave
    from tmath.pydux.ismultispec import ismultispec
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
    if (not ismultispec(inputfile)):
        print('Not a multispec file!')
        return hop
    if (len(fitsdata) > 1):
        print('Multi-Extension File!  Bailing out!')
        return hop
    flux=fitsdata[0].data
    flux=flux.astype(float)
    fhdr=fitsdata[0].header
    objectname=fhdr['OBJECT']
    if (not isinstance(objectname, str)):
        objectname=inputfile
    print('Object is {}'.format(objectname))
    #FIX for aps with different wavelength solutions, you need to look
    # at other keywords, but I don't have an example file right now
    # this is a long-term fix
    num_apertures=flux.shape[1]
    print('\nThere are {} apertures in this file.'.format(num_apertures))
    ap_choice=0
    if (num_apertures > 1):
        while (ap_choice < 1) or (ap_choice > num_apertures):
            ap_choice=inputter('Which aperture do you want? ','int',False)
    else:
        ap_choice=1
    ap_choice=ap_choice - 1
    wave=getmswave(fhdr,ap_choice)
    wdelt=wave[1]-wave[0]
    num_bands=flux.shape[0]
    print('\nThere are {} bands in this aperture.'.format(num_bands))
    # clean up keywords
    delkeylist=['WAT0_001', 'WAT1_001', 'WAT2_001', 'WAT0_002', \
                'WAT1_002', 'WAT2_002', 'WAT3_001', 'WAT2_003', \
                'CTYPE1', 'CTYPE2', 'CTYPE3', 'CD1_1', 'CD2_2', \
                'CD3_3', 'LTM1_1', 'LTM2_2', 'LTM3_3', 'WCSDIM']
    for k in delkeylist:
        try:
            fhdr.remove(k)
        except KeyError:
            pass
    fhdr.set('CRPIX1', 1)
    fhdr.set('CRVAL1', wave[0])
    fhdr.set('CDELT1', wdelt)
    fhdr.set('CTYPE1', 'LINEAR')
    if (num_bands > 2):
        var=flux[3,ap_choice,:]**2
    for i in range(0,num_bands):
        hopnum=0
        while (hopnum < 1) or (hopnum > HOPSIZE):
            hopnum=inputter('Store band {} in which hopper? '.format(i+1),'int',False)
        hop[hopnum].wave=wave.copy()
        hop[hopnum].flux=flux[i,ap_choice,:].copy()
        hop[hopnum].obname=objectname
        # OK, strictly speaking this var only applies to the optimal extraction
        # in band 1 (0 in python index).  We still have access to 'var' in
        # the hopper we put it in.
        hop[hopnum].var=var.copy()
        hop[hopnum].header=fhdr
    fitsdata.close()
    return hop

