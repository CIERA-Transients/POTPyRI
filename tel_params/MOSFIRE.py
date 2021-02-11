#parameter file for MOSFIRE/Keck

def static_mask():
    return './staticmasks/MF.staticmask.fits'

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.1799

def ref_pix():
    # Sorry this ref pix is so crazy.  Haven't figured out how to set this yet.
    return 500174.481085049, 134532.194062595

def WCS_keywords(): #WCS keywords
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tnx axtype=ra lngcor = "3. 4. 4. 2. 34.85092507693425 35.04227'
    WAT2_001= 'wtype=tnx axtype=dec latcor = "3. 5. 5. 2. 34.85092507693425 35.0422'
    WAT1_002= '412406922 9.103584559568529 9.272725992772539 1074.076305826495 -58.'
    WAT1_003= '66970631133986 1.6799195420292 -0.01602888260183174 -127.80486101539'
    WAT1_004= '61 -0.002026891656217711 0. 13.92965631881465 0. -0.5057845760358613'
    WAT2_002= '7412406922 9.103584559568529 9.272725992772539 1129.122011958748 -61'
    WAT2_003= '.95924736013845 1.776895554106438 -0.01701646832378952 0. -133.42572'
    WAT2_004= '85609035 0.01209004508732243 0. 0. 14.46400328733441 0. 0. -0.524309'
    WAT1_005= ' "      '
    WAT2_005= '9455921529 0. 0. "'
    return WAT0_001, WAT1_001, WAT1_002, WAT1_003, WAT1_004, WAT1_005, WAT2_001, WAT2_002, WAT2_003, WAT2_004, WAT2_005
