#parameter file for MOSFIRE/Keck

def run_wcs():
    return True

def wcs_extension():
    return 0

def pixscale():
    return 0.1799

def ref_pix():
    return 1033.0, 1035.0

def WCS_keywords(): #WCS keywords
    WAT0_001 = 'system=image'
    WAT1_001 = 'wtype=tnx axtype=ra'
    WAT1_002 = ''
    WAT1_003 = ''
    WAT1_004 = ''
    WAT1_005 = ''
    WAT2_001 = 'wtype=tnx axtype=dec'
    WAT2_002 = ''
    WAT2_003 = ''
    WAT2_004 = ''
    WAT2_005 = ''
    return WAT0_001, WAT1_001, WAT1_002, WAT1_003, WAT1_004, WAT1_005, WAT2_001, WAT2_002, WAT2_003, WAT2_004, WAT2_005