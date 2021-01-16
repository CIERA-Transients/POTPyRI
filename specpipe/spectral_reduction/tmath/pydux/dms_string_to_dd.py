def dms_string_to_dd(dmsst):
    """ takes d/h:m:s in string form and returns decimal degrees/hours 
    needs sign part for -00:01:02 issues
    """
    sign=1.
    dmsst=dmsst.strip()
    if ('-' in dmsst):
        sign=-1.
    dmsstlist=dmsst.split(':')
    dd=(sign*float(dmsstlist[0])+(float(dmsstlist[1])+float(dmsstlist[2])/60.)/60.)
    dd=dd*sign
    return dd
