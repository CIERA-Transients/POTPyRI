def parsehourangle(hainput):
    """ return hour angle in decimal hours
    some headers have HA as hours in decimal,
    some have HA in hours as sexigesimal
    """
    from tmath.pydux.dms_string_to_dd import dms_string_to_dd

    if isinstance(hainput, float):
        return hainput
    elif isinstance(hainput, str):
        # hasign to avoid any -0:01:02 problems
        #hasign=1.
        #if (hainput[0] == '-'):
        #    hasign=-1.
        #halist=hainput.split(':')
        #ha=(hasign*float(halist[0])+(float(halist[1])+float(halist[2])/60.)/60.)
        #ha=hasign*ha
        ha=dms_string_to_dd(hainput)
        return ha
    else:
        return 0
