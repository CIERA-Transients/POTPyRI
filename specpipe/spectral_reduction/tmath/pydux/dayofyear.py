def dayofyear(year,month,day):
    """
    Calculates day of the year from calendar date

    :param year: int 
    :param month: int
    :param day: float
    :returns: numd
    :rtype: int

    """
    from tmath.pydux.leapyear import leapyear
    islyr=leapyear(year)
    numd=int((275.0*month)/9.0)-(2-islyr)*int((month+9)/12.0) + day -30
    return numd
