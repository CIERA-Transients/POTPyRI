def calday(jd):
    """
    Computes calendar day given Julian day, day can be fractional

    :param jd: Julian day
    :returns: year, month, day
    :rtype: int, int, float

    """
    ljd=jd+0.5
    z=int(ljd)
    f=ljd-z
    if z < 2299161:
        a=z
    else:
        alpha=int(((z-1867216.25)/(36524.25)))
        a=z+1.0+alpha-int(alpha/4.0)

    b=a+1524.0
    c=int((b-122.1)/365.25)
    d=int(365.25*c)
    e=int((b-d)/30.6001)

    day = b-d-int(30.6001*e)+f
    if e < 13.5: 
        month=int(e)-1
    else:
        month=int(e)-13
      
    if month > 2.5:
        year=int(c)-4716
    else:
        year=int(c)-4715     

    return year, month, day
