def julday(year,month,day):
    """
    Compute Julian day

    Parameters
    ----------
    year : int
    month : int
    day : float

    Returns
    -------
    jd : float
        Julian day at year, month, day

    Notes
    computes julian day from year, month, day day can be a float, if day
    is integer, returns JD of 0h UT, so the value of jd ends in 0.5 as
    the day is defined to start at 12h UT
    """
    a=0.0
    b=0.0

    if month <= 2:
        year=year-1
        month=month+12

    if year > 1582:
        a=int(year/100.0)
        b=2-a+int(a/4.0)

    jd=int(365.25*(year+4716)) + int(30.6001*(month+1)) + day+b-1524.5
    return jd
