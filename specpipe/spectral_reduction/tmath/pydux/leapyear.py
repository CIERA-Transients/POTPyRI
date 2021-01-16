def leapyear(year):
    """
    Evaluates if a given year is a leap year

    :param year: int
    :returns: islyr
    :rtype: int

    Notes ----- Could be boolean, but I use the value in some
    calculations--this still might be ok, but needs checking

    """
    is_lyr = 0
    if ((year/4)*4. == year) and (((year/100)*100.0 != year) or \
            ((year/400)*400.0 == year)):
        is_lyr = 1

    return is_lyr
