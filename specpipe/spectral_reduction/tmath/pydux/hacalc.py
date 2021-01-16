def hacalc(ra,longitude,gsid):
    """
    Calculate hour angle

    :param ra: float
    :param longitude: float
    :param gsid: float
    :returns: ha
    :rtype: float

    """
    ha=gsid-longitude-ra
    return ha
