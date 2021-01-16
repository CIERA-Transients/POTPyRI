def timeparser(ut):
    """parse hour, minute, second from UT string"""
    if ('T' in ut):
        ut=ut.split('T')[1]
    # some UT/UTMIDDLE values are ISO format
    hour=int(ut.split(':')[0])
    minute=int(ut.split(':')[1])
    second=float(ut.split(':')[2])
    return hour, minute, second
