def ddaytodhms(dday):
    """convert decimal day to day, hour, minute,second"""
    aday=int(dday)
    hours=(dday-aday)*24.
    ahour=int(hours)
    minutes=(hours-ahour)*60.
    aminute=int(minutes)
    asecond=(minutes-aminute)*60.
    return aday, ahour, aminute, asecond
