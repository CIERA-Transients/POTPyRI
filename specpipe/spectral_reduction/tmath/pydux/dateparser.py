def dateparser(date):
    """parse year, month, day from DATE-OBS keyword"""

    date=date.strip()
    # get year,month,day from ISO Y2K format YYYY-MM-DDThh:mm:ss.ssss
    # the time following T is optional, so we get T from UTMIDDLE
    if (date[4] == '-'):
        year=int(date.split('-')[0])
        month=int(date.split('-')[1])
        day=int(date.split('-')[2].split('T')[0])
    else:
        # Old date format DD/MM/YY
        year=int(date[6:8])
        month=int(date[3:5])
        day=int(date[0:2])
        # try to catch old format written after 2000
        # good until 2050, by which time no one should
        # be doing this
        if (year > 50):
            year=year+1900
        else:
            year=year+2000
    return year, month, day
