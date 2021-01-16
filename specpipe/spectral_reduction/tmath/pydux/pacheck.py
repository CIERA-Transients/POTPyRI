def pacheck(head):
    from tmath.wombat.getch import getch
    from astropy.io import fits
    observat=head['OBSERVAT'].strip().lower()
    airmass=float(head['AIRMASS'])
    opt_pa=float(head['OPT_PA'])
    date=head['DATE-OBS'].strip()
    # modern date-obs is YYYY-MM-DD, old version used DD/MM/YY
    if (date[4] == '-'):
        year=int(date[0:4])
        month=int(date[5:7])
    else:
        year=int(date[6:8])
        month=int(date[3:5])
        #last of the Y2K bugs, good until 2050
        if (year > 50):
            year=year+1900
        else:
            year=year+2000
    # all this date stuff is to check for the point when FLWO changed
    # header keys for position angle
    if (year >= 2003) or ((year == 2002) and (month >= 12)):
        flwokey='ROTANGLE'
    else:
        flwokey='POSANGLE'
    posangkey={'keck': ('ROTPOSN',90.),
               'lick': ('TUB', 0.0),
               'palomar': ('CASSPA', 0.0),
               'flwo': (flwokey, 0.0),
               'mmto': ('POSANG', 0.0),
	       'kpno' : ('ROTANGLE', 0.0),
               'sso': ('POSANG', 0.0),
               'vlt': ('POSANG', 0.0),
               'lco': ('ROTANGLE', 60.3),
               'gemini-north': ('PA', 0.0),
               'gemini-south': ('PA', 0.0),
               'gemini-n': ('PA', 0.0),
               'gemini-s': ('PA', 0.0),
               'ctio':  ('PA', 0.0),
               'lapalma': ('FIELD',-90.),
               'soar': ('POSANGLE', 0.0)
               }
    if (observat in posangkey):
        pa=float(head[posangkey[observat][0]])+posangkey[observat][1]
    elif (observat == 'mcdonald'):
        pa=float(head['PARANGLE'])-float(head['RHO_OFFS'])
    else:
        pa=1000.0
    diffpa=abs(opt_pa-pa)
    if (pa >= 999):
        diffpa=0
    diffpa=diffpa % 180
    if (airmass > 1.1) and (((diffpa > 10) and (diffpa < 170))):
        print('************WARNING***************')
        print('Observed position angle: {}'.format(pa))
        print('Optimal parallactic angle: {}'.format(opt_pa))
        print('Airmass: {}'.format(airmass))
        print(' ')
        print('Relative flux may be compromised')
        print('Hit any key to indemnify me against any and all')
        print('problems that may arise from this')
        any=getch()
    head.set('OBS_PA',pa,'observed position angle')
    return head

