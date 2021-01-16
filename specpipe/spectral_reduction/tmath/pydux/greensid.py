def greensid(jd, hour, minute, second):
    import numpy as np
    degrad=180./np.pi
    from tmath.pydux.nutation import nutation
    from tmath.pydux.obliquity import obliquity
#    needs JD at 0h UT on given date (end in 0.5), hour/minute second
#    is time past 0h first finds sidereal time at 0h, then adds 
#    this is UT!, not ET/TT
    t=(jd-2451545.0)/36525.0
    sid0=100.46061837+36000.770053608*t+0.000387933*t*t-(t*t*t)/38710000.0
    sid0=np.fmod(sid0,360.0)
    if (sid0 < 0):
        sid0=sid0+360.0
    sid0=sid0/15.0
#     add hours here
    dech=hour+(minute+(second/60.0))/60.0
    sidcor=dech*1.00273790935
    gsid=sid0 + sidcor
    
    dpsi, deps=nutation(jd)
    epsilon=obliquity(jd)
    sidcor=dpsi*np.cos(epsilon/degrad)/60.0/60.0/15.0
    gsid=gsid+sidcor
    while gsid > 24.0:
        gsid=gsid-24.0

    return gsid
