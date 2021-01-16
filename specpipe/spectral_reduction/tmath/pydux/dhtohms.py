def dhtohms(dh):
    rah=int(dh)
    ram=int(abs((dh-rah))*60.0)
    ras=((((abs(dh-rah))*60.0)-ram)*60.0)
    return rah, ram, ras
