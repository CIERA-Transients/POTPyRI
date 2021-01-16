def jdcnv(year, month, day, hour):
    if (month < 3):
        leap=-1
    else:
        leap=0

    julian=day - 32075 + 1461 * (year +4800 + leap)//4 + \
            367*(month -2 -leap*12)//12 - \
            3*((year + 4900 + leap)//100)//4
    julian=float(julian)+hour/24.0 -0.5
    return julian

