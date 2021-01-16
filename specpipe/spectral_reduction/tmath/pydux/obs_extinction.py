def obs_extinction(observat):
    import numpy as np
    observatory_heights={'keck':4160,
                         'gemini-north': 4213.4,
                         'gemini-south': 2737.0,
                         'gemini-n': 4213.4,
                         'gemini-s': 2737.0,
                         'soar': 2737.0,
                         'kpno': 2064.0,
                         'lick': 1285.0,
                         'palomar': 1706.0 ,
                         'mcdonald': 2075.0,
                         'flwo': 2320.0,
                         'mmto': 2600.0,
                         'sso': 1149.0,
                         'vlt': 2635.0,
                         'lco': 2282.0,
                         'lco-imacs': 2282.0,
                         'ctio': 2399.0
                         }
    if (observat in observatory_heights):
        height=observatory_heights[observat]
    else:
        print('OBSERVATORY UNKNOWN!!!!!!!!')
        height=0.0
    # 8300 meters is scale height of troposphere
    sitefactor=np.exp(-1.0*height/8300.0)
    return sitefactor

