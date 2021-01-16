def wave_telluric(wave, mode):
    import numpy as np
    loc1=np.logical_and(wave >= 3216.,wave <= 3420.)
    loc2=np.logical_and(wave >= 5600.,wave <= 6050.)
    loc3=np.logical_and(wave >= 6250.,wave <= 6360.)
    loc4=np.logical_and(wave >= 6450.,wave <= 6530.)
    loc5=np.logical_and(wave >= 6840.,wave <= 7410.)
    loc6=np.logical_and(wave >= 7560.,wave <= 8410.)
    loc7=np.logical_and(wave >= 8925.,wave <= 9900.)
    if (mode == 'high'):
        loc=loc1+loc2+loc3+loc4+loc5+loc6+loc7
    else:
        loc=loc1+loc3+loc5+loc6+loc7
    return loc

