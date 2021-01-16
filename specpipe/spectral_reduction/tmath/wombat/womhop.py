def womhop(hop):
    import copy
    from tmath.wombat.inputter import inputter
    from tmath.wombat import HOPSIZE
    hopnum=0
    while (hopnum < 1) or (hopnum > HOPSIZE):
        hopnum=inputter('Store in which hopper: ','int',False)
    hop[hopnum]=copy.deepcopy(hop[0])
    return hop

