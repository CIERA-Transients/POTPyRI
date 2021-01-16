def womrdhop(hop):
    import copy
    from tmath.wombat import HOPSIZE
    from tmath.wombat.inputter import inputter
    hopnum=0
    while (hopnum < 1) or (hopnum > HOPSIZE):
        hopnum=inputter('Read from which hopper: ','int',False)
    if (len(hop[hopnum].flux) == 0):
        print('Nothing in hopper {}'.format(hopnum))
        print('Active spectrum still {}'.format(hop[0].obname))
    else:
        hop[0]=copy.deepcopy(hop[hopnum])
    print('Object is {}'.format(hop[0].obname))
    return hop
    
