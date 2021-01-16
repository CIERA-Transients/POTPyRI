def womcommanyselector(hop):
    """select combine many equally"""
    from tmath.wombat.womcom import womcom
    from tmath.wombat import HOPSIZE
    hop=womcom(hop,HOPSIZE+1,False)
    return hop

