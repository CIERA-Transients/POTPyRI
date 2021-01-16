def womcommanywselector(hop):
    """select combine many with weights"""
    from tmath.wombat.womcom import womcom
    hop=womcom(hop,HOPSIZE+1,True)
    return hop

