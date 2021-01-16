def womscale(hop):
    """scale a spectrum by multiplicative value"""
    from tmath.wombat.inputter import inputter
    factor=inputter('Enter multiplicative scale factor: ','float',False)
    hop[0].flux=hop[0].flux*factor
    return hop

