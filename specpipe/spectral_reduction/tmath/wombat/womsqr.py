def womsqr(hop):
    """square flux"""
    hop[0].flux=hop[0].flux * hop[0].flux
    print('\nSpectrum has been squared\n')
    return hop

