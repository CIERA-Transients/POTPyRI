def wombluen(hop):
    """bluen spectrum with scattering law"""
    from tmath.wombat.inputter import inputter
    print('\nThis routine will bluen a spectrum with a power law')
    print('of lambda^{-a}.  You will supply the value of a.\n')

    factor=inputter('Enter exponential factor: ','float',False)

    wavefac=hop[0].wave**(-1.0*factor)
    wavefac=wavefac/wavefac[-1]
    hop[0].flux=hop[0].flux*wavefac

    return hop

