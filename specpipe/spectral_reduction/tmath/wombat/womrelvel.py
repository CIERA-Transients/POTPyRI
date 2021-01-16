def womrelvel(hop):
    """relativistic velocity calculation"""
    from tmath.wombat.inputter import inputter
    light_speed=2.99792458e5
    print('Relativistic velocity calculation\n')
    
    lambda1=inputter('Observed wavelength ','float',False)
    print(' ')
    lambda0=inputter('Rest wavelength ','float',False)
    if (lambda0 <= 0):
        print('Invalid rest wavelength')
        return hop
    z=(lambda1-lambda0)/lambda0
    sq=(z+1.0)**2
    vel=z*light_speed
    relvel=((sq-1.0)/(sq+1.0))*light_speed

    print('\nz: {}'.format(z))
    print('Velocity: {}'.format(vel))
    print('Relativistic velocity {}'.format(relvel))
    
    return hop

