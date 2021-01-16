def womms(hop):
    """arithmetic on spectra"""
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat import HOPSIZE
    print('\nThis routine will perfrom arithmatic operations with two hoppers\n')
    print('Do you want to (a)dd, (s)ubtract, (m)ultiply, or (d)ivide spectra?\n')
    choice=inputter_single('(a/s/m/d) ','asmd')
    if (choice == 's'):
        print('Second spectrum will be subtracted from the first.')
    if (choice == 'd'):
        print('First spectrum will be divided by the first.')
    done = False
    while (not done):
        hopnum1=inputter('\nEnter the first hopper: ','int',False)
        print(' ')
        hopnum2=inputter('Enter the second hopper: ','int',False)
        if (hopnum1 < 1) or (hopnum1 > HOPSIZE) or (hopnum2 < 1) or \
           (hopnum2 > HOPSIZE):
            print('\nHopper must be in range 1-{}'.format(HOPSIZE))
        else:
            done = True
    if (hop[hopnum1].wave[0] != hop[hopnum2].wave[0])  \
       or (hop[hopnum1].wave[1] != hop[hopnum2].wave[1]) \
       or (hop[hopnum1].wave[-1] != hop[hopnum2].wave[-1]):
        print('Hoppers to not have the same wavelength scale!')
        return hop
    if (choice == 'a'):
        hop[0].flux=hop[hopnum1].flux + hop[hopnum2].flux
        hop[0].var=hop[hopnum1].var + hop[hopnum2].var
    elif (choice == 's'):
        hop[0].flux=hop[hopnum1].flux - hop[hopnum2].flux
        hop[0].var=hop[hopnum1].var + hop[hopnum2].var       
    elif (choice == 'm'):
        hop[0].flux=hop[hopnum1].flux * hop[hopnum2].flux
        hop[0].var=(hop[hopnum2].flux)**2*hop[hopnum1].var + \
                    (hop[hopnum1].flux)**2*hop[hopnum2].var
    elif (choice == 'd'):
        hop[0].flux=hop[hopnum1].flux / hop[hopnum2].flux
        hop[0].var=(hop[hopnum2].flux)**2*hop[hopnum1].var + \
                    (hop[hopnum1].flux)**2*hop[hopnum2].var
    hop[0].wave=hop[hopnum1].wave
    hop[0].obname=hop[hopnum1].obname
    hop[0].header=hop[hopnum1].header
    return hop

