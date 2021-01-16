def womreddening(hop):
    """redden or deredden with various reddening laws"""
    import matplotlib.pyplot as plt
    import extinction
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.inputter import inputter
    from tmath.wombat.yesno import yesno
    r_v=3.1
    print('Redden or deredden a spectrum')
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid',color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    flux=hop[0].flux.copy()
    action=inputter_single('(r)edden or (d)eredden the spectrum? (r/d) ', 'rd')
    print(' ')
    type=inputter_single('Do you want to enter the (c)olor excess, or (v)isual extinction? ','cv')
    print(' ')
    if (type == 'v'):
        av=inputter('Enter A_V in magnitudes: ','float',False)
    else:
        ebv=inputter('Enter E(B-V) in magnitudes: ','float',False)
        av=r_v*ebv
    print(' ')
    print('Do you want to use: ')
    print('(c)ardelli, Clayton, Mathis 1989')
    print("(o)'donnell 1994")
    print('(f)itzpatrick 1999\n')
    
    method=inputter_single('(c/o/f) ','cof')
    if (action == 'r'):
        if (method == 'c'):
            newflux=extinction.apply(extinction.ccm89(hop[0].wave,av,r_v),flux)
        elif (method == 'o'):
            newflux=extinction.apply(extinction.odonnell94(hop[0].wave,av,r_v),flux)
        else:
            newflux=extinction.apply(extinction.fitzpatrick99(hop[0].wave,av,r_v),flux)
    else:
        if (method == 'c'):
            ext=extinction.ccm89(hop[0].wave,av,r_v)
        elif (method == 'o'):
            ext=extinction.odonnell94(hop[0].wave,av,r_v)
        else:
            ext=extinction.fitzpatrick99(hop[0].wave,av,r_v)
        newflux=flux*10**(0.4*ext)
    plt.plot(hop[0].wave,newflux,drawstyle='steps-mid',color='r')
    print('\nOriginal spectrum in black, red/dered in red\n')
    print('Is this OK?\n')
    answer=yesno('y')
    if (answer == 'y'):
        hop[0].flux=newflux.copy()
        print('\nActive spectrum now changed')
    else:
        print('\nSorry to disappoint you, active spectrum unchanged')
    return hop

