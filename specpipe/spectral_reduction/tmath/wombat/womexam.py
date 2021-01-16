def womexam(hop):
    """print spectrum bin values"""
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.yesno import yesno
    done = False
    while (not done):
        answer='y'
        wave,flux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
        indexblue=womget_element(hop[0].wave,wave[0])
        indexred=womget_element(hop[0].wave,wave[-1])
        var=hop[0].var[indexblue:indexred+1]
        if (len(wave) > 30):
            print('This will print out {} values'.format(len(wave)))
            print('Do you really want to do that?')
            answer=yesno('n')
        if (answer == 'y'):
            done=True
    print('\nWave    Flux')
    for i,_ in enumerate(wave):
        print('{} {} {}'.format(wave[i],flux[i], var[i]))
    return hop
            
