def waveparse(wave,oldwaveb,oldwaver):
    import re
    done=False
    while (not done):
        waveparse=input('Enter wavelength range: ')
        if (waveparse != ''):
            waveparse=waveparse.strip()
            wavesplit=re.split('\W+',waveparse)
            if (len(wavesplit) > 1):
                try:
                    waveb=float(wavesplit[0])
                    waver=float(wavesplit[1])
                except ValueError:
                    print('Please enter numbers, #### ####,\n')
                    print('####,####,or ####-####\n')
                    print('You entered {}\n'.format(waveparse))
                else:
                    done=True
            else:
                print('Enter more than one number\n')
            if (waveb <= 0) or (waver <= 0)  or \
               (waveb < wave[0]) or (waver < wave[0]) or \
               (waveb > wave[-1]) or (waver > wave[-1]) or \
               (waveb > waver):
                print('Wavelengths out of bounds.  Try again')
                done = False
        else:
            waveb=oldwaveb
            waver=oldwaver
            done = True
    return waveb,waver

