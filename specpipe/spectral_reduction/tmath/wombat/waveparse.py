def waveparse():
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
    return waveb,waver

