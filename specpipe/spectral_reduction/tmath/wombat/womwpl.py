def womwpl(hop):
    import numpy as np
    print('\nObject is {}\n'.format(hop[0].obname))
    spectxt=input('Enter the name for the output file: ')
    if (spectxt == ''):
        return hop
    spectxt=spectxt.strip()
    np.savetxt(spectxt,np.transpose([hop[0].wave,hop[0].flux,np.sqrt(hop[0].var)]))
    return hop
