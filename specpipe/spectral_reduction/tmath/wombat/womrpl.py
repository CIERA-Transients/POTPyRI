def womrpl(hop):
    import numpy as np
    done=False
    while (not done):
        inputfile=input('Name of file to be read? ')
        inputfile=inputfile.strip()
        if (inputfile == ''):
            return hop
        try:
            data=np.loadtxt(inputfile)
        except OSError:
            print('File {} cannot be opened.'.format(inputfile))
        else:
            done=True
    wave=data[:,0]
    flux=data[:,1]
    if (data.shape[1] > 2):
        var=data[:,2]
    else:
        var=np.ones(wave.shape)

    var=var**2
    print('Found {} bins.'.format(len(wave)))

    hop[0].wave=wave.copy()
    hop[0].flux=flux.copy()
    hop[0].var=var.copy()
    hop[0].obname=inputfile
    hop[0].header=''
    return hop

