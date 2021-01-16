def womspl(hop,fig):
    """fit spline to spectrum"""
    import matplotlib.pyplot as plt
    import numpy as np
    import copy
    from tmath.wombat.womplot import womplot
    from tmath.wombat.onclick import onclick
    from scipy.interpolate import splrep,splev
    from tmath.wombat.inputter import inputter
    from tmath.wombat.yesno import yesno
    from tmath.wombat import HOPSIZE
    import tmath.wombat.womconfig as womconfig
#    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    print('\nObject is {}\n'.format(hop[0].obname))
    womplot(hop)
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    womconfig.nsplinepoints=0
    womconfig.tmpsplptsx=[]
    womconfig.tmpsplptsy=[]

    done=False
    while (not done):
        plt.cla()
        plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
        if (len(womconfig.tmpsplptsx) > 0):
            plt.plot(womconfig.tmpsplptsx,womconfig.tmpsplptsy,'ro')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title(hop[0].obname)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        print('\nClick on continuum points for spline fit.')
        print('Left button    = add point')
        print('Middle button  = delete point')
        print('Right button   = done\n')
        womconfig.pflag=''
        while (womconfig.pflag != 'done'):
            plt.pause(0.01)
        fig.canvas.mpl_disconnect(cid)

        splptsy=[z for _,z in sorted(zip(womconfig.tmpsplptsx,womconfig.tmpsplptsy))]
        splptsx=sorted(womconfig.tmpsplptsx)
        spline=splrep(splptsx,splptsy,k=3)
        splineresult=splev(hop[0].wave,spline)
        plt.plot(hop[0].wave,splineresult,drawstyle='steps-mid')
        plt.pause(0.01)
        print('Is this fit OK? ')
        answer=yesno('y')
        if (answer == 'y'):
            done=True
    print('\nSubtract spline fit from flux?\n')
    sub=yesno('n')
    if (sub == 'y'):
        hop[0].flux=hop[0].flux - splineresult
    print('\nStore spline in hopper?\n')
    store=yesno('y')
    if (store == 'y'):
        hopnum=0
        while (hopnum < 1) or (hopnum > HOPSIZE):
            hopnum=inputter('Store in which hopper: ','int',False)
        hop[hopnum]=copy.deepcopy(hop[0])
        hop[hopnum].flux=splineresult.copy()
        hop[hopnum].obname=hop[hopnum].obname+'spline'
        hop[hopnum].var=np.zeros(len(hop[0].wave))
    return hop

