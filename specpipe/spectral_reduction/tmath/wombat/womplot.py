def womplot(hop):
    import matplotlib.pyplot as plt
    from tmath.wombat.yesno import yesno
    from tmath.wombat.wshow import wshow
    plt.cla()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    if (hop[0].wave[0] < 0):
        plt.xlabel('Velocity')
    plt.ylabel('Flux')
    plt.title(hop[0].obname)
    plt.pause(0.01)
    wshow()
    print('Change scale?')
    answer=yesno('n')
    if (answer == 'y'):
        xmin,xmax=plt.xlim()
        ymin,ymax=plt.ylim()
        done=False
        while (not done):
            print('Click corners of box to change plot scale')
            newlims=plt.ginput(2, timeout=-1)
            xmin=newlims[0][0]
            ymin=newlims[0][1]
            xmax=newlims[1][0]
            ymax=newlims[1][1]
            plt.cla()
            plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
            plt.xlabel('Wavelength')
            if (hop[0].wave[0] < 0):
                plt.xlabel('Velocity')
            plt.ylabel('Flux')
            plt.title(hop[0].obname)
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            print('Is this OK?')
            loopanswer=yesno('y')
            if (loopanswer == 'y'):
                done=True
            
    return hop

