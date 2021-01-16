def onkeypress(event):
#    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
#        event.button, event.x, event.y, event.xdata, event.ydata)
    import matplotlib.pyplot as plt
    import numpy as np
    import tmath.wombat.womconfig as womconfig
#    global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
    if (event.key == 'a'):
        plt.plot(event.xdata,event.ydata,'ro')
        womconfig.tmpsplptsx.append(event.xdata)
        womconfig.tmpsplptsy.append(event.ydata)
        womconfig.nsplinepoints=womconfig.nsplinepoints+1
    if (event.key == 's') and (womconfig.nsplinepoints != 0):
        delindex=np.sqrt((womconfig.tmpsplptsx - event.xdata)**2 + (womconfig.tmpsplptsy - event.ydata)**2).argmin()
        plt.plot(womconfig.tmpsplptsx[delindex],womconfig.tmpsplptsy[delindex],'gx',markersize=12)
        del womconfig.tmpsplptsx[delindex]
        del womconfig.tmpsplptsy[delindex]
        womconfig.nsplinepoints=womconfig.nsplinepoints-1
    if (event.key == 'd') and (womconfig.nsplinepoints <= 3):
        tpl=plt.title('Spline requires FOUR points!',size=24)
        plt.pause(1.0)
        tpl.set_visible(False)
    if (event.key == 'd') and (womconfig.nsplinepoints > 3):
        womconfig.pflag='done'
        return
