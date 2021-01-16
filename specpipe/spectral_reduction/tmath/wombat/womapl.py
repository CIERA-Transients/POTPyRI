def womapl(hop):
    import matplotlib.pyplot as plt
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.plot(hop[0].wave,hop[0].flux,drawstyle='steps-mid')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    return hop


