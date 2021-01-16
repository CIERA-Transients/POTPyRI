def finalscaler(flux):
    import numpy as np
    
#    sortflux=np.sort(flux)
#    der=np.gradient(sortflux)
#    locs=np.where(np.abs(der) < 1.*np.abs(np.mean(der)))
#    sortfluxclip=sortflux[locs]
#    ymin=min(sortfluxclip)
#    ymax=max(sortfluxclip)
#    abandoning gradient method for clip/percentile method

    half=5
    nsig=3.
    boxsize=11
    for i in range(half,len(flux)-half):
        sample=flux[i-half:i+half+1]
        medval=np.median(sample)
        varpix=(sample[half]-medval)**2
        sum=-1.0*varpix
        for k in range(0,boxsize):
            diff=sample[k]-medval
            sum=sum+diff*diff
        sum=sum/(boxsize-1)
        if (varpix > nsig*nsig*sum):
            flux[i]=medval

    ymin=np.percentile(flux,.1)
    ymax=np.percentile(flux,99.9)

    return ymin, ymax

