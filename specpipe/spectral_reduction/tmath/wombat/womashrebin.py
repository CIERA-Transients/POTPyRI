def womashrebin(wave,flux,nwave):
    """moo"""
    import numpy as np
    from tmath.wombat.womrebin_tform import womrebin_tform
    from tmath.wombat.womply import womply
    ACCBIN=0.0001
    arcone=np.zeros(2)
    arctwo=np.zeros(2)
    darcone=np.zeros(2)
    npix=len(wave)
    nrbin=len(nwave)
    rbin=np.zeros(nrbin)
    arc1=np.array([wave[0]-(wave[-1]-wave[0])/(npix-1), \
                   ((wave[-1]-wave[0])*npix)/(npix-1)])
    arc2=np.array([nwave[0]-(nwave[-1]-nwave[0])/(nrbin-1), \
                   ((nwave[-1]-nwave[0])*nrbin)/(nrbin-1)])
    inarc1=2
    inarc2=2
    in11=inarc1-1
    in12=inarc1-2
    in21=inarc2-1
    lstart=0
    acc=ACCBIN/npix
    rnpix=1./npix
    rnrbin=1./nrbin
    rx2=0.5*rnrbin
    for i in range(0,inarc1):
        arcone[inarc1-1-i]=arc1[i]
        darcone[inarc1-1-i]=i*arc1[i]
    for i in range(0, inarc2):
        arctwo[inarc2-1-i]=arc2[i]
    rslope=1./womply(darcone,in12,0.2)
    x1in=0.2
    x1=womrebin_tform(rx2,x1in,arctwo,darcone,in21,acc,rslope,in11,in12,arcone)
    x1=x1*npix
    dx=0
    nsgn=1
    if((womply(arctwo,in21,1)-arc2[1])*rslope < 0):
        nsgn=-1
    nstop=nrbin
    j1=np.round(x1)-1
    for k in range(0, nstop):
        rx2=rx2+rnrbin
        x2in=(x1+dx)*rnpix
        x2=womrebin_tform(rx2,x2in,arctwo,darcone,in21,acc,rslope,in11,in12,arcone)
        x2=x2*npix
        dx=x2-x1
        j2=np.round(x2)-1
        d=0
        if (lstart == 0):
            lstart=1
            m1=int(max([min([j1-1,npix-1]),1]))
            m2=int(max([min([j1,npix-1]),1]))
            m3=int(max([min([j1+1,npix-1]),1]))
            a=(flux[m1]+flux[m3])*0.5
            b=(a-flux[m1])*0.5
            c=(13./12.)*flux[m2]-a/12.0
            a=(a-flux[m2])/3.0
            y=x1-j1-1
            dd=nsgn*((((a*y)+b)*y+c)*y-b*0.25)+a*0.125+c*0.5

        m1=int(max([min([j2-1,npix-1]),1]))
        m2=int(max([min([j2,npix-1]),1]))
        m3=int(max([min([j2+1,npix-1]),1]))
        a=(flux[m1]+flux[m3])*0.5
        b=(a-flux[m1])*0.5
        c=1.0833333333333333*flux[m2]-a*0.0833333333333333
        a=(a-flux[m2])*0.333333333333333
        y=x2-j2-1
        d=d-dd
        dd=nsgn*((((a*y)+b)*y+c)*y-b*0.25)
        ddd=a*0.125+c*0.5
        d=d+dd-ddd
        dd=dd+ddd
      #  print(d,dd,ddd)
        for kk in range(int(j1), int(j2)+1*int(nsgn), int(nsgn)):
            d=d+flux[int(max([min([kk,npix-1]),1]))]
       # print(d)

        rbin[k]=d/np.abs(dx)
        x1=x2
        j1=j2
    nflux=rbin.copy()
    return nflux

