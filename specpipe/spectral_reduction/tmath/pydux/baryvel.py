def baryvel(dje):
#       Based on Fortran routine of P. STUMPFF (IBM-VERSION 1979):
#                 ASTRON. ASTROPHYS.  SUPPL. SER. 41, 1 (1980)
#                 M. H. SLOVAK (VAX 11/780 IMPLEMENTATION 1986)
#                 G. TORRES (1989)
    import numpy as np

    au=1.4959787e8
    dct0=2415020.0
    dcjul=36525.0
    dcbes=0.313
    dctrop=365.24219572
    dc1900=1900.0
    dc2pi=2.0*np.pi
    forbel=[None]*7
    sorbel=[None]*17
    sn=[None]*4
    dvelh=[None]*3
    dvelb=[None]*3

#  constants dcfel(i,k) of fast-changing elements

#                  i = 1           i = 2         i = 3
    dcfel=[[1.7400353e+00, 6.2833195099091e+02, 5.2796e-06], \
                 [6.2565836e+00, 6.2830194572674e+02,-2.6180e-06], \
                 [4.7199666e+00, 8.3997091449254e+03,-1.9780e-05], \
                 [1.9636505e-01, 8.4334662911720e+03,-5.6044e-05], \
                 [4.1547339e+00, 5.2993466764997e+01, 5.8845e-06], \
                 [4.6524223e+00, 2.1354275911213e+01, 5.6797e-06], \
                 [4.2620486e+00, 7.5025342197656e+00, 5.5317e-06], \
                 [1.4740694e+00, 3.8377331909193e+00, 5.6093e-06]]

#  constants dceps and ccsel(i,k) of slowly changing elements

#                  i = 1      i = 2       i = 3

    dceps=[4.093198e-01,-2.271110e-04,-2.860401e-08]

    ccsel=[[1.675104e-02,-4.179579e-05,-1.260516e-07], \
                 [2.220221e-01, 2.809917e-02, 1.852532e-05], \
                 [1.589963e+00, 3.418075e-02, 1.430200e-05], \
                 [2.994089e+00, 2.590824e-02, 4.155840e-06], \
                 [8.155457e-01, 2.486352e-02, 6.836840e-06], \
                 [1.735614e+00, 1.763719e-02, 6.370440e-06], \
                 [1.968564e+00, 1.524020e-02,-2.517152e-06], \
                 [1.282417e+00, 8.703393e-03, 2.289292e-05], \
                 [2.280820e+00, 1.918010e-02, 4.484520e-06], \
                 [4.833473e-02, 1.641773e-04,-4.654200e-07], \
                 [5.589232e-02,-3.455092e-04,-7.388560e-07], \
                 [4.634443e-02,-2.658234e-05, 7.757000e-08], \
                 [8.997041e-03, 6.329728e-06,-1.939256e-09], \
                 [2.284178e-02,-9.941590e-05, 6.787400e-08], \
                 [4.350267e-02,-6.839749e-05,-2.714956e-07], \
                 [1.348204e-02, 1.091504e-05, 6.903760e-07], \
                 [3.106570e-02,-1.665665e-04,-1.590188e-07]]

#  constants of the arguments of the short-period perturbations by
#  the planets:  dcargs(i,k)

#                 i = 1           i = 2
    
    dcargs=[[5.0974222e+00,-7.8604195454652e+02], \
                  [3.9584962e+00,-5.7533848094674e+02], \
                  [1.6338070e+00,-1.1506769618935e+03], \
                  [2.5487111e+00,-3.9302097727326e+02], \
                  [4.9255514e+00,-5.8849265665348e+02], \
                  [1.3363463e+00,-5.5076098609303e+02], \
                  [1.6072053e+00,-5.2237501616674e+02], \
                  [1.3629480e+00,-1.1790629318198e+03], \
                  [5.5657014e+00,-1.0977134971135e+03], \
                  [5.0708205e+00,-1.5774000881978e+02], \
                  [3.9318944e+00, 5.2963464780000e+01], \
                  [4.8989497e+00, 3.9809289073258e+01], \
                  [1.3097446e+00, 7.7540959633708e+01], \
                  [3.5147141e+00, 7.9618578146517e+01], \
                  [3.5413158e+00,-5.4868336758022e+02]]
    
#  amplitudes ccamps(n,k) of the short-period perturbations
    
#       n = 1      n = 2      n = 3      n = 4      n = 5
    
    ccamps= \
        [[-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5,-2.490817e-7], \
               [-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5,-1.823138e-7], \
               [ 6.593466e-7, 1.322572e-5, 9.258695e-6,-4.674248e-7,-3.646275e-7], \
               [ 1.140767e-5,-2.049792e-5,-4.747930e-6,-2.638763e-6,-1.245408e-7], \
               [ 9.516893e-6,-2.748894e-6,-1.319381e-6,-4.549908e-6,-1.864821e-7], \
               [ 7.310990e-6,-1.924710e-6,-8.772849e-7,-3.334143e-6,-1.745256e-7], \
               [-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6,-1.655307e-7], \
               [-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6,-3.736225e-7], \
               [ 3.442177e-7, 2.671323e-6, 1.832858e-6,-2.394688e-7,-3.478444e-7], \
               [ 8.702406e-6,-8.421214e-6,-1.372341e-6,-1.455234e-6,-4.998479e-8], \
               [-1.488378e-6,-1.251789e-5, 5.226868e-7,-2.049301e-7, 0.0], \
               [-8.043059e-6,-2.991300e-6, 1.473654e-7,-3.154542e-7, 0.0], \
               [ 3.699128e-6,-3.316126e-6, 2.901257e-7, 3.407826e-7, 0.0], \
               [ 2.550120e-6,-1.241123e-6, 9.901116e-8, 2.210482e-7, 0.0], \
               [-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.0]]

#  constants of the secular perturbations in longitude ccsec3 and
#  ccsec(n,k)
    
#                 n = 1      n = 2        n = 3

    ccsec3=-7.757020e-08

    ccsec=[[1.289600e-06, 5.550147e-01, 2.076942e+00], \
                 [3.102810e-05, 4.035027e+00, 3.525565e-01], \
                 [9.124190e-06, 9.990265e-01, 2.622706e+00], \
                 [9.793240e-07, 5.508259e+00, 1.559103e+01]]
    
#  sidereal rate dcsld in longitude, rate ccsgd in mean anomaly

    dcsld=1.990987e-07
    ccsgd=1.990969e-07

#  some constants used in the calculation of the lunar contribution

    cckm=3.122140e-05
    ccmld=2.661699e-06
    ccfdi=2.399485e-07

#  constants dcargm(i,k) of the arguments of the perturbations of the
#  motion of the moon

#                  i = 1           i = 2

    dcargm=[[5.1679830e+00, 8.3286911095275e+03], \
       [5.4913150e+00,-7.2140632838100e+03], \
       [5.9598530e+00, 1.5542754389685e+04]]
    
#  amplitudes ccampm(n,k) of the perturbations of the moon
    
#        n = 1       n = 2       n = 3       n = 4
    
    ccampm = \
        [[1.097594e-01, 2.896773e-07, 5.450474e-02, 1.438491e-07], \
               [-2.223581e-02, 5.083103e-08, 1.002548e-02,-2.291823e-08], \
               [1.148966e-02, 5.658888e-08, 8.249439e-03, 4.063015e-08]]

#  ccpamv = a*m*dl/dt (planets); dc1mme = 1 - mass(earth+moon)
    
    ccpamv=[8.326827e-11,1.843484e-11,1.988712e-12,1.881276e-12]
    dc1mme=0.99999696
    
#  program execution begins
    
#  time-arguments
    
    dt = (dje - dct0)/dcjul
    t = dt
    dtsq = dt * dt
    tsq = dtsq
#  values of all elements for the instant dje

    for k in range(0,8):
        dlocal=np.fmod(dcfel[k][0]+dt*dcfel[k][1]+dtsq*dcfel[k][2], dc2pi)
        if (k == 0):
            dml = dlocal
        if (k != 0):
            forbel[k-1]=dlocal

    g=forbel[0]
    deps=np.fmod(dceps[0]+dt*dceps[1]+dtsq*dceps[2], dc2pi)

    for k in range(0,17):
        sorbel[k] = np.fmod(ccsel[k][0]+t*ccsel[k][1]+tsq*ccsel[k][2], dc2pi)
    e=sorbel[0]
#  secular perturbations in longitude
    for k in range(0,4):
        a =np.fmod(ccsec[k][1]+t*ccsec[k][2], dc2pi)
        sn[k] = np.sin(a)

#  periodic perturbations of the emb (earth-moon barycenter)

    pertl = ccsec[0][0]*sn[0] +ccsec[1][0]*sn[1] \
        +(ccsec[2][0]+t*ccsec3)*sn[2] +ccsec[3][0]*sn[3]

    pertld = 0.0
    pertr =  0.0
    pertrd = 0.0
    
    for k in range(0,15):
        a = np.fmod(dcargs[k][0]+dt*dcargs[k][1], dc2pi)
        cosa = np.cos(a)
        sina = np.sin(a)
        pertl = pertl + ccamps[k][0]*cosa + ccamps[k][1]*sina
        pertr = pertr + ccamps[k][2]*cosa + ccamps[k][3]*sina
        if (k < 11):
            pertld=pertld+(ccamps[k][1]*cosa-ccamps[k][0]*sina)*ccamps[k][4]
            pertrd=pertrd+(ccamps[k][3]*cosa-ccamps[k][2]*sina)*ccamps[k][4]


#  elliptic part of the motion of the emb
    esq = e * e
    dparam = 1.0 - esq
    param = dparam
    twoe = e + e
    twog = g + g
    phi = twoe*((1.0 - esq*(1.0/8.0))*np.sin(g) + e*(5.0/8.0)*np.sin(twog) \
                      + esq*0.5416667*np.sin(g+twog) )
    f = g + phi
    sinf = np.sin(f)
    cosf = np.cos(f)
    dpsi = dparam/(1.0 + e*cosf)
    phid = twoe*ccsgd*((1.0+esq*1.5)*cosf+e*(1.25-sinf*sinf*0.5))
    psid = ccsgd*e*sinf/np.sqrt(param)

#  perturbed heliocentric motion of the emb

    d1pdro = (1.0 + pertr)
    drd = d1pdro*(psid+dpsi*pertrd)
    drld = d1pdro*dpsi*(dcsld+phid+pertld)
    dtl = np.fmod(dml+phi+pertl, dc2pi)
    dsinls = np.sin(dtl)
    dcosls = np.cos(dtl)
    dxhd = drd*dcosls - drld*dsinls
    dyhd = drd*dsinls + drld*dcosls

#  influence of eccentricity, evection and variation of the geocentric
#  motion of the moon

    pertl =  0.0
    pertld = 0.0
    pertp =  0.0
    pertpd = 0.0
    
    for k in range(0,3):
        a = np.fmod(dcargm[k][0] + dt*dcargm[k][1], dc2pi)
        sina = np.sin(a)
        cosa = np.cos(a)
        pertl   = pertl  + ccampm[k][0]*sina
        pertld  = pertld + ccampm[k][1]*cosa
        pertp   = pertp  + ccampm[k][2]*cosa
        pertpd  = pertpd - ccampm[k][3]*sina


#  heliocentric motion of the earth

    tl =  forbel[1] + pertl
    sinlm = np.sin(tl)
    coslm = np.cos(tl)
    sigma = cckm/(1.0 + pertp)
    a = sigma*(ccmld+pertld)
    b = sigma*pertpd
    dxhd = dxhd + a*sinlm + b*coslm
    dyhd = dyhd - a*coslm + b*sinlm
    dzhd =      - sigma*ccfdi*np.cos(forbel[2])

#  barycentric motion of the earth

    dxbd = dxhd*dc1mme
    dybd = dyhd*dc1mme
    dzbd = dzhd*dc1mme

    for k in range(0,4):
        plon = forbel[k+3]
        pomg = sorbel[k+1]
        pecc = sorbel[k+9]
        tl = np.fmod(plon+2.0*pecc*np.sin(plon-pomg), dc2pi)
#        sinlp[k] = sin(tl)
#        coslp[k] = np.cos(tl)
        dxbd = dxbd + ccpamv[k]*(np.sin(tl) + pecc*np.sin(pomg))
        dybd = dybd - ccpamv[k]*(np.cos(tl) + pecc*np.cos(pomg))
        dzbd = dzbd - ccpamv[k]*sorbel[k+13]*np.cos(plon-sorbel[k+5])

#  transition to mean equator of date

    dcosep = np.cos(deps)
    dsinep = np.sin(deps)
    dyahd = dcosep*dyhd - dsinep*dzhd
    dzahd = dsinep*dyhd + dcosep*dzhd
    dyabd = dcosep*dybd - dsinep*dzbd
    dzabd = dsinep*dybd + dcosep*dzbd
#
    dvelh[0] = dxhd*au
    dvelh[1] = dyahd*au
    dvelh[2] = dzahd*au
    dvelb[0] = dxbd*au
    dvelb[1] = dyabd*au
    dvelb[2] = dzabd*au

    return dvelh,dvelb

