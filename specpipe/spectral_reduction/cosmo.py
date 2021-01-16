from __future__ import print_function
#!/usr/bin/env python
# 2009-07-23 S.Rodney
# adapted from James Schombert's python version 
# of Ned Wright's cosmology calculator 
#  (www.astro.ucla.edu/~wright/CosmoCalc.html)
#
#  WORK IN PROGRESS : 
#    STILL NEED TO EXTRACT SOME OTHER FUNCTIONS FROM CALCULATE

helpstring = '''
Cosmology calculator ala Ned Wright (www.astro.ucla.edu/~wright)
input values = redshift, Ho, Omega_m, Omega_vac
ouput values = age at z, distance in Mpc, kpc/arcsec, apparent to abs mag conversion

Options:   -h for this message 
           -v for verbose response 
'''
from math import *

c = 299792.458 # velocity of light in km/sec
Tyr = 977.8    # coefficent for converting 1/H into Gyr

# Omega(radiation)
# includes 3 massless neutrino species, T0 = 2.72528
Or = lambda h : 4.165E-5/(h*h)
#Or = lambda h : 0     

def main(): 
    import sys, getopt
    
    # Default Parameter Values
    H0 = 70        # Hubble constant
    Om = 0.3       # Omega(matter)
    Ode = None      # Omega(dark energy) [sometimes Omega(lambda)
    w = -1         # DE equation of state
    n  = 1000      # number of points in integrals
    verbose = 0

    try:  opt,arg = getopt.getopt(
        sys.argv[1:],'hv',
        longopts=["help", "verbose","H0=","Om=","Ode=","z=","cz=","nsteps="] )
    except getopt.GetoptError:
        # print help info and exit:
        print("%s  \nfor help use --help"%getopt.GetoptError)
        return(-1)

    for o, a in opt:
        if o in ["-h", "--help"]:
            print(helpstring)
            sys.exit()
        elif o in ["-v", "--verbose"]:
            verbose=1
        elif o == "--H0" :
            H0 = float(a)
        elif o == "--Ode" :
            WL = float(a)
        elif o == "--Om" :
            Om = float(a)
        elif o == "--w" :
            w = float(a)
        elif o == "--z" :
            z = float(a)
        elif o == "--cz" :
            z = float(a) / c
        elif o == "--nsteps" :
            n = int(a)

    # make the universe flat unless Ode is provided
    if not Ode : Ode = 1.0 - Om

    # RUN THE CALCULATIONS HERE 
    age_Gyr, zage_Gyr, DTT_Gyr, DCMR_Mpc, V_Gpc, DA_Mpc, kpc_DA, DL_Mpc, mu = calculate( z, H0=70, Om=0.3, Ode=Ode, n=1000 )

    if verbose == 1:
        print('For H_o = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % Om + ', Omega_vac = ',)
        print('%1.2f' % Ode + ', z = ' + '%1.3f' % z)
        print('It is now ' + '%1.1f' % age_Gyr + ' Gyr since the Big Bang.')
        print('The age at redshift z was ' + '%1.1f' % zage_Gyr + ' Gyr.')
        print('The light travel time was ' + '%1.1f' % DTT_Gyr + ' Gyr.')
        print('The comoving radial distance, which goes into Hubbles law, is',)
        print('%1.1f' % DCMR_Mpc + ' Mpc or ' + '%1.1f' % (Tyr/H0)*DCMR + ' Gly.')
        print('The comoving volume within redshift z is ' + '%1.1f' % V_Gpc + ' Gpc^3.')
        print('The angular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc or',)
        print('%1.1f' % DA_Gyr + ' Gly.')
        print('This gives a scale of ' + '%.2f' % kpc_DA + ' kpc/".')
        print('The luminosity distance D_L is ' + '%1.1f' % DL_Mpc + ' Mpc or ' + '%1.1f' % DL_Gyr + ' Gly.')
        print('The distance modulus, m-M, is '+'%1.2f' % (5*log10(DL_Mpc*1e6)-5))
    else:
        print('age = %1.2f Gyr' % zage_Gyr)
        print('DC = %1.2f Mpc' % DCMR_Mpc)
        print('DL = %1.2f Mpc' % DL_Mpc)
        print('DA = %1.2f kpc' % kpc_DA)
        print('mu = %1.2f' % (5*log10(DL_Mpc*1e6)-5))


def calculate( z, H0=70, Om=0.3, Ode=0.7, n=1000 ):
    """ all the cosmo calculations together """
    # initialize constants
    Ok = 0.        # Omega curvaturve = 1-Omega(total)
    DTT = 0.5      # time from z to now in units of 1/H0
    DTT_Gyr = 0.0  # value of DTT in Gyr
    age = 0.5      # age of Universe in units of 1/H0
    age_Gyr = 0.0  # value of age in Gyr
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0 
    DCMR_Gyr = 0.0
    DA = 0.0       # angular size distance
    DA_Mpc = 0.0
    DA_Gyr = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    DL_Gyr = 0.0   # DL in units of billions of light years
    V_Gpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))

    h = H0/100.
    Ok = 1-Om-Ode-Or(h)
    az = 1.0/(1+1.0*z)
    age = 0.
    for i in range(n):
        a = az*(i+0.5)/n
        adot = sqrt(Ok+(Om/a)+(Or(h)/(a*a))+(Ode*a*a))
        age = age + 1./adot

    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(Ok+(Om/a)+(Or(h)/(a*a))+(Ode*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR

    # tangential comoving distance
    ratio = 1.00
    x = sqrt(abs(Ok))*DCMR
    if x > 0.1:
        if Ok > 0:
            ratio =  0.5*(exp(x)-exp(-x))/x 
        else:
            ratio = sin(x)/x
    else:
        y = x*x
        if Ok < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DA_Gyr = (Tyr/H0)*DA
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL

    mu = (5*log10(DL_Mpc*1e6)-5)  # Distance modulus
    
    return( age_Gyr, zage_Gyr, DTT_Gyr, DCMR_Mpc, V_Gpc, DA_Mpc, kpc_DA, DL_Mpc, mu )
    
def zfromt( t, H0=70, Om=0.3, unit='Gyr', debug=False ):
    """ The redshift z when the universe has an age t.
    set unit=None for ages in units of 1/H0
    From Peebles:1993 p.315  
    This is crude, and only applicable for a flat 
    lambdaCDM universe.  Errors approach 1% at z=25.
    """
    from numpy import sqrt, sinh
    if unit=='Gyr' :  tx = t * (H0 / Tyr)
    if unit=='Myr' :  tx = t * (H0 / Tyr) * 1e3
    elif unit=='yr' : tx = t * (H0 / Tyr) * 1e9
    A = sqrt( Om/(1-Om)) 
    B = (3*tx/2) * sqrt(1-Om)
    z = ( A * sinh( B ) )**(-2/3.) - 1
    return( z ) 

def agez( z, H0=70, Om=0.3, Ode=0.7, n=1000, unit='Gyr', debug=False ):
    """ The age of the universe at redshift z
    set unit=None to get the age in units of 1/H0"""
    # initialize constants
    from scipy import integrate as scint
    from numpy import iterable, array, sqrt, append, abs, sin, sinh, ndarray
    if debug : import pdb; pdb.set_trace()

    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    h = H0/100.
    Ok = 1-Om-Ode-Or(h)

    if not iterable(z) : z = array([ z ])
    elif not isinstance( z, ndarray) :  z = array( z )

    # do integral over a=1/(1+z) from az to 1, Quadrature integration
    integrand = lambda a : 1./sqrt(Ok+(Om/a)+(Or(h)/(a*a))+(Ode*a*a))

    zage = []
    for zz in z : 
        az = 1.0/(1+1.0*zz)
        zage.append( scint.quad( integrand, 1e-6, az )[0] )
    zage = array(zage)
    if len(zage)==1: zage = zage[0]
    zage_Gyr = (Tyr/H0)*zage

    if unit=='Gyr' : return( zage_Gyr )
    elif unit=='yr' : return( zage_Gyr*1e9 )
    elif unit==None : return( zage )


    

def DC( z, H0=70, Om=0.3, Ode=0.7, n=1000, unit=None ):
    """ Comoving radial distance out to redshift z,
    which goes into the Hubble law. 
    Use unit=None for distance in units of c/H0"""
    from numpy import sqrt, iterable, array

    h = H0/100.
    Ok = 1-Om-Ode-Or(h)
    az = 1.0/(1+1.0*z)
    if not iterable( az ) : az = array( [az] ) 

    try : 
        # Faster with SciPy
        # do integral over a=1/(1+z) from az to 1, 
        # using quadrature integration
        from scipy import integrate as scint
        integrand = lambda a : 1 / ( a * sqrt( 
                Ok + (Om/a) + (Or(h)/(a*a)) + (Ode*a*a) ) )
        DCMR = []
        for azi in az : 
            DCMR.append( scint.quad( integrand, azi, 1)[0] )
    except  : 
        # If no SciPy, do the slow way
        # do integral over a=1/(1+z) from az to 1 in n steps, 
        # using the midpoint rule
        from numpy import sqrt
        DCMR = 0.0
        for i in range(n):
            a = az+(1-az)*(i+0.5)/n
            adot = sqrt(Ok+(Om/a)+(Or(h)/(a*a))+(Ode*a*a))
            DCMR = DCMR + 1./(a*adot)
        DCMR = (1.-az)*DCMR/n
    if len(DCMR) == 1 : DCMR = DCMR[0]
    else : DCMR = array( DCMR )
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR

    if unit=="Mpc" : return( DCMR_Mpc )
    elif unit=="Gpc" : return( DCMR_Gpc )
    elif unit==None :  return( DCMR )

def DL_Wright( z, H0=70, Om=0.3, Ode=0.7, n=1000, unit=None ):
    """ 
    Ned Wright's version. No w.
    Luminosity distance out to redshift z.
    Units may be Mpc or Gpc, or use unit=None 
    for distance in units of c/H0"""
    from numpy import iterable, array
    h = H0/100.
    Ok = 1-Om-Ode-Or(h)
    az = 1.0/(1+1.0*z)

    # radial comoving distance
    DCMR = DC( z, H0=H0, Om=Om, Ode=Ode, n=n, unit=None )

    # tangential comoving distance
    ratio = 1.00
    x = sqrt(abs(Ok))*DCMR

    if not iterable(x) : x = array([x])
    ratio = []
    for xi in x :
        if xi > 0.1:
            if Ok > 0:
                ratio.append( 0.5*(exp(xi)-exp(-xi))/xi )
            else:
                ratio.append( sin(xi)/xi )
        else:
            y = xi*xi
            if Ok < 0: y = -y
            ratio.append(  1. + y/6. + y*y/120. )
    if len(ratio) == 1 : ratio = ratio[0]
    else : ratio = array(ratio)
    DCMT = ratio*DCMR
    DA = az*DCMT
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL

    if unit=="Mpc" : return( DL_Mpc )
    elif unit=="Gpc" : return( DL_Gpc )
    elif unit==None :  return( DL )


def DLFw( z, w=-1, Om=0.3, H0=70, unit=None, debug=False ):
    """ luminosity distance for constant w in a flat universe """
    return( DL( z, w0=w, wa=0, Om=Om, H0=H0, unit=unit, Flat=True ) )


def DL( z, Om=0.3, Ode=0.7, w0=-1, wa=0,
        H0=70, unit=None, Flat=False, Lambda=False,
        debug=False ):
    """ luminosity distance calculation allowing w!=-1 
    and allowing w to vary with time using a linear parameterization:
       w(a) = w0 + wa(1-a)
    To get a constant w, just set wa=0 and use w0 as w. 
    Set Flat==True to enforce a flat universe (Ode=1-Om)
    Lambda=True to force a constant w=-1
    """
    from scipy import integrate as scint
    from numpy import iterable, array, sqrt, append, \
        abs, sin, sinh, ndarray, exp, zeros, isfinite
    if debug : import pdb; pdb.set_trace()
    
    if Flat : Ode = 1-Om  # Omega_{DarkEnergy} for flat universe
    Ok = 1-Om-Ode  # Omega_curvature
    if Lambda : w0,wa = -1, 0  # DE is Cosmological constant

    if not iterable(z) : z = array([ z ])
    elif not isinstance( z, ndarray) :  z = array( z )
    
    # the integrand is 1/E(z)  where E(z) = H(z)/H0 is 
    # the dimensionless expansion rate
    invE = lambda zz : 1 / sqrt( Om*(1+zz)**3  + Ok*(1+zz)**2 + Ode*(1+zz)**(3*(1+w0+wa))/exp(3*wa*zz/(1+zz)) )
    isorted = z.argsort()
    integrated = zeros( len(z) )
    lastz, lastint = 0, 0 
    # slightly quicker to handle the redshifts in order 
    for i in isorted :
        zz = z[i] 
        #if isfinite( invE( zz ) ) : 
        #    integresult=0
        #    continue
        integresult,err = scint.quad( invE, lastz, zz )
        integrated[ i ] = integresult + lastint
        lastz = zz
        lastint = integrated[i]
    
    if Ok<0 : 
        DL = ( (1+z) / sqrt(abs(Ok)) ) * sin( sqrt(abs(Ok)) * integrated  ) 
    if Ok>0 : 
        DL = ( (1+z) / sqrt(abs(Ok)) ) * sinh( sqrt(abs(Ok)) * integrated  ) 
    elif Ok==0: 
        DL = (1+z) * integrated 
    if len(DL)==1 : DL = DL[0]

    DL_Mpc = (c/H0) * DL
    DL_Gyr = (Tyr/H0) * DL
    if unit=="Mpc" : return( DL_Mpc )
    elif unit=="Gpc" : return( DL_Gpc )
    elif unit==None :  return( DL )


def mu( z, H0=70, Om=0.3, Ode=0.7, w0=-1, wa=0):
    """ Distance modulus to redshift z, 
    for the given cosmology """
    from numpy import log10
    DL_Mpc = DL(  z, H0=H0, Om=Om, Ode=Ode, 
                  w0=w0, wa=wa, unit='Mpc' )
    return( 5*log10( DL_Mpc) + 25 )

def mue( z, H0=70, Om=0.3, Ode=0.7, w0=-1, wa=0, 
         Flat=False, Lambda=False):
    """ Distance modulus to redshift z, 
    relative to an empty universe.
    Set Flat=True to enforce a flat universe (Ok=0)
    by overriding the Ode parameter. """
    from numpy import log10
    if Flat : Ode = 1 - Om
    DL_Mpc = DL(  z, H0=H0, Om=Om, Ode=Ode, w0=w0, wa=wa, unit='Mpc' )
    mu_empty = 5 * (log10( DL_Mpc ) - log10( c*z/H0*(1+z/2)) )

    return( mu_empty )


def zfromd(d, H0=70, Om=0.3, Ode=0.7, w0=-1, wa=0 ):
    """ 
    given a distance in Mpc, compute the redshift 
    for the assumed cosmology
    """
    from scipy.optimize import fmin

    zstart = d*H0/c
    mini = lambda z : (d -  DL(z,H0=H0,Om=Om,Ode=Ode,
                               w0=w0, wa=wa, unit='Mpc'))**2

    zmin = fmin(mini,zstart,disp=0)
    return(zmin[0])

    
def volume( z, H0=70, Om=0.3, Ode=0.7, n=1000, unit="Mpc"):
    """ comoving volume out to redshift z.
    unit may by 'Mpc' or 'Gpc'"""
    from numpy import iterable,where,zeros,exp,sqrt,sin,array
    if unit not in ["Mpc","Gpc"]:
        raise exceptions.RuntimeError(
            "Unrecognized unit %s.   Use 'Mpc' or 'Gpc'"%unit)

    h = H0/100.
    Ok = 1-Om-Ode-Or(h)

    # compute the comoving radial distance
    DCMR = DC( z, H0=H0, Om=Om, Ode=Ode, n=n, unit=None )
    ratio = 1.00
    x = sqrt(abs(Ok))*DCMR

    # divide the distance vector at a comoving distance of 0.1
    if not iterable(x) : x = array([x])
    i1 = where(x>0.1)
    i2 = where(x<=0.1)[0]
    x1 = x[ i1 ]
    x2 = x[ i2 ]
    ratio = zeros( len(x) )

    # fill the ratio array for x>0.1
    if len(i1)>0: 
        if Ok > 0:
            ratio[i1] = (0.125*(exp(2.*x1)-exp(-2.*x1))-x1/2.)/(x1*x1*x1/3.)
        else:
            ratio[i1] = (x1/2. - sin(2.*x1)/4.)/(x1*x1*x1/3.)

    # fill the ratio array for x<=0.1
    if len(i2)>0:
        y = x2*x2
        if Ok < 0: y = -y
        ratio[i2] = 1. + y/5. + (2./105.)*y*y

    # combine ratio and distance to get volume
    VCM = ratio*DCMR*DCMR*DCMR/3.
    V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM

    if iterable( V_Gpc ) :
        if len(V_Gpc)==1: 
            V_Gpc = V_Gpc[0]

    if unit=="Mpc" : return( V_Gpc * 1e9 )
    elif unit=="Gpc" : return( V_Gpc )

def E( z, Om=0.3, Ode=0.7, w0=-1, wa=0, H0=70, 
       debug=False ):
    """ The dimensionless expansion rate :
    E(z) = H(z) / H0 
    (squared, this is he LHS of the Friedmann eqn)
    """
    from numpy import exp, sqrt 
    Ok = 1-Om-Ode  # Omega_curvature

    # the integrand is 1/E(z)  where E(z) = H(z)/H0 is 
    # the dimensionless expansion rate
    M = lambda zz : Om*(1+zz)**3   # Matter 
    K = lambda zz : Ok*(1+zz)**2   # Kurvature
    DE = lambda zz : Ode*(1+zz)**(3*(1+w0+wa))/exp(3*wa*zz/(1+zz)) # Dark Energy
    return(  sqrt( M(z) + K(z) + DE(z) )  )

#  THIS IS MUCH FASTER TO INTEGRATE AS A LAMBDA FUNCTION
#def invE( z, Om=0.3, Ode=0.7, w0=-1, wa=0, H0=70, 
#       debug=False ):
#    """ The inverse of the dimensionless expansion rate :
#     1/E(z) =  H0 / H(z)
#    (this is the integrand for computing A, R, Dl, etc)
#    """
#    return(  1 / E(z, Om=Om, Ode=Ode, w0=w0, wa=wa, H0=H0) )

def DA( z, Om=0.3, Ode=0.7, w0=-1, wa=0, H0=70, 
        unit=None, debug=False ):
    """ angular diameter distance """
    az = 1/(1+1.0*z) # scale factor of the universe
    DA = az*az * DL(z,Om=Om,Ode=Ode,w0=w0,wa=wa,H0=H0,unit=None)
    DA_Mpc = (c/H0) * DA
    DA_Gyr = (Tyr/H0) * DA
    if unit=="Mpc" : return( DA_Mpc )
    elif unit=="Gpc" : return( DA_Gpc )
    elif unit==None :  return( DA )
    
       
def rsDv( z, Om=0.3, Ode=0.7, w0=-1, wa=0, H0=70, 
       debug=False ):
    """ rs / Dv  : BAO constraint from Percival:2010 and Komatsu:2010 
    rs = sound horizon distance
    Dv = baryon drag epoch
    """
   
    h = (H0/100.0)
    Ob = 0.0445
    H = H0* E(z, Om, Ode, w0, wa, H0 ) 
    DAng = DA(z,Om=Om,Ode=Ode,w0=w0,wa=wa,H0=H0,unit='Mpc')
    Dv = (( 1 + z)**2 * DAng**2 * c * z / H )**(1/3.)
    rs = 153.5 * (Ob * h**2 /0.02273 )**(-0.134) * (Om*h**2/0.1326)**(-0.255)
    return( rs / Dv)

def A( z, Om=0.3, Ode=0.7, w0=-1, wa=0, H0=70, 
       debug=False ):
    """ angular diameter distance with Omega_matter,
    from Eisenstein:2005 (see also Kessler:2009) 
    """
    from scipy import integrate as scint
    from numpy import iterable, array, sqrt, append, abs, sin, sinh, ndarray, isfinite
    if debug : import pdb; pdb.set_trace()
    
    Ok = 1-Om-Ode  # Omega_curvature
    if not iterable(z) : z = array([ z ])
    elif not isinstance( z, ndarray) :  z = array( z )
    # the integrand is 1/E(z)  where E(z) = H(z)/H0 is 
    # the dimensionless expansion rate
    invE = lambda zz : 1 / sqrt( Om*(1+zz)**3  + Ok*(1+zz)**2 + Ode*(1+zz)**(3*(1+w0+wa))/exp(3*wa*zz/(1+zz)) )
    integrated = array([])
    for zz in z : 
        #if not isfinite( invE( zz ) ) : 
        #    integresult=0
        #    continue
        integresult,err = scint.quad( invE, 0, zz )
        integrated = append( integrated, integresult )
    if Ok<0 : 
        Sk = 1/sqrt(abs(Ok)) *  sin( sqrt(abs(Ok)) * integrated  ) 
    if Ok>0 : 
        Sk = 1/sqrt(Ok) * sinh( sqrt(abs(Ok)) * integrated  ) 
    elif Ok==0: 
        Sk = integrated 
    A = ( sqrt(Om) /  E(z,Om=Om,Ode=Ode,w0=w0,wa=wa,H0=H0)**(1/3.) ) * ( Sk / z )**(2/3.)
    if len(A)==1 : A = A[0]
    return( A )


def R( z, Om=0.3, Ode=0.7, w0=-1, wa=0, H0=70, 
       debug=False ):
    """ The shift parameter, constrained by CMB
    measurments in Komatsu:2008 (see also Kessler:2009) 
    """
    from scipy import integrate as scint
    from numpy import iterable, array, sqrt, append, abs, \
        sin, sinh, ndarray, zeros, isfinite, exp
    if debug : import pdb; pdb.set_trace()
    
    Ok = 1-Om-Ode  # Omega_curvature
    if not iterable(z) : z = array([ z ])
    elif not isinstance( z, ndarray) :  z = array( z )

    # the integrand is 1/E(z)  where E(z) = H(z)/H0 is 
    # the dimensionless expansion rate
    invE = lambda zz : 1 / sqrt( Om*(1+zz)**3  + Ok*(1+zz)**2 + Ode*(1+zz)**(3*(1+w0+wa))/exp(3*wa*zz/(1+zz)) )
    isorted = z.argsort()
    integrated = zeros( len(z) )
    lastz, lastint = 0, 0 
    # slightly quicker to handle the redshifts in order 
    for i in isorted :
        zz = z[i] 
        #if not isfinite( invE( zz ) ) : 
        #    integresult=0
        #    continue
        # WARNING : hiding error messages b/c R is undefined for some of parameter space
        integresult,err = scint.quad( invE, lastz, zz, full_output=1 )[:2]
        integrated[ i ] = integresult + lastint
        lastz = zz
        lastint = integrated[i]
    if Ok<0 : 
        Sk = 1/sqrt(abs(Ok)) *  sin( sqrt(abs(Ok)) * integrated  ) 
    if Ok>0 : 
        Sk = 1/sqrt(Ok) * sinh( sqrt(abs(Ok)) * integrated  ) 
    elif Ok==0: 
        Sk = integrated 
    R = sqrt(Om) * Sk 
    if len(R)==1 : R = R[0]
    return( R )



if __name__ == "__main__":
    main()
