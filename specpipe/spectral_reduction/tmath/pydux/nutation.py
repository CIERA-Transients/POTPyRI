def nutation(jd):
    import numpy as np
    from tmath.pydux.mod360 import mod360
    degrad = 180./np.pi
# returns dpsi and deps in arcsec
    t = (jd-2451545.0)/36525.0
#     Sun's mean anomaly
    m = 357.52772 + 35999.050340*t - 0.0001603*t*t - (t*t*t)/300000.0
    m = mod360(m)
    m = m/degrad
#     Moon's mean anomaly
    mp = 134.96298 + 477198.867398*t + 0.0086972*t*t + (t*t*t)/56250.0
    mp = mod360(mp)
    mp = mp/degrad
#     Longitude of moon's ascending node
    omega = 125.04452 - 1934.136261*t + 0.0020708*t*t + (t*t*t)/450000.0
    omega = mod360(omega)
    omega = omega/degrad
#     moon's argument of latitude
    f = 93.27191 + 483202.017538*t - 0.0036825*t*t + (t*t*t)/327270.0
    f = mod360(f)
    f = f/degrad
#     mean elongation of the moon from the sun
    d = 297.85036 + 445267.111480*t - 0.0019142*t*t + (t*t*t)/189474.0
    d = mod360(d)
    d = d/degrad
#     Delta(psi) is nutation in longitude
#     Delta(epsilon) is nutation in obliquity

#     In units of 0.0001 arcsec:
    dpsi = (-171996.0 - 174.2*t)*np.sin(omega) \
           +(-13187.0 - 1.6*t)*np.sin(-2.0*d + 2.0*f + 2.0*omega) \
           +(-2274.0 - 0.2*t)*np.sin(2.0*f + 2.0*omega) \
           +(2062.0 + 0.2*t)*np.sin(2.0*omega) \
           +(1426.0 - 3.4*t)*np.sin(m) \
           +(712.0 + 0.1*t)*np.sin(mp) \
           +(-517.0 + 1.2*t)*np.sin(-2.0*d + m + 2.0*f + 2.0*omega) \
           +(-386.0 - 0.4*t)*np.sin(2.0*f + omega) \
           -301.0*np.sin(mp + 2.0*f + 2.0*omega) \
           +(217.0 - 0.5*t)*np.sin(-2.0*d - m + 2.0*f + 2.0*omega) \
           -158.0*np.sin(-2.0*d + mp) \
           +(129.0 + 0.1*t)*np.sin(-2.0*d + 2.0*f + omega) \
           +123.0*np.sin(-mp + 2.0*f + 2.0*omega) \
           +63.0*np.sin(2.0*d) \
           +(63.0 + 0.1*t)*np.sin(mp + omega) \
           -59.0*np.sin(2.0*d - mp + 2.0*f + 2.0*omega) \
           +(-58.0 - 0.1*t)*np.sin(-mp + omega) \
           -51.0*np.sin(mp + 2.0*f + omega) \
           +48.0*np.sin(-2.0*d + 2.0*mp) \
           +46.0*np.sin(-2.0*mp + 2.0*f + omega) \
           -38.0*np.sin(2.0*d + 2.0*f + 2.0*omega) \
           -31.0*np.sin(2.0*mp + 2.0*f + 2.0*omega) \
           +29.0*np.sin(2.0*mp) \
           +29.0*np.sin(-2.0*d + mp + 2.0*f + 2.0*omega) \
           +26.0*np.sin(2.0*f) \
           -22.0*np.sin(-2.0*d + 2.0*f) \
           +21.0*np.sin(-mp + 2.0*f + omega) \
           +(17.0 - 0.1*t)*np.sin(2.0*m) \
           +16.0*np.sin(2.0*d - mp + omega) \
           +(-16.0 + 0.1*t)*np.sin(-2.0*d + 2.0*m + 2.0*f + 2.0*omega) \
           -15.0*np.sin(m + omega) \
           -13.0*np.sin(-2.0*d + mp + omega) \
           -12.0*np.sin(-m + omega) \
           +11.0*np.sin(2.0*mp - 2.0*f) \
           -10.0*np.sin(2.0*d - mp + 2.0*f + omega) \
           -8.0*np.sin(2.0*d + mp + 2.0*f + 2.0*omega) \
           +7.0*np.sin(m + 2.0*f + 2.0*omega) \
           -7.0*np.sin(-2.0*d + m + mp) \
           -7.0*np.sin(-m + 2.0*f + 2.0*omega) \
           -7.0*np.sin(2.0*d + 2.0*f + omega) \
           +6.0*np.sin(2.0*d + mp) \
           +6.0*np.sin(-2.0*d + 2.0*mp + 2.0*f + 2.0*omega) \
           +6.0*np.sin(-2.0*d + mp + 2.0*f + omega) \
           -6.0*np.sin(2.0*d - 2.0*mp + omega) \
           -6.0*np.sin(2.0*d + omega) \
           +5.0*np.sin(-m + mp) \
           -5.0*np.sin(-2.0*d - m + 2.0*f + omega) \
           -5.0*np.sin(-2.0*d + omega) \
           -5.0*np.sin(2.0*mp + 2.0*f + omega) \
           +4.0*np.sin(-2.0*d + 2.0*mp + omega) \
           +4.0*np.sin(-2.0*d + m + 2.0*f + omega) \
           +4.0*np.sin(mp - 2.0*f) \
           -4.0*np.sin(-d + mp) \
           -4.0*np.sin(-2.0*d + m) \
           -4.0*np.sin(d) \
           +3.0*np.sin(mp + 2.0*f) \
           -3.0*np.sin(-2.0*mp + 2.0*f + 2.0*omega) \
           -3.0*np.sin(-d - m + mp) \
           -3.0*np.sin(m + mp) \
           -3.0*np.sin(-m + mp + 2.0*f + 2.0*omega) \
           -3.0*np.sin(2.0*d - m -mp + 2.0*f + 2.0*omega) \
           -3.0*np.sin(3.0*mp + 2.0*f + 2.0*omega) \
           -3.0*np.sin(2.0*d - m + 2.0*f + 2.0*omega)
    dpsi = dpsi*0.0001

    deps = (92025.0 + 8.9*t)*np.cos(omega) \
           +(5736.0 - 3.1*t)*np.cos(-2.0*d + 2.0*f + 2.0*omega) \
           +(977.0 - 0.5*t)*np.cos(2.0*f + 2.0*omega) \
           +(-895.0 + 0.5*t)*np.cos(2.0*omega) \
           +(54.0 - 0.1*t)*np.cos(m) \
           -7.0*np.cos(mp) \
           +(224.0 - 0.6*t)*np.cos(-2.0*d + m + 2.0*f + 2.0*omega) \
           +200.0*np.cos(2.0*f + omega) \
           +(129.0 - 0.1*t)*np.cos(mp + 2.0*f + 2.0*omega) \
           +(-95.0 + 0.3*t)*np.cos(-2.0*d - m + 2.0*f + 2.0*omega) \
           -70.0*np.cos(-2.0*d + 2.0*f + omega) \
           -53.0*np.cos(-mp + 2.0*f + 2.0*omega) \
           -33.0*np.cos(mp + omega) \
           +26.0*np.cos(2.0*d - mp + 2.0*f + 2.0*omega) \
           +32.0*np.cos(-mp + omega) \
           +27.0*np.cos(mp + 2.0*f + omega) \
           -24.0*np.cos(-2.0*mp + 2.0*f + omega) \
           +16.0*np.cos(2.0*d + 2.0*f + 2.0*omega) \
           +13.0*np.cos(2.0*mp + 2.0*f + 2.0*omega) \
           -12.0*np.cos(-2.0*d + mp + 2.0*f + 2.0*omega) \
           -10.0*np.cos(-mp + 2.0*f + omega) \
           -8.0*np.cos(2.0*d - mp + omega) \
           +7.0*np.cos(-2.0*d + 2.0*m + 2.0*f + 2.0*omega) \
           +9.0*np.cos(m + omega) \
           +7.0*np.cos(-2.0*d + mp + omega) \
           +6.0*np.cos(-m + omega) \
           +5.0*np.cos(2.0*d - mp + 2.0*f + omega) \
           +3.0*np.cos(2.0*d + mp + 2.0*f + 2.0*omega) \
           -3.0*np.cos(m + 2.0*f + 2.0*omega) \
           +3.0*np.cos(-m + 2.0*f + 2.0*omega) \
           +3.0*np.cos(2.0*d + 2.0*f + omega) \
           -3.0*np.cos(-2.0*d + 2.0*mp + 2.0*f + 2.0*omega) \
           -3.0*np.cos(-2.0*d + mp + 2.0*f + omega) \
           +3.0*np.cos(2.0*d - 2.0*mp + omega) \
           +3.0*np.cos(2.0*d + omega) \
           +3.0*np.cos(-2.0*d - m + 2.0*f + omega) \
           +3.0*np.cos(-2.0*d + omega) \
           +3.0*np.cos(2.0*mp + 2.0*f + omega)
    deps = deps*0.0001

    return dpsi, deps
