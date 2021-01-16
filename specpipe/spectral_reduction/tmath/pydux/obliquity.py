def obliquity(julian_day):
    """
    Fit of Laskar, 1986, A&A, 157, 68
    this is *mean* obliquity, need nutation correction for *true* obliquity
    calculation in arcsec
    Note that this formula is only valid for values of ten_k_julian_years in
    the range +/-1 (i.e., -8000BC to 12000AD).
    """
    import numpy as np

    julian_century = (julian_day-2451545.0)/36525.0
    ten_k_julian_years = julian_century/100.0
    laskars_polynomial = [2.45000e+00, 5.79000e+00, 2.78700e+01, 7.12000e+00,
                          -3.90500e+01, -2.49670e+02, -5.13800e+01, 1.99925e+03,
                          -1.55000e+00, -4.68093e+03, 0.00000e+00]
#    a=[0,-4680.93,-1.55, 1999.25,-51.38,-249.67, -39.05,7.12,27.87,5.79,2.45]
    corr = np.polyval(laskars_polynomial, ten_k_julian_years)
    epsilon = 23.4392911111+corr/60.0/60.0
    return epsilon
