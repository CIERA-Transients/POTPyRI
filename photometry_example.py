from psf import *
from absphot import *

def main():

    #do_phot('2020jfo.g.ut200621.051434_1.sw.fits')
    zp_cal = absphot()
    zp, zp_err = zp_cal.find_zeropoint('sGRB191031D_z_z_left_wcs.pcmp', 'z', 'PS1',
        plot=True)

main()
