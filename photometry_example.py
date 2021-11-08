from psf import *
from absphot import *

def main():

    #do_phot('2020jfo.g.ut200621.051434_1.sw.fits')
    zp_cal = absphot()
    zp, zp_err = zp_cal.find_zeropoint('GRB130131A_Y_Y_wcs.pcmp', 'Y', '2MASS',
        plot=True)

main()
