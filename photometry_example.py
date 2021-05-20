from psf import *
from absphot import *

def main():

    #do_phot('2020jfo.g.ut200621.051434_1.sw.fits')
    ap = absphot()
    ap.find_zeropoint('sGRB140930B_J_J_wcs.pcmp', 'J', '2MASS')

main()
