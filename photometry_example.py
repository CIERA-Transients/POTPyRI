from psf import *
from absphot import *

def main():

    do_phot('FRB190523_J_J_wcs.fits')
    ap = absphot()
    ap.find_zeropoint('FRB190523_J_J_wcs.pcmp', 'J', '2MASS')

main()
