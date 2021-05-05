from psf import *
from absphot import *

def main():

    do_phot('2020jfo.g.ut200621.051434_1.sw.fits')
    ap = absphot()
    ap.find_zeropoint('2020jfo.g.ut200621.051434_1.sw.pcmp', 'g', 'PS1')

main()
