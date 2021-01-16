#!/usr/bin/env python
#from specclass import Spec
import getpass
import logging
import copy
# need copy.deepcopy() for Spec class
import os
import pdb
import sys
import glob
import datetime
from astropy.io import fits
import tmath.pydux as pydux
from tmath.wombat.yesno import yesno
from tmath.wombat.inputter import inputter
from tmath.pydux.getfitsfile import getfitsfile
# needs extinction, astropy.io, scipy

import matplotlib
matplotlib.use('TkAgg')
matplotlib.rcParams["savefig.directory"] = "."
matplotlib.rcParams["xtick.minor.visible"] = True
matplotlib.rcParams["ytick.minor.visible"] = True
matplotlib.rcParams["lines.linewidth"] = 0.5

from matplotlib.widgets import Cursor

import matplotlib.pyplot as plt

CALVERSION = 0.1



def main():
    import headerfix
    
    secondord = False
    gratcode2 = ''
    # logging straight from docs.python.org cookbook page
    #  INFO level to screen and cal.log
    #  DEBUG level only to cal.log
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename='cal.log',
                        filemode='a')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    user = getpass.getuser()
    print('Hello, {} \n'.format(user))
    logging.debug('CAL starts')
    logging.debug('{} running things'.format(user))

    print(' ')
    print('This program is a driver for the calibration routines.')
    print('It expects files that have been dispersion calibrated,')
    print('usually through IRAF (you know, with a -d- in front).')
    print('(In other words, include the d in your file names.)')
    print(' ')
    print('Do you want to use second-order correction?\n')
    secondord = False
    # answer = yesno('n')
    # if (answer == 'y'):
    #     secondord = True
    #axarr=fig.subplots(2,sharex=True)
    #fig.subplots_adjust(hspace=0)

    print('\nWe now need a grating code, such as opt or ir2.')
    print('This will be used to keep track of the fluxstar and')
    print('bstar as in fluxstaropt.fits or bstarir2.fits\n')
    gratcode = inputter('Enter the grating code (blue/red): ', 'string', False)
    print(' ')
    if (secondord):
        gratcode2 = inputter('Enter the second-order grating code: ', 'string', False)
    print(' ')
    # print('Enter the file containing the list of objects')
    # print('(should be dispersion-corrected)\n')
    done = False
    while (not done):
        # objectlist = inputter('Object list file: ', 'string', False)
        dfiles=glob.glob('d*.fits')
        for d in dfiles:
            if gratcode in d:
                listfile= open(d.split('_ex.fits')[0],"w+")
                listfile.write(d)
                listfile.close()
                objectlist = listfile.name
        if (os.path.isfile(objectlist)):
            done = True
        else:
            print('No such file')
    print('\nDo you want to fit a flux star?\n')
    answer_flux = yesno('y')
    if (answer_flux == 'y'):
        plt.close()
        fluxfile = getfitsfile('flux star', '.fits', gratcode=gratcode)
        pydux.mkfluxstar(fluxfile, gratcode)
        if (secondord):
            print(' ')
            fluxfile2 = getfitsfile('second flux star', '.fits')
            pydux.mkfluxstar(fluxfile2, gratcode2)
    else:
        try:
            fluxfits = fits.open('fluxstar'+gratcode+'.fits')
            fluxfile = fluxfits[0].header['FLUXFILE']
            fluxfits.close()
        # except FileNotFoundError:
        except KeyError:
            fluxfile = 'Unknown'
        if (secondord):
            try:
                fluxfits2 = fits.open('fluxstar'+gratcode2+'.fits')
                fluxfile2 = fluxfits2[0].header['FLUXFILE']
                fluxfits2.close()
            # except FileNotFoundError:
            except KeyError:
                fluxfile2 = 'Unknown'

    print('\nDo you want to apply the flux star(s) to the data?\n')
    answer = yesno('y')
    if (answer == 'y'):
        pydux.calibrate(objectlist, gratcode, secondord, gratcode2, answer_flux=answer_flux)


    print('\nDo you want to fit a b-star?\n')
    answer = yesno('y')
    if (answer == 'y'):
        plt.close()
        print('\nDo you want to use the flux star {}'.format(fluxfile))
        print('as the b-star?\n')
        same = yesno('y')
        print(' ')
        if (same == 'n'):
            bfile = getfitsfile('b-star', '.fits')
        else:
            bfile = fluxfile
        bfile = 'c'+gratcode+bfile
        pydux.mkbstar(bfile, gratcode)
        if (secondord):
            print('\nDo you want to use the flux star {}'.format(fluxfile2))
            print('as the b-star?\n')
            same = yesno('y')
            if (same == 'n'):
                bfile2 = getfitsfile('b-star', '.fits')
            else:
                bfile2 = fluxfile2
            bfile2 = 'c'+gratcode+bfile2
            pydux.mkbstar(bfile2, gratcode2)

    print('\nFinal calibration (atmos. removal, sky line wave. adjust, etc.)?')
    answer = yesno('y')
    if (answer == 'y'):
        plt.close()
        pydux.final(objectlist, gratcode, secondord, gratcode2, user)

    if (answer_flux == 'y'):
        print('\nAre these master calibrations?\n')
        answer=yesno('y')
        if (answer == 'y'):
            os.system('mv ' + 'fluxstar' + gratcode + '.fits' + ' ../../master_files/')
            os.system('mv ' + 'bstar' + gratcode + '.fits'+ ' ../../master_files/')
        print('\nThere, was that so hard?')

if __name__ == '__main__':
    main()
