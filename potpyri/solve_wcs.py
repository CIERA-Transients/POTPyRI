#!/usr/bin/env python

"Python script for WCS solution."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "3.12" #last updated 01/10/2021

import numpy as np
import time
import astropy
import argparse
import subprocess
import os
import astropy.units as u
import astropy.wcs as wcs

from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from scipy.optimize import curve_fit
from utilities import util

#turn Astropy warnings off
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

def run_sextractor(input_file, cat_name, tel, sex_config_dir='./Config', 
    log=None):

    if not os.path.exists(os.path.join(sex_config_dir,'config')):
        e = 'ERROR: Could not find source extractor config file!'
        if log:
            log.error(e)
        raise Exception(e)
    elif not os.path.exists(os.path.join(sex_config_dir,'params')):
        e = 'ERROR: Could not find source extractor param file!'
        if log:
            log.error(e)
        raise Exception(e)

    config = os.path.join(sex_config_dir,'config')
    param_file = os.path.join(sex_config_dir,'params')

    cmd=['sex','-c',config,input_file,'-CATALOG_NAME',cat_name]
    subprocess.call(cmd)

    params = []
    with open(param_file) as f:
        for line in f:
            params.append(line.split()[0].strip())

    if os.path.exists(cat_name):
        table = ascii.read(cat_name, names=params, comment='#')
        if log:
            log.info('SExtractor succesfully run, %d sources extracted.'%len(table))
        return(table)
    else:
        e = 'ERROR: source extractor did not successfully produce {0}'
        if log:
            log.error(e)
        raise Exception(e.format(cat_name))

def clean_up_astrometry(directory, file, exten):
    filelist = [file.replace(exten,'.axy'),
                file.replace(exten,'.corr'),
                file.replace(exten,'-indx.xyls'),
                file.replace(exten,'.match'),
                file.replace(exten,'.rdls'),
                file.replace(exten,'.solved'),
                file.replace(exten,'.wcs')]

    filelist = [os.path.join(directory, f) for f in filelist]

    for f in filelist:
        if os.path.exists(f):
            os.remove(f)

def solve_astrometry(directory, file, tel, radius=0.5, replace=True, 
    shift_only=True, log=None):

    # Starting solve, print file and directory for reference
    if os.path.basename(file)==file:
        fullfile = os.path.join(directory, file)
    else:
        fullfile = file

    if log: 
        log.info(f'Trying to solve file: {file}')
    else:
        print(f'Trying to solve file: {file}')

    if not os.path.exists(fullfile):
        return(False)

    hdu = fits.open(fullfile)
    data = hdu[0].data
    header = hdu[0].header
    hkeys = list(header.keys())

    exten = '.'+file.split('.')[-1]
    if not replace:
        if os.path.exists(fullfile.replace(exten,'.solved.fits')):
            if log: 
                log.info(f'SUCCESS: solved {fullfile}')
            else:
                print(f'SUCCESS: solved {fullfile}')
            return(True)

    exten = '.'+file.split('.')[-1]

    check_pairs = [('RA','DEC'),('CRVAL1','CRVAL2'),('OBJCTRA','OBJCTDEC')]
    coord = None

    for pair in check_pairs:
        if pair[0] in header.keys() and pair[1] in header.keys():
            ra = header[pair[0]]
            dec = header[pair[1]]
            coord = util.parse_coord(ra, dec)
            if coord:
                break

    if not coord:
        if log: log.error(f'Could not parse RA/DEC from header of {file}')
        return(False)

    if 'solved.fits' in fullfile:
        newfile = 'tmp.fits'
    else:
        newfile = fullfile.replace(exten,'.solved.fits')

    # Handle pixel scale guess
    scale = tel.pixscale()
    scale_high = float('%.4f'%(scale * 1.2))
    scale_low = float('%.4f'%(scale * 0.8))

    # Get RA and Dec
    ra = float('%.6f'%coord.ra.degree)
    dec = float('%.6f'%coord.dec.degree)

    cmd = 'solve-field'
    args = '--scale-units arcsecperpix '
    args += f'--scale-low {scale_low} --scale-high {scale_high} '
    args += f'--ra {ra} --dec {dec} '
    args += f' --radius {radius} --no-plots '
    args += f'--overwrite -N {newfile} --dir {directory} '

    extra_opts = '--downsample 2 --no-verify --odds-to-tune-up 1e4 --objs 15'

    tries = 1
    good = False
    while tries < 4 and not good:
        input_args = args + extra_opts

        if log: 
            log.info(f'Try #{tries} with astrometry.net...')
        else:
            print(f'Try #{tries} with astrometry.net...')
        if log: 
            log.info(input_args)
        else:
            print(input_args)

        process = [cmd,fullfile]+input_args.split()

        FNULL = open(os.devnull, 'w')

        p = subprocess.Popen(process, stdout=FNULL, stderr=subprocess.STDOUT)
        try:
            p.wait(30)
        except subprocess.TimeoutExpired:
            p.kill()

        if os.path.exists(newfile):
            good = True
        else:
            tries += 1
            if 'downsample' in extra_opts:
                extra_opts='--objs 15'
            else:
                extra_opts=''

    file_exists=os.path.exists(newfile)
    if log: 
        log.info(f'{newfile} exists: {file_exists}')
    else:
        print(f'{newfile} exists: {file_exists}')

    if os.path.exists(newfile):

        clean_up_astrometry(directory, file, exten)
        if log: 
            log.info(f'SUCCESS: solved {fullfile}')
        else:
            print(f'SUCCESS: solved {fullfile}')

        if replace or newfile=='tmp.fits':
            output_file = fullfile

            # astrometry.net replaces extra extensions, so instead of renaming
            # copy the new header into the old file
            newhdu = fits.open(newfile)
            hdu = fits.open(output_file)

            if shift_only:
                hdu[0].header['CRPIX1'] = newhdu[0].header['CRPIX1']
                hdu[0].header['CRPIX2'] = newhdu[0].header['CRPIX2']
                hdu[0].header['CRVAL1'] = newhdu[0].header['CRVAL1']
                hdu[0].header['CRVAL2'] = newhdu[0].header['CRVAL2']
            else:
                hdu[0].header = newhdu[0].header
                
            hdu.writeto(output_file, overwrite=True, output_verify='silentfix')
            os.remove(newfile)
        else:
            output_file = fullfile.replace(exten,'.solved.fits')

        hdu = fits.open(output_file)
        for i,h in enumerate(hdu):

            if 'COMMENT' in hdu[i].header.keys():
                del hdu[i].header['COMMENT']
            if 'HISTORY' in hdu[i].header.keys():
                del hdu[i].header['HISTORY']

            # Delete other WCS header keys
            wcs_key = ['CSYER1','CSYER2','CRDER1','CRDER2','CD1_1','CD1_2',
                'CD2_1','CD2_2','CRPIX1','CRPIX2','CUNIT1','CUNIT2','EQUINOX',
                'RADESYS','CNAME1','CNAME2','CTYPE1','CTYPE2','WCSNAME',
                'CRVAL1','CRVAL2']

            for key in [w + 'C' for w in wcs_key]:
                if key in list(hdu[i].header.keys()):
                    del hdu[i].header[key]

            for key in [w + 'S' for w in wcs_key]:
                if key in list(hdu[i].header.keys()):
                    del hdu[i].header[key]

            # Delete old WCS keywords starting with '_'
            for key in list(hdu[i].header.keys()):
                if key.startswith('_'):
                    del hdu[i].header[key]

        try:
            hdu.writeto(output_file, overwrite=True, output_verify='silentfix')
        except TypeError:
            if log: 
                log.error(f'FAILURE: could not save file {fullfile}')
            else:
                print(f'FAILURE: could not save file {fullfile}')
            return(False)

        return(True)

    else:
        if log: 
            log.error(f'FAILURE: did not solve {fullfile}')
        else:
            print(f'FAILURE: did not solve {fullfile}')
        clean_up_astrometry(directory, file, exten)
        return(False)

if __name__ == "__main__":
    pass
