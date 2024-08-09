#!/usr/bin/env python

"""Automatic generalized pipeline for imaging reduction. Creates median coadded 
images for each target. Individual images should be checked and removed from 
sci path.

Author: Kerry Paterson
This project was funded by AST 
If you use this code for your work, please consider citing."""

# Refactor on 2024-07-19
__version__ = "1.22"

import sys
import os
import time
import shutil
import logging
import importlib
import numpy as np

from astropy.io import fits

# Internal dependencies
from custom_logger import get_log
from utilities import util
import Sort_files
import solve_wcs
import photometry
import absphot
import calibration
import image_procs

def main_pipeline(telescope:str,
                  data_path:str,
                  input_target:list=None,
                  proc:str=None,
                  incl_bad:bool=None,
                  no_redo_sort:bool=None,
                  phot_sn_min:float=None,
                  phot_sn_max:float=None,
                  fwhm_init:float=None,
                  file_list_name:str=None)->None:

    # start time
    t1 = time.time()
    
    # import telescope parameter file
    tel = importlib.import_module(f'params.{telescope}')

    # Get path to code directory
    paths = {'data': os.path.abspath(data_path)}
    paths['abspath']=os.path.abspath(sys.argv[0])
    paths['code']=os.path.split(paths['abspath'])[0]
    paths['config']=os.path.join(paths['code'], 'config')
    paths['raw']=os.path.join(data_path, 'raw') #path containing the raw data
    paths['bad']=os.path.join(data_path, 'bad') #path containing the unused data
    paths['red']=os.path.join(data_path, 'red') #path to write the reduced files
    paths['log']=os.path.join(data_path, 'red', 'log')
    paths['cal']=os.path.join(data_path, 'red', 'cals')
    paths['work']=os.path.join(data_path, 'red', 'workspace')
    for key in paths.keys():
        if not os.path.exists(paths[key]): os.makedirs(paths[key])

    # Copy config directory to data path
    if not os.path.exists(os.path.join(paths['data'], 'config')):
        shutil.copytree(paths['config'], os.path.join(paths['data'], 'config'))

    # Generate logfile
    log = get_log(paths['log'])

    log.info(f'Running main pipeline version {__version__}')
    log.info(f'Running telescope paramater file version {tel.__version__}')

    # CDK - added editing for headers with LRIS raw data files (*[b,r]*.fits)
    # TODO: Refactor so this is a method that happens for every telescope,
    # rename edit_headers with proc as a variable, then call every time
    if ((telescope=='LRIS' and proc and str(proc)=='raw') or
        (telescope=='DEIMOS' and proc and str(proc)=='raw')):
        log.info('Edit raw headers')
        tel.edit_raw_headers(paths['data'])
        tel.edit_raw_headers(paths['raw'])

    file_list = os.path.join(paths['data'], file_list_name)

    # This contains all of the file data
    file_table = Sort_files.handle_files(file_list, tel, paths, 
        incl_bad=incl_bad, proc=proc, no_redo=no_redo_sort, log=log)

    if not file_table or len(file_table)==0:
        log.critical('Could not generate file table.  Check file paths.')
        logging.shutdown()
        sys.exit(-1)

    # Calibration images
    ####################
    # Master bias, dark, and flat creation (will skip if unnecessary)
    bias_files = file_table[file_table['Type']=='BIAS']
    flat_files = file_table[file_table['Type']=='FLAT']
    dark_files = file_table[file_table['Type']=='DARK']
    calibration.do_bias(bias_files, tel, paths['cal'], log=log)
    calibration.do_dark(dark_files, tel, paths['cal'], log=log)
    calibration.do_flat(flat_files, tel, paths['cal'], log=log)

    # Science images
    ################
    # Basic processing
    # Does it even need to run?
    science_data = file_table[file_table['Type']=='SCIENCE']
    if len(science_data)==0:
        if log: log.critical('No science files, check data before re-running.')
        logging.shutdown()
        sys.exit(-1)

    # Begin processing
    for tar in np.unique(science_data['CalType']):
        mask = science_data['CalType']==tar
        target_table = science_data[mask]

        if input_target is not None:
            if target_table['Target'][0]!=input_target: continue
        else:
            if log: log.info(f'User input target for reduction: {input_target}')

        files = target_table['File']
        stack = image_procs.image_proc(target_table, tel, paths, proc=proc, 
            log=log)

        # Auto-WCS solution
        if log: log.info('Solving WCS.')
        solve_wcs.solve_astrometry(paths['red'], stack, tel, log=log)

        # Photometry and flux calibration
        if log: log.info('Running psf photometry loop for zero point.')
        epsf, fwhm=photometry.photloop(stack, phot_sn_min=phot_sn_min,
            fwhm_init=fwhm_init, phot_sn_max=phot_sn_max, log=log)

        log.info('Calculating zeropint.')
        photfile = stack.replace('.fits','.pcmp')
        if os.path.exists(photfile):
            hdu = fits.open(stack)
            cat = tel.catalog_zp(hdu[0].header)
            cal = absphot.absphot()
            cal.find_zeropoint(photfile, target_table['Filter'][0], cat, 
                log=log)
        else:
            if log: log.error(f'Could not generate PSF photometry for {stack}')
        
    t2 = time.time()
    if log: log.info('Pipeline finshed.')
    if log: log.info(f'Total runtime: {t2-t1} sec')

def main():
    import argparse
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('instrument', 
        default=None, 
        choices=['Binospec','DEIMOS','GMOS','LRIS','MMIRS','MOSFIRE','TEST'],
        help='''Name of instrument (must be in params folder) of data to 
        reduce. Required to run pipeline. Use TEST for pipeline test.''')
    params.add_argument('data_path', 
        default=None, 
        help='''Path of data to reduce. See manual for specific details. 
        Required to run pipeline.''')
    params.add_argument('--target', 
        type=str, 
        default=None, 
        help='''Option to only reduce a specific target. String used here must 
        be contained within the target name in file headers. Optional 
        parameter.''')
    params.add_argument('--proc', 
        type=str, 
        default=True, 
        help='''Option to use the _proc data from MMT. Optional parameter. 
        Default is False.''')
    params.add_argument('--include-bad','--incl-bad', 
        default=False,
        action='store_true', 
        help='''Include files flagged as bad in the file list.''')
    params.add_argument('--no-redo-sort',
        default=False,
        action='store_true', 
        help='''Do not resort files and generate new file list.''')
    params.add_argument('--file-list-name', 
        type=str, 
        default='file_list.txt', 
        help='''Change the name of the archive file list.''')
    params.add_argument('--phot-sn-min',
        type=float,
        default=3.0,
        help='''Minimum signal-to-noise to try in photometry loop.''')
    params.add_argument('--phot-sn-max',
        type=float,
        default=40.0,
        help='''Maximum signal-to-noise to try in photometry loop.''')
    params.add_argument('--fwhm-init',
        type=float,
        default=5.0,
        help='''Initial FWHM (in pixels) for photometry loop.''')
    args = params.parse_args()

    main_pipeline(args.instrument, args.data_path, input_target=args.target, 
        proc=args.proc, no_redo_sort=args.no_redo_sort, 
        file_list_name=args.file_list_name, phot_sn_min=args.phot_sn_min, 
        phot_sn_max=args.phot_sn_max, fwhm_init=args.fwhm_init)

if __name__ == "__main__":
    main()
