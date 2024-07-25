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

# Internal dependencies
from custom_logger import get_log
from utilities import util
import Sort_files
import solve_wcs
import psf
import absphot
import cal_procs
import image_procs
import Find_target_phot as tp

def main_pipeline(telescope:str,
                  data_path:str,
                  input_target:list=None,
                  skip_red:bool=None,
                  proc:str=None,
                  reset:bool=None,
                  incl_bad:bool=None,
                  no_redo_sort:bool=None,
                  file_list_name:str='file_list.txt')->None:

    # start time
    t1 = time.time()
    # import telescope parameter file
    global tel
    try:
        tel = importlib.import_module(f'params.{telescope}')
    except ImportError:
        telfile=f'params.{telescope}'
        print(f'''No such telescope file {telfile}, please check that you have 
            entered the correct name or this telescope is available.''')
        sys.exit(-1)

    # Get path to code directory
    paths = {'data': data_path}
    paths['abspath']=os.path.abspath(sys.argv[0])
    paths['code']=os.path.split(paths['abspath'])[0]
    paths['config']=os.path.join(paths['code'], 'config')
    paths['raw']=os.path.join(data_path, 'raw') #path containing the raw data
    paths['bad']=os.path.join(data_path, 'bad') #path containing the unused data
    paths['red']=os.path.join(data_path, 'red') #path to write the reduced files
    for key in paths.keys():
        if not os.path.exists(paths[key]): os.makedirs(paths[key])

    # Copy config directory to data path
    if not os.path.exists(os.path.join(paths['data'], 'config')):
        shutil.copytree(paths['config'], os.path.join(paths['data'], 'config'))

    # Generate logfile
    log = get_log(paths['red'])

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

    # Master bias, dark, and flat creation (will skip if unnecessary)
    bias_files = file_table[file_table['Type']=='BIAS']
    flat_files = file_table[file_table['Type']=='FLAT']
    dark_files = file_table[file_table['Type']=='DARK']
    cal_procs.do_bias(bias_files, tel, paths['red'], skip_red=skip_red, log=log)
    cal_procs.do_dark(dark_files, tel, paths['red'], skip_red=skip_red, log=log)
    cal_procs.do_flat(flat_files, tel, paths['red'], skip_red=skip_red, log=log)

    # Science files
    ##################
    # Basic processing
    # Does it even need to run?
    science_data = file_table[file_table['Type']=='SCIENCE']
    if len(science_data)==0:
        log.critical('No science files, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)

    if input_target is not None:
        log.info(f'User input target for reduction: {input_target}')

    # Begin processing
    for tar in np.unique(science_data['CalType']):
        mask = science_data['CalType']==tar
        target_table = science_data[mask]

        if input_target is not None:
            if target_table['Target'][0]!=input_target: continue

        files = target_table['File']
        stack = image_procs.image_proc(target_table, tel, paths['red'], 
            proc=proc, log=log)

        # Auto-WCS solution
        if tel.run_wcs():
            log.info('Solving WCS.')
            solve_wcs.solve_astrometry(paths['red'], stack, tel, log=log)

        # Photometry and flux calibration
        if tel.run_phot():
            log.info('Running psf photometry.')
            for snthresh_final in [5.0, 10.0, 20.0]:
                try:
                    log.info(f'Trying photometry with final S/N={snthresh_final}')
                    star_param={'snthresh_psf': snthresh_final*2.0, 
                                'fwhm_init': 5.0, 
                                'snthresh_final': snthresh_final}
                    psf.do_phot(stack, star_param=star_param, log=log)
                    break
                except:
                    pass

            log.info('Calculating zeropint.')
            photfile = stack.replace('.fits','.pcmp')
            if os.path.exists(photfile):
                cat = tel.catalog_zp()
                cal = absphot.absphot()
                cal.find_zeropoint(photfile, target_table['Filter'][0], cat, 
                    log=log)
            else:
                log.error(f'Could not generate PSF photometry for {stack}')
        
    t2 = time.time()
    log.info('Pipeline finshed.')
    log.info(f'Total runtime: {t2-t1} sec')


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
    params.add_argument('--skip_red', 
        type=str, 
        default=None, 
        help='''Option to skip reduction i.e. the creation of master files used 
        for calibration or the stacked image for science targets. See manual 
        for specific details. Optional parameter.''')
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
    params.add_argument('--reset', 
        type=str, 
        default=None, 
        help='''Option to reset data files sorted previously by pipeline. 
        Optional parameter.''')
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
    args = params.parse_args()

    main_pipeline(args.instrument, args.data_path, input_target=args.target, 
        skip_red=args.skip_red, proc=args.proc, reset=args.reset,
        no_redo_sort=args.no_redo_sort,
        file_list_name=args.file_list_name)

if __name__ == "__main__":
    main()
