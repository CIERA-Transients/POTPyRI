"""
Automatic generalized pipeline for imaging reduction. Creates coadded images 
for each target.  Files are automatically sorted and categorized based on 
expected header keywords.  Pixel-level processing is performed with ccdproc, 
WCS alignment is performed with astrometry.net and astropy/photutils, 
photometry is performed with photutils, and flux calibration with astroquery
and numpy.

Authors: Kerry Paterson, Charlie Kilpatrick
This project was funded with support by the National Science Foundation under
grant Nos. AST-1814782, AST-1909358 and CAREER grant No. AST-2047919.
If you use this code for your work, please consider citing the package release.
"""

# Last updated 26/06/2025
__version__ = "2.2"

import sys
import os
import time
import logging
import importlib
import re
import numpy as np

from astropy.io import fits

# Internal dependenciess
from potpyri.utils import logger
from potpyri.utils import options
from potpyri.stages import sort_files
from potpyri.stages import solve_wcs
from potpyri.stages import photometry
from potpyri.stages import absphot
from potpyri.stages import calibration
from potpyri.stages import image_procs

# Check options.py - all named parameters need to correspond to options that are
# provided via argparse.
def main_pipeline(instrument:str,
                  data_path:str,
                  target:list=None,
                  proc:str=None,
                  incl_bad:bool=None,
                  no_redo_sort:bool=None,
                  phot_sn_min:float=None,
                  phot_sn_max:float=None,
                  fwhm_init:float=None,
                  skip_skysub:bool=None,
                  file_list_name:str=None,
                  fieldcenter:list=None,
                  out_size:int=None,
                  skip_flatten:bool=None,
                  skip_cr:bool=None,
                  skip_gaia:bool=None,
                  keep_all_astro:bool=None,
                  **kwargs)->None:

    # start time
    t1 = time.time()
    
    # import telescope parameter file
    module = importlib.import_module(f'potpyri.instruments.{instrument.upper()}')
    tel = getattr(module, instrument.upper())()

    if skip_flatten: tel.flat=False

    # Generate code and data paths based on input path
    paths = options.add_paths(data_path, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    if log: log.info(f'Running main pipeline version {__version__}')
    if log: log.info(f'Running instrument paramater file version {tel.version}')

    # This contains all of the file data
    file_list = os.path.join(paths['data'], file_list_name)
    file_table = sort_files.handle_files(file_list, tel, paths, 
        incl_bad=incl_bad, proc=proc, no_redo=no_redo_sort, log=log)

    if not file_table or len(file_table)==0:
        log.critical('Could not generate file table.  Check file paths.')
        logging.shutdown()
        sys.exit(-1)

    # Calibration images
    ####################
    # Master bias, dark, and flat creation (will skip if unnecessary)
    kwds = tel.filetype_keywords
    bias_match = tel.match_type_keywords(kwds['BIAS'], file_table)
    flat_match = tel.match_type_keywords(kwds['FLAT'], file_table)
    dark_match = tel.match_type_keywords(kwds['DARK'], file_table)
    science_match = tel.match_type_keywords(kwds['SCIENCE'], file_table)
    bias_files = file_table[bias_match]
    flat_files = file_table[flat_match]
    dark_files = file_table[dark_match]
    science_data = file_table[science_match]

    calibration.do_bias(bias_files, tel, paths, log=log)
    calibration.do_dark(dark_files, tel, paths, log=log)

    # If there are no flats, use caldb flat instead of nightly flat, otherwise
    # make flat as normal
    if tel.flat and len(flat_files)==0:
        paths['cal'] = paths['caldb']
    else:
        calibration.do_flat(flat_files, tel, paths, log=log)

    # Science images
    ################
    # Basic processing
    # Does it even need to run?
    if len(science_data)==0:
        if log: log.critical('No science files, check data before re-running.')
        logging.shutdown()
        sys.exit(-1)

    # Begin processing
    for tar in np.unique(science_data['TargType']):
        mask = science_data['TargType']==tar
        target_table = science_data[mask]

        if target is not None:
            if target_table['Target'][0]!=target: continue
        else:
            if log: log.info(f'User input target for reduction: {target}')

        files = target_table['File']
        cosmic_ray = not skip_cr
        stack = image_procs.image_proc(target_table, tel, paths,
            skip_skysub=skip_skysub, fieldcenter=fieldcenter, out_size=out_size,
            cosmic_ray=cosmic_ray, skip_gaia=skip_gaia, keep_all_astro=keep_all_astro,
            log=log)

        if stack is None:
            if log: log.error(f'Could not generate a stack for {tar}')
            continue

        # Photometry and flux calibration
        if log: log.info('Running psf photometry loop for zero point.')
        photometry.photloop(stack, phot_sn_min=phot_sn_min,
            fwhm_init=fwhm_init, phot_sn_max=phot_sn_max, log=log)

        if log: log.info('Calculating zeropint.')
        hdu = fits.open(stack)
        cat = tel.get_catalog(hdu[0].header)
        cal = absphot.absphot()
        try:
            cal.find_zeropoint(stack, target_table['Filter'][0], cat, log=log)
        except:
            if log: log.error(f'Could not calibrate photometry for {stack}')
        
    if log: 
        t2 = time.time()
        log.info('Pipeline finshed.')
        log.info(f'Total runtime: {t2-t1} sec')

def main():
    options.test_for_dependencies()
    args = options.add_options()
    main_pipeline(**vars(args))

if __name__ == "__main__":
    main()
