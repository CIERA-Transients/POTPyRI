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

__version__ = "2.2" # Last updated 03/10/2025

import sys
import time
import numpy as np

# Internal dependenciess
from potpyri.instruments import instrument_getter
from potpyri.utils import logger
from potpyri.utils import options
from potpyri.primitives import sort_files
from potpyri.primitives import solve_wcs
from potpyri.primitives import photometry
from potpyri.primitives import absphot
from potpyri.primitives import calibration
from potpyri.primitives import image_procs

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
    tel = instrument_getter(instrument)

    if skip_flatten: tel.flat=False

    # Generate code and data paths based on input path
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])
    log.info(f'Running main pipeline version {__version__}')
    log.info(f'Running instrument paramater file version {tel.version}')

    # This contains all of the file data
    file_table = sort_files.handle_files(paths['filelist'], paths, tel,
        incl_bad=incl_bad, proc=proc, no_redo=no_redo_sort, log=log)

    # Calibrations
    ##############
    calibration.do_bias(file_table, tel, paths, log=log)
    calibration.do_dark(file_table, tel, paths, log=log)
    calibration.do_flat(file_table, tel, paths, log=log)

    # Begin processing
    ##################
    for tar in np.unique(file_table[tel.match_type_keywords(tel.filetype_keywords['SCIENCE'], file_table)]['TargType']):

        # Image pixel calibration, WCS, and stacking procedures
        #######################################################
        log.info(f'Generating stack for {tar}')
        stack = image_procs.image_proc(file_table[file_table['TargType']==tar], tel, paths,
            skip_skysub=skip_skysub, fieldcenter=fieldcenter, out_size=out_size,
            cosmic_ray=not skip_cr, skip_gaia=skip_gaia, keep_all_astro=keep_all_astro,
            log=log)

        # Photometry step
        #################
        log.info('Running PSF photometry.')
        photometry.photloop(stack, phot_sn_min=phot_sn_min,
            fwhm_init=fwhm_init, phot_sn_max=phot_sn_max, log=log)

        # Zero point/flux calibration step
        ##################################
        log.info('Calculating zeropint.')
        absphot.find_zeropoint(stack, tel, log=log)
        
    t2 = time.time()
    log.info('Pipeline finshed.')
    log.info(f'Total runtime: {t2-t1} sec')
    log.shutdown()

def main():
    options.test_for_dependencies()
    args = options.add_options()
    main_pipeline(**vars(args))

if __name__ == "__main__":
    main()
