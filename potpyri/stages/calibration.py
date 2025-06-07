"Wrapper methods for creating bias, dark, and flat images from pipeline data."
"Authors: Charlie Kilpatrick"

# Initial version tracking on 09/29/2024
__version__ = "1.1"

import os
import time
import logging
import numpy as np
import sys

def do_bias(bias_table, tel, paths, log=None):

    # Exit if telescope does not require bias
    if not tel.bias:
        if log:
            log.info('No bias is required.')
        else:
            print('No bias is required.')
        return(None)

    # Exit if bias_table has no images
    if len(bias_table)==0:
        if log:
            log.info('No bias images were provided for this setup.')
        else:
            print('No bias images were provided for this setup.')
        return(None)

    bias_num = 0
    for cal_type in np.unique(bias_table['CalType']):
        mask = bias_table['CalType']==cal_type
        cal_table = bias_table[mask]
        
        bias_num += 1
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        bias_name = tel.get_mbias_name(paths, amp, binn)
        
        if os.path.exists(bias_name):
            if log: log.info(f'Master bias {bias_name} exists.')
        else:
            t1 = time.time()
            if log: log.info('Processing bias files.')
            tel.create_bias(cal_table['File'], amp, binn, paths, 
                log=log)
            t2 = time.time()

            if log: log.info(f'Master bias creation completed in {t2-t1} sec')

    if bias_num==0:
        if log: log.critical('No bias present, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)

def do_dark(dark_table, tel, paths, log=None):

    # Exit if telescope does not require dark
    if not tel.dark:
        return(None)

    # Exit if dark_table has no images
    if len(dark_table)==0:
        if log:
            log.info('No dark images were provided for this setup.')
        else:
            print('No dark images were provided for this setup.')
        return(None)

    for cal_type in np.unique(dark_table['CalType']):
        mask = dark_table['CalType']==cal_type
        cal_table = dark_table[mask]
            
        exp = cal_table['Exp'][0]
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        dark_name = tel.get_mdark_name(paths, amp, binn)

        if os.path.exists(dark_name):
            if log: log.info(f'Master dark {dark_name} exists.')
        else:
            mbias = None
            if tel.bias:
                if log: log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(paths, amp, binn)
                except:
                    if log: log.error(f''''No master bias found for this 
                        configuration, skipping master dark creation for 
                        exposure {exp}, {amp} amps and {binn} binning.''')
                    continue

            t1 = time.time()
            tel.create_dark(cal_table['File'], amp, binn,
                paths, mbias=mbias, log=log)
            t2 = time.time()
            if log: log.info(f'Master dark creation completed in {t2-t1} sec.')

def do_flat(flat_table, tel, paths, log=None):

    # Exit if telescope does not require dark
    if not tel.flat:
        return(None)

    # Exit if flat_table has no images
    if len(flat_table)==0:
        if log:
            log.info('No flat images were provided for this setup.')
        else:
            print('No flat images were provided for this setup.')
        return(None)

    for cal_type in np.unique(flat_table['CalType']):
        mask = flat_table['CalType']==cal_type
        cal_table = flat_table[mask]

        fil = cal_table['Filter'][0]
        amp = cal_table['Amp'][0]
        binn = cal_table['Binning'][0]

        is_science = np.any([f=='SCIENCE' for f in cal_table['Type']])

        flat_name = tel.get_mflat_name(paths, fil, amp, binn)

        if os.path.exists(flat_name):
            if log: log.info(f'Master flat {flat_name} exists.')
        else:
            mbias = None
            mdark = None
            if tel.bias:
                if log: log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(paths, amp, binn)
                except:
                    if log: log.error(f'''No master bias found for this 
                        configuration, skipping master flat creation for 
                        filter {fil}, {amp} amps, {binn} binning.''')
                    continue

            if tel.dark:
                if log: log.info('Loading master dark.')
                try:
                    mdark = tel.load_dark(paths, amp, binn)
                except:
                    if log: log.error(f''''No master dark found for this 
                        configuration, skipping master flat creation for 
                        filter {fil}, {amp} amps and {binn} binning.''')
                    continue

            t1 = time.time()
            tel.create_flat(cal_table['File'], fil, amp, binn,
                paths, mbias=mbias, mdark=mdark, is_science=is_science, 
                log=log)
            t2 = time.time()            
            if log: log.info(f'Master flat creation completed in {t2-t1} sec')

