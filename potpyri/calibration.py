import os
import time
import logging
import numpy as np
import sys

def do_bias(bias_table, tel, red_path, staticmask=None, log=None):

    # Exit if telescope does not require bias
    if not tel.bias:
        return(None)

    bias_num = 0
    for cal_type in np.unique(bias_table['CalType']):
        mask = bias_table['CalType']==cal_type
        cal_table = bias_table[mask]
        
        bias_num += 1
        amp = cal_type.split('_')[0]
        binn = cal_type.split('_')[1]

        bias_name = tel.get_mbias_name(red_path, amp, binn)
        
        if os.path.exists(bias_name):
            if log: log.info(f'Master bias {bias_name} exists.')
        else:
            t1 = time.time()
            if log: log.info('Processing bias files.')
            tel.create_bias(cal_table['File'], amp, binn, red_path, 
                staticmask=staticmask, log=log)
            t2 = time.time()

            if log: log.info(f'Master bias creation completed in {t2-t1} sec')

    if bias_num==0:
        if log: log.critical('No bias present, check data before rerunning.')
        logging.shutdown()
        sys.exit(-1)

def do_dark(dark_table, tel, red_path, staticmask=None, log=None):

    # Exit if telescope does not require dark
    if not tel.dark:
        return(None)

    for cal_type in np.unique(dark_table['CalType']):
        mask = dark_table['CalType']==cal_type
        cal_table = dark_table[mask]
            
        exp = cal_type.split('_')[0]
        amp = cal_type.split('_')[1]
        binn = cal_type.split('_')[2]

        dark_name = tel.get_mdark_name(red_path, amp, binn)

        if os.path.exists(dark_name):
            if log: log.info(f'Master dark {dark_name} exists.')
        else:
            mbias = None
            if tel.bias:
                if log: log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(red_path, amp, binn)
                except:
                    if log: log.error(f''''No master bias found for this 
                        configuration, skipping master dark creation for 
                        exposure {exp}, {amp} amps and {binn} binning.''')
                    continue

            t1 = time.time()
            tel.create_dark(cal_table['File'], amp, binn,
                red_path, mbias=mbias, staticmask=staticmask, log=log)
            t2 = time.time()
            if log: log.info(f'Master dark creation completed in {t2-t1} sec.')

def do_flat(flat_table, tel, red_path, staticmask=None, log=None):

    # Exit if telescope does not require dark
    if not tel.flat:
        return(None)

    for cal_type in np.unique(flat_table['CalType']):
        mask = flat_table['CalType']==cal_type
        cal_table = flat_table[mask]

        fil = flat_table[mask]['Filter'][0]
        amp = flat_table[mask]['Amp'][0]
        binn = flat_table[mask]['Binning'][0]
        is_science = flat_table[mask]['Type'][0]=='SCIENCE'

        flat_name = tel.get_mflat_name(red_path, fil, amp, binn)

        if os.path.exists(flat_name):
            if log: log.info(f'Master flat {flat_name} exists.')
        else:
            mbias = None
            mdark = None
            if tel.bias:
                if log: log.info('Loading master bias.')
                try:
                    mbias = tel.load_bias(red_path, amp, binn)
                except:
                    if log: log.error(f'''No master bias found for this 
                        configuration, skipping master flat creation for 
                        filter {fil}, {amp} amps, {binn} binning.''')
                    continue

            if tel.dark:
                if log: log.info('Loading master dark.')
                try:
                    mdark = tel.load_dark(red_path, amp, binn)
                except:
                    if log: log.error(f''''No master dark found for this 
                        configuration, skipping master flat creation for 
                        filter {fil}, {amp} amps and {binn} binning.''')
                    continue

            t1 = time.time()
            tel.create_flat(cal_table['File'], fil, amp, binn,
                red_path, mbias=mbias, mdark=mdark, is_science=is_science, 
                staticmask=staticmask, log=log)
            t2 = time.time()            
            if log: log.info(f'Master flat creation completed in {t2-t1} sec')

