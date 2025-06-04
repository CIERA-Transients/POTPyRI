from potpyri import main_pipeline
from potpyri import options
from potpyri import logger
from potpyri import sort_files

import numpy as np
import astropy.utils.data
import importlib
import os

from tests.utils import download_github_file

def test_potpyri_sort(tmp_path):

    # Science file (just tellurics)
    cache_dir = astropy.utils.data._get_download_cache_loc()
    git_file_path = 'Binospec/raw/sci_img_2024.0812.034220_proc.fits.fz'
    file_path = download_github_file(git_file_path, use_cached=True)

    data_path, basefile = os.path.split(file_path)

    instrument = 'BINOSPEC'

    module = importlib.import_module(f'potpyri.instruments.{instrument.upper()}')
    tel = getattr(module, instrument.upper())()

    # Generate code and data paths based on input path
    paths = options.add_paths(data_path, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    # This contains all of the file data
    file_list = os.path.join(paths['data'], 'file_list.txt')
    file_table = sort_files.handle_files(file_list, tel, paths, 
        incl_bad=True, proc='proc', no_redo=False, log=log)

    # Run validation checks on file_table
    assert len(file_table)==1
    assert file_table[0]['Target']=='GRB240809A_r_ep1'
    assert file_table[0]['Filter']=='r'
    assert file_table[0]['Type']=='SCIENCE'
    assert file_table[0]['TargType']=='GRB240809A_r_ep1_r_2_11'
    assert file_table[0]['CalType']=='r_2_11'
    assert file_table[0]['Exp']=="120.0"
    assert file_table[0]['Binning']=="11"
    assert file_table[0]['Amp']=="2"

    # Test the no_redo flag
    new_file_table = sort_files.handle_files(file_list, tel, paths, 
        incl_bad=True, proc='proc', no_redo=True, log=log)

    assert new_file_table==file_table
