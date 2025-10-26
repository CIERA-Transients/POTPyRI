from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import sort_files
from potpyri.instruments import instrument_getter

import os

from tests.utils import download_gdrive_file

def test_sort(tmp_path):

    instrument = 'BINOSPEC'
    file_list_name = 'files.txt'

    # Raw science file
    file_path = download_gdrive_file('Binospec/raw/sci_img_2024.0812.034220_proc.fits.fz', use_cached=True)

    data_path, basefile = os.path.split(file_path)
    data_path, _ = os.path.split(data_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    # This contains all of the file data
    file_table = sort_files.handle_files(paths['filelist'], paths, tel,
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
    new_file_table = sort_files.handle_files(paths['filelist'], paths, tel,
        incl_bad=True, proc='proc', no_redo=True, log=log)

    assert len(new_file_table)==1
    assert file_table[0]['Target']==new_file_table[0]['Target']
    assert file_table[0]['Filter']==new_file_table[0]['Filter']
    assert file_table[0]['Type']==new_file_table[0]['Type']
    assert file_table[0]['TargType']==new_file_table[0]['TargType']
    assert file_table[0]['CalType']==new_file_table[0]['CalType']
    assert file_table[0]['Exp']==new_file_table[0]['Exp']
    assert file_table[0]['Binning']==new_file_table[0]['Binning']
    assert file_table[0]['Amp']==new_file_table[0]['Amp']
