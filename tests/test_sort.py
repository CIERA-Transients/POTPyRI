"""Tests for sort_files.handle_files and file table content (BINOSPEC), and classifiers (is_bad, is_spec, is_flat, is_dark, is_bias, is_science)."""
from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import sort_files
from potpyri.instruments import instrument_getter

import os

from astropy.io import fits

import pytest
from tests.utils import download_gdrive_file


@pytest.mark.integration
def test_sort(tmp_path):
    """Run handle_files on BINOSPEC proc file; check file table rows and no_redo reuse."""
    instrument = 'BINOSPEC'
    file_list_name = 'files.txt'

    # Raw science file (download to tmp_path so log is writable)
    file_path = download_gdrive_file('Binospec/raw/sci_img_2024.0812.034220_proc.fits.fz', output_dir=str(tmp_path), use_cached=True)

    data_path, basefile = os.path.split(file_path)
    data_path, _ = os.path.split(data_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    # This contains all of the file data
    try:
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
    finally:
        log.close()

    assert len(new_file_table)==1
    assert file_table[0]['Target']==new_file_table[0]['Target']
    assert file_table[0]['Filter']==new_file_table[0]['Filter']
    assert file_table[0]['Type']==new_file_table[0]['Type']
    assert file_table[0]['TargType']==new_file_table[0]['TargType']
    assert file_table[0]['CalType']==new_file_table[0]['CalType']
    assert file_table[0]['Exp']==new_file_table[0]['Exp']
    assert file_table[0]['Binning']==new_file_table[0]['Binning']
    assert file_table[0]['Amp']==new_file_table[0]['Amp']


def test_is_bad_binospec():
    """is_bad True when header matches bad_keywords/bad_values (BINOSPEC: MASK=mira)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'mira'
    # No CCDSUM so binning check does not overwrite; keyword match gives bad=True
    assert sort_files.is_bad(hdr, tel) == True
    hdr['MASK'] = 'imaging'
    hdr['CCDSUM'] = '1,1'  # valid binning so keyword result (False) is kept
    assert sort_files.is_bad(hdr, tel) == False


def test_is_spec_binospec():
    """is_spec True when header matches spec keywords (BINOSPEC: MASK=spectroscopy)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'spectroscopy'
    assert sort_files.is_spec(hdr, tel) == True
    hdr['MASK'] = 'imaging'
    assert sort_files.is_spec(hdr, tel) == False


def test_is_flat_binospec():
    """is_flat True when MASK=imaging, SCRN=deployed (BINOSPEC)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'imaging'
    hdr['SCRN'] = 'deployed'
    hdr['CCDSUM'] = '1,1'
    assert sort_files.is_flat(hdr, tel) == True
    hdr['SCRN'] = 'stowed'
    assert sort_files.is_flat(hdr, tel) == False


def test_is_science_binospec():
    """is_science True when MASK=imaging, SCRN=stowed and exptime >= min (BINOSPEC)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'imaging'
    hdr['SCRN'] = 'stowed'
    hdr['CCDSUM'] = '1,1'
    hdr['EXPTIME'] = 120.0
    assert sort_files.is_science(hdr, tel) == True
    hdr['SCRN'] = 'deployed'
    assert sort_files.is_science(hdr, tel) == False


def test_is_bias_binospec():
    """is_bias False for BINOSPEC (no bias keywords)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    assert sort_files.is_bias(hdr, tel) is False


def test_is_dark_binospec():
    """is_dark False for BINOSPEC (no dark keywords)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    assert sort_files.is_dark(hdr, tel) is False


def test_is_bias_gmos():
    """is_bias True when SHUTTER=closed, OBSCLASS=daycal, OBSTYPE=bias (GMOS)."""
    tel = instrument_getter('GMOS')
    hdr = fits.Header()
    hdr['SHUTTER'] = 'closed'
    hdr['OBSCLASS'] = 'daycal'
    hdr['OBSTYPE'] = 'bias'
    hdr['CCDSUM'] = '2,2'
    assert sort_files.is_bias(hdr, tel) == True
    hdr['OBSTYPE'] = 'object'
    assert sort_files.is_bias(hdr, tel) == False
