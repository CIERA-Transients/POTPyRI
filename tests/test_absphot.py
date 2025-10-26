from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import absphot
from potpyri.instruments import instrument_getter

from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord

import os
import numpy as np

from astropy.io import fits

from tests.utils import download_gdrive_file

def test_absphot(tmp_path):

    instrument = 'LRIS'
    file_list_name = 'files.txt'

    # Raw science file
    file_path = download_gdrive_file('LRISr/red/R155_host.R.ut240629.1R.11.stk.fits.fz', use_cached=True)
    cat_path = download_gdrive_file('LRISr/red/test_catalog.dat', use_cached=True)
    input_catalog = ascii.read(cat_path)

    data_path, basefile = os.path.split(file_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    hdu = fits.open(file_path)

    absphot.find_zeropoint(file_path, tel, log=log)

    hdu = fits.open(file_path)
    header = hdu['PRIMARY'].header

    assert header['ZPTCAT']=='PS1'
    assert header['FILTER']=='R'
    assert np.abs(27.576871-header['ZPTMAG'])<0.01
    assert header['ZPTMUCER']<0.01
