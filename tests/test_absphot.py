from potpyri.utils import options
from potpyri.utils import logger
from potpyri.stages import absphot

from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord

import os
import numpy as np

from astropy.io import fits

from tests.utils import download_gdrive_file

def test_absphot(tmp_path):

    instrument = 'LRIS'

    # Raw science file
    file_path = download_gdrive_file('LRISr/red/R155_host.R.ut240629.1R.11.stk.fits.fz', use_cached=True)
    cat_path = download_gdrive_file('LRISr/red/test_catalog.dat', use_cached=True)
    input_catalog = ascii.read(cat_path)

    data_path, basefile = os.path.split(file_path)
    paths, tel = options.initialize_telescope(instrument, data_path)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    hdu = fits.open(file_path)

    cal = absphot.absphot()
    cat = tel.get_catalog(hdu['SCI'].header)
    filt = hdu['SCI'].header['FILTER']

    zpt, zpterr = cal.find_zeropoint(file_path, hdu['SCI'].header['FILTER'], 
        cat, input_catalog=input_catalog, log=log)

    assert cat=='PS1'
    assert filt=='R'
    assert np.abs(zpt-27.576871)<0.01
    assert zpterr<0.01
