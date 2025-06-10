from potpyri.utils import options
from potpyri.utils import logger
from potpyri.stages import absphot

import os
import numpy as np

from astropy.io import fits

from tests.utils import download_gdrive_file

def test_absphot(tmp_path):

    instrument = 'GMOS'

    # Raw science file
    file_path = download_gdrive_file('GMOS/red/sGRB240615A-GRB.i.ut240618.12.22.stk.fits.fz', use_cached=True)

    # Strip out extensions added by photometry loop if they exist
    hdu = fits.open(file_path)
    for key in ['APPHOT','PSFPHOT','PSFSTARS','RESIDUAL','PSF']:
        if key in [h.name for h in hdu]:
            del hdu[key]

    hdu.writeto(file_path, overwrite=True)

    data_path, basefile = os.path.split(file_path)
    paths, tel = options.initialize_telescope(instrument, data_path)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    hdu = fits.open(file_path)
    cat = tel.get_catalog(hdu['SCI'].header)
    filt = hdu['SCI'].header['FILTER']
    cal = absphot.absphot()
    zpt, zpterr = cal.find_zeropoint(file_path, hdu['SCI'].header['FILTER'], cat, log=log)

    hdu = fits.open(file_path)
    m3sigma = hdu['PRIMARY'].header['M3SIGMA']

    assert cat=='PS1'
    assert filt=='i'
    assert np.abs(zpt-33.649132)<0.01
    assert zpterr<0.01
    assert np.abs(m3sigma-25.808789)<0.01
