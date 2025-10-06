from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import solve_wcs

import os
import numpy as np

from astropy.io import fits

from tests.utils import download_gdrive_file

def test_wcs(tmp_path):

    instrument = 'GMOS'

    # Download science file
    file_path = download_gdrive_file('GMOS/red/workspace/sGRB240615A-GRB.i.ut240618.12.22.1403185114.fits.fz', use_cached=True)
    astm_path = download_gdrive_file('astrometry/index-52m1-14.fits')

    data_path, basefile = os.path.split(file_path)
    paths, tel = options.initialize_telescope(instrument, data_path[0:-14])

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    # Write out basic file to solve WCS on
    hdu = fits.open(file_path)
    binn = tel.get_binning(hdu[1].header)
    file_path = file_path.replace('.fz','')
    del hdu[1].header['RADISP']
    del hdu[1].header['DEDISP']
    fits.writeto(file_path,hdu[1].data,hdu[1].header,overwrite=True)

    # Solve WCS
    solve_wcs.solve_astrometry(file_path, tel, binn, paths, index=astm_path, log=log)
    solve_wcs.align_to_gaia(file_path, tel, radius=0.5, log=log)

    # Test dispersion
    hdu = fits.open(file_path)
    ra_disp = hdu[0].header['RADISP']
    dec_disp = hdu[0].header['DEDISP']

    assert ra_disp < 0.5
    assert dec_disp < 0.5

    assert ra_disp/tel.pixscale<1
    assert dec_disp/tel.pixscale<1

if __name__=="__main__":
    test_wcs('.')
