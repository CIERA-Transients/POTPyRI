from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import image_procs
from potpyri.instruments import instrument_getter

import numpy as np
import os

from ccdproc import CCDData

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

from . import utils

def test_error(tmp_path):

    instrument = 'MOSFIRE'
    file_list_name = 'files.txt'

    # Science file (just tellurics)
    tel = instrument_getter(instrument)
    paths = options.add_paths(tmp_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    aligned_data = []
    masks = []
    errors = []

    rdnoise = 10.
    imgval = 1000.

    imdata = np.ones(shape=(2048, 2048), dtype=np.float32) * imgval
    maskim = np.zeros(shape=(2048, 2048), dtype=np.uint16)
    maskhdu = fits.PrimaryHDU(maskim)

    outfile = os.path.join(tmp_path, 'test.fits')
    hdu = fits.PrimaryHDU(imdata)
    hdu.header['SATURATE'] = 35000.0
    hdu.writeto(outfile, overwrite=True, output_verify='silentfix')

    error_hdu = image_procs.create_error(outfile, maskhdu, rdnoise)

    rms = 0.5 * (
            np.percentile(imdata[~maskim.astype(bool)], 84.13)
            - np.percentile(imdata[~maskim.astype(bool)], 15.86)
    )

    theoretical_error = np.sqrt(imdata + rms**2 + rdnoise)

    np.testing.assert_array_equal(error_hdu.data, theoretical_error)
