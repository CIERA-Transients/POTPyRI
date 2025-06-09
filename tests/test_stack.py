from potpyri.utils import options
from potpyri.utils import logger
from potpyri.stages import image_procs

import numpy as np
import os

from ccdproc import CCDData

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

from . import utils

def test_stack(tmp_path):

    instrument = 'MMIRS'
    nims = 5

    # Science file (just tellurics)
    paths, tel = options.initialize_telescope(instrument, tmp_path)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    aligned_data = []
    masks = []
    errors = []
    for i in range(nims):
        imdata = np.ones(shape=(2048, 2048), dtype=np.float32) * (i+1)
        errimg = np.sqrt(imdata)
        maskim = np.zeros(shape=(2048, 2048), dtype=np.uint16)
        header = {'EXPTIME': 1.0}

        wcshead = utils.make_dummy_wcs()
        header.update(wcshead)

        wcs = WCS(header)

        image = CCDData(imdata, meta=header, wcs=wcs, unit=u.electron)
        aligned_data.append(image)
        masks.append(fits.ImageHDU(maskim))
        errors.append(fits.ImageHDU(errimg))

    sci_med = image_procs.stack_data(aligned_data, tel, masks, errors, log=log)

    np.testing.assert_array_equal(sci_med[0].data, np.median(np.arange(nims)+1))
