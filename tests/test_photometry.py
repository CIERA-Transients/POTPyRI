from potpyri.utils import options
from potpyri.utils import logger
from potpyri.stages import photometry

import os
import numpy as np

from astropy.io import fits

from tests.utils import download_gdrive_file

def test_photometry(tmp_path):

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

    phot_sn_min = 20.0
    phot_sn_max = 100.0
    fwhm_init = 5.0

    photometry.photloop(file_path, phot_sn_min=phot_sn_min,
        fwhm_init=fwhm_init, phot_sn_max=phot_sn_max, log=log)

    hdu = fits.open(file_path)

    assert np.any(['PSFPHOT' in h.name for h in hdu])
    assert np.any(['PSFSTARS' in h.name for h in hdu])
    assert np.any(['APPPHOT' in h.name for h in hdu])
    assert np.any(['PSF' in h.name for h in hdu])

    assert len(hdu['APPPHOT'].data)>100
