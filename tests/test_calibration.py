from potpyri.utils import options
from potpyri.utils import logger
from potpyri.stages import image_procs

import os
import numpy as np

from astropy.io import fits

from tests.utils import download_gdrive_file

def test_cal(tmp_path):

    instrument = 'GMOS'

    # Raw science file plus calibration files
    file_path = download_gdrive_file('GMOS/raw/N20240618S0015.fits.bz2', use_cached=True)
    download_gdrive_file('GMOS/red/cals/mbias_12_22.fits.fz', use_cached=True)
    download_gdrive_file('GMOS/red/cals/mflat_i_12_22.fits.fz', use_cached=True)

    data_path, basefile = os.path.split(file_path)
    paths, tel = options.initialize_telescope(instrument, data_path[0:-4])

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    mbias = tel.load_bias(paths, '12', '22')
    mflat = tel.load_flat(paths, 'i', '12', '22')

    sci_full = tel.import_image(file_path, '12', log=log)

    data = (sci_full-mbias.data)/(mflat.data/np.nanmean(mflat.data))
    data = data.astype(np.float32)

    processed = tel.process_science([file_path], 'i', '12', '22', paths,
        mbias=mbias, mflat=mflat, mdark=None, skip_skysub=True, log=log)
    
    np.testing.assert_array_equal(processed[0].data[~np.isnan(processed[0].data)], data[~np.isnan(processed[0].data)])
