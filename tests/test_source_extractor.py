from potpyri.stages import photometry
from tests.utils import download_gdrive_file
import numpy as np

def test_source_extractor(tmp_path):

    # Stacked science file
    file_path = download_gdrive_file('LRISr/red/R155_host.R.ut240629.1R.11.stk.fits.fz', use_cached=True)
    sky_cat = photometry.run_sextractor(file_path, log=None)

    assert len(sky_cat)>1000
    assert np.all([key in sky_cat.keys() for key in ['X_IMAGE','Y_IMAGE',
        'MAG_AUTO','MAGERR_AUTO','FWHM_IMAGE']])
