from potpyri.stages import photometry

import os

from tests.utils import download_gdrive_file

def test_source_extractor(tmp_path):

    # Stacked science file
    file_path = download_gdrive_file('LRISr/red/R155_host.R.ut240629.1R.11.stk.fits.fz', use_cached=True)
    sky_cat = photometry.run_sextractor(file_path, log=None)

    breakpoint()
