from potpyri.scripts import main_pipeline
from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import sort_files

import numpy as np
import astropy.utils.data
import importlib
import os

def test_staticmask(tmp_path):

    instruments = {
        'F2': {'INST': 'F2', 'INSTRUME': 'F2'},
        'GMOS-N': {'INST': 'GMOS', 'INSTRUME': 'GMOS-N', 'CCDSUM': '2,2'},
        'GMOS-S': {'INST': 'GMOS', 'INSTRUME': 'GMOS-S', 'CCDSUM': '2,2'},
        'MOSFIRE': {'INST': 'MOSFIRE', 'INSTRUME': 'MOSFIRE'},
        'BINOSPEC': {'INST': 'BINOSPEC', 'INSTRUME': 'BINOSPEC', 'CCDSUM': '1,1'},
    }

    # Test initialization of each POTPyRI instrument
    for instrument in instruments.keys():
        hdr = instruments[instrument]
        instname = hdr['INST']
        paths, tel = options.initialize_telescope(instname, tmp_path)
        assert tel.name.upper()==instname.upper()

        staticmask = tel.load_staticmask(hdr, paths)

        assert staticmask is not None
