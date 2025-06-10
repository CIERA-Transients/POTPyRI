from potpyri.scripts import main_pipeline
from potpyri.utils import options
from potpyri.utils import logger
from potpyri.stages import sort_files

import numpy as np
import astropy.utils.data
import importlib
import os

def test_init(tmp_path):

    # Parse allowed instruments directly from options so this test matches
    # to values that the argument parser allows as input
    instruments = []
    params = options.init_options()
    for pos in params._get_positional_actions():
        if pos.dest=='instrument':
            instruments = pos.choices

    # Test initialization of each POTPyRI instrument
    for instrument in instruments:
        paths, tel = options.initialize_telescope(instrument, tmp_path)
        assert tel.name.upper()==instrument.upper()
