from potpyri.utils import options
from potpyri.instruments import instrument_getter

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
        tel = instrument_getter(instrument)
        assert tel.name.upper()==instrument.upper()
