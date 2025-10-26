from potpyri.utils import options
from potpyri.instruments import instrument_getter

def test_staticmask(tmp_path):

    file_list_name = 'files.txt'

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

        tel = instrument_getter(instname)
        paths = options.add_paths(tmp_path, file_list_name, tel)

        assert tel.name.upper()==instname.upper()

        staticmask = tel.load_staticmask(hdr, paths)

        assert staticmask is not None
