"""Tests for sort_files.handle_files and file table content (BINOSPEC), and classifiers (is_bad, is_spec, is_flat, is_dark, is_bias, is_science)."""
from potpyri.utils import options
from potpyri.utils import logger
from potpyri.primitives import sort_files
from potpyri.instruments import instrument_getter

import os

from astropy.io import fits

import pytest
from tests.utils import download_gdrive_file


@pytest.mark.integration
def test_sort(tmp_path):
    """Run handle_files on BINOSPEC proc file; check file table rows and no_redo reuse."""
    instrument = 'BINOSPEC'
    file_list_name = 'files.txt'

    # Raw science file (download to tmp_path so log is writable)
    file_path = download_gdrive_file('Binospec/raw/sci_img_2024.0812.034220_proc.fits.fz', output_dir=str(tmp_path), use_cached=True)

    data_path, basefile = os.path.split(file_path)
    data_path, _ = os.path.split(data_path)
    tel = instrument_getter(instrument)
    paths = options.add_paths(data_path, file_list_name, tel)

    # Generate log file in corresponding directory for log
    log = logger.get_log(paths['log'])

    # This contains all of the file data
    try:
        file_table = sort_files.handle_files(paths['filelist'], paths, tel,
            incl_bad=True, proc='proc', no_redo=False, log=log)

        # Run validation checks on file_table
        assert len(file_table)==1
        assert file_table[0]['Target']=='GRB240809A_r_ep1'
        assert file_table[0]['Filter']=='r'
        assert file_table[0]['Type']=='SCIENCE'
        assert file_table[0]['TargType']=='GRB240809A_r_ep1_r_2_11'
        assert file_table[0]['CalType']=='r_2_11'
        assert file_table[0]['Exp']=="120.0"
        assert file_table[0]['Binning']=="11"
        assert file_table[0]['Amp']=="2"

        # Test the no_redo flag
        new_file_table = sort_files.handle_files(paths['filelist'], paths, tel,
            incl_bad=True, proc='proc', no_redo=True, log=log)
    finally:
        log.close()

    assert len(new_file_table)==1
    assert file_table[0]['Target']==new_file_table[0]['Target']
    assert file_table[0]['Filter']==new_file_table[0]['Filter']
    assert file_table[0]['Type']==new_file_table[0]['Type']
    assert file_table[0]['TargType']==new_file_table[0]['TargType']
    assert file_table[0]['CalType']==new_file_table[0]['CalType']
    assert file_table[0]['Exp']==new_file_table[0]['Exp']
    assert file_table[0]['Binning']==new_file_table[0]['Binning']
    assert file_table[0]['Amp']==new_file_table[0]['Amp']


def test_is_bad_binospec():
    """is_bad True when header matches bad_keywords/bad_values (BINOSPEC: MASK=mira)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'mira'
    # No CCDSUM so binning check does not overwrite; keyword match gives bad=True
    assert sort_files.is_bad(hdr, tel) == True
    hdr['MASK'] = 'imaging'
    hdr['CCDSUM'] = '1,1'  # valid binning so keyword result (False) is kept
    assert sort_files.is_bad(hdr, tel) == False


def test_is_spec_binospec():
    """is_spec True when header matches spec keywords (BINOSPEC: MASK=spectroscopy)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'spectroscopy'
    assert sort_files.is_spec(hdr, tel) == True
    hdr['MASK'] = 'imaging'
    assert sort_files.is_spec(hdr, tel) == False


def test_is_flat_binospec():
    """is_flat True when MASK=imaging, SCRN=deployed (BINOSPEC)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'imaging'
    hdr['SCRN'] = 'deployed'
    hdr['CCDSUM'] = '1,1'
    assert sort_files.is_flat(hdr, tel) == True
    hdr['SCRN'] = 'stowed'
    assert sort_files.is_flat(hdr, tel) == False


def test_is_science_binospec():
    """is_science True when MASK=imaging, SCRN=stowed and exptime >= min (BINOSPEC)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    hdr['MASK'] = 'imaging'
    hdr['SCRN'] = 'stowed'
    hdr['CCDSUM'] = '1,1'
    hdr['EXPTIME'] = 120.0
    assert sort_files.is_science(hdr, tel) == True
    hdr['SCRN'] = 'deployed'
    assert sort_files.is_science(hdr, tel) == False


def test_is_bias_binospec():
    """is_bias False for BINOSPEC (no bias keywords)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    assert sort_files.is_bias(hdr, tel) is False


def test_is_dark_binospec():
    """is_dark False for BINOSPEC (no dark keywords)."""
    tel = instrument_getter('BINOSPEC')
    hdr = fits.Header()
    assert sort_files.is_dark(hdr, tel) is False


def test_is_bias_gmos():
    """is_bias True when SHUTTER=closed, OBSCLASS=daycal, OBSTYPE=bias (GMOS)."""
    tel = instrument_getter('GMOS')
    hdr = fits.Header()
    hdr['SHUTTER'] = 'closed'
    hdr['OBSCLASS'] = 'daycal'
    hdr['OBSTYPE'] = 'bias'
    hdr['CCDSUM'] = '2,2'
    assert sort_files.is_bias(hdr, tel) == True
    hdr['OBSTYPE'] = 'object'
    assert sort_files.is_bias(hdr, tel) == False


def _lris_imaging_hdr(**kwargs):
    """Minimal raw LRIS imaging header (no KOAIMTYP)."""
    hdr = fits.Header()
    hdr['SLITNAME'] = 'direct'
    hdr['GRANAME'] = 'mirror'
    hdr['TRAPDOOR'] = 'open'
    hdr['GRISTRAN'] = 'stowed'
    hdr['BINNING'] = '2,2'
    hdr['ELAPTIME'] = 240
    hdr['OBJECT'] = 'FRB20260326A'
    hdr.update(kwargs)
    return hdr


def test_is_science_lris_raw_without_koaimtyp():
    """Raw LRIS science matches slit/grating/door; missing KOAIMTYP does not raise."""
    tel = instrument_getter('LRIS')
    hdr = _lris_imaging_hdr()
    assert 'KOAIMTYP' not in hdr
    assert sort_files.is_science(hdr, tel) == True
    assert sort_files.is_flat(hdr, tel) == False
    assert sort_files.is_bias(hdr, tel) == False

    hdr['ELAPTIME'] = 5
    assert sort_files.is_science(hdr, tel) == False


def test_is_flat_lris_object_name():
    """LRIS flats match OBJECT containing 'flat' (dome flat direct)."""
    tel = instrument_getter('LRIS')
    hdr = _lris_imaging_hdr(OBJECT='dome flat direct', ELAPTIME=60)
    assert sort_files.is_flat(hdr, tel) == True
    assert sort_files.is_bias(hdr, tel) == False


def test_is_flat_lris_object_twilight_name():
    """LRIS flats also match OBJECT containing twilight or skyflat."""
    tel = instrument_getter('LRIS')
    hdr = _lris_imaging_hdr(OBJECT='twilight', ELAPTIME=20)
    assert sort_files.is_flat(hdr, tel) == True
    hdr['OBJECT'] = 'skyflat R'
    assert sort_files.is_flat(hdr, tel) == True


def test_is_flat_lris_twilight_sun_altitude():
    """Unlabeled imaging frames in twilight with 15k-50k sky are flats."""
    tel = instrument_getter('LRIS')
    from astropy.time import Time
    import numpy as np

    twilight = Time('2026-07-16T06:00:00')
    night = Time('2026-07-16T12:00:00')
    twilight_alt = tel.get_sun_altitude(
        _lris_imaging_hdr(MJD=twilight.mjd, OBJECT='FRB20251015A', ELAPTIME=30)
    )
    night_alt = tel.get_sun_altitude(
        _lris_imaging_hdr(MJD=night.mjd, OBJECT='FRB20251015A', ELAPTIME=180)
    )
    assert twilight_alt is not None
    assert tel.twilight_sunalt_min <= twilight_alt <= tel.twilight_sunalt_max
    assert night_alt < tel.twilight_sunalt_min

    sky = np.full((32, 32), 20000.0)
    hdr = _lris_imaging_hdr(
        MJD=twilight.mjd, OBJECT='FRB20251015A', ELAPTIME=30, FLAMP1='off', FLAMP2='off',
    )
    assert tel.is_twilight_setup(hdr) == True
    assert tel.is_twilight_flat(hdr) == False
    assert tel.is_twilight_flat(hdr, data=sky) == True
    assert sort_files.is_flat(hdr, tel, data=sky) == True
    assert sort_files.is_flat(hdr, tel, data=np.full((32, 32), 1000.0)) == False
    assert sort_files.is_flat(hdr, tel, data=np.full((32, 32), 60000.0)) == False

    hdr['MJD'] = night.mjd
    hdr['ELAPTIME'] = 180
    assert tel.is_twilight_setup(hdr) == False
    assert tel.is_twilight_flat(hdr, data=sky) == False
    assert sort_files.is_flat(hdr, tel, data=sky) == False


def test_lris_blue_twilight_uses_vidinp2_and_vidinp3():
    """LRISblue sky level ignores unilluminated VidInp1 and VidInp4."""
    import numpy as np
    from astropy.time import Time

    tel = instrument_getter('LRIS')
    primary = fits.PrimaryHDU(header=_lris_imaging_hdr(
        INSTRUME='LRISBLUE', MJD=Time('2026-07-16T06:00:00').mjd,
        OBJECT='ZTF26abnuyoi', ELAPTIME=13, FLAMP1='off', FLAMP2='off',
        BLUFILT='G',
    ))
    hdus = [primary]
    for name, level in (
        ('VidInp1', 1300.0),
        ('VidInp2', 22000.0),
        ('VidInp3', 23000.0),
        ('VidInp4', 1500.0),
    ):
        hdus.append(fits.ImageHDU(np.full((16, 16), level, dtype=np.float32), name=name))
    hdul = fits.HDUList(hdus)

    data = tel.get_twilight_image_data(hdul, 0)
    assert data is not None
    assert data.size == 16 * 16 * 2
    level = tel.twilight_sky_level(data)
    assert 20000.0 <= level <= 24000.0
    assert tel.is_twilight_flat(hdul[0].header, data=data) == True
    assert sort_files.is_flat(hdul[0].header, tel, data=data) == True


def test_is_flat_lris_focus_loop_flatlamp():
    """KOA dome flats keep OBJECT='Focus loop' and KOAIMTYP=flatlamp."""
    tel = instrument_getter('LRIS')
    hdr = _lris_imaging_hdr(
        OBJECT='Focus loop', TARGNAME='unknown', KOAIMTYP='flatlamp',
        ELAPTIME=2, FLAMP1='off', FLAMP2='off',
    )
    assert sort_files.is_flat(hdr, tel) == True
    assert sort_files.is_science(hdr, tel) == False


def test_is_flat_lris_lamp_on():
    """LRIS lamp-on frames are flats even when OBJECT is not named flat."""
    tel = instrument_getter('LRIS')
    hdr = _lris_imaging_hdr(OBJECT='HORIZON STOW', ELAPTIME=15, FLAMP1='on', FLAMP2='on')
    assert sort_files.is_flat(hdr, tel) == True
    assert sort_files.is_science(hdr, tel) == True

    hdr['FLAMP1'] = 'off'
    hdr['FLAMP2'] = 'off'
    assert sort_files.is_flat(hdr, tel) == False


def test_is_bias_lris_object_and_horizon_stow():
    """LRIS bias from OBJECT=bias or closed-door zero-second HORIZON STOW."""
    tel = instrument_getter('LRIS')
    hdr = _lris_imaging_hdr(
        OBJECT='bias', TARGNAME='HORIZON STOW', TRAPDOOR='closed',
        SLITNAME='long_1.0', GRANAME='600/7500', GRISTRAN='deployed',
        ELAPTIME=0,
    )
    assert sort_files.is_bias(hdr, tel) == True

    hdr['OBJECT'] = 'HORIZON STOW'
    assert sort_files.is_bias(hdr, tel) == True

    hdr['ELAPTIME'] = 120
    assert sort_files.is_bias(hdr, tel) == False


def test_lris_missing_keywords_do_not_raise():
    """Missing sort keywords are a non-match, not a KeyError."""
    tel = instrument_getter('LRIS')
    hdr = fits.Header()
    hdr['BINNING'] = '2,2'
    hdr['ELAPTIME'] = 240
    assert sort_files.is_flat(hdr, tel) == False
    assert sort_files.is_bias(hdr, tel) == False
    assert sort_files.is_spec(hdr, tel) == False
    assert sort_files.is_science(hdr, tel) == False


def test_handle_files_logs_discovery_when_no_files(tmp_path, capsys):
    """handle_files prints search paths and glob pattern before exiting on empty input."""
    tel = instrument_getter('F2')
    raw_dir = tmp_path / 'raw'
    data_dir = tmp_path
    bad_dir = tmp_path / 'bad'
    for d in (raw_dir, bad_dir):
        d.mkdir()
    paths = {
        'raw': str(raw_dir),
        'data': str(data_dir),
        'bad': str(bad_dir),
        'filelist': str(tmp_path / 'file_list.txt'),
    }
    with pytest.raises(SystemExit):
        sort_files.handle_files(
            paths['filelist'], paths, tel, proc='fits', log=None,
        )
    out = capsys.readouterr().out
    assert 'Glob pattern' in out
    assert '*.fits' in out
    assert str(raw_dir) in out
    assert 'No files matched' in out 
