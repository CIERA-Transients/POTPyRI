"""Regression: main_pipeline must not call absphot when photometry raises."""
from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from potpyri.primitives import photometry
from potpyri.scripts import main_pipeline as mp


class _FakeTel:
    version = 'test'
    flat = True
    filetype_keywords = {'SCIENCE': 'SCIENCE'}

    def match_type_keywords(self, kwds, file_table):
        return file_table['Type'] == 'SCIENCE'


def _minimal_file_table():
    return Table(
        {
            'TargType': ['science'],
            'Type': ['SCIENCE'],
            'Filename': ['dummy.fits'],
        },
    )


class _FakeLog:
    def info(self, *args, **kwargs):
        pass

    def error(self, *args, **kwargs):
        pass

    def exception(self, *args, **kwargs):
        pass

    def close(self):
        pass

    def shutdown(self):
        pass


def test_main_pipeline_propagates_photometry_error_without_calling_absphot(
    monkeypatch, tmp_path,
):
    stack_path = tmp_path / 'stack.fits'
    data = np.ones((8, 8), dtype=np.float32)
    hdr = fits.Header({'FILTER': 'r', 'EXTNAME': 'SCI'})
    fits.HDUList(
        [
            fits.PrimaryHDU(data=data, header=hdr),
            fits.ImageHDU(data=np.zeros((8, 8), np.uint8), name='MASK'),
            fits.ImageHDU(data=np.ones((8, 8), np.float32), name='ERROR'),
        ],
    ).writeto(stack_path, overwrite=True)

    absphot_called = []

    monkeypatch.setattr(mp, 'instrument_getter', lambda name: _FakeTel())
    monkeypatch.setattr(
        mp.sort_files,
        'handle_files',
        lambda *a, **k: _minimal_file_table(),
    )
    monkeypatch.setattr(mp.calibration, 'do_bias', lambda *a, **k: None)
    monkeypatch.setattr(mp.calibration, 'do_dark', lambda *a, **k: None)
    monkeypatch.setattr(mp.calibration, 'do_flat', lambda *a, **k: None)
    monkeypatch.setattr(mp.image_procs, 'image_proc', lambda *a, **k: str(stack_path))
    monkeypatch.setattr(
        mp.photometry,
        'photloop',
        lambda *a, **k: (_ for _ in ()).throw(
            photometry.PhotometryError('photometry failed'),
        ),
    )
    monkeypatch.setattr(
        mp.absphot,
        'find_zeropoint',
        lambda *a, **k: absphot_called.append(True),
    )
    monkeypatch.setattr(
        mp.options,
        'add_paths',
        lambda *a, **k: {
            'log': str(tmp_path / 'log'),
            'work': str(tmp_path),
            'filelist': str(tmp_path / 'files.txt'),
        },
    )
    monkeypatch.setattr(mp.logger, 'get_log', lambda *a, **k: _FakeLog())
    (tmp_path / 'log').mkdir(exist_ok=True)

    with pytest.raises(photometry.PhotometryError, match='photometry failed'):
        mp.main_pipeline(
            instrument='LRIS',
            data_path=str(tmp_path),
            target=None,
            file_list_name='files.txt',
        )

    assert absphot_called == []
