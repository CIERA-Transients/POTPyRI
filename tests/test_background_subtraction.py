"""Tests for per-frame optical background subtraction modes."""
import numpy as np
import pytest
import astropy.units as u
from astropy.nddata import CCDData

from potpyri.instruments import instrument_getter
from potpyri.instruments.instrument import (
    BKG_SUB_CONSTANT,
    BKG_SUB_LOCAL,
    BKG_SUB_NONE,
)


def _flat_frame(level=1000.0, shape=(128, 128)):
    data = np.full(shape, level, dtype=np.float64)
    mask = np.zeros(shape, dtype=bool)
    return CCDData(data, unit=u.electron, mask=mask)


def _gradient_frame(shape=(128, 128)):
    y, x = np.indices(shape)
    data = 500.0 + 2.0 * x + 3.0 * y
    mask = np.zeros(shape, dtype=bool)
    return CCDData(data, unit=u.electron, mask=mask)


def test_constant_background_subtracts_uniform_level():
    tel = instrument_getter('GMOS')
    frame = _flat_frame(level=1000.0)
    result = tel.apply_optical_background_subtraction(
        frame, bkg_sub=BKG_SUB_CONSTANT)
    assert result.header['BKGSUB'] == BKG_SUB_CONSTANT
    assert result.header['SKYBKG'] == pytest.approx(1000.0)
    finite = result.data[np.isfinite(result.data)]
    assert np.nanmedian(finite) == pytest.approx(0.0, abs=1e-5)


def test_constant_background_preserves_spatial_gradient():
    tel = instrument_getter('GMOS')
    frame = _gradient_frame()
    result = tel.apply_optical_background_subtraction(
        frame, bkg_sub=BKG_SUB_CONSTANT)
    # Constant subtraction removes only the median level, not the gradient.
    assert np.std(result.data) > np.std(frame.data) * 0.5


def test_none_background_leaves_data_unchanged():
    tel = instrument_getter('GMOS')
    frame = _flat_frame(level=1000.0)
    result = tel.apply_optical_background_subtraction(frame, bkg_sub=BKG_SUB_NONE)
    assert result.header['BKGSUB'] == BKG_SUB_NONE
    assert result.header['SKYBKG'] == 0.0
    np.testing.assert_array_equal(result.data, frame.data)


def test_local_background_reduces_spatial_gradient_more_than_constant():
    tel = instrument_getter('GMOS')
    frame = _gradient_frame(shape=(256, 256))
    constant = tel.apply_optical_background_subtraction(
        frame, bkg_sub=BKG_SUB_CONSTANT)
    local = tel.apply_optical_background_subtraction(frame, bkg_sub=BKG_SUB_LOCAL)
    assert local.header['BKGSUB'] == BKG_SUB_LOCAL
    assert np.nanstd(local.data) < np.nanstd(constant.data)


def test_invalid_bkg_sub_raises():
    tel = instrument_getter('GMOS')
    frame = _flat_frame()
    with pytest.raises(ValueError, match='bkg_sub must be one of'):
        tel.apply_optical_background_subtraction(frame, bkg_sub='polynomial')
