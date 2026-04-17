"""Catalog-based absolute zeropoint for stacked science images."""
from __future__ import annotations

from potpyri._version import __version__

from .catalog_zeropoint import (
    DEFAULT_MAGSYS,
    GF18_TWOMASS_TO_VISTA_Y_SLOPE,
    GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR,
    GF18_VISTA_Y_VEGA_TO_AB,
    VIZIER_MIRRORS,
    ZeropointFitter,
)
from .zeropoint_primitive import ZeropointCalibrationPrimitive, find_zeropoint

__all__ = [
    'DEFAULT_MAGSYS',
    'GF18_TWOMASS_TO_VISTA_Y_SLOPE',
    'GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR',
    'GF18_VISTA_Y_VEGA_TO_AB',
    'VIZIER_MIRRORS',
    'ZeropointCalibrationPrimitive',
    'ZeropointFitter',
    '__version__',
    'find_zeropoint',
]
