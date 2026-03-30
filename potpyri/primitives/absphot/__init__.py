"""Absolute photometry zeropoint calibration using catalog magnitudes."""
from __future__ import annotations

from potpyri._version import __version__

from .primitive import AbsPhotZeropointPrimitive, find_zeropoint
from .zeropoint import (
    DEFAULT_MAGSYS,
    GF18_TWOMASS_TO_VISTA_Y_SLOPE,
    GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR,
    GF18_VISTA_Y_VEGA_TO_AB,
    VIZIER_MIRRORS,
    absphot,
)

__all__ = [
    'AbsPhotZeropointPrimitive',
    'DEFAULT_MAGSYS',
    'GF18_TWOMASS_TO_VISTA_Y_SLOPE',
    'GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR',
    'GF18_VISTA_Y_VEGA_TO_AB',
    'VIZIER_MIRRORS',
    '__version__',
    'absphot',
    'find_zeropoint',
]
