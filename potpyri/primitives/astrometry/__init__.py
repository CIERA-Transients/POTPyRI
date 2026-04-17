"""WCS solution using astrometry.net and Gaia DR3 alignment.

Provides coarse solution via solve-field and fine alignment by matching
detected sources to Gaia. Writes RADISP/DEDISP and updated WCS to FITS.
Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from __future__ import annotations

from .astrometry_net_primitive import AstrometryNetPrimitive
from .wcs_solution import (
    _log_gaia,
    _validate_refined_wcs,
    _write_gaia_fallback_header,
    align_to_gaia,
    clean_up_astrometry,
    get_gaia_catalog,
    solve_astrometry,
)

__all__ = [
    'AstrometryNetPrimitive',
    '_log_gaia',
    '_validate_refined_wcs',
    '_write_gaia_fallback_header',
    'align_to_gaia',
    'clean_up_astrometry',
    'get_gaia_catalog',
    'solve_astrometry',
]
