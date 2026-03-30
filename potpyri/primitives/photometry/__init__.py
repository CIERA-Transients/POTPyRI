"""Aperture and PSF photometry for pipeline stacks.

Uses Source Extractor for initial catalogs, photutils for ePSF building
and PSF fitting. Writes APPPHOT, PSFPHOT, PSFSTARS, PSF, and RESIDUAL
extensions to the stack FITS. Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from __future__ import annotations

from .core import (
    create_conv,
    create_params,
    do_phot,
    extract_aperture_stats,
    extract_fwhm_from_epsf,
    generate_epsf,
    get_star_catalog,
    run_photometry,
    run_sextractor,
)
from .loop import PhotometryPrimitive, _photloop_worker, photloop

__all__ = [
    'PhotometryPrimitive',
    '_photloop_worker',
    'create_conv',
    'create_params',
    'do_phot',
    'extract_aperture_stats',
    'extract_fwhm_from_epsf',
    'generate_epsf',
    'get_star_catalog',
    'photloop',
    'run_photometry',
    'run_sextractor',
]
