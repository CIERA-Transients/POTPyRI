"""Aperture and PSF photometry on pipeline stacks.

Uses Source Extractor for initial catalogs, photutils for ePSF building
and PSF fitting. Writes APPPHOT, PSFPHOT, PSFSTARS, PSF, and RESIDUAL
extensions to the stack FITS. Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from __future__ import annotations

from .catalog_psf import (
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
from .photometry_loop import (
    StackPhotometryPrimitive,
    _stack_photometry_worker,
    run_stack_photometry,
)

__all__ = [
    'StackPhotometryPrimitive',
    '_stack_photometry_worker',
    'create_conv',
    'create_params',
    'do_phot',
    'extract_aperture_stats',
    'extract_fwhm_from_epsf',
    'generate_epsf',
    'get_star_catalog',
    'run_photometry',
    'run_sextractor',
    'run_stack_photometry',
]
