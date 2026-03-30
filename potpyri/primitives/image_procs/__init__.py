"""Image calibration, masking, error arrays, alignment, and stacking.

Provides WCS handling, alignment to a common grid, cosmic-ray rejection,
satellite masking, master calibration combination, and stacked science
output. Authors: Charlie Kilpatrick.
"""
from __future__ import annotations

from potpyri._version import __version__

from .align import align_images
from .mask_error import create_error, create_mask
from .primitive import ImageProcPrimitive, image_proc
from .satellites import mask_satellites
from .stack_ops import (
    add_stack_mask,
    compute_relative_scales,
    detrend_stack,
    stack_data,
)
from .wcs_geom import generate_wcs, get_fieldcenter, remove_pv_distortion
from .worker import _image_proc_worker

__all__ = [
    'ImageProcPrimitive',
    '__version__',
    '_image_proc_worker',
    'add_stack_mask',
    'align_images',
    'compute_relative_scales',
    'create_error',
    'create_mask',
    'detrend_stack',
    'generate_wcs',
    'get_fieldcenter',
    'image_proc',
    'mask_satellites',
    'remove_pv_distortion',
    'stack_data',
]
