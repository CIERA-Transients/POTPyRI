"""Science image calibration, masking, alignment, and stacking.

Provides WCS handling, alignment to a common grid, cosmic-ray rejection,
satellite masking, and stacked science output. Authors: Charlie Kilpatrick.
"""
from __future__ import annotations

from potpyri._version import __version__

from .align_reproject import align_images
from .combine_and_detrend import (
    add_stack_mask,
    compute_relative_scales,
    detrend_stack,
    stack_data,
)
from .field_geometry import generate_wcs, get_fieldcenter, remove_pv_distortion
from .mask_and_error import create_error, create_mask
from .satellite_trails import mask_satellites
from .stacking_pipeline_worker import _stacking_pipeline_worker
from .stacking_primitive import ScienceStackingPrimitive, stack_science_frames

__all__ = [
    'ScienceStackingPrimitive',
    '__version__',
    '_stacking_pipeline_worker',
    'add_stack_mask',
    'align_images',
    'compute_relative_scales',
    'create_error',
    'create_mask',
    'detrend_stack',
    'generate_wcs',
    'get_fieldcenter',
    'mask_satellites',
    'remove_pv_distortion',
    'stack_data',
    'stack_science_frames',
]
