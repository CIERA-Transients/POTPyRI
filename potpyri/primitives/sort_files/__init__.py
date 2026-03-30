"""File sorting and file-list generation for the main pipeline.

Classifies raw files (science, flat, bias, dark, bad) from header keywords
and writes a fixed-width file list used by calibration and reduction.
Authors: Owen Eskandari, Kerry Paterson, Charlie Kilpatrick.
"""
from __future__ import annotations

from .classify import (
    is_bad,
    is_bias,
    is_dark,
    is_flat,
    is_science,
    is_spec,
)
from .handle_files_fn import handle_files
from .handle_files_impl import _handle_files_worker
from .primitive import SortFilesPrimitive
from .sort_table import sort_files

__all__ = [
    'SortFilesPrimitive',
    '_handle_files_worker',
    'handle_files',
    'is_bad',
    'is_bias',
    'is_dark',
    'is_flat',
    'is_science',
    'is_spec',
    'sort_files',
]
