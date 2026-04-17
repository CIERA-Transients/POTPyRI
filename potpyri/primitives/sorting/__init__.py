"""Raw file discovery and classification for the main pipeline.

Classifies FITS files (science, flat, bias, dark, bad) from header keywords
and writes a fixed-width file list for calibration and stacking.
Authors: Owen Eskandari, Kerry Paterson, Charlie Kilpatrick.
"""
from __future__ import annotations

from .collect_file_table import collect_file_table
from .collect_file_table_worker import _collect_file_table_worker
from .header_classification import (
    is_bad,
    is_bias,
    is_dark,
    is_flat,
    is_science,
    is_spec,
)
from .raw_file_table import build_raw_file_table
from .sorting_primitive import FileSortingPrimitive

__all__ = [
    'FileSortingPrimitive',
    '_collect_file_table_worker',
    'build_raw_file_table',
    'collect_file_table',
    'is_bad',
    'is_bias',
    'is_dark',
    'is_flat',
    'is_science',
    'is_spec',
]
