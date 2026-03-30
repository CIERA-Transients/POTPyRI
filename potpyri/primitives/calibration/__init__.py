"""Master bias, dark, and flat creation from pipeline file tables.

Submodules :mod:`bias`, :mod:`dark`, and :mod:`flat` each define a primitive
class and ``do_*`` entry point. This package re-exports the public API for
backward compatibility with ``from potpyri.primitives import calibration``.

Authors: Charlie Kilpatrick.
"""
from __future__ import annotations

from .bias import BiasPrimitive
from .bias import bias_calibration_worker
from .bias import do_bias
from .dark import DarkPrimitive
from .dark import dark_calibration_worker
from .dark import do_dark
from .flat import FlatPrimitive
from .flat import flat_calibration_worker
from .flat import do_flat

# Deprecated aliases for tests / external code that used private worker names
_bias_calibration_worker = bias_calibration_worker
_dark_calibration_worker = dark_calibration_worker
_flat_calibration_worker = flat_calibration_worker

__all__ = [
    'BiasPrimitive',
    'DarkPrimitive',
    'FlatPrimitive',
    'bias_calibration_worker',
    'dark_calibration_worker',
    'flat_calibration_worker',
    'do_bias',
    'do_dark',
    'do_flat',
    '_bias_calibration_worker',
    '_dark_calibration_worker',
    '_flat_calibration_worker',
]
