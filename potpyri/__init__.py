"""POTPyRI: Data reduction pipeline for imaging from large aperture telescopes.

This package provides instrument-specific reduction workflows for bias, dark,
flat calibration, image alignment, stacking, WCS solving, and photometry.
"""

from ._version import version as __version__

from .instruments import *
from .primitives import *
