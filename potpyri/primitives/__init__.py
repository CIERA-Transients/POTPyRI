"""Pipeline primitives: sorting, calibration, stacking, astrometry, photometry, zeropoint."""

from .base_primitive import BasePrimitive

from .astrometry import AstrometryNetPrimitive
from .calibration import BiasPrimitive, DarkPrimitive, FlatPrimitive
from .sorting import FileSortingPrimitive
from .stacking import ScienceStackingPrimitive
from .photometry import StackPhotometryPrimitive
from .zeropoint import ZeropointCalibrationPrimitive

from . import astrometry
from . import calibration
from . import sorting
from . import stacking
from . import photometry
from . import zeropoint

__all__ = [
    'BasePrimitive',
    'AstrometryNetPrimitive',
    'BiasPrimitive',
    'DarkPrimitive',
    'FlatPrimitive',
    'FileSortingPrimitive',
    'ScienceStackingPrimitive',
    'StackPhotometryPrimitive',
    'ZeropointCalibrationPrimitive',
    'astrometry',
    'calibration',
    'sorting',
    'stacking',
    'photometry',
    'zeropoint',
]
