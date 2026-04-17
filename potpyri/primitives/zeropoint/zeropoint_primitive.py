"""Zeropoint calibration primitive and :func:`find_zeropoint` entry point."""
from __future__ import annotations

from potpyri.primitives.base_primitive import BasePrimitive

from .catalog_zeropoint import ZeropointFitter


class ZeropointCalibrationPrimitive(BasePrimitive):
    """Fit and write catalog-based zeropoint (see :func:`find_zeropoint`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'magsys': None,
    }

    def _perform(self):
        cal = ZeropointFitter(magsys=self.magsys)
        cal.find_zeropoint(self.input, self.tel, log=self.log)
        return {'output': None}


def find_zeropoint(stack, tel, log=None, magsys=None):
    """Compute and write zeropoint to stack FITS using the instrument catalog (e.g. PS1).

    Delegates to :class:`ZeropointCalibrationPrimitive`.

    Parameters
    ----------
    stack : str
        Path to stacked science FITS file (SCI + APPPHOT extensions).
    tel : Instrument
        Instrument instance (for catalog and filter).
    log : ColoredLogger, optional
        Logger for progress and errors.
    magsys : str, optional
        Magnitude system for ``ZPTMAG`` / ``MAGSYS``; default :data:`DEFAULT_MAGSYS`.

    Returns
    -------
    None
        Stack FITS is updated in place with ZPTMAG, MAGSYS, ZPTNSTAR, etc.
    """
    return ZeropointCalibrationPrimitive(magsys=magsys).apply(
        input=stack, tel=tel, log=log)
