"""Astrometry.net WCS primitive."""
from __future__ import annotations

from potpyri.primitives.base_primitive import BasePrimitive

from .core import solve_astrometry


class AstrometryNetPrimitive(BasePrimitive):
    """Coarse WCS via astrometry.net ``solve-field`` (see :func:`solve_astrometry`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'radius': 0.5,
        'replace': True,
        'shift_only': False,
    }

    def _perform(self):
        ok = solve_astrometry(
            self.input,
            self.tel,
            self.binn,
            self.paths,
            radius=self.radius,
            replace=self.replace,
            shift_only=self.shift_only,
            log=self.log,
        )
        return {'output': ok}
