"""Science stacking primitive and :func:`stack_science_frames` entry point."""
from __future__ import annotations

from potpyri.primitives.base_primitive import BasePrimitive

from .stacking_pipeline_worker import _stacking_pipeline_worker


class ScienceStackingPrimitive(BasePrimitive):
    """Stack and calibrate science frames for one target (see :func:`stack_science_frames`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'skip_skysub': False,
        'fieldcenter': None,
        'out_size': None,
        'satellites': True,
        'cosmic_ray': True,
        'skip_gaia': False,
        'keep_all_astro': False,
        'relative_calibration': False,
    }

    def _perform(self):
        out = _stacking_pipeline_worker(
            self.input,
            self.tel,
            self.paths,
            skip_skysub=self.skip_skysub,
            fieldcenter=self.fieldcenter,
            out_size=self.out_size,
            satellites=self.satellites,
            cosmic_ray=self.cosmic_ray,
            skip_gaia=self.skip_gaia,
            keep_all_astro=self.keep_all_astro,
            relative_calibration=self.relative_calibration,
            log=self.log,
        )
        return {'output': out}


def stack_science_frames(image_data, tel, paths, skip_skysub=False,
    fieldcenter=None, out_size=None, satellites=True, cosmic_ray=True,
    skip_gaia=False, keep_all_astro=False, relative_calibration=False,
    log=None):
    """Align, mask, stack, and optionally detrend science frames for one table row group.

    Delegates to :class:`ScienceStackingPrimitive`.

    Parameters
    ----------
    image_data : astropy.table.Table
        File table subset (e.g. one TargType) with 'File' column.
    tel : Instrument
        Instrument instance.
    paths : dict
        Paths dict from options.add_paths.
    skip_skysub : bool, optional
        If True, skip sky subtraction. Default is False.
    fieldcenter : sequence, optional
        [ra, dec] for alignment.
    out_size : int, optional
        Output image size.
    satellites : bool, optional
        If True, mask satellite trails. Default is True.
    cosmic_ray : bool, optional
        If True, run cosmic-ray rejection. Default is True.
    skip_gaia : bool, optional
        If True, skip Gaia alignment. Default is False.
    keep_all_astro : bool, optional
        If True, keep all images regardless of astrometric dispersion.
    relative_calibration : bool, optional
        If True, run relative flux calibration before stacking.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    str or None
        Path to stacked output FITS, or None if processing failed.
    """
    return ScienceStackingPrimitive(
        skip_skysub=skip_skysub,
        fieldcenter=fieldcenter,
        out_size=out_size,
        satellites=satellites,
        cosmic_ray=cosmic_ray,
        skip_gaia=skip_gaia,
        keep_all_astro=keep_all_astro,
        relative_calibration=relative_calibration,
    ).apply(input=image_data, tel=tel, paths=paths, log=log)
