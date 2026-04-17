"""
Automatic generalized pipeline for imaging reduction. Creates coadded images
for each target.  Files are automatically sorted and categorized based on
expected header keywords.  Pixel-level processing is performed with ccdproc,
WCS alignment is performed with astrometry.net and astropy/photutils,
photometry is performed with photutils, and flux calibration with astroquery
and numpy.

The reduction flow is implemented by :class:`potpyri.pipelines.ImagingReductionPipeline`
as an ordered sequence of stages (similar in spirit to ``hispecdrp`` pipelines).

Authors: Kerry Paterson, Charlie Kilpatrick
This project was funded with support by the National Science Foundation under
grant Nos. AST-1814782, AST-1909358 and CAREER grant No. AST-2047919.
If you use this code for your work, please consider citing the package release.
"""

from __future__ import annotations

from typing import Any

from potpyri.pipelines import ImagingReductionPipeline
from potpyri.utils import options


def main_pipeline(instrument: str, data_path: str, **kwargs: Any) -> None:
    """Run the full reduction pipeline: sort files, calibrate, stack, photometry, zeropoint.

    Positional ``instrument`` and ``data_path`` match the CLI:
    ``main_pipeline INSTRUMENT DATA_PATH [options]``.

    Other parameters are the same as :func:`options.add_options` (e.g. ``target``,
    ``proc``, ``include_bad``, ``no_redo_sort``, ``file_list_name``, photometry
    and processing flags). Undeclared keys on ``kwargs`` are ignored.

    Delegates to :class:`~potpyri.pipelines.ImagingReductionPipeline`.
    """
    kw = dict(kwargs)
    if "include_bad" in kw and "incl_bad" not in kw:
        kw["incl_bad"] = kw.pop("include_bad")

    allowed = set(ImagingReductionPipeline.ARGUMENTS.keys())
    # Omit None so :class:`ImagingReductionPipeline` keeps defaults from ARGUMENTS;
    # keep False and 0.
    apply_kw = {k: v for k, v in kw.items() if k in allowed and v is not None}

    ImagingReductionPipeline().apply(
        input=data_path,
        instrument=instrument,
        **apply_kw,
    )


def main() -> None:
    """Entry point: check dependencies, parse args, and run :func:`main_pipeline`."""
    options.test_for_dependencies()
    args = options.add_options()
    main_pipeline(**vars(args))


if __name__ == "__main__":
    main()
