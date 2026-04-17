"""Named stages for the imaging reduction pipeline (modular, composable)."""
from __future__ import annotations

import time
from collections.abc import Callable, Sequence
import numpy as np

from potpyri._version import __version__
from potpyri.instruments import instrument_getter
from potpyri.primitives import calibration
from potpyri.primitives import photometry
from potpyri.primitives import sorting
from potpyri.primitives import stacking
from potpyri.primitives import zeropoint
from potpyri.utils import logger as logger_mod
from potpyri.utils import options

from .context import PipelineContext

StageFn = Callable[[PipelineContext], None]


def stage_setup(ctx: PipelineContext) -> None:
    """Load instrument, paths, and log; honor ``skip_flatten`` on flat field."""
    p = ctx.pipeline
    p.tel = instrument_getter(p.instrument)
    if p.skip_flatten:
        p.tel.flat = False
    flist = p.file_list_name or "file_list.txt"
    p.paths = options.add_paths(p.input, flist, p.tel)
    p.log = logger_mod.get_log(p.paths["log"])
    p.log.info("Running main pipeline version %s", __version__)
    p.log.info("Running instrument parameter file version %s", p.tel.version)


def stage_file_sorting(ctx: PipelineContext) -> None:
    """Build or reload the pipeline file table."""
    p = ctx.pipeline
    p.file_table = sorting.collect_file_table(
        p.paths["filelist"],
        p.paths,
        p.tel,
        incl_bad=bool(p.incl_bad),
        proc=p.proc,
        no_redo=bool(p.no_redo_sort),
        log=p.log,
    )


def stage_calibration(ctx: PipelineContext) -> None:
    """Combine bias, dark, and flat masters."""
    p = ctx.pipeline
    calibration.do_bias(p.file_table, p.tel, p.paths, log=p.log)
    calibration.do_dark(p.file_table, p.tel, p.paths, log=p.log)
    calibration.do_flat(p.file_table, p.tel, p.paths, log=p.log)


def stage_science_targets(ctx: PipelineContext) -> None:
    """Stack each science ``TargType``, then run photometry and zeropoint."""
    p = ctx.pipeline
    log = p.log
    ft = p.file_table
    tel = p.tel
    paths = p.paths

    sci_mask = tel.match_type_keywords(tel.filetype_keywords["SCIENCE"], ft)
    science_rows = ft[sci_mask]
    targtypes = np.unique(science_rows["TargType"])

    if p.target:
        needle = str(p.target)
        targtypes = np.array(
            [t for t in targtypes if needle in str(t)],
            dtype=targtypes.dtype,
        )
        if len(targtypes) == 0:
            log.warning("No TargType matched --target %r; skipping science loop.", needle)
            return

    skip_cr = bool(p.skip_cr)
    rel_cal = bool(p.relative_calibration or False)

    for tar in targtypes:
        log.info("Generating stack for %s", tar)
        subset = ft[ft["TargType"] == tar]
        stack_path = stacking.stack_science_frames(
            subset,
            tel,
            paths,
            skip_skysub=bool(p.skip_skysub),
            fieldcenter=p.fieldcenter,
            out_size=p.out_size,
            cosmic_ray=not skip_cr,
            skip_gaia=bool(p.skip_gaia),
            keep_all_astro=bool(p.keep_all_astro),
            relative_calibration=rel_cal,
            log=log,
        )
        log.info("Running PSF photometry.")
        photometry.run_stack_photometry(
            stack_path,
            phot_sn_min=p.phot_sn_min,
            fwhm_init=p.fwhm_init,
            phot_sn_max=p.phot_sn_max,
            log=log,
        )
        log.info("Calculating zeropoint.")
        zeropoint.find_zeropoint(stack_path, tel, log=log)


def stage_finalize(ctx: PipelineContext) -> None:
    """Log total runtime and shut down the reduction logger."""
    p = ctx.pipeline
    t_end = time.time()
    p.log.info("Pipeline finished.")
    p.log.info("Total runtime: %s sec", t_end - ctx.t_start)
    p.log.shutdown()


# Default stage order (HISPEC-style explicit sequence; subclasses may replace).
DEFAULT_STAGE_SEQUENCE: Sequence[tuple[str, StageFn]] = (
    ("setup", stage_setup),
    ("file_sorting", stage_file_sorting),
    ("calibration", stage_calibration),
    ("science_targets", stage_science_targets),
    ("finalize", stage_finalize),
)
