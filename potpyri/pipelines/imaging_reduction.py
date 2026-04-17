"""Full imaging reduction pipeline as a modular :class:`BasePrimitive` (HISPEC-style stages)."""
from __future__ import annotations

import logging
import time
from collections.abc import Sequence
from typing import Any, ClassVar

from potpyri.primitives.base_primitive import BasePrimitive

from .context import PipelineContext
from .stages import DEFAULT_STAGE_SEQUENCE, StageFn

_stage_log = logging.getLogger(__name__)


class ImagingReductionPipeline(BasePrimitive):
    """Ordered stages (setup → sorting → calibration → science targets → finalize).

    Mirrors the structure of ``hispecdrp.pipelines.BasePipeline``: a declared
    sequence of named steps, each a callable taking a :class:`PipelineContext`.
    Call :meth:`apply` with ``input=data_path`` and ``instrument=...`` plus
    processing keyword arguments (same names as the CLI).

    Attributes set during :meth:`_perform`
    --------------------------------------
    tel, paths, log, file_table
        Populated by stages; useful for tests or follow-on code.
    """

    ARGUMENTS: ClassVar[dict[str, Any]] = {
        **BasePrimitive.ARGUMENTS,
        "instrument": None,
        "target": None,
        "proc": None,
        "incl_bad": False,
        "no_redo_sort": False,
        "phot_sn_min": 3.0,
        "phot_sn_max": 20.0,
        "fwhm_init": 5.0,
        "skip_skysub": False,
        "file_list_name": None,
        "fieldcenter": None,
        "out_size": None,
        "skip_flatten": False,
        "skip_cr": False,
        "skip_gaia": False,
        "keep_all_astro": False,
        "relative_calibration": False,
    }

    #: Override in subclasses to replace or reorder stages.
    STAGE_SEQUENCE: ClassVar[Sequence[tuple[str, StageFn]]] = DEFAULT_STAGE_SEQUENCE

    def get_stage_sequence(self) -> Sequence[tuple[str, StageFn]]:
        """Return the ordered ``(name, callable)`` stages for this run."""
        return self.STAGE_SEQUENCE

    def set_args(self, input: Any = None, **kwargs: Any) -> None:
        """Bind ``input`` (data path) and pipeline kwargs; accept ``include_bad`` alias."""
        if "include_bad" in kwargs and "incl_bad" not in kwargs:
            kwargs = {**kwargs, "incl_bad": kwargs.pop("include_bad")}
        super().set_args(input=input, **kwargs)

    def _perform(self) -> dict[str, Any]:
        if not self.instrument:
            raise ValueError("instrument is required")
        if self.input is None:
            raise ValueError("input (data_path) is required")

        ctx = PipelineContext(pipeline=self, t_start=time.time())

        for name, stage_fn in self.get_stage_sequence():
            t0 = time.time()
            _stage_log.info("Pipeline stage: %s", name)
            stage_fn(ctx)
            ctx.stage_runtimes[name] = time.time() - t0
            log = getattr(self, "log", None)
            if log is not None and name != "finalize":
                log.info("Finished stage %s in %.2f s", name, ctx.stage_runtimes[name])

        return {"output": self.file_table}

    def __repr__(self) -> str:
        lines = [f"{self.__class__.__name__}:", "  Stages:"]
        for name, fn in self.get_stage_sequence():
            lines.append(f"    {name}: {fn.__name__}")
        return "\n".join(lines)
