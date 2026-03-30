"""SortFilesPrimitive."""
from __future__ import annotations

from potpyri.primitives.base_primitive import BasePrimitive

from .handle_files_impl import _handle_files_worker


class SortFilesPrimitive(BasePrimitive):
    """Discover raw files, classify, and build the pipeline file list table."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'incl_bad': False,
        'proc': None,
        'no_redo': False,
    }

    def _perform(self):
        out = _handle_files_worker(
            self.input,
            self.paths,
            self.tel,
            incl_bad=self.incl_bad,
            proc=self.proc,
            no_redo=self.no_redo,
            log=self.log,
        )
        return {'output': out}
