"""Primitive that builds the pipeline raw-file table."""
from __future__ import annotations

from potpyri.primitives.base_primitive import BasePrimitive

from .collect_file_table_worker import _collect_file_table_worker


class FileSortingPrimitive(BasePrimitive):
    """Discover raw files, classify by header, and build the pipeline file list table."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'incl_bad': False,
        'proc': None,
        'no_redo': False,
    }

    def _perform(self):
        out = _collect_file_table_worker(
            self.input,
            self.paths,
            self.tel,
            incl_bad=self.incl_bad,
            proc=self.proc,
            no_redo=self.no_redo,
            log=self.log,
        )
        return {'output': out}
