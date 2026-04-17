"""Mutable context passed through imaging pipeline stages."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class PipelineContext:
    """Carries the live pipeline instance and per-stage timing.

    Stages read and write attributes on ``pipeline`` (e.g. ``tel``, ``paths``,
    ``log``, ``file_table``).
    """

    pipeline: Any
    t_start: float = 0.0
    stage_runtimes: dict[str, float] = field(default_factory=dict)
