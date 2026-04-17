"""High-level imaging pipelines (modular stage sequences)."""

from .context import PipelineContext
from .imaging_reduction import ImagingReductionPipeline
from .stages import DEFAULT_STAGE_SEQUENCE

__all__ = [
    "DEFAULT_STAGE_SEQUENCE",
    "ImagingReductionPipeline",
    "PipelineContext",
]
