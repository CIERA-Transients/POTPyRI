"""Base class for POTPyRI pipeline primitives.

Follows the same structural pattern as ``hispecdrp.primitives.BasePrimitive``:
configuration, default arguments, ``apply`` orchestration, and hooks
(``on_skip``, ``on_exception``, ``on_pass``).  This implementation is
standalone (no datamodel or calibration-store dependencies).

Authors: Charlie Kilpatrick.
"""
from __future__ import annotations

import copy
import datetime
import logging
from abc import ABC, abstractmethod
from typing import Any

from potpyri._version import __version__

logger = logging.getLogger(__name__)

__all__ = ["BasePrimitive"]


class BasePrimitive(ABC):
    """Configurable pipeline step with ``apply``, configuration, and lifecycle hooks.

    Subclasses define :attr:`ARGUMENTS` defaults, optional :attr:`CALIBRATIONS`,
    and implement :meth:`_perform`.

    Parameters
    ----------
    config : dict or None
        Optional configuration dict. Merged with defaults from :attr:`ARGUMENTS`
        and :attr:`LOGGING`.
    **kwargs
        Overrides for keys in :attr:`ARGUMENTS` (and passed to :meth:`init_args`).
    """

    ARGUMENTS: dict[str, Any] = {
        "run": True,
        "save_result": False,
        "output_dir": None,
    }

    LOGGING: dict[str, Any] = {
        "level": "INFO",
    }

    CALIBRATIONS: dict[str, Any] = {}

    def __init__(self, config: dict | None = None, **kwargs: Any) -> None:
        self._output: dict[str, Any] | None = None
        self._status: str | None = None
        self._datetime_start: str | None = None
        self._datetime_end: str | None = None
        self.input: Any = None

        self.init_config(config)
        self.init_args(**kwargs)

    def init_config(self, config: dict | None = None) -> dict[str, Any]:
        """Build ``self.config`` from defaults and optional user ``config`` dict."""
        class_name = f"{self.__class__.__module__}.{self.__class__.__name__}"
        if config is None:
            self.config = {
                "class": class_name,
                "arguments": copy.deepcopy(self.ARGUMENTS),
                "logging": copy.deepcopy(self.LOGGING),
            }
            return self.config

        if not isinstance(config, dict):
            raise TypeError(f"config must be dict or None, got {type(config)!r}")

        self.config = copy.deepcopy(config)
        self.config["class"] = class_name
        merged_args = copy.deepcopy(self.ARGUMENTS)
        merged_args.update(self.config.get("arguments", {}))
        self.config["arguments"] = merged_args
        merged_log = copy.deepcopy(self.LOGGING)
        merged_log.update(self.config.get("logging", {}))
        self.config["logging"] = merged_log
        return self.config

    def init_args(self, **kwargs: Any) -> None:
        """Set attributes from ``self.config['arguments']`` and ``kwargs``."""
        args = copy.deepcopy(self.config["arguments"])
        args.update(kwargs)
        for key, value in args.items():
            setattr(self, key, value)

    def set_args(self, input: Any = None, **kwargs: Any) -> None:
        """Called at the start of :meth:`apply` to bind ``input`` and call-time kwargs."""
        self.input = input
        for key, value in kwargs.items():
            setattr(self, key, value)

    def _pre_condition(self) -> bool:
        """Return False to skip :meth:`_perform` (sets status SKIPPED)."""
        return bool(getattr(self, "run", True))

    def apply(self, input: Any = None, **kwargs: Any) -> Any:
        """Run the primitive: :meth:`set_args`, optional skip, :meth:`_perform`, hooks.

        Returns
        -------
        Any
            The main result, normally ``self._output['output']`` from :meth:`_perform`.
        """
        self._output = None
        self._status = None
        self._datetime_start = datetime.datetime.now().isoformat(timespec="milliseconds")

        try:
            self.set_args(input=input, **kwargs)

            if not self._pre_condition():
                logger.info("Skipping %s with args %s", self.__class__.__name__, self._repr_args())
                self._status = "SKIPPED"
                self.on_skip()
                self._output = {"output": self.input}
                self.finalize_output()
                return self._output.get("output")

            logger.info("Applying %s with args %s", self.__class__.__name__, self._repr_args())
            self._output = self._perform()
            if not isinstance(self._output, dict):
                self._output = {"output": self._output}
            if "output" not in self._output:
                self._output = {"output": self._output}
            if self._status is None:
                self._status = "COMPLETE"
            out = self._output.get("output")
            self.on_pass(out)
            self.finalize_output()
            return out

        except Exception as exc:
            logger.error(
                "Failed executing primitive %s: %s",
                self.__class__.__name__,
                exc,
                exc_info=True,
            )
            self._status = "FAILED"
            self._output = {"output": input, "error": exc}
            self.on_exception(exc)
            self.finalize_output()
            raise

    def on_exception(self, exc: BaseException) -> None:
        """Hook when :meth:`apply` raises (after status is set to FAILED)."""
        del exc  # default: no-op

    def on_skip(self) -> None:
        """Hook when :meth:`_pre_condition` is False."""
        pass

    def on_pass(self, output: Any) -> None:
        """Hook after successful :meth:`_perform` (receives main output)."""
        del output

    def finalize_output(self) -> None:
        """Record end time; subclasses may add provenance or I/O."""
        self._datetime_end = datetime.datetime.now().isoformat(timespec="milliseconds")

    def _repr_args(self) -> str:
        parts = [f"input={self.input!r}"]
        for key in self.ARGUMENTS:
            if hasattr(self, key):
                parts.append(f"{key}={getattr(self, key)!r}")
        return ", ".join(parts)

    def __repr__(self) -> str:
        s = self.__class__.__name__
        for key in self.ARGUMENTS:
            if hasattr(self, key):
                s += f"\n  {key}: {getattr(self, key)}"
        return s

    @abstractmethod
    def _perform(self) -> Any:
        """Execute the primitive; return the main result or ``dict`` with ``'output'`` key."""

    def add_receipt_metadata(self, header: Any) -> None:
        """Optional: write processing provenance to a FITS header (subclasses may use)."""
        if header is None:
            return
        try:
            header["POTPYVER"] = (__version__, "POTPyRI version")
            header["POTPRIM"] = (self.config.get("class", ""), "Primitive class")
            if self._datetime_start:
                header["POTPTIME"] = (self._datetime_start, "Primitive start time")
        except Exception:
            logger.debug("Could not add receipt metadata to header", exc_info=True)
