"""Instrument implementations and factory for POTPyRI.

Each instrument module defines a subclass of :class:`~potpyri.instruments.instrument.Instrument`
with detector keywords, calibration behavior, and file-sorting rules. To add a
new instrument, create its module and add the canonical name to ``__all__`` below.
"""
from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .instrument import Instrument

# Canonical instrument names (single source of truth). Each name must match a
# submodule ``potpyri.instruments.<NAME>`` with class ``<NAME>``.
__all__ = [
    'BINOSPEC',
    'DEIMOS',
    'F2',
    'FOURSTAR',
    'GMOS',
    'IMACS',
    'LRIS',
    'MMIRS',
    'MOSFIRE',
]

# Optional CLI / shorthand aliases (not duplicated in __all__).
_INSTRUMENT_ALIASES: dict[str, str] = {
    'BINO': 'BINOSPEC',
    'MMIR': 'MMIRS',
}

_MODULES = {name: import_module(f'.{name}', __name__) for name in __all__}

# Re-export instrument submodules: ``from potpyri.instruments import GMOS``
for _name, _mod in _MODULES.items():
    globals()[_name] = _mod


class UnknownInstrumentError(ValueError):
    """Raised when an instrument name is not supported."""

    def __init__(self, name: str, *, resolved: str | None = None) -> None:
        self.name = name
        self.resolved = resolved
        supported = ', '.join(__all__)
        aliases = ', '.join(f'{k}->{v}' for k, v in sorted(_INSTRUMENT_ALIASES.items()))
        hint = f' Supported instruments: {supported}.'
        if aliases:
            hint += f' Aliases: {aliases}.'
        if resolved and resolved != name.strip().upper():
            detail = f'{name!r} (resolved to {resolved!r})'
        else:
            detail = repr(name)
        super().__init__(f'Instrument {detail} is not supported by POTPyRI.{hint}')


def supported_instruments() -> tuple[str, ...]:
    """Return supported canonical instrument names."""
    return tuple(__all__)


def resolve_instrument_name(name: str) -> str:
    """Normalize *name* to a canonical instrument key (uppercase, aliases applied).

    Parameters
    ----------
    name : str
        User-supplied instrument name or alias (e.g. ``'gmos'``, ``'BINO'``).

    Returns
    -------
    str
        Canonical name from :data:`__all__` when recognized, otherwise the
        uppercased input (may still be unknown to :func:`instrument_getter`).
    """
    key = name.strip().upper()
    if key in __all__:
        return key
    if key in _INSTRUMENT_ALIASES:
        return _INSTRUMENT_ALIASES[key]
    # Legacy substring aliases used by older scripts / partial names.
    if 'BINO' in key:
        return 'BINOSPEC'
    if 'MMIR' in key:
        return 'MMIRS'
    return key


def instrument_getter(instname: str, log=None) -> Instrument | None:
    """Return an :class:`~potpyri.instruments.instrument.Instrument` for *instname*.

    Parameters
    ----------
    instname : str
        Instrument name or alias (e.g. ``'GMOS'``, ``'BINO'``, ``'MMIR'``).
    log : ColoredLogger, optional
        If provided and the instrument is unknown, log an error and return
        ``None`` instead of raising.

    Returns
    -------
    Instrument or None
        Configured instrument instance, or ``None`` if unknown and *log* is set.

    Raises
    ------
    UnknownInstrumentError
        If *instname* is not supported and *log* is ``None``.
    TypeError
        If *instname* is missing or not a string.
    """
    if not instname or not isinstance(instname, str):
        raise TypeError(
            f'instrument name must be a non-empty string, got {instname!r}'
        )

    canonical = resolve_instrument_name(instname)
    if canonical not in __all__:
        exc = UnknownInstrumentError(instname, resolved=canonical)
        if log is not None:
            log.error(str(exc))
            return None
        raise exc

    cls = getattr(_MODULES[canonical], canonical)
    return cls()
