"""Instrument implementations and factory for POTPyRI.

Each instrument module defines a subclass of Instrument with detector keywords,
calibration behavior, and file-sorting rules. New instruments must be added
here and to __all__.
"""
from . import BINOSPEC
from . import DEIMOS
from . import F2
from . import FOURSTAR
from . import GMOS
from . import IMACS
from . import LRIS
from . import MMIRS
from . import MOSFIRE

# New instruments need to be added here and to the function below
__all__ = ["BINOSPEC","DEIMOS","F2","FOURSTAR","GMOS","IMACS","LRIS","MMIRS","MOSFIRE"]


def instrument_getter(instname, log=None):
    """Return an Instrument instance for the given instrument name.

    Parameters
    ----------
    instname : str
        Instrument name (e.g. 'GMOS', 'LRIS', 'BINOSPEC').
    log : ColoredLogger, optional
        Logger for error messages. If None, raises Exception on unsupported
        instrument.

    Returns
    -------
    Instrument
        Subclass instance for the requested instrument.

    Raises
    ------
    Exception
        If instname is not in __all__ and log is None.
    """
    instname = instname.upper()

    if instname not in __all__:
        if log:
            log.error(f'Instrument {instname} not supported by POTPyRI')
        else:
            raise Exception(f'Instrument {instname} not supported by POTPyRI')

    tel = None

    if instname=="BINOSPEC":
        tel = BINOSPEC.BINOSPEC()
    elif instname=="DEIMOS":
        tel = DEIMOS.DEIMOS()
    elif instname=="F2":
        tel = F2.F2()
    elif instname=="FOURSTAR":
        tel = FOURSTAR.FOURSTAR()
    elif instname=="GMOS":
        tel = GMOS.GMOS()
    elif instname=="IMACS":
        tel = IMACS.IMACS()
    elif instname=="LRIS":
        tel = LRIS.LRIS()
    elif instname=="MMIRS":
        tel = MMIRS.MMIRS()
    elif instname=="MOSFIRE":
        tel = MOSFIRE.MOSFIRE()

    return(tel)
