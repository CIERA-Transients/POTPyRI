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
