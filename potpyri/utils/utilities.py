"""General utilities for coordinates and numeric parsing.

Vizier catalog IDs and :func:`find_catalog` live in :mod:`potpyri.utils.catalogs`;
they are re-exported here for backward compatibility.
"""
from astropy import units as u
from astropy.coordinates import SkyCoord
import warnings
import csv
import sys
csv.field_size_limit(sys.maxsize)

warnings.filterwarnings('ignore')

from . import catalogs as _catalogs

viziercat = _catalogs.viziercat
find_catalog = _catalogs.find_catalog
POINT_SOURCE_CALIBRATION_CATALOGS = _catalogs.POINT_SOURCE_CALIBRATION_CATALOGS

def is_number(num):
    """Return True if the value can be interpreted as a number.

    Parameters
    ----------
    num : str or number
        Value to test.

    Returns
    -------
    bool
        True if float(num) succeeds, False otherwise.
    """
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)


def parse_coord(ra, dec):
    """Parse RA and Dec strings into an astropy SkyCoord.

    Accepts decimal degrees or sexagesimal (e.g. '12:30:00' or '12.5').

    Parameters
    ----------
    ra : str or float
        Right ascension.
    dec : str or float
        Declination.

    Returns
    -------
    SkyCoord or None
        ICRS coordinate, or None if parsing fails.
    """
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in str(ra) and ':' in str(dec)):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)


