from astropy import units as u
from astropy.coordinates import SkyCoord
import warnings
import csv
import sys
csv.field_size_limit(sys.maxsize)

warnings.filterwarnings('ignore')

viziercat = {
    'sdssdr12': {'name':'V/147',
        'columns': ['RA_ICRS', 'DE_ICRS','class','umag', 'e_umag',
            'gmag','e_gmag', 'rmag','e_rmag', 'imag','i_mag', 'zmag',
            'e_zmag', 'zph']
    },
    '2mass': {'name':'II/246',
        'columns':['RAJ2000', 'DEJ2000', 'Jmag','e_Jmag','Hmag','e_Hmag',
            'Kmag', 'e_Kmag']
    },
    'unwise': {'name':'II/363',
        'columns': ['RAJ2000', 'DEJ2000', 'FW1','e_FW1', 'FW2','e_FW2']
    },
    'glade': {'name':'VII/281',
        'columns': ['RAJ2000', 'DEJ2000', 'Dist', 'e_Dist', 'Bmag', 'Jmag',
            'Hmag', 'Kmag', 'z']
    },
    'des': {'name':'II/357',
        'columns': ['RAJ2000', 'DEJ2000', 'S/Gg', 'S/Gr', 'S/Gi', 'S/Gz',
            'gmag','e_gmag', 'rmag','e_rmag', 'imag','e_imag', 'zmag','e_zmag']
    }
}

# Sort the calibration files:
def find_catalog(catalog,fil):

    '''

    Function to return catalog ID for Vizier query.

    Parameters
    ----------

    :param catalog: str
        Name of catalog.

    Returns
    -------

    :return: list (str)
        Catalog ID and column information to from Vizier.

    '''

    catalog_ID, ra, dec, mag, err = None, None, None, None, None

    # If these catalogs are to be updated in the future, select mag columns that correspond to the
    # PSF mags.
    if catalog.upper() == 'SDSS':
        if fil.lower() not in ['u','g','r','i','z']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'V/154', 'RA_ICRS', 'DE_ICRS', fil.lower()+'mag', 'e_'+fil.lower()+'mag'
    elif catalog.upper() == '2MASS':
        if fil.upper() not in ['J','H','K']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/246', 'RAJ2000', 'DEJ2000', fil.upper()+'mag', 'e_'+fil.upper()+'mag'
    elif catalog.upper() == 'UKIRT':
        if fil.upper() not in ['Y','J','H','K']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/319', 'ra', 'dec', fil.upper() + 'mag', 'e_'+fil.upper()+'mag'
    elif catalog.upper() == 'PS1':
        if fil.lower() not in ['g','r','i','z','y']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/349', 'RAJ2000', 'DEJ2000', fil.lower()+'mag', 'e_'+fil.lower()+'mag'
    elif catalog.upper() == 'SKYMAPPER':
        if fil.lower() not in ['u','v','g','r','i','z']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/379', 'RAICRS', 'DEICRS', fil.lower()+'PSF', 'e_'+fil.lower()+'PSF'

    return(catalog_ID, ra, dec, mag, err)

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def parse_coord(ra, dec):
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


