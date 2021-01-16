#!/usr/bin/env python

"Python function to return commonly used catalog ID for zp calculation."
"Author: Kerry Paterson."

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

    if catalog == 'SDSS':
        catalog_ID, ra, dec, mag = 'V/147', 'RA_ICRS', 'DE_ICRS', fil.upper()+'mag'
    elif catalog == '2MASS':
        catalog_ID, ra, dec, mag = 'II/246', 'RAJ2000', 'DEJ2000', fil.upper()+'mag'
    elif catalog_ID == 'UKIRT':
        catalog_ID, ra, dec, mag = 'II/319', 'ra', 'dec', fil.lower() + '_m'

    return catalog_ID, ra, dec, mag
