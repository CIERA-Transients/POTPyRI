#!/usr/bin/env python

"Python function to return commonly used catalog ID for zp calculation."
"Author: Kerry Paterson."

# Updated by CDK - 13 Mar 2021

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

    if catalog == 'SDSS':
        if fil.lower() not in ['u','g','r','i','z']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'V/147', 'RA_ICRS', 'DE_ICRS', fil.lower()+'mag', 'e_'+fil.lower()+'mag'
    elif catalog == '2MASS':
        if fil.upper() not in ['J','H','K']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/246', 'RAJ2000', 'DEJ2000', fil.upper()+'mag', 'e_'+fil.upper()+'mag'
    elif catalog == 'UKIRT':
        if fil.upper() not in ['Y','J','H','K']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/319', 'ra', 'dec', fil.upper() + 'mag', 'e_'+fil.upper()+'mag'
    elif catalog == 'PS1':
        if fil.lower() not in ['g','r','i','z','y']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/349', 'RAJ2000', 'DEJ2000', fil.lower()+'mag', 'e_'+fil.lower()+'mag'
    elif catalog == 'SkyMapper':
        if fil.lower() not in ['u','v','g','r','i','z']: return(catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = 'II/358', 'RAICRS', 'DEICRS', fil.lower()+'PSF', 'e_'+fil.lower()+'PSF'

    return(catalog_ID, ra, dec, mag, err)
