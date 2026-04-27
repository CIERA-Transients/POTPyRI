"""VizieR and catalog metadata for astrometry and photometric calibration.

Photometric reference catalogs (PS1, 2MASS, SkyMapper, SDSS, etc.) and
astrometric queries (e.g. Gaia DR3) live here so primitives such as
``solve_wcs`` and ``absphot`` can share mirror fallback and column metadata
without duplicating ``astroquery`` usage.

Point-source catalogs (stars and compact objects) suitable for calibration are
documented in :data:`POINT_SOURCE_CALIBRATION_CATALOGS` with VizieR IDs and
suggested column lists for ``Vizier(columns=...)`` queries. Extended-source
catalogs (e.g. galaxy lists) are intentionally excluded.

Authors: Kerry Paterson, Charlie Kilpatrick.
"""
import numpy as np
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.table import Table

# Vizier mirror servers (hostnames only); tried in order when a query fails.
VIZIER_MIRRORS = [
    'vizier.cds.unistra.fr',   # CDS / Strasbourg, France
    'vizier.nao.ac.jp',        # ADAC / Tokyo, Japan
    'vizier.cfa.harvard.edu',  # CfA / Harvard, USA
    'vizier.inasan.ru',        # INASAN / Moscow, Russia
    'vizier.iucaa.in',         # IUCAA / Pune, India
    'vizier.idia.ac.za',       # IDIA / South Africa
]

# Photometric / multi-band catalogs: VizieR IDs and default column sets.
# Keys are internal registry names (``sdss`` replaces the former ``sdssdr12``).
viziercat = {
    'sdss': {'name': 'V/147',
        'columns': ['RA_ICRS', 'DE_ICRS', 'class', 'umag', 'e_umag',
            'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'i_mag', 'zmag',
            'e_zmag', 'zph']
    },
    '2mass': {'name': 'II/246',
        'columns': ['RAJ2000', 'DEJ2000', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag',
            'Kmag', 'e_Kmag']
    },
    'unwise': {'name': 'II/363',
        'columns': ['RAJ2000', 'DEJ2000', 'FW1', 'e_FW1', 'FW2', 'e_FW2']
    },
    'des': {'name': 'II/357',
        'columns': ['RAJ2000', 'DEJ2000', 'S/Gg', 'S/Gr', 'S/Gi', 'S/Gz',
            'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag']
    },
    'skymapper': {'name': 'II/379/smssdr4',
        'columns': ['RAICRS', 'DEICRS', 'uPSF', 'e_uPSF', 'gPSF', 'e_gPSF',
            'rPSF', 'e_rPSF', 'iPSF', 'e_iPSF', 'zPSF', 'e_zPSF']
    },
}

#: Point-source VizieR catalogs useful for WCS fine alignment and/or photometric
#: zeropoint work. ``vizier_id`` is passed to ``astroquery`` as ``catalog=``;
#: ``default_vizier_columns`` are typical inputs to ``Vizier(columns=...)``.
#: Always confirm column names against the current VizieR ReadMe before adding
#: a new code path—CDS tables evolve.
POINT_SOURCE_CALIBRATION_CATALOGS = {
    'gaia_dr3': {
        'description': 'Gaia DR3 (primary astrometry; G/BP/RP photometry)',
        'vizier_id': 'I/355/gaiadr3',
        'ra_column': 'RA_ICRS',
        'dec_column': 'DE_ICRS',
        'default_vizier_columns': [
            'RA_ICRS', 'DE_ICRS', 'Source', 'Gmag', 'BPmag', 'RPmag',
            'e_Gmag', 'e_BPmag', 'e_RPmag', 'Plx', 'PM', 'PSS',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'sdss': {
        'description': 'SDSS DR12 photometric catalog (ugriz)',
        'vizier_id': 'V/147',
        'ra_column': 'RA_ICRS',
        'dec_column': 'DE_ICRS',
        'default_vizier_columns': [
            'RA_ICRS', 'DE_ICRS', 'class', 'umag', 'e_umag', 'gmag', 'e_gmag',
            'rmag', 'e_rmag', 'imag', 'i_mag', 'zmag', 'e_zmag', 'zph',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'ps1': {
        'description': 'Pan-STARRS1 DR1 stacked photometry',
        'vizier_id': 'II/349',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'gmag', 'e_gmag', 'rmag', 'e_rmag',
            'imag', 'e_imag', 'zmag', 'e_zmag', 'ymag', 'e_ymag',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'twomass': {
        'description': '2MASS point source catalog (JHK)',
        'vizier_id': 'II/246',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag',
            'Kmag', 'e_Kmag',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'skymapper': {
        'description': 'SkyMapper Southern Survey DR4 (uvgriz PSF mags)',
        'vizier_id': 'II/379/smssdr4',
        'ra_column': 'RAICRS',
        'dec_column': 'DEICRS',
        'default_vizier_columns': [
            'RAICRS', 'DEICRS', 'uPSF', 'e_uPSF', 'gPSF', 'e_gPSF',
            'rPSF', 'e_rPSF', 'iPSF', 'e_iPSF', 'zPSF', 'e_zPSF',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'des': {
        'description': 'DES DR1 wide-field griz (star–galaxy flags per band)',
        'vizier_id': 'II/357',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'S/Gg', 'S/Gr', 'S/Gi', 'S/Gz',
            'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'unwise': {
        'description': 'unWISE forced photometry at W1/W2 (compact sources)',
        'vizier_id': 'II/363',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'FW1', 'e_FW1', 'FW2', 'e_FW2',
        ],
        'roles': ('photometry',),
    },
    'apass9': {
        'description': 'APASS DR9 optical gri (homogeneous all-sky)',
        'vizier_id': 'II/336',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'g_mag', 'e_g_mag', 'r_mag', 'e_r_mag',
            'i_mag', 'e_i_mag',
        ],
        'roles': ('photometry',),
        'notes': 'If column names differ for your VizieR table version, check II/336 ReadMe.',
    },
    'ucac4': {
        'description': 'UCAC4 (sub-100 mas positions; JHK from 2MASS; optical magnitudes)',
        'vizier_id': 'I/322',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'pmRA', 'pmDE', 'Jmag', 'Hmag', 'Kmag', 'Apmag',
        ],
        'roles': ('astrometry', 'photometry'),
        'notes': 'Apmag is UCAC aperture magnitude between Landolt R and I; verify columns.',
    },
    'usno_b1': {
        'description': 'USNO-B1.0 (photographic BRI; dense reference)',
        'vizier_id': 'I/284',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'B1mag', 'R1mag', 'Imag',
        ],
        'roles': ('astrometry', 'photometry'),
    },
    'allwise': {
        'description': 'AllWISE source catalog (W1–W4)',
        'vizier_id': 'II/328/allwise',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'W1mag', 'W1sigmag', 'W2mag', 'W2sigmag',
            'W3mag', 'W3sigmag', 'W4mag', 'W4sigmag',
        ],
        'roles': ('photometry',),
    },
    'galex_ais': {
        'description': 'GALEX AIS (NUV/FUV; point sources and shallow galaxies)',
        'vizier_id': 'II/312/ais',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'nuv_mag', 'nuv_magerr', 'fuv_mag', 'fuv_magerr',
        ],
        'roles': ('photometry',),
    },
    'gsc23': {
        'description': 'Guide Star Catalog II v2.3 (plate-based optical)',
        'vizier_id': 'I/305',
        'ra_column': 'RAJ2000',
        'dec_column': 'DEJ2000',
        'default_vizier_columns': [
            'RAJ2000', 'DEJ2000', 'Jmag', 'Fmag', 'Class',
        ],
        'roles': ('astrometry', 'photometry'),
        'notes': 'Confirm I/305 table and column names against VizieR (GSC2.3 release).',
    },
}

# VizieR catalog ID for Gaia DR3 astrometry/photometry (used by WCS alignment).
GAIA_DR3_VIZIER_ID = 'I/355/gaiadr3'
GAIA_DR3_ASTROMETRY_COLUMNS = ['RA_ICRS', 'DE_ICRS', 'Plx', 'PSS', 'PM']


def find_catalog(catalog, fil, coord_ra, coord_dec):
    """Return Vizier catalog ID and column names for the given catalog and filter.

    Supports SDSS, 2MASS, UKIRT, PS1, SKYMAPPER. For southern u-band, uses
    SkyMapper automatically.

    Parameters
    ----------
    catalog : str
        Catalog name (e.g. 'PS1', 'SDSS', '2MASS').
    fil : str
        Filter band (e.g. 'r', 'g', 'J').
    coord_ra : float
        Right ascension (used for catalog selection).
    coord_dec : float
        Declination (used for catalog selection; <0 can trigger SkyMapper for u-band).

    Returns
    -------
    tuple
        (catalog, catalog_ID, ra_col, dec_col, mag_col, err_col) for use in
        Vizier queries. catalog_ID/ra/dec/mag/err may be None if filter not supported.
    """
    catalog_ID, ra, dec, mag, err = None, None, None, None, None

    if coord_dec < 0 and fil.lower() == 'u':
        catalog = 'skymapper'

    if catalog.upper() == 'SDSS':
        if fil.lower() not in ['u', 'g', 'r', 'i', 'z']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'V/154', 'RA_ICRS', 'DE_ICRS', fil.lower() + 'mag', 'e_' + fil.lower() + 'mag')
    elif catalog.upper() == '2MASS':
        fil_2mass = fil.upper()
        if fil_2mass in ('KS', 'KSPEC'):
            fil_2mass = 'K'
        if fil_2mass not in ['J', 'H', 'K']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/246', 'RAJ2000', 'DEJ2000', fil_2mass + 'mag', 'e_' + fil_2mass + 'mag')
    elif catalog.upper() == 'UKIRT':
        if fil.upper() not in ['Y', 'J', 'H', 'K']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/319', 'ra', 'dec', fil.upper() + 'mag', 'e_' + fil.upper() + 'mag')
    elif catalog.upper() == 'PS1':
        if fil.lower() not in ['g', 'r', 'i', 'z', 'y']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/349', 'RAJ2000', 'DEJ2000', fil.lower() + 'mag', 'e_' + fil.lower() + 'mag')
    elif catalog.upper() == 'SKYMAPPER':
        if fil.lower() not in ['u', 'v', 'g', 'r', 'i', 'z']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/379/smssdr4', 'RAICRS', 'DEICRS',
            fil.lower() + 'PSF', 'e_' + fil.lower() + 'PSF')

    return (catalog, catalog_ID, ra, dec, mag, err)


def query_vizier_region(center, width, catalog_id, columns, log=None):
    """Query one VizieR catalog in a sky region, trying mirror servers in order.

    Parameters
    ----------
    center : astropy.coordinates.SkyCoord
        Region center.
    width : astropy.units.Quantity
        Angular width (e.g. ``0.5 * u.deg`` or ``20 * u.arcmin``).
    catalog_id : str
        VizieR catalog identifier (e.g. ``'II/349'`` for PS1, ``'I/355/gaiadr3'``).
    columns : list of str
        Column names to request.
    log : ColoredLogger, optional
        Logger for progress and warnings.

    Returns
    -------
    astropy.table.Table or None
        First table of the query result, or None if all mirrors fail or the
        catalog returns no table.
    """
    Vizier.clear_cache()
    last_error = None
    for server in VIZIER_MIRRORS:
        try:
            vizier = Vizier(columns=columns, vizier_server=server)
            vizier.ROW_LIMIT = -1
            result = vizier.query_region(center, width=width, catalog=catalog_id)
            if result is not None and len(result) > 0:
                if log:
                    log.info(f'Vizier query succeeded via {server} ({catalog_id})')
                return result[0]
        except Exception as e:
            last_error = e
            if log:
                log.warning(f'Vizier mirror {server} failed for {catalog_id}: {e}')
            else:
                print(f'Vizier mirror {server} failed for {catalog_id}: {e}')
    if log and last_error is not None:
        log.warning(
            f'All Vizier mirrors failed for {catalog_id}; last error: {last_error}')
    return None


def query_gaia_dr3_region(coord, width=20 * u.arcmin, log=None, max_rounds=4):
    """Query Gaia DR3 in a region (VizieR ``I/355/gaiadr3``), with retries.

    Used for astrometric alignment; photometry primitives can use
    :func:`query_vizier_region` with other catalog IDs.

    Parameters
    ----------
    coord : astropy.coordinates.SkyCoord
        Field center (ICRS).
    width : astropy.units.Quantity, optional
        Search box width. Default 20 arcmin.
    log : ColoredLogger, optional
        Logger.
    max_rounds : int, optional
        Number of full retry rounds if no table is returned (timeouts, etc.).

    Returns
    -------
    astropy.table.Table
        Unfiltered source table from VizieR.

    Raises
    ------
    Exception
        If no usable catalog response is obtained after all rounds.
    """
    for tries in range(max_rounds):
        tab = query_vizier_region(
            coord, width, GAIA_DR3_VIZIER_ID, GAIA_DR3_ASTROMETRY_COLUMNS, log=log)
        if tab is not None:
            return tab
        if log:
            log.error(f'Gaia did not return catalog. Try #{tries + 1}')
    raise Exception('ERROR: could not get Gaia catalog')


# Fine WCS alignment (after astrometry.net): supported reference catalogs (CLI / API).
FINE_ALIGN_CATALOG_CHOICES = (
    'gaia', 'panstarrs', 'sdss', 'legacy', 'twomass', 'skymapper')


def normalize_fine_align_catalog(catalog):
    """Normalize user/catalog string to a key in :data:`FINE_ALIGN_CATALOG_CHOICES`.

    Parameters
    ----------
    catalog : str
        e.g. ``'gaia'``, ``'2mass'``, ``'PS1'``.

    Returns
    -------
    str
        One of ``FINE_ALIGN_CATALOG_CHOICES``.

    Raises
    ------
    ValueError
        If *catalog* is not recognized.
    """
    if catalog is None:
        return 'gaia'
    n = str(catalog).strip().lower()
    aliases = {
        'ps1': 'panstarrs',
        'pan-starrs': 'panstarrs',
        'panstarrs1': 'panstarrs',
        '2mass': 'twomass',
        '2masspsc': 'twomass',
        'sdssdr12': 'sdss',
    }
    n = aliases.get(n, n)
    if n not in FINE_ALIGN_CATALOG_CHOICES:
        raise ValueError(
            f'Unknown fine-alignment catalog {catalog!r}; '
            f'expected one of {FINE_ALIGN_CATALOG_CHOICES}')
    return n


def _table_with_icrs_radec(tab, ra_col, dec_col):
    """Return *tab* with ``RA_ICRS`` and ``DE_ICRS`` float columns (degrees)."""
    out = Table(tab)
    out['RA_ICRS'] = np.asarray(out[ra_col], dtype=float)
    out['DE_ICRS'] = np.asarray(out[dec_col], dtype=float)
    return out


def _fetch_sdss_dr12_v147(coord, field_width, log=None):
    """SDSS DR12 photometry on VizieR (V/147); shared by ``sdss`` and ``legacy`` keys."""
    tab = query_vizier_region(
        coord, field_width, 'V/147',
        ['RA_ICRS', 'DE_ICRS', 'gmag'], log=log)
    if tab is None or len(tab) == 0:
        return None
    g = np.asarray(tab['gmag'], dtype=float)
    mask = np.isfinite(g) & (g < 22.0)
    tab = tab[mask]
    if len(tab) == 0:
        return None
    tab['RA_ICRS'] = np.asarray(tab['RA_ICRS'], dtype=float)
    tab['DE_ICRS'] = np.asarray(tab['DE_ICRS'], dtype=float)
    return tab


def fetch_astrometry_reference_table(coord, catalog, field_width, log=None):
    """Query a reference catalog for fine WCS alignment; return ICRS positions.

    All rows include ``RA_ICRS`` and ``DE_ICRS`` in degrees. Catalog-specific
    quality cuts reduce crowding for dense surveys.

    Parameters
    ----------
    coord : astropy.coordinates.SkyCoord
        Field center (ICRS).
    catalog : str
        One of :data:`FINE_ALIGN_CATALOG_CHOICES` (or accepted aliases).
    field_width : astropy.units.Quantity
        VizieR box width (e.g. ``0.5 * u.deg``).
    log : ColoredLogger, optional
        Logger.

    Returns
    -------
    astropy.table.Table or None
        Table with at least ``RA_ICRS``, ``DE_ICRS``, or None if the query
        failed or no sources remain after cuts.
    """
    key = normalize_fine_align_catalog(catalog)

    if key == 'gaia':
        try:
            w = field_width.to(u.arcmin)
            if w < 20 * u.arcmin:
                w = 20 * u.arcmin
            tab = query_gaia_dr3_region(coord, width=w, log=log)
        except Exception:
            return None
        if tab is None or len(tab) == 0:
            return None
        mask = (tab['PSS'] > 0.99) & (tab['Plx'] < 20) & (tab['PM'] < 10)
        tab = tab[mask]
        if len(tab) == 0:
            return None
        tab['RA_ICRS'] = np.asarray(tab['RA_ICRS'], dtype=float)
        tab['DE_ICRS'] = np.asarray(tab['DE_ICRS'], dtype=float)
        return tab

    if key == 'panstarrs':
        tab = query_vizier_region(
            coord, field_width, 'II/349',
            ['RAJ2000', 'DEJ2000', 'gmag'], log=log)
        if tab is None or len(tab) == 0:
            return None
        g = np.asarray(tab['gmag'], dtype=float)
        mask = np.isfinite(g) & (g < 21.5)
        tab = tab[mask]
        if len(tab) == 0:
            return None
        return _table_with_icrs_radec(tab, 'RAJ2000', 'DEJ2000')

    if key == 'sdss':
        return _fetch_sdss_dr12_v147(coord, field_width, log=log)

    # ``legacy`` is kept as a separate fine-align option from ``sdss``; both
    # currently use the same VizieR table (V/147) until ``legacy`` is retargeted.
    if key == 'legacy':
        return _fetch_sdss_dr12_v147(coord, field_width, log=log)

    if key == 'twomass':
        tab = query_vizier_region(
            coord, field_width, 'II/246',
            ['RAJ2000', 'DEJ2000', 'Jmag'], log=log)
        if tab is None or len(tab) == 0:
            return None
        j = np.asarray(tab['Jmag'], dtype=float)
        mask = np.isfinite(j) & (j > 7.0) & (j < 16.5)
        tab = tab[mask]
        if len(tab) == 0:
            return None
        return _table_with_icrs_radec(tab, 'RAJ2000', 'DEJ2000')

    if key == 'skymapper':
        tab = query_vizier_region(
            coord, field_width, 'II/379/smssdr4',
            ['RAICRS', 'DEICRS', 'gPSF'], log=log)
        if tab is None or len(tab) == 0:
            return None
        g = np.asarray(tab['gPSF'], dtype=float)
        mask = np.isfinite(g) & (g < 21.0)
        tab = tab[mask]
        if len(tab) == 0:
            return None
        return _table_with_icrs_radec(tab, 'RAICRS', 'DEICRS')

    return None
