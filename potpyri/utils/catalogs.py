"""VizieR and catalog metadata for astrometry and photometric calibration.

Photometric reference catalogs (PS1, 2MASS, SkyMapper, SDSS, DES, DECaLS, etc.)
and astrometric queries (e.g. Gaia DR3) live here so primitives such as
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
from astropy.coordinates import SkyCoord
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
    'decals': {
        'description': (
            'DECaLS / DESI Legacy Surveys Tractor PSF photometry '
            '(NOIRLab Data Lab TAP: ls_dr10.tractor)'
        ),
        'vizier_id': None,
        'datalab_table': 'ls_dr10.tractor',
        'ra_column': 'ra',
        'dec_column': 'dec',
        'default_vizier_columns': [
            'ra', 'dec', 'type', 'dered_mag_g', 'dered_mag_r', 'dered_mag_i',
            'dered_mag_z', 'flux_g', 'flux_r', 'flux_i', 'flux_z',
            'flux_ivar_g', 'flux_ivar_r', 'flux_ivar_i', 'flux_ivar_z',
        ],
        'roles': ('astrometry', 'photometry'),
        'notes': 'Queried via Data Lab TAP, not VizieR.',
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

# NOIRLab Astro Data Lab TAP endpoint for Legacy Surveys / DECaLS Tractor catalogs.
DATALAB_TAP_URL = 'https://datalab.noirlab.edu/tap'
DECALS_TRACTOR_TABLE = 'ls_dr10.tractor'
DECALS_TRACTOR_FALLBACK_TABLE = 'ls_dr9.tractor'

#: Catalogs used for photometric zeropoint calibration (instrument ``catalog_zp``).
FLUX_CALIBRATION_CATALOGS = (
    'PS1', 'SDSS', '2MASS', 'UKIRT', 'SKYMAPPER', 'DES', 'DECALS',
)

#: Point-source selection rules keyed by catalog ID returned from :func:`find_catalog`.
#: Each rule is applied by :func:`apply_point_source_cut`.
POINT_SOURCE_CUTS = {
    'II/349': {  # Pan-STARRS1
        'method': 'psf_kron',
        'kron_template': '{filt}Kmag',
        'threshold': 0.1,
        'description': 'PSF - Kron < 0.1 mag',
    },
    'V/154': {  # SDSS (find_catalog ID)
        'method': 'class_equals',
        'column': 'class',
        'column_aliases': ('class', 'cl', 'Class'),
        'value': 6,
        'description': 'SDSS photometric class == 6 (star)',
    },
    'V/147': {  # SDSS DR12 (legacy / fine-align ID)
        'method': 'class_equals',
        'column': 'class',
        'column_aliases': ('class', 'cl', 'Class'),
        'value': 6,
        'description': 'SDSS photometric class == 6 (star)',
    },
    'II/246': {  # 2MASS PSC
        'method': 'class_equals',
        'column': 'Xflg',
        'column_aliases': ('Xflg', 'X', 'ext_key'),
        'value': 0,
        'description': '2MASS extended-source flag Xflg == 0',
    },
    'II/319': {  # UKIRT / UKIDSS
        'method': 'class_equals',
        'column': 'mergedClass',
        'column_aliases': ('mergedClass', 'mergedclass', 'Class'),
        'value': -1,
        'description': 'UKIDSS mergedClass == -1 (star)',
    },
    'II/379/smssdr4': {  # SkyMapper DR4
        'method': 'score_above',
        'column': 'ClassStar',
        'column_aliases': ('ClassStar', 'class_star', 'Class_Star'),
        'threshold': 0.9,
        'description': 'SkyMapper ClassStar > 0.9',
    },
    'II/357': {  # DES DR1
        'method': 'des_sg',
        'column_template': 'S/G{filt}',
        'threshold': 0.5,
        'description': 'DES S/G flag > 0.5 (star-like)',
    },
    'ls_dr10.tractor': {  # DECaLS / Legacy Surveys
        'method': 'type_equals',
        'column': 'type',
        'value': 'PSF',
        'description': "Tractor morphological type == 'PSF'",
    },
    'ls_dr9.tractor': {
        'method': 'type_equals',
        'column': 'type',
        'value': 'PSF',
        'description': "Tractor morphological type == 'PSF'",
    },
}


def normalize_flux_catalog_name(catalog):
    """Normalize a flux-calibration catalog name to a canonical key."""
    if catalog is None:
        return None
    key = str(catalog).strip().upper()
    aliases = {
        'PANSTARRS': 'PS1',
        'PAN-STARRS': 'PS1',
        'PANSTARRS1': 'PS1',
        'PS1DR1': 'PS1',
        'TWOMASS': '2MASS',
        '2MASSPSC': '2MASS',
        'SMSS': 'SKYMAPPER',
        'SMSSDR4': 'SKYMAPPER',
        'DECALS_DR10': 'DECALS',
        'DECALS_DR9': 'DECALS',
        'LEGACY': 'DECALS',
        'LEGACYSURVEY': 'DECALS',
        'LEGACY_SURVEY': 'DECALS',
        'DESDR1': 'DES',
    }
    return aliases.get(key, key)


def find_catalog(catalog, fil, coord_ra, coord_dec):
    """Return catalog ID and column names for the given catalog and filter.

    Supports PS1, SDSS, 2MASS, UKIRT, SKYMAPPER, DES, and DECALS. For southern
    u-band, uses SkyMapper automatically when the requested catalog is PS1.

    Parameters
    ----------
    catalog : str
        Catalog name (e.g. 'PS1', 'SDSS', '2MASS', 'DECALS').
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
        Vizier / Data Lab queries. catalog_ID/ra/dec/mag/err may be None if
        the filter is not supported.
    """
    catalog_ID, ra, dec, mag, err = None, None, None, None, None
    catalog = normalize_flux_catalog_name(catalog)

    if coord_dec < 0 and fil.lower() == 'u' and catalog in ('PS1', 'DES', 'DECALS'):
        catalog = 'SKYMAPPER'

    if catalog == 'SDSS':
        if fil.lower() not in ['u', 'g', 'r', 'i', 'z']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'V/154', 'RA_ICRS', 'DE_ICRS',
            fil.lower() + 'mag', 'e_' + fil.lower() + 'mag')
    elif catalog == '2MASS':
        fil_2mass = fil.upper()
        if fil_2mass in ('KS', 'KSPEC'):
            fil_2mass = 'K'
        if fil_2mass not in ['J', 'H', 'K']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/246', 'RAJ2000', 'DEJ2000',
            fil_2mass + 'mag', 'e_' + fil_2mass + 'mag')
    elif catalog == 'UKIRT':
        if fil.upper() not in ['Y', 'J', 'H', 'K']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/319', 'ra', 'dec',
            fil.upper() + 'mag', 'e_' + fil.upper() + 'mag')
    elif catalog == 'PS1':
        if fil.lower() not in ['g', 'r', 'i', 'z', 'y']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/349', 'RAJ2000', 'DEJ2000',
            fil.lower() + 'mag', 'e_' + fil.lower() + 'mag')
    elif catalog == 'SKYMAPPER':
        if fil.lower() not in ['u', 'v', 'g', 'r', 'i', 'z']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        catalog_ID, ra, dec, mag, err = (
            'II/379/smssdr4', 'RAICRS', 'DEICRS',
            fil.lower() + 'PSF', 'e_' + fil.lower() + 'PSF')
    elif catalog == 'DES':
        if fil.lower() not in ['g', 'r', 'i', 'z', 'y']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        # DES DR1 on Vizier uses y as Y in some releases; prefer lowercase ymag.
        f = fil.lower()
        mag_col = f + 'mag' if f != 'y' else 'Ymag'
        err_col = 'e_' + mag_col if f != 'y' else 'e_Ymag'
        # Prefer standard lowercase columns used in viziercat for griz.
        if f in ('g', 'r', 'i', 'z'):
            mag_col = f + 'mag'
            err_col = 'e_' + f + 'mag'
        catalog_ID, ra, dec, mag, err = (
            'II/357', 'RAJ2000', 'DEJ2000', mag_col, err_col)
    elif catalog == 'DECALS':
        if fil.lower() not in ['g', 'r', 'i', 'z']:
            return (catalog, catalog_ID, ra, dec, mag, err)
        f = fil.lower()
        catalog_ID, ra, dec, mag, err = (
            DECALS_TRACTOR_TABLE, 'ra', 'dec',
            f'dered_mag_{f}', f'e_dered_mag_{f}')

    return (catalog, catalog_ID, ra, dec, mag, err)


def point_source_extra_columns(catalog_id, filt):
    """Return extra VizieR/TAP columns needed for the point-source cut."""
    rule = POINT_SOURCE_CUTS.get(catalog_id)
    if rule is None:
        return []
    method = rule['method']
    if method == 'psf_kron':
        return [rule['kron_template'].format(filt=filt.lower())]
    if method in ('class_equals', 'score_above', 'type_equals'):
        cols = [rule['column']]
        for alias in rule.get('column_aliases', ()):
            if alias not in cols:
                cols.append(alias)
        return cols
    if method == 'des_sg':
        f = filt.lower()
        # DES S/G columns are S/Gg, S/Gr, ...
        return [f'S/G{f}', f'S/G{f.upper()}']
    return []


def _resolve_cut_column(table, names):
    """Return the first column name in *names* present in *table*, else None."""
    for name in names:
        if name in table.colnames:
            return name
    return None


def apply_point_source_cut(table, catalog_id, filt, mag_col=None, log=None):
    """Apply the catalog-specific point-source selection to *table*.

    Parameters
    ----------
    table : astropy.table.Table
        Catalog rows (may include cut-helper columns).
    catalog_id : str
        Catalog ID from :func:`find_catalog` (e.g. ``'II/349'``, ``'ls_dr10.tractor'``).
    filt : str
        Filter used for PSF–Kron or DES S/G cuts.
    mag_col : str, optional
        Magnitude column for PSF–Kron (defaults to ``{filt}mag`` / catalog default).
    log : ColoredLogger, optional
        Logger.

    Returns
    -------
    astropy.table.Table
        Filtered table. If the cut columns are missing, returns *table* unchanged
        (with a warning) so network/schema drift does not hard-fail calibration.
    """
    if table is None or len(table) == 0:
        return table

    rule = POINT_SOURCE_CUTS.get(catalog_id)
    if rule is None:
        if log:
            log.warning(f'No point-source cut defined for catalog_id={catalog_id!r}')
        return table

    n0 = len(table)
    method = rule['method']
    filt_l = str(filt).lower()

    if method == 'psf_kron':
        kron_col = rule['kron_template'].format(filt=filt_l)
        if mag_col is None:
            mag_col = f'{filt_l}mag'
        if mag_col not in table.colnames or kron_col not in table.colnames:
            msg = (
                f'Point-source cut skipped for {catalog_id}: missing '
                f'{mag_col!r} or {kron_col!r}'
            )
            if log:
                log.warning(msg)
            else:
                print(msg)
            return table
        mag = np.asarray(table[mag_col], dtype=float)
        kron = np.asarray(table[kron_col], dtype=float)
        mask = np.isfinite(mag) & np.isfinite(kron) & (
            (mag - kron) < float(rule['threshold']))
        out = table[mask]
    elif method == 'class_equals':
        col = _resolve_cut_column(
            table, (rule['column'],) + tuple(rule.get('column_aliases', ())))
        if col is None:
            msg = (
                f'Point-source cut skipped for {catalog_id}: missing '
                f'{rule["column"]!r}'
            )
            if log:
                log.warning(msg)
            else:
                print(msg)
            return table
        vals = np.asarray(table[col])
        # Allow numeric or string-encoded integers (e.g. Vizier).
        try:
            vals_num = np.asarray(vals, dtype=float)
            mask = vals_num == float(rule['value'])
        except (TypeError, ValueError):
            mask = np.asarray([str(v).strip() == str(rule['value']) for v in vals])
        out = table[mask]
    elif method == 'score_above':
        col = _resolve_cut_column(
            table, (rule['column'],) + tuple(rule.get('column_aliases', ())))
        if col is None:
            msg = (
                f'Point-source cut skipped for {catalog_id}: missing '
                f'{rule["column"]!r}'
            )
            if log:
                log.warning(msg)
            else:
                print(msg)
            return table
        score = np.asarray(table[col], dtype=float)
        mask = np.isfinite(score) & (score > float(rule['threshold']))
        out = table[mask]
    elif method == 'type_equals':
        col = rule['column']
        if col not in table.colnames:
            msg = f'Point-source cut skipped for {catalog_id}: missing {col!r}'
            if log:
                log.warning(msg)
            else:
                print(msg)
            return table
        types = np.asarray([str(v).strip().upper() for v in table[col]])
        mask = types == str(rule['value']).strip().upper()
        out = table[mask]
    elif method == 'des_sg':
        candidates = [
            f'S/G{filt_l}',
            f'S/G{filt_l.upper()}',
            rule.get('column_template', 'S/G{filt}').format(filt=filt_l),
        ]
        col = _resolve_cut_column(table, candidates)
        if col is None:
            msg = f'Point-source cut skipped for {catalog_id}: missing DES S/G column'
            if log:
                log.warning(msg)
            else:
                print(msg)
            return table
        sg = np.asarray(table[col], dtype=float)
        mask = np.isfinite(sg) & (sg > float(rule['threshold']))
        out = table[mask]
    else:
        if log:
            log.warning(f'Unknown point-source cut method {method!r}')
        return table

    msg = (
        f'Point-source cut ({rule["description"]}): {n0} -> {len(out)} sources'
    )
    if log:
        log.info(msg)
    else:
        print(msg)
    return out


def nanomaggy_to_ab_mag(flux, flux_ivar=None):
    """Convert Legacy Surveys nanomaggy flux (+ optional ivar) to AB mag / magerr."""
    flux = np.asarray(flux, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = 22.5 - 2.5 * np.log10(flux)
    magerr = np.full(flux.shape, np.nan, dtype=float)
    if flux_ivar is not None:
        flux_ivar = np.asarray(flux_ivar, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            snr = flux * np.sqrt(flux_ivar)
            magerr = (2.5 / np.log(10.0)) / snr
            magerr[~np.isfinite(snr) | (snr <= 0)] = np.nan
    mag[~np.isfinite(flux) | (flux <= 0)] = np.nan
    return mag, magerr


def query_decals_region(center, width, filt, log=None, table=None):
    """Query DECaLS / Legacy Surveys Tractor PSF sources via Data Lab TAP.

    Selects ``type='PSF'`` rows and returns a table with columns ``ra``, ``dec``,
    ``mag``, ``mag_err``, and ``type`` — the same product schema used by
    :meth:`potpyri.primitives.absphot.absphot.get_catalog` after renaming.

    Parameters
    ----------
    center : astropy.coordinates.SkyCoord
        Field center.
    width : astropy.units.Quantity
        Search box width (cone radius uses half-width).
    filt : str
        Optical band ``g``, ``r``, ``i``, or ``z``.
    log : ColoredLogger, optional
        Logger.
    table : str, optional
        Tractor TAP table (default :data:`DECALS_TRACTOR_TABLE`).

    Returns
    -------
    astropy.table.Table or None
        Point-source photometry table, or None on failure / empty result.
    """
    f = str(filt).lower()
    if f not in ('g', 'r', 'i', 'z'):
        raise ValueError(f'DECaLS filter must be g/r/i/z, got {filt!r}')

    if not isinstance(center, SkyCoord):
        center = SkyCoord(center)

    radius_deg = 0.5 * float(width.to(u.deg).value)
    radius_deg = max(radius_deg, 0.1)
    ra = float(center.ra.degree)
    dec = float(center.dec.degree)
    # Data Lab Tractor tables do not expose q3c_radial_query in all TAP setups;
    # use a conservative RA/Dec box (no cos(dec) shrink — slightly larger area).
    ra_lo, ra_hi = ra - radius_deg, ra + radius_deg
    dec_lo, dec_hi = dec - radius_deg, dec + radius_deg
    tables = [table or DECALS_TRACTOR_TABLE, DECALS_TRACTOR_FALLBACK_TABLE]
    # Deduplicate while preserving order.
    seen = set()
    tables = [t for t in tables if not (t in seen or seen.add(t))]

    try:
        from astroquery.utils.tap.core import TapPlus
    except Exception as e:
        if log:
            log.error(f'Data Lab TAP client unavailable: {e}')
        return None

    last_error = None
    for tap_table in tables:
        query = f"""
SELECT ra, dec, type,
       dered_mag_{f} AS mag,
       flux_{f}, flux_ivar_{f}
FROM {tap_table}
WHERE ra BETWEEN {ra_lo} AND {ra_hi}
  AND dec BETWEEN {dec_lo} AND {dec_hi}
  AND type = 'PSF'
  AND flux_{f} > 0
  AND flux_ivar_{f} > 0
  AND dered_mag_{f} > 0
  AND dered_mag_{f} < 30
""".strip()
        try:
            if log:
                log.info(
                    f'Querying DECaLS Tractor via Data Lab ({tap_table}), '
                    f'filter={f}, box half-width={radius_deg:.3f} deg'
                )
            tap = TapPlus(url=DATALAB_TAP_URL)
            # Default TAP maxrec is often 2000; raise it for typical ZP fields.
            job = tap.launch_job(query, maxrec=100000)
            result = job.get_results()
            if result is None or len(result) == 0:
                continue
            out = Table(result)
            # Prefer dereddened mag; recompute magerr from flux ivar.
            mag = np.asarray(out['mag'], dtype=float)
            _, magerr = nanomaggy_to_ab_mag(out['flux_' + f], out['flux_ivar_' + f])
            # If dered mag is missing/NaN, fall back to nanomaggy conversion.
            bad = ~np.isfinite(mag)
            if np.any(bad):
                mag_fb, _ = nanomaggy_to_ab_mag(
                    out['flux_' + f], out['flux_ivar_' + f])
                mag = np.where(bad, mag_fb, mag)
            out['mag'] = mag
            out['mag_err'] = magerr
            keep = (
                np.isfinite(out['mag']) & np.isfinite(out['mag_err'])
                & (out['mag_err'] > 0)
            )
            out = out[keep]
            # Drop intermediate flux columns for a stable product schema.
            for col in list(out.colnames):
                if col.startswith('flux'):
                    out.remove_column(col)
            if log:
                log.info(f'DECaLS query returned {len(out)} PSF sources')
            return out if len(out) else None
        except Exception as e:
            last_error = e
            if log:
                log.warning(f'DECaLS TAP query failed for {tap_table}: {e}')

    if log and last_error is not None:
        log.warning(f'All DECaLS TAP attempts failed; last error: {last_error}')
    return None


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
