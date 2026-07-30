"""Absolute photometry zeropoint calibration using catalog magnitudes.

Queries VizieR (PS1, SDSS, 2MASS, UKIRT, SkyMapper, DES) or DECaLS / Legacy
Surveys Tractor PSF photometry via NOIRLab Data Lab TAP
(:mod:`potpyri.utils.catalogs`), applies point-source cuts, matches sources,
and fits zeropoint via iterative ODR. Writes ZPTMAG and related keywords to
the stack header.

All zeropoints and limiting magnitudes are reported in the **AB** system.
Native-Vega catalogs (notably 2MASS) are converted to AB before the fit; see
:data:`TWOMASS_VEGA_TO_AB` and header keyword ``MAGSYS``.

Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from potpyri._version import __version__

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits

# Internal dependencies
from potpyri.utils import catalogs

# 2MASS PSC is published on the Vega system. Offsets below are
# (ABmag - Vegamag) = (mag_zero_AB - mag_zero_Vega) from the Tokunaga & Vacca
# (2005) / Bessell style zeropoints used historically in this pipeline:
#   J: 4.56 - 3.65 = 1.91
#   H: 4.71 - 3.32 = 1.39
#   K/Ks: 5.14 - 3.29 = 1.85
TWOMASS_VEGA_TO_AB = {
    'J': 4.56 - 3.65,
    'H': 4.71 - 3.32,
    'K': 5.14 - 3.29,
    'Ks': 5.14 - 3.29,
}

# FITS MAGSYS value written with ZPTMAG / M*SIGMA (always AB after conversion).
MAGSYS_AB = 'AB'


def _fits_extension_names(hdulist):
    """Return EXTNAME values for an HDUList (empty string if unset)."""
    return [h.name if h.name else '' for h in hdulist]


def _log_or_print(msg, log, level='info'):
    """Emit *msg* via *log* at *level*, or print if *log* is None."""
    if log is None:
        print(msg, flush=True)
        return
    getattr(log, level)(msg)


def catalog_native_magsys(catalog):
    """Return the native magnitude system of a photometric reference catalog.

    Parameters
    ----------
    catalog : str
        Catalog name (e.g. ``'PS1'``, ``'2MASS'``).

    Returns
    -------
    str
        ``'VEGA'`` for 2MASS; ``'AB'`` for PS1 / SDSS / SkyMapper and others
        treated as AB in this pipeline.
    """
    if str(catalog).upper() in ('2MASS', 'TWOMASS'):
        return 'VEGA'
    return 'AB'


def apply_catalog_to_ab(cat, catalog, filt, log=None):
    """Convert catalog magnitudes in-place to AB if they are native Vega.

    Parameters
    ----------
    cat : astropy.table.Table
        Table with a ``mag`` column (modified in place).
    catalog : str
        Catalog name.
    filt : str
        Catalog filter (``J``, ``H``, ``K``, ``Ks``, ...).
    log : ColoredLogger, optional
        Logger.

    Returns
    -------
    str
        Magnitude system of the returned photometry (``'AB'``).
    """
    native = catalog_native_magsys(catalog)
    if native == 'VEGA' and filt in TWOMASS_VEGA_TO_AB:
        offset = TWOMASS_VEGA_TO_AB[filt]
        cat['mag'] = cat['mag'] + offset
        _log_or_print(
            f'Converted {catalog} {filt} magnitudes Vega→AB '
            f'(+{offset:.2f} mag); MAGSYS will be {MAGSYS_AB!r}',
            log,
        )
    elif native == 'VEGA':
        _log_or_print(
            f'WARNING: {catalog} is Vega-native but no AB offset is defined '
            f'for filter {filt!r}; magnitudes left unchanged',
            log,
            level='warning',
        )
    return MAGSYS_AB


class absphot(object):
    """Zeropoint fitter using catalog magnitudes and iterative sigma clipping."""

    def __init__(self, iterations=5, sigma=5):
        """Initialize zeropoint fitter.

        Parameters
        ----------
        iterations : int, optional
            Number of sigma-clip iterations. Default is 5.
        sigma : float, optional
            Sigma threshold for clipping. Default is 5.
        """
        self.iterations = iterations
        self.sigma = sigma

        self.odr_iter = 2

    def get_zeropoint(self, flux, fluxerr, mag, magerr):
        """Fit zeropoint (and error) from flux/mag arrays using ODR.

        Parameters
        ----------
        flux : array-like
            Object fluxes (e.g. from aperture photometry).
        fluxerr : array-like
            Flux errors.
        mag : array-like
            Catalog magnitudes.
        magerr : array-like
            Catalog magnitude errors.

        Returns
        -------
        tuple of float
            (zeropoint, zeropoint_error).
        """
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message=r'.*scipy\.odr.*', category=DeprecationWarning)
            from scipy.odr import ODR
            from scipy.odr import Model
            from scipy.odr import RealData
        
        def magnitude(zpt, flux):
            return(zpt - 2.5*np.log10(flux))

        zpt_guess = np.nanmedian(mag + 2.5*np.log10(flux))

        data = RealData(np.array(flux, dtype=float),
                        np.array(mag, dtype=float), 
                        sx=np.array(fluxerr, dtype=float),
                        sy=np.array(magerr, dtype=float))

        model = Model(magnitude)
        odr = ODR(data, model, beta0=[zpt_guess], maxit=self.odr_iter)

        output = odr.run()

        zpt = output.beta[0]
        zpterr = output.sd_beta[0]

        return(zpt, zpterr)

    def zpt_iteration(self, flux, fluxerr, mag, magerr, log=None):
        """Iteratively sigma-clip and fit zeropoint via ODR.

        Parameters
        ----------
        flux : array-like
            Object fluxes.
        fluxerr : array-like
            Flux errors.
        mag : array-like
            Catalog magnitudes.
        magerr : array-like
            Catalog magnitude errors.
        log : ColoredLogger, optional
            Logger for progress.

        Returns
        -------
        tuple
            (zeropoint, zeropoint_error, master_mask). master_mask is boolean
            array of sources used in final fit.
        """
        flux = np.array(flux)
        fluxerr = np.array(fluxerr)
        mag = np.array(mag)
        magerr = np.array(magerr)

        nobj_orig = len(mag)
        idx = np.arange(nobj_orig)

        for i in np.arange(self.iterations):

            nobj = len(flux)
            zpt, zpterr = self.get_zeropoint(flux, fluxerr, mag, magerr)

            mag_deriv = -2.5*np.log10(flux)+zpt
            magerr_deriv = np.array(2.5/np.log(10) * fluxerr/flux)

            total_err = np.sqrt(magerr**2+magerr_deriv**2+zpterr**2)

            mask = np.abs(mag - mag_deriv) < self.sigma * total_err
            idx = idx[mask]

            flux=flux[mask] ; fluxerr=fluxerr[mask]
            mag=mag[mask] ; magerr=magerr[mask]

            zpt = float('%.6f'%zpt)
            zpterr = float('%.6f'%zpterr)

            message = f'Iteration {i+1}/{self.iterations}: {nobj} obj, '
            message += f'zpt={zpt}+/-{zpterr}, '
            message += f'{self.sigma}-sigma clip to {len(flux)} obj'
            if log: 
                log.info(message)
            else:
                print(message)
        
        if log:
            log.info(f'{len(flux)} stars used to calculate final zeropoint.')
        else:
            print(f'{len(flux)} stars used to calculate final zeropoint.')
        
        master_mask = np.array([i in idx for i in np.arange(nobj_orig)])

        return(zpt, zpterr, master_mask)

    def get_catalog(self, coords, catalog, filt, log=None):
        """Query a reference catalog for magnitudes around given coordinates.

        Supports VizieR catalogs (PS1, SDSS, 2MASS, UKIRT, SkyMapper, DES) and
        DECaLS / Legacy Surveys Tractor PSF photometry via NOIRLab Data Lab TAP.
        Applies catalog-specific point-source cuts and returns a uniform product
        with columns ``ra``, ``dec``, ``mag``, ``mag_err``.

        Parameters
        ----------
        coords : astropy.coordinates.SkyCoord
            Target coordinates (used for region size).
        catalog : str
            Catalog name (e.g. 'PS1', '2MASS', 'DECALS').
        filt : str
            Filter name (e.g. 'g', 'r', 'J').
        log : ColoredLogger, optional
            Logger for progress.

        Returns
        -------
        tuple or (None, None, None)
            (astropy.table.Table, catalog_name, catalog_ID) or (None, None, None)
            if query fails.
        """
        if log: log.info(f'Searching for catalog {catalog}')

        coord_ra = np.median([c.ra.degree for c in coords])
        coord_dec = np.median([c.dec.degree for c in coords])

        catalog, cat_ID, cat_ra, cat_dec, cat_mag, cat_err = catalogs.find_catalog(
            catalog, filt, coord_ra, coord_dec)

        if cat_ID is None:
            m = (f'ERROR: catalog {catalog!r} does not support filter {filt!r}')
            if log:
                log.error(m)
            else:
                print(m)
            return (None, None, None)

        med_coord = SkyCoord(coord_ra, coord_dec, unit='deg')
        seps = med_coord.separation(coords)
        max_sep = np.max(seps.to(u.deg).value)
        width = np.max([2.0 * max_sep, 0.5]) * u.degree

        # DECaLS / Legacy Surveys: Data Lab TAP (already PSF-selected).
        if catalog == 'DECALS' or (
                cat_ID in (catalogs.DECALS_TRACTOR_TABLE,
                           catalogs.DECALS_TRACTOR_FALLBACK_TABLE)):
            if log:
                log.info(
                    f'Getting DECaLS Tractor catalog ({cat_ID}) in filt {filt}'
                )
                log.info(f'Querying around {coord_ra}, {coord_dec} deg')
            cat = catalogs.query_decals_region(
                med_coord, width, filt, log=log, table=cat_ID)
            if cat is None or len(cat) == 0:
                m = ('ERROR: cat {0}, ra {1}, dec {2} did not return a catalog'
                     ).format(catalog, coord_ra, coord_dec)
                if log:
                    log.error(m)
                else:
                    print(m)
                return (None, None, None)
            # Product already has ra/dec/mag/mag_err; re-apply cut for consistency.
            cat = catalogs.apply_point_source_cut(
                cat, cat_ID, filt, mag_col='mag', log=log)
            if cat is None or len(cat) == 0:
                return (None, None, None)
            # Ensure exact product columns.
            for required in ('ra', 'dec', 'mag', 'mag_err'):
                if required not in cat.colnames:
                    return (None, None, None)
            return (cat['ra', 'dec', 'mag', 'mag_err'], catalog, cat_ID)

        cols = [cat_ra, cat_dec, cat_mag, cat_err]
        cols.extend(catalogs.point_source_extra_columns(cat_ID, filt))
        # Preserve order while dropping duplicates.
        seen = set()
        cols = [c for c in cols if not (c in seen or seen.add(c))]

        if log:
            log.info(f'Getting {catalog} catalog with ID {cat_ID} in filt {filt}')
            log.info(f'Querying around {coord_ra}, {coord_dec} deg')
        cat = catalogs.query_vizier_region(
            med_coord, width, cat_ID, cols, log=log)

        if cat is not None:
            cat = cat[~np.isnan(cat[cat_mag])]
            cat = cat[cat[cat_err] > 0.]

            cat = catalogs.apply_point_source_cut(
                cat, cat_ID, filt, mag_col=cat_mag, log=log)
            if cat is None or len(cat) == 0:
                m = (f'ERROR: no point sources remain after cut for {catalog}')
                if log:
                    log.error(m)
                else:
                    print(m)
                return (None, None, None)

            cat.rename_column(cat_ra, 'ra')
            cat.rename_column(cat_dec, 'dec')
            cat.rename_column(cat_mag, 'mag')
            cat.rename_column(cat_err, 'mag_err')

            # Convert native-Vega catalogs (2MASS) to AB; PS1/SDSS already AB.
            apply_catalog_to_ab(cat, catalog, filt, log=log)

            if catalog == '2MASS' and filt == 'Y':
                cat = cat[~np.isnan(cat['Kmag'])]
                cat['mag'], cat['mag_err'] = self.Y_band(
                    cat['mag'], cat['mag_err'], cat['Kmag'], cat['e_Kmag'])

            return (cat, catalog, cat_ID)

        else:
            m = 'ERROR: cat {0}, ra {1}, dec {2} did not return a catalog'
            m = m.format(catalog, coord_ra, coord_dec)
            if log:
                log.error(m)
            else:
                print(m)
            return (None, None, None)

    def find_zeropoint(self, cmpfile, tel, match_radius=2.5*u.arcsec,
        phottable='APPPHOT', input_catalog=None, log=None):
        """Compute zeropoint from cmpfile photometry and catalog; write to FITS header.

        Matches sources to catalog (e.g. PS1, 2MASS), runs iterative ODR fit, and
        updates ZPTMAG, ZPTNSTAR, ZPTCAT, MAGSYS, etc. in the stack FITS.
        Magnitudes are always stored in the AB system (``MAGSYS='AB'``); 2MASS
        Vega values are converted before the fit.

        Parameters
        ----------
        cmpfile : str
            Path to stacked/comparison FITS with SCI and phottable extensions.
        tel : Instrument
            Instrument instance (for get_catalog).
        match_radius : astropy.units.Quantity, optional
            Matching radius for catalog. Default is 2.5 arcsec.
        phottable : str, optional
            FITS extension with photometry table. Default is 'APPPHOT'.
        input_catalog : astropy.table.Table, optional
            Pre-loaded catalog; if None, catalog is queried via get_catalog.
        log : ColoredLogger, optional
            Logger for progress.

        Returns
        -------
        bool
            True if ZPT keywords were written; False if photometry or catalog step
            did not complete (cmpfile may be unchanged).
        """
        _log_or_print(f'Zeropoint: reading stack {cmpfile!r}', log)

        with fits.open(cmpfile) as hdu:
            extnames = _fits_extension_names(hdu)
            if phottable not in hdu:
                _log_or_print(
                    f'Zeropoint aborted: FITS extension {phottable!r} not found in '
                    f'{cmpfile!r}. Present extensions (EXTNAME): {extnames!r}. '
                    'Run photometry first (photometry.photloop / pipeline photometry '
                    'step) so APPPHOT is written; if photometry already ran, check logs '
                    'for PhotometryError or earlier tracebacks.',
                    log,
                    level='error',
                )
                return False

            header = hdu['SCI'].header
            filtorig = header['FILTER']
            catalog = tel.get_catalog(header)

            table = Table(hdu[phottable].data, meta=hdu[phottable].header)
            required = ('RA', 'Dec', 'flux', 'flux_err')
            missing = [c for c in required if c not in table.colnames]
            if missing:
                _log_or_print(
                    f'Zeropoint aborted: extension {phottable!r} is missing columns '
                    f'{missing!r}; have {list(table.colnames)!r}. Photometry output '
                    'may be corrupt or from an incompatible version.',
                    log,
                    level='error',
                )
                return False

            coords = SkyCoord(table['RA'], table['Dec'], unit='deg')

            # New metadata to update
            metadata = {}

            cat = None
            cat_ID = None

            filt = self.convert_filter_name(filtorig)
            if input_catalog is not None:
                cat = input_catalog
            else:
                _log_or_print(
                    f'Zeropoint: querying {catalog} catalog for filter {filt}', log)

                cat, catalog, cat_ID = self.get_catalog(coords, catalog, filt, log=log)

            if cat is not None and len(cat) > 0:
                min_mag = self.get_minmag(filt)
                cat = cat[cat['mag'] > min_mag]
                if len(cat) == 0:
                    _log_or_print(
                        f'Zeropoint not written: no catalog sources fainter than '
                        f'bright limit min_mag={min_mag} for filter {filt!r}.',
                        log,
                        level='warning',
                    )
                    return False

                coords_cat = SkyCoord(cat['ra'], cat['dec'], unit='deg')

                idx, d2, d3 = coords_cat.match_to_catalog_sky(coords)

                # Get matches from calibration catalog and cmpfile
                mask = d2 < match_radius
                cat = cat[mask]
                idx = idx[mask]

                if len(cat) == 0:
                    _log_or_print(
                        'Zeropoint not written: no catalog stars match photometry '
                        f'within {match_radius} (try a larger match radius or check '
                        'WCS/field).',
                        log,
                        level='warning',
                    )
                    return False

                match_table = None
                for i in idx:
                    if not match_table:
                        match_table = Table(table[i])
                    else:
                        match_table.add_row(table[i])

                # Sort by flux
                flux_idx = np.argsort(match_table['flux'])
                match_table = match_table[flux_idx]
                cat = cat[flux_idx]

                # Do basic cuts on flux, fluxerr, catalog magnitude, cat magerr
                flux = match_table['flux'].data.astype('float32')
                fluxerr = match_table['flux_err'].data.astype('float32')
                cat_mag = cat['mag'].data.astype('float32')
                cat_magerr = cat['mag_err'].data.astype('float32')

                if len(flux) == 0:
                    _log_or_print(
                        'Zeropoint not written: no rows left after flux/magnitude cuts.',
                        log,
                        level='warning',
                    )
                    return False

                zpt, zpterr, master_mask = self.zpt_iteration(flux, fluxerr,
                    cat_mag, cat_magerr, log=log)

                # Set header variables
                metadata['ZPTNSTAR'] = len(flux)
                metadata['ZPTMAG'] = zpt
                metadata['ZPTMUCER'] = zpterr
                metadata['ZPTCAT'] = catalog
                metadata['ZPTCATID'] = cat_ID
                metadata['ZPTPHOT'] = phottable
                metadata['FILTER'] = filtorig
                # Zeropoint and limiting mags are always on the AB system
                # (2MASS Vega→AB applied in get_catalog / apply_catalog_to_ab).
                metadata['MAGSYS'] = (
                    MAGSYS_AB,
                    'Magnitude system for ZPTMAG and M*SIGMA (AB)',
                )

                # Add limiting magnitudes
                if 'FWHM' in header.keys() and 'SKYSIG' in header.keys():
                    fwhm = header['FWHM']
                    sky = header['SKYSIG']
                    Npix_per_FWHM_Area = 2.5 * 2.5 * fwhm * fwhm
                    skysig_per_FWHM_Area = np.sqrt(Npix_per_FWHM_Area * (sky * sky))
                    m3sigma = -2.5 * np.log10(3.0 * skysig_per_FWHM_Area) + zpt
                    m5sigma = -2.5 * np.log10(5.0 * skysig_per_FWHM_Area) + zpt
                    m10sigma = -2.5 * np.log10(10.0 * skysig_per_FWHM_Area) + zpt
                    m3sigma = float('%.6f' % m3sigma)
                    m5sigma = float('%.6f' % m5sigma)
                    m10sigma = float('%.6f' % m10sigma)
                    metadata['M3SIGMA'] = m3sigma
                    metadata['M5SIGMA'] = m5sigma
                    metadata['M10SIGMA'] = m10sigma
                    _log_or_print(
                        f'3-sigma limiting mag of image is {m3sigma} ({MAGSYS_AB})',
                        log,
                    )

                hdu['PRIMARY'].header.update(metadata)
                hdu['SCI'].header.update(metadata)
                hdu[phottable].header.update(metadata)

                hdu.writeto(cmpfile, overwrite=True)
                _log_or_print(
                    f'Zeropoint written to {cmpfile!r} (ZPTMAG={zpt:.4f}, N={len(flux)})',
                    log,
                )
                return True

            if input_catalog is not None:
                _log_or_print(
                    'Zeropoint not written: input_catalog is missing or has zero rows.',
                    log,
                    level='warning',
                )
            else:
                _log_or_print(
                    'Zeropoint not written: catalog query returned no usable table '
                    f'(catalog_name={catalog!r}, filter={filt!r}).',
                    log,
                    level='warning',
                )
            return False

    def Y_band(self, J, J_err, K, K_err):
        """Compute Y-band mag and error from J and K (2MASS relation).

        Parameters
        ----------
        J, J_err : array-like or float
            J magnitude(s) and error(s).
        K, K_err : array-like or float
            K magnitude(s) and error(s).

        Returns
        -------
        tuple
            (Y_mag, Y_err).
        """
        Y = J+0.46*(J-K)
        JK_err = np.sqrt(J_err**2+K_err**2)
        JKY_err = 0.46*(J-K)*np.sqrt((0.02/0.46)**2+(JK_err/(J-K))**2)
        Y_err = np.sqrt(J_err**2+JKY_err**2)
        return Y, Y_err

    def convert_filter_name(self, filt):
        """Map instrument filter names to catalog filter names (e.g. PS1).

        Parameters
        ----------
        filt : str
            Instrument filter keyword (e.g. 'gG0301', 'RG850').

        Returns
        -------
        str
            Catalog filter name ('u', 'g', 'r', 'i', 'z', 'J', etc.).
        """
        if filt=='uG0308' or filt=='uG0332' or filt=='U':
            return 'u'
        if filt=='gG0301' or filt=='gG0325' or filt=='G' or filt=='V' or filt=='B':
            return 'g'
        if filt=='rG0303' or filt=='rG0326' or filt=='R' or filt=='Rs':
            return 'r'
        if filt=='iG0302' or filt=='iG0327' or filt=='I':
            return 'i'
        if filt=='zG0304' or filt=='zG0328' or filt=='Z':
            return 'z'
        if filt=='RG850':
            return 'z'
        if filt=='Y':
            return 'J'
        # K, Ks, Kspec all use 2MASS K-band for calibration
        if filt in ('K', 'Ks', 'Kspec'):
            return 'K'
        else:
            return filt
    
    def get_minmag(self, filt):
        """Return minimum catalog magnitude to use for zeropoint (bright limit).

        Parameters
        ----------
        filt : str
            Filter name.

        Returns
        -------
        float
            Minimum magnitude (brighter limit).
        """
        if filt=='J':
            return 15.5
        if filt=='K':
            return 13.0
        if filt=='Y':
            return 15.0
        else:
            return 16.0

def find_zeropoint(stack, tel, log=None):
    """Compute and write zeropoint to stack FITS using instrument catalog (e.g. PS1).

    Parameters
    ----------
    stack : str
        Path to stacked science FITS file (SCI + APPPHOT extensions).
    tel : Instrument
        Instrument instance (for catalog and filter).
    log : ColoredLogger, optional
        Logger for progress and errors.

    Returns
    -------
    bool
        True if ZPT keywords were written; False otherwise (see log).
    """
    cal = absphot()
    return cal.find_zeropoint(stack, tel, log=log)
