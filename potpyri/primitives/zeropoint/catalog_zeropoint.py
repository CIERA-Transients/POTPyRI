"""Absolute photometry zeropoint calibration using catalog magnitudes.

Queries Vizier (e.g. PS1), matches sources, and fits zeropoint via
iterative ODR. Writes ZPTMAG, MAGSYS, and related keywords to the stack header.
Synthetic 2MASS-based :math:`Y` magnitudes follow González-Fernández et al.
(2018, MNRAS 474, 5459–5478; doi:10.1093/mnras/stx3073): equation (6) for
VISTA :math:`Y` from 2MASS :math:`J` and :math:`K_s`, and Appendix D for
VISTA-to-AB offsets when ``MAGSYS`` is AB.

Uses multiple Vizier mirrors (CDS, Tokyo, CfA, INASAN, IUCAA, IDIA) with
fallback when the first attempt fails.
Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from potpyri._version import __version__

import numpy as np

from astroquery.vizier import Vizier
from astropy import units as u

# Vizier mirror servers (hostnames only); tried in order when a query fails.
VIZIER_MIRRORS = [
    'vizier.cds.unistra.fr',   # CDS / Strasbourg, France
    'vizier.nao.ac.jp',        # ADAC / Tokyo, Japan
    'vizier.cfa.harvard.edu',  # CfA / Harvard, USA
    'vizier.inasan.ru',        # INASAN / Moscow, Russia
    'vizier.iucaa.in',         # IUCAA / Pune, India
    'vizier.idia.ac.za',       # IDIA / South Africa
]
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits

# Internal dependency
from potpyri.utils import utilities

# FITS header value for AB magnitudes (zeropoint fit uses AB catalog mags by default).
DEFAULT_MAGSYS = "ABMAG"

# -----------------------------------------------------------------------------
# 2MASS / VISTA Y-band: González-Fernández et al. (2018, MNRAS 474, 5459–5478;
# doi:10.1093/mnras/stx3073, hereafter GF18).  Equation (6) gives VISTA Y from
# 2MASS Vega J and (J−Ks).  Appendix D equation (D3) converts VISTA Y (Vega-like)
# to AB magnitudes (via Hewett et al. 2006 zero-points as summarized in GF18).
# -----------------------------------------------------------------------------
# C_Y in GF18 eq. (6); formal uncertainty GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR.
GF18_TWOMASS_TO_VISTA_Y_SLOPE = 0.46
GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR = 0.02
# Y(AB) − Y(VISTA) from GF18 Appendix D eq. (D3); see Hewett et al. 2006.
GF18_VISTA_Y_VEGA_TO_AB = 0.600


class ZeropointFitter(object):
    """Zeropoint fitter using catalog magnitudes and iterative sigma clipping."""

    def __init__(self, iterations=5, sigma=5, magsys=None):
        """Initialize zeropoint fitter.

        Parameters
        ----------
        iterations : int, optional
            Number of sigma-clip iterations. Default is 5.
        sigma : float, optional
            Sigma threshold for clipping. Default is 5.
        magsys : str, optional
            Magnitude system for ``ZPTMAG`` and catalog mags used in the fit,
            recorded as FITS keyword ``MAGSYS``. Default is :data:`DEFAULT_MAGSYS`
            (``\"ABMAG\"``). Non-AB systems may be supported in future by
            extending :meth:`apply_catalog_magnitude_system`.
        """
        self.iterations = iterations
        self.sigma = sigma
        self.magsys = magsys if magsys is not None else DEFAULT_MAGSYS

        self.odr_iter = 2

    def get_magsys(self):
        """Return the magnitude system name written to ``MAGSYS`` and used for the fit.

        Returns
        -------
        str
            Magnitude system identifier (e.g. ``DEFAULT_MAGSYS``).
        """
        return self.magsys

    def set_magsys(self, magsys):
        """Set the magnitude system for subsequent catalog conversion and headers.

        Parameters
        ----------
        magsys : str
            Identifier stored in ``MAGSYS`` (e.g. ``ABMAG``).
        """
        self.magsys = magsys

    def apply_catalog_magnitude_system(self, cat, catalog, filt, log=None):
        """Convert catalog magnitudes into the system used for zeropoint fitting.

        For ``MAGSYS`` = :data:`DEFAULT_MAGSYS`, Vega-native catalogs (e.g. 2MASS)
        are converted to AB using the standard offsets; see
        :meth:`twomass_vega_to_ab`.  Synthetic :math:`Y` from 2MASS uses
        :meth:`twomass_y_mag_to_ab` (GF18 eq. 6 and Appendix D eq. D3).

        Parameters
        ----------
        cat : astropy.table.Table
            Table with columns ``mag`` and ``mag_err`` (and 2MASS Y-band extras
            as needed).
        catalog : str
            Catalog name (e.g. ``'2MASS'``, ``'PS1'``).
        filt : str
            Catalog filter name.
        log : ColoredLogger, optional
            Logger for progress.

        Returns
        -------
        astropy.table.Table
            ``cat`` with magnitudes transformed in place when applicable.
        """
        if self.get_magsys() != DEFAULT_MAGSYS:
            if log:
                log.warning(
                    'Catalog magnitude conversion is only implemented for '
                    f'{DEFAULT_MAGSYS!r}; magsys={self.get_magsys()!r} unchanged.'
                )
            else:
                print(
                    'WARNING: non-AB MAGSYS; catalog mags not converted '
                    f'({self.get_magsys()}).'
                )
            return cat

        if catalog == '2MASS':
            cat['mag'] = self.twomass_vega_to_ab(cat['mag'], filt)

        if catalog == '2MASS' and filt == 'Y':
            cat = cat[~np.isnan(cat['Kmag'])]
            cat['mag'], cat['mag_err'] = self.twomass_y_mag_to_ab(
                cat['mag'], cat['mag_err'], cat['Kmag'], cat['e_Kmag'])

        return cat

    @staticmethod
    def twomass_vega_to_ab(mag, filt):
        """Convert 2MASS Vega magnitudes to AB for the given band.

        Parameters
        ----------
        mag : array-like
            Vega magnitude(s).
        filt : str
            Band: ``J``, ``H``, ``K``, or ``Ks``.

        Returns
        -------
        ndarray
            AB magnitude(s).
        """
        mag = np.asarray(mag, dtype=float)
        out = mag.copy()
        if filt == 'J':
            out = out + (4.56 - 3.65)
        elif filt == 'H':
            out = out + (4.71 - 3.32)
        elif filt in ('K', 'Ks'):
            out = out + (5.14 - 3.29)
        return out

    def twomass_y_mag_to_ab(self, J, J_err, Ks, Ks_err):
        """Synthetic Y-band magnitude in AB from 2MASS Vega :math:`J` and :math:`K_s`.

        Uses GF18 equation (6) for VISTA :math:`Y` from 2MASS, then GF18
        Appendix D equation (D3) to convert that VISTA :math:`Y` to AB.  The
        2MASS point-source catalog ``Kmag`` column is the :math:`K_s` band.

        Parameters
        ----------
        J, J_err : array-like or float
            2MASS :math:`J` magnitude and error (Vega).
        Ks, Ks_err : array-like or float
            2MASS :math:`K_s` magnitude and error (Vega); e.g. ``Kmag``, ``e_Kmag``.

        Returns
        -------
        tuple
            ``(Y_AB, Y_err)`` for use when ``MAGSYS`` is AB.
        """
        y_vega, y_err = self.twomass_j_ks_to_vista_y_vega(J, J_err, Ks, Ks_err)
        y_ab = self.vista_y_vega_to_ab(y_vega)
        return y_ab, y_err

    def get_zeropoint(self, flux, fluxerr, mag, magerr):
        """Fit zeropoint (and error) from flux/mag arrays using ODR.

        Parameters
        ----------
        flux : array-like
            Object fluxes (e.g. from aperture photometry).
        fluxerr : array-like
            Flux errors.
        mag : array-like
            Catalog magnitudes in the system given by :meth:`get_magsys` (e.g. AB
            when ``MAGSYS`` is :data:`DEFAULT_MAGSYS`).
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
        """Query VizieR for catalog magnitudes in a filter around given coordinates.

        Parameters
        ----------
        coords : astropy.coordinates.SkyCoord
            Target coordinates (used for region size).
        catalog : str
            Catalog name (e.g. 'PS1', '2MASS').
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

        catalog, cat_ID, cat_ra, cat_dec, cat_mag, cat_err = utilities.find_catalog(catalog, filt, coord_ra, coord_dec)
        
        med_coord = SkyCoord(coord_ra, coord_dec, unit='deg')

        seps = med_coord.separation(coords)
        max_sep = np.max(seps.to(u.deg).value)

        cols = [cat_ra, cat_dec, cat_mag, cat_err]
        # Add Kron mag if the catalog is PS1
        if cat_ID=='II/349':
            cols.append(f'{filt}Kmag')
        
        if log:
            log.info(f'Getting {catalog} catalog with ID {cat_ID} in filt {filt}')
            log.info(f'Querying around {coord_ra}, {coord_dec} deg')
        Vizier.clear_cache()
        width = np.max([2.0 * max_sep, 0.5])
        cat = None
        last_error = None
        for server in VIZIER_MIRRORS:
            try:
                vizier = Vizier(columns=cols, vizier_server=server)
                vizier.ROW_LIMIT = -1
                cat = vizier.query_region(med_coord, width=width*u.degree,
                    catalog=cat_ID)
                if cat is not None and len(cat) > 0:
                    if log:
                        log.info(f'Vizier query succeeded via {server}')
                    break
                cat = None
            except Exception as e:
                last_error = e
                if log:
                    log.warning(f'Vizier mirror {server} failed: {e}')
                else:
                    print(f'Vizier mirror {server} failed: {e}')
                cat = None
        if cat is None and last_error is not None and log:
            log.warning('All Vizier mirrors failed; last error: {}'.format(last_error))

        if cat is not None and len(cat)>0:
            cat = cat[0]
            cat = cat[~np.isnan(cat[cat_mag])]
            cat = cat[cat[cat_err]>0.]

            if cat_ID=='II/349':
                if log:
                    log.info('Cutting on Kron magnitudes')
                else:
                    print('Cutting on Kron magnitudes')

                nsources = len(cat)
                cat_kron = f'{filt}Kmag'
                mask = cat[cat_mag]-cat[cat_kron] < 0.1
                cat = cat[mask]
                nkron = len(cat)

                if log:
                    log.info(f'Cut catalog from {nsources} to {nkron}')
                else:
                    print(f'Cut catalog from {nsources} to {nkron}')

            cat.rename_column(cat_ra, 'ra')
            cat.rename_column(cat_dec, 'dec')
            cat.rename_column(cat_mag, 'mag')
            cat.rename_column(cat_err, 'mag_err')

            cat = self.apply_catalog_magnitude_system(cat, catalog, filt, log=log)

            return(cat, catalog, cat_ID)

        else:
            m='ERROR: cat {0}, ra {1}, dec {2} did not return a catalog'
            m=m.format(catalog, coord_ra, coord_dec)
            if log:
                log.error(m)
            else:
                print(m)
            return(None, None, None)

    def find_zeropoint(self, cmpfile, tel, match_radius=2.5*u.arcsec,
        phottable='APPPHOT', input_catalog=None, log=None):
        """Compute zeropoint from cmpfile photometry and catalog; write to FITS header.

        Matches sources to catalog (e.g. PS1), runs iterative ODR fit, and
        updates ZPTMAG, MAGSYS, ZPTNSTAR, ZPTCAT, etc. in the stack FITS.

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
        None
            cmpfile is updated in place.
        """
        if log:
            log.info(f'Importing catalog from file: {cmpfile}')
        else:
            print(f'Importing catalog from file: {cmpfile}')

        hdu = fits.open(cmpfile)
        header = hdu['SCI'].header
        filtorig = header['FILTER']
        catalog = tel.get_catalog(header)


        table = Table(hdu[phottable].data, meta=hdu[phottable].header)
        coords = SkyCoord(table['RA'], table['Dec'], unit='deg')

        # New metadata to update
        metadata = {}

        cat = None ; cat_ID = None

        filt = self.convert_filter_name(filtorig)
        if input_catalog is not None:
            cat = input_catalog
        else:
            if log:
                log.info(f'Downloading {catalog} catalog in {filt}')
            else:
                print(f'Downloading {catalog} catalog in {filt}')

            cat, catalog, cat_ID = self.get_catalog(coords, catalog, filt, log=log)

        if cat:
            min_mag = self.get_minmag(filt)
            cat = cat[cat['mag']>min_mag]

            coords_cat = SkyCoord(cat['ra'], cat['dec'], unit='deg')

            idx, d2, d3 = coords_cat.match_to_catalog_sky(coords)

            # Get matches from calibration catalog and cmpfile
            mask = d2 < match_radius
            cat = cat[mask]
            idx = idx[mask]

            if len(cat)==0:
                if log:
                    log.info('No star matches within match radius.')
                else:
                    print('No star matches within match radius.')
                return(None, None)

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

            if len(flux)==0:
                if log:
                    log.info('No suitable stars to calculate zeropoint .')
                return(None, None)

            zpt, zpterr, master_mask = self.zpt_iteration(flux, fluxerr,
                cat_mag, cat_magerr, log=log)

            # Set header variables
            metadata['ZPTNSTAR']=len(flux)
            metadata['ZPTMAG']=zpt
            metadata['ZPTMUCER']=zpterr
            metadata['MAGSYS']=self.get_magsys()
            metadata['ZPTCAT']=catalog
            metadata['ZPTCATID']=cat_ID
            metadata['ZPTPHOT']=phottable
            metadata['FILTER']=filtorig

            # Add limiting magnitudes
            if 'FWHM' in header.keys() and 'SKYSIG' in header.keys():
                fwhm = header['FWHM']
                sky = header['SKYSIG']
                Npix_per_FWHM_Area = 2.5 * 2.5 * fwhm * fwhm
                skysig_per_FWHM_Area = np.sqrt(Npix_per_FWHM_Area * (sky*sky))
                m3sigma = -2.5*np.log10(3.0*skysig_per_FWHM_Area)+zpt
                m5sigma = -2.5*np.log10(5.0*skysig_per_FWHM_Area)+zpt
                m10sigma = -2.5*np.log10(10.0*skysig_per_FWHM_Area)+zpt
                m3sigma = float('%.6f'%m3sigma)
                m5sigma = float('%.6f'%m5sigma)
                m10sigma = float('%.6f'%m10sigma)
                metadata['M3SIGMA']=m3sigma
                metadata['M5SIGMA']=m5sigma
                metadata['M10SIGMA']=m10sigma
                if log:
                    log.info(f'3-sigma limiting mag of image is {m3sigma}')
                else:
                    print(f'3-sigma limiting mag of image is {m3sigma}')

            hdu['PRIMARY'].header.update(metadata)
            hdu['SCI'].header.update(metadata)
            hdu[phottable].header.update(metadata)

            hdu.writeto(cmpfile, overwrite=True)
            
        elif log:
            log.info('No zeropoint calculated.')

    @staticmethod
    def twomass_j_ks_to_vista_y_vega(j, j_err, ks, ks_err,
            c_y=None, sigma_c_y=None):
        """VISTA :math:`Y` magnitude (Vega-like) from 2MASS :math:`J` and :math:`K_s` (Vega).

        Implements GF18 equation (6),

        .. math:: Y_V = J_2 + C_Y \\cdot (J-K_s)_2,

        with :math:`C_Y = 0.46 \\pm 0.02` (their notation: :math:`Y_V`, :math:`J_2`,
        :math:`(J-K_s)_2` for 2MASS).

        Error propagation uses

        .. math:: \\sigma_Y^2 = (1+C_Y)^2\\sigma_J^2 + C_Y^2\\sigma_{K_s}^2
            + (J-K_s)^2\\sigma_{C_Y}^2.

        Parameters
        ----------
        j, j_err : array-like
            2MASS :math:`J` and uncertainty (Vega).
        ks, ks_err : array-like
            2MASS :math:`K_s` and uncertainty (Vega).
        c_y, sigma_c_y : float, optional
            Override slope and its uncertainty (defaults: GF18 eq. 6).

        Returns
        -------
        tuple
            ``(Y_V, sigma_Y)`` in the VISTA Vega-like system (not AB).

        References
        ----------
        González-Fernández, C., Hodgkin, S. T., Irwin, M. J., et al. 2018,
        MNRAS, 474, 5459 (eq. 6); doi:10.1093/mnras/stx3073.
        """
        if c_y is None:
            c_y = GF18_TWOMASS_TO_VISTA_Y_SLOPE
        if sigma_c_y is None:
            sigma_c_y = GF18_TWOMASS_TO_VISTA_Y_SLOPE_ERR
        j = np.asarray(j, dtype=float)
        ks = np.asarray(ks, dtype=float)
        j_err = np.asarray(j_err, dtype=float)
        ks_err = np.asarray(ks_err, dtype=float)
        color_j_ks = j - ks
        y_v = j + c_y * color_j_ks
        var_y = ((1.0 + c_y) ** 2) * (j_err ** 2) + (c_y ** 2) * (ks_err ** 2)
        var_y = var_y + (sigma_c_y * color_j_ks) ** 2
        y_err = np.sqrt(var_y)
        return y_v, y_err

    @staticmethod
    def vista_y_vega_to_ab(y_vega):
        """Convert VISTA :math:`Y` (Vega-like) to AB using GF18 Appendix D eq. (D3).

        Parameters
        ----------
        y_vega : array-like
            VISTA :math:`Y` magnitude in the pipeline Vega-like system.

        Returns
        -------
        ndarray
            :math:`Y` in AB magnitudes.

        References
        ----------
        González-Fernández et al. 2018, Appendix D eq. (D3); see also
        Hewett et al. 2006, MNRAS, 367, 454.
        """
        return np.asarray(y_vega, dtype=float) + GF18_VISTA_Y_VEGA_TO_AB

    def y_band_from_jk(self, J, J_err, K, K_err):
        """VISTA :math:`Y` (Vega-like) from 2MASS :math:`J` and :math:`K_s`; use :meth:`twomass_j_ks_to_vista_y_vega`.

        Parameters ``K``, ``K_err`` are the 2MASS :math:`K_s` magnitude and error
        (catalog names ``Kmag``, ``e_Kmag``).  For AB output and zeropoint
        fitting with ``MAGSYS='ABMAG'``, use :meth:`twomass_y_mag_to_ab` instead,
        which applies GF18 eq. (6) and Appendix D (D3).

        Parameters
        ----------
        J, J_err : array-like or float
            2MASS :math:`J` (Vega).
        K, K_err : array-like or float
            2MASS :math:`K_s` (Vega).

        Returns
        -------
        tuple
            ``(Y_V, sigma_Y)`` as GF18 eq. (6).
        """
        return self.twomass_j_ks_to_vista_y_vega(J, J_err, K, K_err)

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

    Y_band = y_band_from_jk

