import numpy as np
import os
import sys
from Vizier_catalogs import find_catalog
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
from psf import *

class absphot(object):
    def __init__(self, iterations=5, sigma=5):

        self.iterations = iterations
        self.sigma = sigma

        self.odr_iter = 10

    def get_zeropoint(self, flux, fluxerr, mag, magerr):

        from scipy.odr import ODR, Model, Data, RealData
        def magnitude(zpt, flux):
            return(zpt - 2.5*np.log10(flux))

        zpt_guess = np.median(mag + 2.5*np.log10(flux))

        data = RealData(flux, mag, fluxerr, magerr)
        model = Model(magnitude)
        odr = ODR(data, model, [zpt_guess], maxit=self.odr_iter)

        output = odr.run()

        zpt = output.beta[0]
        zpterr = output.sd_beta[0]

        return(zpt, zpterr)

    def zpt_iteration(self, flux, fluxerr, mag, magerr):

        flux = np.array(flux)
        fluxerr = np.array(fluxerr)
        mag = np.array(mag)
        magerr = np.array(magerr)

        for i in np.arange(self.iterations):

            nobj = len(flux)
            zpt, zpterr = self.get_zeropoint(flux, fluxerr, mag, magerr)

            mag_deriv = -2.5*np.log10(flux)+zpt
            magerr_deriv = np.array(2.5/np.log(10) * fluxerr/flux)

            total_err = np.sqrt(magerr**2+magerr_deriv**2+zpterr**2)

            mask = np.abs(mag - mag_deriv) < self.sigma * total_err

            flux=flux[mask] ; fluxerr=fluxerr[mask]
            mag=mag[mask] ; magerr=magerr[mask]

            m='Iteraction {i}: {N} obj, zpt = {zpt}+/-{zpterr}, {s}-sigma '+\
                'clip to {M} obj'
            print(m.format(i=i, N=nobj, zpt='%2.4f'%zpt, zpterr='%2.4f'%zpterr,
                s=self.sigma, M=len(flux)))

        return(zpt, zpterr)

    def get_catalog(self, coords, catalog, filt, size=0.5):

        cat_ID, cat_ra, cat_dec, cat_mag, cat_err = find_catalog(catalog, filt)
        coord_ra = np.median([c.ra.degree for c in coords])
        coord_dec = np.median([c.dec.degree for c in coords])
        med_coord = SkyCoord(coord_ra, coord_dec, unit='deg')
        Vizier.ROW_LIMIT = -1
        cat = Vizier.query_region(med_coord, width=0.5*u.degree, catalog=cat_ID)
        if len(cat)>0:
            cat = cat[0]
        else:
            m='ERROR: cat {0}, ra {1}, dec {2} did not return a catalog'
            m=m.format(catalog, coord_ra, coord_dec)
            raise RuntimeError(m)

        cat.rename_column(cat_ra, 'ra')
        cat.rename_column(cat_dec, 'dec')
        cat.rename_column(cat_mag, 'mag')
        cat.rename_column(cat_err, 'mag_err')

        return(cat)

    def find_zeropoint(self, cmpfile, filt, catalog, match_radius=2.0*u.arcsec):

        header, table = import_catalog(cmpfile)
        coords = SkyCoord(table['RA'], table['Dec'], unit='deg')
        cat = self.get_catalog(coords, catalog, header['FILTER'])

        coords_cat = SkyCoord(cat['ra'], cat['dec'], unit='deg')

        idx, d2, d3 = coords_cat.match_to_catalog_sky(coords)

        # Get matches from calibration catalog and cmpfile
        mask = d2 < match_radius
        cat = cat[mask]
        idx = idx[mask]

        match_table = None
        for i in idx:
            if not match_table:
                match_table = Table(table[i])
            else:
                match_table.add_row(table[i])

        metadata = {}
        metadata['nmatch']=len(match_table)

        zpt, zpterr = self.zpt_iteration(match_table['flux'],
            match_table['flux_err'], cat['mag'], cat['mag_err'])

        metadata['ZPTMAG']=zpt
        metadata['ZPTMUCER']=zpterr

        # Add limiting magnitudes
        if 'FWHM' in header.keys():
            fwhm = header['FWHM']
            Npix_per_FWHM_Area = 2.5 * 2.5 * fwhm * fwhm

        modify_catalog(cmpfile)

        print(zpt, zpterr)
