import numpy as np
import os
import sys

from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
from Vizier_catalogs import find_catalog

from utilities.util import *

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

    def zpt_iteration(self, flux, fluxerr, mag, magerr, log=None):

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

            m='Iteration {i}: {N} obj, zpt = {zpt}+/-{zpterr}, {s}-sigma '+\
                'clip to {M} obj'
            if log:
                log.info(m.format(i=i, N=nobj, zpt='%2.4f'%zpt, zpterr='%2.4f'%zpterr,
                s=self.sigma, M=len(flux)))
            else:
                print(m.format(i=i, N=nobj, zpt='%2.4f'%zpt, zpterr='%2.4f'%zpterr,
                s=self.sigma, M=len(flux)))

        return(zpt, zpterr)

    def get_catalog(self, coords, catalog, filt, size=0.5, log=None):

        if log:
            log.info('Searching for catalog '+catalog)
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
            if log:
                log.error(m)
            raise RuntimeError(m)

        cat.rename_column(cat_ra, 'ra')
        cat.rename_column(cat_dec, 'dec')
        cat.rename_column(cat_mag, 'mag')
        cat.rename_column(cat_err, 'mag_err')

        return(cat)

    def find_zeropoint(self, cmpfile, filt, catalog, match_radius=2.0*u.arcsec, log=None):

        header, table = import_catalog(cmpfile)
        coords = SkyCoord(table['RA'], table['Dec'], unit='deg')
        cat = self.get_catalog(coords, catalog, filt, log=log)

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
        metadata['ZPTNSTAR']=len(match_table)

        zpt, zpterr = self.zpt_iteration(match_table['flux'],
            match_table['flux_err'], cat['mag'], cat['mag_err'], log=log)
        if log:
            log.info('Calculating AB correction.')
        cor = self.AB_conversion(catalog,filt)
        zpt+=cor            

        if log:
            log.info('Added AB correction of %2.3f mag.'%cor)
            log.info('Final zeropoint calculated:')
            log.info('zpt = %2.4f +/- %2.4f AB mag'%(zpt,zpterr))

        metadata['ZPTMAG']=zpt
        metadata['ZPTMUCER']=zpterr

        # Add limiting magnitudes
        if 'FWHM' in header.keys() and 'SKYADU' in header.keys():
            fwhm = header['FWHM']
            sky = header['SKYADU']
            Npix_per_FWHM_Area = 2.5 * 2.5 * fwhm * fwhm
            skysig_per_FWHM_Area = np.sqrt(Npix_per_FWHM_Area * (sky*sky))
            metadata['M3SIGMA']=-2.5*np.log10(3.0*skysig_per_FWHM_Area)+zpt
            metadata['M5SIGMA']=-2.5*np.log10(5.0*skysig_per_FWHM_Area)+zpt
            metadata['M10SIGMA']=-2.5*np.log10(10.0*skysig_per_FWHM_Area)+zpt

        modify_catalog(cmpfile, metadata)
    
    def AB_conversion(self, catalog, fil):
        if catalog=='2MASS':
            if fil == 'K' or 'Ks':
                cor = 1.85
            elif fil == 'J':
                cor = 0.91
            elif fil == 'H':
                cor = 1.39
            elif fil == 'Y':
                cor = 0.634
        else:
            cor = 0
        return cor
