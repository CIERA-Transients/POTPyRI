"Function for calculating zero points during flux calibration."
"Authors: Kerry Paterson, Charlie Kilpatrick"

# Initial version tracking on 09/21/2024
__version__ = "1.0"

import numpy as np

from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits

# Internal dependency
from potpyri.utils import utilities

class absphot(object):
    def __init__(self, iterations=5, sigma=5):

        self.iterations = iterations
        self.sigma = sigma

        self.odr_iter = 2

    def get_zeropoint(self, flux, fluxerr, mag, magerr):

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

        if log: log.info(f'Searching for catalog {catalog}')
        
        cat_ID, cat_ra, cat_dec, cat_mag, cat_err = utilities.find_catalog(catalog, filt)
        
        coord_ra = np.median([c.ra.degree for c in coords])
        coord_dec = np.median([c.dec.degree for c in coords])
        
        med_coord = SkyCoord(coord_ra, coord_dec, unit='deg')

        seps = med_coord.separation(coords)
        max_sep = np.max(seps.to(u.deg).value)
        
        vizier = Vizier(columns=[cat_ra, cat_dec, cat_mag, cat_err])
        vizier.ROW_LIMIT = -1
        if log: 
            log.info(f'Getting {catalog} catalog with ID {cat_ID} in filt {filt}')
            log.info(f'Querying around {coord_ra}, {coord_dec} deg')
        cat = vizier.query_region(med_coord, width=1.2*max_sep*u.degree, 
            catalog=cat_ID)

        if len(cat)>0:
            cat = cat[0]
            cat = cat[~np.isnan(cat[cat_mag])]
            cat = cat[cat[cat_err]>0.]

            cat.rename_column(cat_ra, 'ra')
            cat.rename_column(cat_dec, 'dec')
            cat.rename_column(cat_mag, 'mag')
            cat.rename_column(cat_err, 'mag_err')

            # Convert to AB magnitudes
            if catalog=='2MASS':
                if filt=='J': cat['mag'] = cat['mag'] + (4.56-3.65)
                if filt=='H': cat['mag'] = cat['mag'] + (4.71-3.32)
                if filt=='K': cat['mag'] = cat['mag'] + (5.14-3.29)
                if filt=='Ks': cat['mag'] = cat['mag'] + (5.14-3.29)

            if catalog=='2MASS' and filt=='Y':
                cat = cat[~np.isnan(cat['Kmag'])]
                cat['mag'], cat['mag_err'] = self.Y_band(cat['mag'], 
                    cat['mag_err'], cat['Kmag'], cat['e_Kmag'])

            return(cat)

        else:
            m='ERROR: cat {0}, ra {1}, dec {2} did not return a catalog'
            m=m.format(catalog, coord_ra, coord_dec)
            if log:
                log.error(m)
            else:
                print(m)
            return(None)

    def find_zeropoint(self, cmpfile, filt, catalog, match_radius=2.5*u.arcsec,
        phottable='APPPHOT', input_catalog=None, log=None):

        if log:
            log.info(f'Importing catalog from file: {cmpfile}')
        else:
            print(f'Importing catalog from file: {cmpfile}')

        hdu = fits.open(cmpfile)
        header = hdu['PRIMARY'].header
        table = Table(hdu[phottable].data, meta=hdu[phottable].header)
        coords = SkyCoord(table['RA'], table['Dec'], unit='deg')

        # New metadata to update
        metadata = {}

        filt = self.convert_filter_name(filt)
        if input_catalog is not None:
            cat = input_catalog
        else:
            if log:
                log.info(f'Downloading {catalog} catalog in {filt}')
            else:
                print(f'Downloading {catalog} catalog in {filt}')

            cat = self.get_catalog(coords, catalog, filt, log=log)

        min_mag = self.get_minmag(filt)
        cat = cat[cat['mag']>min_mag]

        if cat:
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
            metadata['ZPTCAT']=catalog
            metadata['ZPTPHOT']=phottable

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
            hdu[phottable].header.update(metadata)

            hdu.writeto(cmpfile, overwrite=True)
            
            return(zpt, zpterr)
        else:
            if log:
                log.info('No zeropoint calculated.')
            return(None, None)

    def Y_band(self, J, J_err, K, K_err):
        Y = J+0.46*(J-K)
        JK_err = np.sqrt(J_err**2+K_err**2)
        JKY_err = 0.46*(J-K)*np.sqrt((0.02/0.46)**2+(JK_err/(J-K))**2)
        Y_err = np.sqrt(J_err**2+JKY_err**2)
        return Y, Y_err

    def convert_filter_name(self, filt):
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
            return 'r'
        if filt=='Y':
            return 'J'
        else:
            return filt
    
    def get_minmag(self, filt):
        if filt=='J':
            return 15.5
        if filt=='K':
            return 13.0
        if filt=='Y':
            return 15.0
        else:
            return 16.0

if __name__=="__main__":
    zp_cal = absphot()

    zp, zp_err = zp_cal.find_zeropoint(
        '/Users/ckilpatrick/Dropbox/Data/POTPyRI/test/MMIRS/red/GRB231117A_J.J.ut231204.32.11.stk.fits', 
        'J', '2MASS', log=None)
