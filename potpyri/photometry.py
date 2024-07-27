from photutils.aperture import ApertureStats
from photutils.aperture import CircularAperture
from photutils.background import MMMBackground
from photutils.background import MADStdBackgroundRMS
from photutils.background import Background2D
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFModel
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder
from photutils.psf import PSFPhotometry
from photutils.psf import SourceGrouper

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.nddata import NDData
from astropy.table import Table, unique, Column, hstack, vstack
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.modeling import models, fitting, functional_models

import matplotlib.pyplot as plt

from scipy.ndimage import rotate
from scipy.stats import sigmaclip
import numpy as np
import sys
import glob
import copy
import os

from utilities import util

import warnings
warnings.filterwarnings('ignore')


def extract_aperture_stats(img_file, stars, aperture_radius=10.0):

    apertable = Table([[0.],[0.],[0.],[0.],[0.],[0.]], 
        names=('fwhm','semimajor_sigma','semiminor_sigma',
        'orientation','eccentricity','signal_to_noise')).copy()[:0]

    img_hdu = fits.open(img_file)
    data = img_hdu[0].data 
    mask = img_hdu[1].data.astype(bool)
    error = img_hdu[2].data

    for star in stars:

        aper = CircularAperture((star['xcentroid'], star['ycentroid']),
            aperture_radius)

        aperstats = ApertureStats(data, aper, mask=mask, error=error)

        apertable.add_row([aperstats.fwhm.value, aperstats.semimajor_sigma.value,
            aperstats.semiminor_sigma.value, aperstats.orientation.value,
            aperstats.eccentricity,aperstats.sum/aperstats.sum_err])
    
    return(apertable)


def generate_epsf(img_file, x, y, size=11, oversampling=2, maxiters=5,
    log=None):
    # Construct stars table from bright
    stars_tbl = Table()
    stars_tbl['x'] = x
    stars_tbl['y'] = y

    img_hdu = fits.open(img_file)
    ndimage = NDData(data=img_hdu[0].data)

    stars = extract_stars(ndimage, stars_tbl, size=size)

    if log:
        log.info('Extracted {0} stars.  Building EPSF...'.format(len(stars)))
    else:
        print('Extracted {0} stars.  Building EPSF...'.format(len(stars)))

    epsf_builder = EPSFBuilder(oversampling=oversampling,
        maxiters=maxiters, progress_bar=True,
        sigma_clip=SigmaClip(sigma=3, sigma_lower=3, sigma_upper=3, 
            maxiters=10, cenfunc='median', stdfunc='std', grow=False))

    epsf, fitted_stars = epsf_builder(stars)

    return(epsf)

def extract_fwhm_from_epsf(epsf, fwhm_init):

    data = epsf.normalized_data
    x = np.arange(data.shape[0])
    y = np.arange(data.shape[1])
    xx, yy = np.meshgrid(x, y) 

    p_init = functional_models.Moffat2D(amplitude=0.5, x_0=data.shape[0]/2.,
        y_0=data.shape[1]/2., gamma=fwhm_init, alpha=1.)
    fit_p = fitting.LevMarLSQFitter()

    p = fit_p(p_init, xx, yy, data)

    return(p.fwhm)

def run_photometry(img_file, epsf, fwhm, threshold, shape, stars):

    img_hdu = fits.open(img_file)
    image = img_hdu[0].data
    ndimage = NDData(data=img_hdu[0].data)
    mask = img_hdu[1].data.astype(bool)
    error = img_hdu[2].data

    psf = copy.copy(epsf)

    stars_tbl = Table()
    stars_tbl['x'] = stars['xcentroid']
    stars_tbl['y'] = stars['ycentroid']
    size = int(np.round(2.5*fwhm))
    if size%2==0: size = size+1
    stars = extract_stars(ndimage, stars_tbl, size=size)

    stars_tbl['flux'] = np.array([stars[0].estimate_flux()])

    targets = Table()
    targets['x_0'] = stars_tbl['x']
    targets['y_0'] = stars_tbl['y']
    targets['flux_0'] = stars_tbl['flux']

    grouper = SourceGrouper(min_separation=shape)
    finder = DAOStarFinder(threshold, fwhm)

    photometry = PSFPhotometry(psf_model=psf, fit_shape=(shape,shape),
        grouper=grouper, aperture_radius=int(shape*1.5), finder=finder,
        progress_bar=True)

    result_tab = photometry(image, mask=mask, error=error,
        init_params=targets)

    mask = ~np.isnan(result_tab['flux_err'])
    result_tab = result_tab[mask]

    mask = result_tab['flux_fit'] > 0.
    result_tab = result_tab[mask]

    mask = ~np.isnan(result_tab['x_err'])
    result_tab = result_tab[mask]

    mask = ~np.isnan(result_tab['y_err'])
    result_tab = result_tab[mask]

    return(result_tab)

def rotate_psf_to_img(psf_file, img_file, pa=0.0):
    psf_hdu = fits.open(psf_file)
    img_hdu = fits.open(img_file)

    mask = img_hdu[0].data == img_hdu[0].data[0,0]
    img_hdu[0].data[mask]=np.nan

    # Rotate PSF to the angle of the input image
    psf_rot = rotate(psf_hdu[0].data, pa)

    return(psf_rot)

def get_star_catalog(img_file, fwhm_init=5.0, threshold=3.5, log=None):

    img_hdu = fits.open(img_file)
    data = img_hdu[0].data
    mask = img_hdu[1].data.astype(bool)
    err = img_hdu[2].data

    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(data)
    daofind = DAOStarFinder(fwhm=fwhm_init, threshold=threshold*std,
        exclude_border=True, min_separation=fwhm_init)

    # Get initial set of stars in image with iraffind
    if log:
        log.info('Finding stars...')
    else:
        print('Finding stars...')

    stars = daofind(data, mask=mask)

    mask = stars['peak'] > 0.
    stars = stars[mask]

    stars.sort('flux')
    
    return(stars)

def write_out_catalog(catalog, img_file, columns, sigfig, outfile, metadata):
    # zfill the column numbers to 2 for standard format with 01, 02, 03, etc.
    zfil = 2

    # Get image header
    header = fits.open(img_file)[0].header

    # Add metadata from photometry
    for key in metadata.keys():
        header[key]=metadata[key]

    header['NCOL']=len(columns)
    for i,col in enumerate(columns):
        header['COLUM{0}'.format(str(i+1).zfill(zfil))]=col

    catdata = []
    for row in catalog:
        outdata=[]
        outfmt=''
        for i,col,sig in zip(np.arange(len(columns)),columns, sigfig):
            fmt = '%7.{0}f'.format(sig)
            data = fmt%row[col]
            outdata.append(fmt%row[col])
            outfmt+='{'+'{0}:>{1}'.format(i,7+int(sig))+'} '
        outline = outfmt.format(*outdata)
        outline = str(outline)
        catdata.append(outline)

    # Finally write out catalog
    # Delete the outfile if it already exists:
    if os.path.exists(outfile):
        os.remove(outfile)
        
    util.write_catalog(outfile, header, catdata)

def do_phot(img_file, write_out_back=False,
    write_out_epsf_img=True, write_out_epsf_file=True, write_out_psf_stars=True,
    fwhm_scale_psf=3.5, log=None, oversampling=2,
    star_param={'snthresh_psf': 20.0, 'fwhm_init': 8.0, 'snthresh_final': 10.0}):

    stars = get_star_catalog(img_file, fwhm_init=star_param['fwhm_init'],
        threshold=star_param['snthresh_final']/2.0, log=log)

    stats = extract_aperture_stats(img_file, stars, 
        aperture_radius=2.5*star_param['fwhm_init'])

    stars = hstack([stars, stats])

    img_hdu = fits.open(img_file)

    # Get image statistics
    data = img_hdu[0].data
    mask = img_hdu[1].data.astype(bool)
    mean, median, std_sky = sigma_clipped_stats(data[~mask], sigma=5.0)

    metadata={'SKYADU':median, 'SKYSIG': std_sky}

    if log:
        log.info(f'Found {len(stars)} stars')
    else:
        print(f'Found {len(stars)} stars')

    # Estimate the uncertainty from sky and flux values
    stars['flux_err'] = stars['flux']/stars['signal_to_noise']

    med_sharp = np.median(stars['sharpness'])
    med_round = np.median(stars['roundness1'])
    std_sharp = np.std(stars['sharpness'])
    std_round = np.std(stars['roundness1'])

    mask = ((stars['sharpness'] < med_sharp + std_sharp) &\
            (stars['roundness1'] < med_round + 3*std_round) &\
            (stars['roundness1'] > med_round - 3*std_round))

    fwhm_stars = stars[mask]

    mask = ~np.isnan(fwhm_stars['fwhm'])
    fwhm_stars = fwhm_stars[mask]

    # Clip fwhm_stars by fwhm
    fwhm_clipped, _, _ = sigmaclip(fwhm_stars['fwhm'])
    fwhm = np.median(fwhm_clipped)
    std_fwhm = np.std(fwhm_clipped)
    mask = (fwhm_stars['fwhm'] > fwhm-3*std_fwhm) &\
        (fwhm_stars['fwhm'] < fwhm+3*std_fwhm)
    fwhm_stars = fwhm_stars[mask]
    fwhm = np.median(fwhm_stars['fwhm'])

    fwhm = float('%.4f'%fwhm)
    std_fwhm = float('%.4f'%std_fwhm)
    
    if log:
        log.info(f'''Masked to {len(fwhm_stars)} stars''')
    else:
        print(f'''Masked to {len(fwhm_stars)} stars''')

    mask = (fwhm_stars['signal_to_noise'] > star_param['snthresh_psf'])
    bright = fwhm_stars[mask]
    if len(bright) < 40:
        if log:
            log.info('Bright stars for PSF are <90, lowering S/N thresh to 5.')
        else:
            print('Bright stars for PSF are <90, lowering S/N tresh to 5.')

        mask = (fwhm_stars['signal_to_noise'] > 5)
        bright = fwhm_stars[mask]

    if log:
        log.info(f'Masked to {len(bright)} PSF stars based on flux.')
    else:
        print(f'Masked to {len(bright)} PSF stars based on flux.')

    metadata['NPSFSTAR']=len(bright)

    # Instantiate EPSF
    size=int(fwhm*fwhm_scale_psf)
    if size%2==0: size=size+1

    if log:
        log.info(f'EPSF size will be {size} pixels')
    else:
        print(f'EPSF size will be {size} pixels')

    epsf = generate_epsf(img_file, bright['xcentroid'], bright['ycentroid'], 
        size=size, oversampling=oversampling, maxiters=5, log=log)

    fwhm = extract_fwhm_from_epsf(epsf, fwhm*oversampling)
    # Scale by oversampling
    fwhm = fwhm/oversampling

    if log:
        log.info(f'FWHM={fwhm} pixels')
    else:
        print(f'FWHM={fwhm} pixels')

    metadata['FWHM']=fwhm

    mask = (stars['signal_to_noise'] > star_param['snthresh_final'])
    final_stars = stars[mask]
    if log:
        log.info(f'Final catalog is {len(final_stars)} stars')
        log.info('Getting final photometry...')
    else:
        print(f'Final catalog is {len(final_stars)} stars')
        print('Getting final photometry...')

    photometry = run_photometry(img_file, epsf, fwhm, 
        star_param['snthresh_final'], size, final_stars)

    # Get RA/Dec from final positions
    w = WCS(img_hdu[0].header)
    coords = w.pixel_to_world(photometry['x_fit'], photometry['y_fit'])

    # Join the photometry and all star catalogs and rename columns
    photometry['SN'] = photometry['flux_fit']/photometry['flux_err']
    photometry['mag'] = -2.5*np.log10(photometry['flux_fit'])
    photometry['mag_err'] = 2.5/np.log(10) * 1./photometry['SN']
    photometry['RA'] = [c.ra.degree for c in coords]
    photometry['Dec'] = [c.dec.degree for c in coords]
    photometry.add_column(Column([fwhm]*len(photometry), name='FWHM'))

    photometry.rename_column('x_fit', 'Xpos')
    photometry.rename_column('y_fit', 'Ypos')
    photometry.rename_column('x_err', 'Xpos_err')
    photometry.rename_column('y_err', 'Ypos_err')
    photometry.rename_column('flux_fit','flux')
    photometry.rename_column('local_bkg','SKY')
    # Sort from brightest to faintest
    photometry.sort('flux', reverse=True)

    # Get final list of column names we want in catalog in order and the number
    # of significant figures they should all have
    colnames=['Xpos','Ypos','Xpos_err','Ypos_err','mag','mag_err','flux',
        'flux_err','SN','SKY','RA','Dec']
    sigfig=[4,4,4,4,4,4,4,4,4,4,7,7]

    if log:
        log.info(f'Got {len(photometry)} stars for final photometry')
    else:
        print(f'Got {len(photometry)} stars for final photometry')

    metadata['NOBJECT']=len(photometry)

    # We should always write out catalog for stars
    phot_file = img_file.replace('.fits','.pcmp')
    write_out_catalog(photometry, img_file, colnames, sigfig, phot_file,
        metadata)

    if write_out_psf_stars:
        bright.sort('flux')
        outname = img_file.replace('.fits','.psf.stars')
        bright.write(outname, format='ascii', overwrite=True)

    if write_out_epsf_img:
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        outname = img_file.replace('.fits','.epsf.png')
        plt.savefig(outname)
        plt.clf()

    if write_out_epsf_file:
        hdu = fits.PrimaryHDU(epsf.data)
        hdu.header['FWHM']=fwhm
        outname = img_file.replace('.fits','.psf.fits')
        hdu.writeto(outname, overwrite=True)

    # Finally return EPSF
    return(epsf, fwhm)

if __name__=="__main__":

    snthresh_final = 10.0
    star_param={'snthresh_psf': snthresh_final*2.0, 
                                'fwhm_init': 5.0, 
                                'snthresh_final': snthresh_final}

    do_phot('/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/red/R155_host.R.ut240629.1R.11.stk.fits',
        star_param=star_param)
