from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.background import Background2D
from photutils.detection import IRAFStarFinder, DAOStarFinder
from photutils.psf import EPSFModel
from photutils.psf import extract_stars, EPSFBuilder
from photutils.psf import PSFPhotometry
from photutils.psf import SourceGrouper
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry

from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
from astropy.stats import gaussian_sigma_to_fwhm, sigma_clipped_stats
from astropy.table import Table, unique, Column
from astropy.visualization import simple_norm
from astropy.wcs import WCS

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.ndimage import rotate
from scipy.stats import sigmaclip
import numpy as np
import sys
import glob
import copy
import os

from utilities.Gaussians import fix_x0, twoD_Gaussian
from utilities import util

import warnings
warnings.filterwarnings('ignore')

bkgrms = MADStdBackgroundRMS()
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()

def generate_epsf(img_file, x, y, size=11, oversampling=3, maxiters=5,
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
        maxiters=maxiters, progress_bar=True)

    epsf, fitted_stars = epsf_builder(stars)

    return(epsf)

def construct_error_array(img_data):

    rms = 0.5 * (
        np.percentile(img_data[img_data != 0.0], 84.13)
        - np.percentile(img_data[img_data != 0.0], 15.86)
    )

    poisson = img_data[img_data < 0.]=0.
    error = np.sqrt(img_data + rms**2)

    return(error)

def run_photometry(img_file, epsf, fwhm, x, y, subtract_back=False,
    forced=False):

    img_hdu = fits.open(img_file)
    if subtract_back:
        bkg = Background2D(img_hdu[0].data, (21,21), filter_size=(3,3))
        image = img_hdu[0].data - bkg.background
        ndimage = NDData(data=image)
    else:
        image = img_hdu[0].data
        ndimage = NDData(data=img_hdu[0].data)

    psf = copy.copy(epsf)

    stars_tbl = Table()
    stars_tbl['x'] = x
    stars_tbl['y'] = y
    stars = extract_stars(ndimage, stars_tbl, size=11)

    stars_tbl['flux'] = np.array([stars[0].estimate_flux()])

    targets = Table()
    targets['x_0'] = stars_tbl['x']
    targets['y_0'] = stars_tbl['y']
    targets['flux_0'] = stars_tbl['flux']

    if forced:
        psf.x_0.fixed = True
        psf.y_0.fixed = True

    grouper = SourceGrouper(min_separation=21)
    finder = DAOStarFinder(6.0, 2.0)

    photometry = PSFPhotometry(psf_model=psf, fit_shape=(11,11),
        grouper=grouper, aperture_radius=10, finder=finder)

    error = construct_error_array(image)
    result_tab = photometry(image, error=error)

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

def get_star_catalog(img_file, fwhm_init=5.0, threshold=3.5):

    img_hdu = fits.open(img_file)
    data = img_hdu[0].data

    std = bkgrms(data)
    iraffind = IRAFStarFinder(threshold=threshold*std,
        fwhm=fwhm_init, minsep_fwhm=0.5,
        roundhi=0.2, roundlo=0,sharplo=0.4, sharphi=1.0)

    # Get initial set of stars in image with iraffind
    print('Finding stars...')
    stars = iraffind(data)
    mask = ((stars['xcentroid'] > 50) & (stars['xcentroid'] < (data.shape[1] -51))
                & (stars['ycentroid'] > 50) & (stars['ycentroid'] < (data.shape[0] -51)))
    stars = stars[mask]
    
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
    outdir='', subtract_back=False, fwhm_scale_psf=2.5, log=None,
    star_param={'snthresh_psf': 20.0, 'fwhm_init': 5.0, 'snthresh_final': 10.0}):

    stars = get_star_catalog(img_file, fwhm_init=star_param['fwhm_init'])
    img_hdu = fits.open(img_file)

    # Get image statistics
    mask = img_hdu[0].data!=0.0
    mean, median, std_sky = sigma_clipped_stats(img_hdu[0].data[mask],
        sigma=5.0)

    if log:
        log.info('Found {0} stars'.format(len(stars)))
    else:
        print('Found {0} stars'.format(len(stars)))

    stars = stars['xcentroid','ycentroid','fwhm','sharpness','roundness',
        'npix','pa','flux','sky']

    # Estimate the uncertainty from sky and flux values
    stars['flux_err'] = np.sqrt(stars['flux']+stars['npix']*stars['sky'])

    mask = ((stars['sharpness'] < np.median(stars['sharpness'])+np.std(stars['sharpness'])) &\
                (stars['roundness'] < np.median(stars['roundness'])+3*np.std(stars['roundness'])) &\
                (stars['roundness'] > np.median(stars['roundness'])-3*np.std(stars['roundness'])))

    fwhm_stars = stars[mask]

    m='Masked to {0} stars based on sharpness, roundness'
    
    if log:
        log.info(m.format(len(fwhm_stars)))
    else:
        print(m.format(len(fwhm_stars)))

    fwhm_clipped, _, _ = sigmaclip(fwhm_stars['fwhm'])
    fwhm = np.median(fwhm_clipped)
    std_fwhm = np.std(fwhm_clipped)
    mask = (fwhm_stars['fwhm'] > fwhm-3*std_fwhm) &\
        (fwhm_stars['fwhm'] < fwhm+3*std_fwhm)
    fwhm_stars = fwhm_stars[mask]
    fwhm = np.median(fwhm_stars['fwhm'])

    if log:
        log.info('Masked to {0} stars based on FWHM'.format(len(fwhm_stars)))
    else:
        print('Masked to {0} stars based on FWHM'.format(len(fwhm_stars)))

    step_size = 0.5
    radii = np.arange(step_size, 2.5*fwhm, step_size)
    x_data = np.arange(0, 2.5*fwhm, step_size)
    apers_area = [np.pi*(step_size**2)]
    for r in radii:
        apers_area.append(np.pi*((r+step_size)**2 - r**2))
    coords = [(fwhm_stars['xcentroid'][i],fwhm_stars['ycentroid'][i]) for i in range(len(fwhm_stars))]
    apertures = [CircularAperture(coords, r=step_size)]     # For circle aperture around center
    for r in radii:
        apertures.append(CircularAnnulus(coords, r_in=r, r_out=r+step_size))  # Annuli apertures
    phot_table = aperture_photometry(img_hdu[0].data, apertures)
    sigmas = []
    for i, source in enumerate(phot_table):
        sums = np.array([source[i] / apers_area[i - 3] for i in range(3, len(source))])
        if median+10*std_sky <= sums[0] <= 60000 and median-3*std_sky < sums[-1] < median+3*std_sky:
            w_fit = x_data[~np.isnan(sums)]
            f_fit = sums[~np.isnan(sums)]
            g, _ = curve_fit(fix_x0, w_fit, f_fit)
            sigmas.append(np.abs(g[1]))
    sigmas_post_clipping, _, _ = sigmaclip(sigmas, low=3, high=3)
    fwhm = 2.35482*np.median(sigmas_post_clipping)
    std_fwhm = np.std(sigmas_post_clipping)
    
    if np.isnan(fwhm):
        if log:
            log.error('Could not refine FWHM, using first value.')
        else:
            print('Could not refine FWHM, using first value.')
        fwhm = np.median(fwhm_stars['fwhm'])
        std_fwhm = np.std(fwhm_stars['fwhm'])
    
    if log:
        log.info('FWHM={0}+/-{1}'.format('%2.4f'%fwhm,'%2.4f'%std_fwhm))
    else:
        print('FWHM={0}+/-{1}'.format('%2.4f'%fwhm,'%2.4f'%std_fwhm))

    metadata={'FWHM':fwhm, 'EFWHM':std_fwhm, 'SKYADU': std_sky}

    if subtract_back:
        bkg = Background2D(img_hdu[0].data, (21,21), filter_size=(3,3))
        backsub = img_hdu[0].data - bkg.background
        ndimage = NDData(data=backsub)
        backhdu = fits.PrimaryHDU(bkg.background)
        backsubhdu = fits.PrimaryHDU(backsub)
    else:
        ndimage = NDData(data=img_hdu[0].data)

    if write_out_back:
        if log:
            log.info('Writing out background and background-subtracted file...')
        else:
            print('Writing out background and background-subtracted file...')
        back_file = img_file.replace('.fits','.back.fits')
        backsub_file = img_file.replace('.fits','.backsub.fits')

        if log:
            log.info('Background file:',back_file)
            log.info('Background-subtracted file:',backsub_file)
        else:
            print('Background file:',back_file)
            print('Background-subtracted file:',backsub_file)

        backhdu.writeto(back_file, overwrite=True)
        backsubhdu.writeto(backsub_file, overwrite=True)

    mask = (fwhm_stars['flux']/fwhm_stars['flux_err'] > star_param['snthresh_psf'])
    bright = fwhm_stars[mask]
    if len(bright) < 90:
        m='The number of bright stars for PSF determination is < 90, lowering flux threshold to 5.'
        if log:
            log.info(m)
        else:
            print(m)        
        mask = (fwhm_stars['flux']/fwhm_stars['flux_err'] > 5)
        bright = fwhm_stars[mask]

    m='Masked to {0} PSF stars based on flux.'
    if log:
        log.info(m.format(len(bright)))
    else:
        print(m.format(len(bright)))

    metadata['NPSFSTAR']=len(bright)

    # Instantiate EPSF
    size=int(fwhm*fwhm_scale_psf)
    if size%2==0: size=size+1

    if log:
        log.info(f'EPSF size will be {size} pixels')
    else:
        print(f'EPSF size will be {size} pixels')

    epsf = generate_epsf(img_file, bright['xcentroid'],
        bright['ycentroid'], size=size, oversampling=2, maxiters=5,
        log=log)

    #Check shape of PSF compared to 2d Gaussian
    m='Checking the 2D shape of the PSF.'
    if log:
        log.error(m)
    else:
        print(m)  

    d = epsf.data
    x_mesh, y_mesh = np.meshgrid(np.linspace(0, np.shape(d)[1] - 1, np.shape(d)[1]),
                                        np.linspace(0, np.shape(d)[0] - 1, np.shape(d)[0]))
    try:
        popt, pcov = curve_fit(twoD_Gaussian, (x_mesh, y_mesh), d.ravel(),p0=(np.max(d),
                                                 25, 25, fwhm/2.35482, fwhm/2.35482, 0, 0))
        sigmax, sigmay = popt[3], popt[4]
        if log:
            log.info('Sigma (x, y) of PSF: {0}, {1}'.format(sigmax, sigmay))
        else:
            print(log.info('Sigma (x, y) of PSF: {0}, {1}'.format(sigmax, sigmay)))
    except:
        m='Error calculating the PSF. Check PSF png for quaility control.'
        if log:
            log.error(m)
        else:
            print(m)          

    mask = (stars['flux']/stars['flux_err'] > star_param['snthresh_final'])
    all_stars = stars[mask]
    if log:
        log.info('Final catalog is {0} stars'.format(len(all_stars)))
        log.info('Getting final photometry...')
    else:
        print('Final catalog is {0} stars'.format(len(all_stars)))
        print('Getting final photometry...')

    photometry = run_photometry(img_file, epsf, fwhm, all_stars['xcentroid'],
        all_stars['ycentroid'], subtract_back=subtract_back)

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
        log.info('Got {0} stars for final photometry'.format(len(photometry)))
    else:
        print('Got {0} stars for final photometry'.format(len(photometry)))
    metadata['NOBJECT']=len(photometry)

    # We should always write out catalog for stars
    phot_file = os.path.join(outdir,img_file.replace('.fits','.pcmp'))
    write_out_catalog(photometry, img_file, colnames, sigfig, phot_file,
        metadata)

    if write_out_psf_stars:
        all_stars.sort('flux')
        outname = os.path.join(outdir,img_file.replace('.fits','.psf.stars'))
        all_stars.write(outname, format='ascii.no_header', overwrite=True)

    if write_out_epsf_img:
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        outname = os.path.join(outdir, img_file.replace('.fits','.epsf.png'))
        plt.savefig(outname)
        plt.clf()

    if write_out_epsf_file:
        hdu = fits.PrimaryHDU(epsf.data)
        hdu.header['FWHM']=fwhm
        outname = os.path.join(outdir, img_file.replace('.fits','.psf.fits'))
        hdu.writeto(outname, overwrite=True)

    # Finally return EPSF
    return(epsf, fwhm)

if __name__=="__main__":

    do_phot('/Users/ckilpatrick/Dropbox/Data/FRB/FRB240619/R_band/red/meertrap_wf_hos_R_1R_11.fits')
