"Functions for performing aperture and PSF photometry on pipeline images."
"Authors: Kerry Paterson, Charlie Kilpatrick"

# Initial version tracking on 09/21/2024
__version__ = "1.0"

from photutils.aperture import ApertureStats
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder
from photutils.psf import extract_stars
from photutils.psf import PSFPhotometry

from astropy.io import fits
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.nddata import NDData
from astropy.table import Table
from astropy.table import Column
from astropy.table import hstack
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.modeling import fitting
from astropy.modeling import functional_models

import matplotlib.pyplot as plt

from scipy.stats import sigmaclip
import numpy as np
import copy
import os

import warnings
warnings.filterwarnings('ignore')


def create_conv(outfile):
    with open(outfile, 'w') as f:
        f.write('CONV NORM \n')
        f.write('1 2 1 \n')
        f.write('2 4 2 \n')
        f.write('1 2 1 \n')

def create_params(outfile):
    with open(outfile, 'w') as f:
        f.write('NUMBER \n')
        f.write('X_IMAGE \n')
        f.write('Y_IMAGE \n')
        f.write('ALPHA_J2000 \n')
        f.write('DELTA_J2000 \n')
        f.write('MAG_AUTO \n')
        f.write('MAGERR_AUTO \n')
        f.write('FLUX_AUTO \n')
        f.write('FLUXERR_AUTO \n')
        f.write('FWHM_IMAGE \n')

def run_sextractor(img_file, log=None):

    hdu = fits.open(img_file)

    data = hdu['SCI'].data

    saturate = hdu['SCI'].header['SATURATE']

    paramfile = img_file.replace('.fits','.param')
    convfile = img_file.replace('.fits','.conv')
    tmpfile = img_file.replace('.fits','.tmp.fits')
    catfile = img_file.replace('.fits','.cat')

    create_params(paramfile)
    create_conv(convfile)

    datahdu = fits.PrimaryHDU()
    datahdu.data = data
    datahdu.header = hdu['SCI'].header
    datahdu.writeto(tmpfile, overwrite=True)

    if log:
        log.info(f'Running source extractor on {img_file}')
    else:
        print(f'Running source extractor on {img_file}')

    cmd = f'sex {tmpfile} -CATALOG_NAME {catfile} -CATALOG_TYPE ASCII_HEAD '
    cmd += f'-PARAMETERS_NAME {paramfile} -FILTER_NAME {convfile} '
    cmd += f'-SATUR_LEVEL {saturate} > /dev/null 2> /dev/null'

    os.system(cmd)

    if os.path.exists(catfile):
        if log:
            log.info(f'Reading {catfile}')
        else:
            print(f'Reading {catfile}')
        table = ascii.read(catfile)
    else:
        if log:
            log.error(f'Reading {catfile}')
        else:
            print(f'Reading {catfile}')
        table = None

    for file in [tmpfile, catfile, paramfile, convfile]:
        if os.path.exists(file):
            os.remove(file)

    return(table)

def extract_aperture_stats(img_data, img_mask, img_error, stars, 
    aperture_radius=10.0, log=None):

    apertable = Table([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],
        [0.],[0.]], 
        names=('fwhm','semimajor_sigma','semiminor_sigma',
        'orientation','eccentricity','signal_to_noise',
        'flux_best', 'flux_best_err','Xpos','Ypos',
        'Xpos_err','Ypos_err')).copy()[:0]

    # Estimate a reasonable aperture radius and centroid for sources
    fwhms=[]
    for i,star in enumerate(stars):
        aper = CircularAperture((star['xcentroid'], star['ycentroid']),
            aperture_radius)
        aperstats = ApertureStats(img_data, aper, mask=img_mask, 
            error=img_error)

        fwhms.append(aperstats.fwhm.value)
        stars[i]['xcentroid']=aperstats.xcentroid
        stars[i]['ycentroid']=aperstats.ycentroid

    if aperture_radius<2.5*np.nanmean(fwhms):
        aperture_radius=2.5*np.nanmean(fwhms)

    if log:
        log.info(f'New aperture radius={aperture_radius}')
    else:
        print(f'New aperture radius={aperture_radius}')

    for star in stars:
        aper = CircularAperture((star['xcentroid'], star['ycentroid']),
            aperture_radius)

        aperstats = ApertureStats(img_data, aper, mask=img_mask, 
            error=img_error)

        apertable.add_row([aperstats.fwhm.value, aperstats.semimajor_sigma.value,
            aperstats.semiminor_sigma.value, aperstats.orientation.value,
            aperstats.eccentricity, aperstats.sum/aperstats.sum_err,
            aperstats.sum, aperstats.sum_err, aperstats.xcentroid,
            aperstats.ycentroid, np.sqrt(aperstats.covar_sigx2.value),
            np.sqrt(aperstats.covar_sigy2.value)])
    
    return(apertable)


def generate_epsf(img_file, x, y, size=11, oversampling=2, maxiters=11,
    log=None):
    # Construct stars table from bright
    stars_tbl = Table()
    stars_tbl['x'] = x
    stars_tbl['y'] = y

    img_hdu = fits.open(img_file)
    ndimage = NDData(data=img_hdu['SCI'].data)

    stars = extract_stars(ndimage, stars_tbl, size=size)

    if log:
        log.info(f'Extracted {len(stars)} stars.  Building EPSF...')
    else:
        print(f'Extracted {len(stars)} stars.  Building EPSF...')

    epsf_builder = EPSFBuilder(oversampling=oversampling,
        maxiters=maxiters, progress_bar=True, smoothing_kernel='quadratic',
        sigma_clip=SigmaClip(sigma=5, sigma_lower=5, sigma_upper=5, 
            maxiters=20, cenfunc='median', stdfunc='std', grow=False))

    epsf, fitted_stars = epsf_builder(stars)

    return(epsf)

def extract_fwhm_from_epsf(epsf, fwhm_init):

    # Get the raw data for the FWHM and size in x and y
    data = epsf.data
    x = np.arange(data.shape[0])
    y = np.arange(data.shape[1])
    xx, yy = np.meshgrid(x, y) 

    # Fit to Moffat2D model in astropy, initial guess is amplitude,
    # centroid in x and y, core width of
    # Moffat model and power index scaling of model
    p_init = functional_models.Moffat2D(amplitude=0.5, x_0=data.shape[0]/2.,
        y_0=data.shape[1]/2., gamma=fwhm_init, alpha=1.)
    
    fit_p = fitting.LevMarLSQFitter()

    # Fit functional model to the data
    p = fit_p(p_init, xx, yy, data)

    # Extract and round FWHM
    fwhm = float('%.4f'%p.fwhm)

    return(p.fwhm)

def run_photometry(img_file, epsf, fwhm, threshold, shape, stars):

    img_hdu = fits.open(img_file)
    image = img_hdu['SCI'].data
    ndimage = NDData(data=img_hdu['SCI'].data)
    mask = img_hdu['MASK'].data.astype(bool)
    error = img_hdu['ERROR'].data

    psf = copy.copy(epsf)

    # Generate initial guesses for star centroid and flux from aperture table
    stars_tbl = Table()
    stars_tbl['x_0'] = stars['xcentroid']
    stars_tbl['y_0'] = stars['ycentroid']
    stars_tbl['flux_0'] = stars['flux_best']
    stars_tbl['local_bkg'] = np.array([0.]*len(stars))

    photometry = PSFPhotometry(psf_model=psf, fit_shape=(shape,shape),
        aperture_radius=int(shape*1.5), progress_bar=True)

    result_tab = photometry(image, mask=mask, error=error,
        init_params=stars_tbl)

    # Also generate a residual image for quality control
    residual_image = photometry.make_residual_image(image, 
        psf_shape=(shape, shape), include_localbkg=True)

    # Mask results table for sources with bad flux error, fit flux, or centroid
    mask = ~np.isnan(result_tab['flux_err'])
    result_tab = result_tab[mask]
    mask = result_tab['flux_fit'] > 0.
    result_tab = result_tab[mask]
    mask = ~np.isnan(result_tab['x_err'])
    result_tab = result_tab[mask]
    mask = ~np.isnan(result_tab['y_err'])
    result_tab = result_tab[mask]

    return(result_tab, residual_image)

# Identifies stars and gets aperture photometry and statistics using 
# photutils.aperture.ApertureStats
def get_star_catalog(img_data, img_mask, img_error, fwhm_init=5.0, 
    threshold=50.0, log=None):

    # Construct the finder with input FWHM and threshold
    daofind = DAOStarFinder(fwhm=fwhm_init, threshold=threshold,
        exclude_border=True, min_separation=fwhm_init)

    # Get initial set of stars in image with iraffind
    if log:
        log.info('Finding stars...')
    else:
        print('Finding stars...')

    # Do the finding...
    stars = daofind(img_data, mask=img_mask)

    # Ignore stars whose peak is below the background-subtracted level
    mask = stars['peak'] > 0.
    stars = stars[mask]

    stars.sort('flux')

    # Extract the aperture stats from each star and append to the output catalog
    stats = extract_aperture_stats(img_data, img_mask, img_error, stars, 
        aperture_radius=2.5*fwhm_init, log=log)
    stars = hstack([stars, stats])
    
    return(stars)

def do_phot(img_file,
    fwhm_scale_psf=4.5, oversampling=1,
    star_param={'snthresh_psf': 20.0, 'fwhm_init': 8.0, 'snthresh_final': 10.0},
    save_psf_img=False,
    save_residual_hdu=False,
    log=None):

    img_hdu = fits.open(img_file)

    # Get image statistics
    data = img_hdu['SCI'].data
    mask = img_hdu['MASK'].data.astype(bool)
    error = img_hdu['ERROR'].data

    # Get sky statistics
    mean, median, std_sky = sigma_clipped_stats(data[~mask], sigma=5.0,
        maxiters=21, grow=1)

    # Threshold for constructing star catalog
    threshold = std_sky*star_param['snthresh_final']/2.0

    stars = get_star_catalog(data, mask, error, 
        fwhm_init=star_param['fwhm_init'], threshold=threshold, log=log)

    metadata={'SKYADU':median, 'SKYSIG': std_sky}

    if log:
        log.info(f'Found {len(stars)} stars')
    else:
        print(f'Found {len(stars)} stars')

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
        log.info(f'Masked to {len(fwhm_stars)} stars based on sharpness, roundness, FWHM')
    else:
        print(f'Masked to {len(fwhm_stars)} stars based on sharpness, roundness, FWHM')

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
        size=size, oversampling=oversampling, maxiters=11, log=log)

    fwhm = extract_fwhm_from_epsf(epsf, fwhm*oversampling)
    # Scale by oversampling
    fwhm = fwhm/oversampling
    fwhm = float('%.4f'%fwhm)

    if log:
        log.info(f'FWHM={fwhm} pixels')
    else:
        print(f'FWHM={fwhm} pixels')

    metadata['FWHM']=fwhm

    # Also update img_hdu['PRIMARY'] with FWHM
    img_hdu['PRIMARY'].header['FWHM']=(fwhm, 
        'Full-width at half-maximum [pixels]')

    mask = (stars['signal_to_noise'] > star_param['snthresh_final'])
    final_stars = stars[mask]

    photometry, residual_image = run_photometry(img_file, epsf, fwhm, 
        star_param['snthresh_final'], size, final_stars)

    # Format final stars table and add as APPPHOT to img_hdu
    w = WCS(img_hdu['SCI'].header)
    coords = w.pixel_to_world(final_stars['Xpos'], final_stars['Ypos'])
    final_stars['flux'] = final_stars['flux_best']
    final_stars['flux_err'] = final_stars['flux_best_err']
    final_stars['SN'] = final_stars['flux']/final_stars['flux_best_err']
    final_stars['mag'] = -2.5*np.log10(final_stars['flux'])
    final_stars['mag_err'] = 2.5/np.log(10) * 1./final_stars['SN']
    final_stars['RA'] = [c.ra.degree for c in coords]
    final_stars['Dec'] = [c.dec.degree for c in coords]
    final_stars['sky'] = np.array([0.]*len(final_stars))
    final_stars['FWHM'] = final_stars['fwhm']
    # Sort from brightest to faintest
    final_stars.sort('flux', reverse=True)

    # Get final list of column names we want in catalog in order and the number
    # of significant figures they should all have
    colnames=['Xpos','Ypos','Xpos_err','Ypos_err','mag','mag_err','flux',
        'flux_err','SN','FWHM','sky','RA','Dec']
    sigfig=[4,4,4,4,4,4,4,4,4,4,4,7,7]

    final_stars = final_stars[*colnames]
    for col,sig in zip(colnames, sigfig):
        final_stars[col] = np.array([float(f'%.{sig}f'%val) 
            for val in final_stars[col].data])

    # Create a new HDU for the photometry table
    newhdu = fits.BinTableHDU(final_stars)
    newhdu.header.update(metadata)
    newhdu.name='APPPHOT'
    if newhdu.name in [h.name for h in img_hdu]:
        img_hdu[newhdu.name] = newhdu
    else:
        img_hdu.append(newhdu)

    # Record final number of objects to write out
    metadata['NOBJECT']=len(photometry)

    if log:
        log.info(f'Final catalog is {len(photometry)} stars')
    else:
        print(f'Final catalog is {len(photometry)} stars')

    # Get RA/Dec from final positions
    w = WCS(img_hdu['SCI'].header)
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
    photometry.rename_column('local_bkg','sky')
    # Sort from brightest to faintest
    photometry.sort('flux', reverse=True)

    # Get final list of column names we want in catalog in order and the number
    # of significant figures they should all have
    colnames=['Xpos','Ypos','Xpos_err','Ypos_err','mag','mag_err','flux',
        'flux_err','SN','FWHM','sky','RA','Dec']
    sigfig=[4,4,4,4,4,4,4,4,4,4,4,7,7]

    photometry = photometry[*colnames]
    for col,sig in zip(colnames, sigfig):
        photometry[col] = np.array([float(f'%.{sig}f'%val) 
            for val in photometry[col].data])

    # Update primary header with metadata
    img_hdu['PRIMARY'].header.update(metadata)

    # Create a new HDU for the photometry table
    newhdu = fits.BinTableHDU(photometry)
    newhdu.header.update(metadata)
    newhdu.name='PSFPHOT'
    if newhdu.name in [h.name for h in img_hdu]:
        img_hdu[newhdu.name] = newhdu
    else:
        img_hdu.append(newhdu)

    # Add PSF stars to img_hdu
    bright.sort('flux')
    newhdu = fits.BinTableHDU(bright)
    newhdu.name='PSFSTARS'
    if newhdu.name in [h.name for h in img_hdu]:
        img_hdu[newhdu.name] = newhdu
    else:
        img_hdu.append(newhdu)

    # Write out a EPSF quality image and add raw data to img_hdu
    # Useful for validating the photometry methods
    if save_psf_img:
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        outname = img_file.replace('.fits','.epsf.png')
        plt.savefig(outname)
        plt.clf()

    # Add the residual image to img_hdu
    # Can be used to validate the quality of PSF fitting to point sources
    if save_residual_hdu:
        newhdu = fits.ImageHDU(residual_image)
        newhdu.name = 'RESIDUAL'
        if newhdu.name in [h.name for h in img_hdu]:
            img_hdu[newhdu.name] = newhdu
        else:
            img_hdu.append(newhdu)

    # Always add the PSF to the image as an extension
    newhdu = fits.ImageHDU(epsf.data)
    newhdu.name='PSF'
    if newhdu.name in [h.name for h in img_hdu]:
        img_hdu[newhdu.name] = newhdu
    else:
        img_hdu.append(newhdu)

    # Finally, write out img_hdu to save all data
    img_hdu.writeto(img_file, overwrite=True)

def photloop(stack, phot_sn_min=3.0, phot_sn_max=40.0, fwhm_init=5.0, log=None):
    signal_to_noise = phot_sn_max

    epsf = None ; fwhm = None
    while signal_to_noise > phot_sn_min:

        if log: log.info(f'Trying photometry with final S/N={signal_to_noise}')
        star_param = {'snthresh_psf': signal_to_noise*2.0,
                      'fwhm_init': fwhm_init,
                      'snthresh_final': signal_to_noise}
        try:
            do_phot(stack, star_param=star_param) 
        except:
            signal_to_noise = signal_to_noise / 2.0
            continue
        break

if __name__=="__main__":

    # Pick image you want to test on
    test = '/Users/ckilpatrick/Dropbox/Data/POTPyRI/test/GMOS/red/sGRB240615A-GRB.i.ut240618.12.22.stk.fits'

    # For testing - give different names so the module doesn't consider them
    # global variables that override the values in methods
    f = 50.0
    s={'snthresh_psf': f*2.0, 'fwhm_init': 5.0, 'snthresh_final': f}

    # Run methods
    do_phot(test, star_param=s)
    #table = run_sextractor(test)
