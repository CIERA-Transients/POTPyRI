"""Aperture and PSF photometry for pipeline stacks.

Uses Source Extractor for initial catalogs, photutils for ePSF building
and PSF fitting. Writes APPPHOT, PSFPHOT, PSFSTARS, PSF, and RESIDUAL
extensions to the stack FITS. Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from potpyri._version import __version__

from photutils.aperture import ApertureStats
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder
from photutils.psf import extract_stars
try:
    # photutils >= 2.x
    from photutils.psf import PSFPhotometry  # type: ignore
    _HAS_PSF_PHOTOMETRY = True
except Exception:  # pragma: no cover (depends on photutils version)
    # photutils 1.x fallback
    from photutils.psf import BasicPSFPhotometry, DAOGroup  # type: ignore
    from photutils.background import MMMBackground  # type: ignore
    _HAS_PSF_PHOTOMETRY = False

from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
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
import traceback

warnings.filterwarnings('ignore')


class PhotometryError(RuntimeError):
    """PSF/aperture photometry did not complete; stack is missing APPPHOT/PSFPHOT."""

    pass


def _normalize_daofind_catalog(stars):
    """Ensure DAOStarFinder columns use POTPyRI's expected names.

    photutils 3.x renamed ``xcentroid``/``ycentroid`` to ``x_centroid``/
    ``y_centroid``; older photutils used the legacy names. This keeps downstream
    code working across versions.
    """
    if stars is None:
        return None
    if len(stars) == 0:
        return stars
    # photutils 3.x wraps DAOStarFinder output in ``DeprecatedColumnQTable``:
    # ``stars['xcentroid']`` resolves to ``x_centroid``, but
    # ``row['xcentroid']`` in ``extract_aperture_stats`` can still raise
    # ``KeyError`` because rows index ``.columns`` without that translation.
    # Materialize a plain ``Table`` so legacy centroid columns are real columns.
    if getattr(stars, 'deprecation_map', None):
        meta = dict(stars.meta) if stars.meta else {}
        stars = Table(
            {name: stars[name] for name in list(stars.colnames)},
            meta=meta,
            copy=True,
        )
    if 'xcentroid' in stars.colnames:
        return stars
    # photutils 3.x naming
    if 'x_centroid' in stars.colnames and 'y_centroid' in stars.colnames:
        stars['xcentroid'] = np.asarray(
            u.Quantity(stars['x_centroid'], copy=False).value, dtype=np.float64)
        stars['ycentroid'] = np.asarray(
            u.Quantity(stars['y_centroid'], copy=False).value, dtype=np.float64)
        return stars
    # Additional centroid aliases seen across photutils / table producers.
    # We prefer to coerce to plain float columns because downstream code expects
    # numeric pixel positions (not Quantity mixins).
    x_aliases = ['xcenter', 'x_center', 'x', 'x_0', 'xpos', 'Xpos']
    y_aliases = ['ycenter', 'y_center', 'y', 'y_0', 'ypos', 'Ypos']
    x_col = next((c for c in x_aliases if c in stars.colnames), None)
    y_col = next((c for c in y_aliases if c in stars.colnames), None)
    if x_col and y_col:
        stars['xcentroid'] = np.asarray(
            u.Quantity(stars[x_col], copy=False).value, dtype=np.float64)
        stars['ycentroid'] = np.asarray(
            u.Quantity(stars[y_col], copy=False).value, dtype=np.float64)
        return stars
    raise PhotometryError(
        'DAOStarFinder output is missing centroid columns '
        '(expected xcentroid/ycentroid, x_centroid/y_centroid, or a known alias). '
        f'Got columns: {list(stars.colnames)}'
    )


def create_conv(outfile):
    """Write a 3x3 CONV NORM filter file for Source Extractor.

    Parameters
    ----------
    outfile : str
        Path to write the convolution kernel file.

    Returns
    -------
    None
    """
    with open(outfile, 'w') as f:
        f.write('CONV NORM \n')
        f.write('1 2 1 \n')
        f.write('2 4 2 \n')
        f.write('1 2 1 \n')

def create_params(outfile):
    """Write Source Extractor parameter file (NUMBER, X_IMAGE, FWHM_IMAGE, etc.).

    Parameters
    ----------
    outfile : str
        Path to write the parameter file.

    Returns
    -------
    None
    """
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

def run_sextractor(img_file, log=None, sextractor_path=None):
    """Run Source Extractor on SCI extension; return catalog table or None.

    Parameters
    ----------
    img_file : str
        Path to FITS with SCI extension (and SATURATE in header).
    log : ColoredLogger, optional
        Logger for progress.
    sextractor_path : str, optional
        Path to Source Extractor binary. If None, uses 'sex' from PATH.

    Returns
    -------
    astropy.table.Table or None
        Source Extractor catalog, or None if run failed.
    """
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

    sex_cmd = (sextractor_path if sextractor_path else 'sex').strip()
    if log:
        log.info(f'Running source extractor on {img_file}')
    else:
        print(f'Running source extractor on {img_file}')

    cmd = f'{sex_cmd} {tmpfile} -CATALOG_NAME {catfile} -CATALOG_TYPE ASCII_HEAD '
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
    """Compute aperture flux and error for a star table; return table with added columns.

    Parameters
    ----------
    img_data : np.ndarray
        Science image data.
    img_mask : np.ndarray
        Boolean or integer mask (True/masked excluded).
    img_error : np.ndarray
        Per-pixel error array.
    stars : astropy.table.Table
        Table with xcentroid, ycentroid (modified in place with refined centroids).
    aperture_radius : float, optional
        Aperture radius in pixels. Default is 10.0.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    astropy.table.Table
        Table with fwhm, flux_best, flux_best_err, Xpos, Ypos, etc.
    """

    apertable = Table([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],
        [0.],[0.]], 
        names=('fwhm','semimajor_sigma','semiminor_sigma',
        'orientation','eccentricity','signal_to_noise',
        'flux_best', 'flux_best_err','Xpos','Ypos',
        'Xpos_err','Ypos_err')).copy()[:0]

    if len(stars) == 0:
        return apertable

    # Be defensive: some call sites may pass tables that didn't come from
    # get_star_catalog() (or DAOStarFinder output with new column names).
    stars = _normalize_daofind_catalog(stars)

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

        covx = np.maximum(aperstats.covar_sigx2.value, 0.0)
        covy = np.maximum(aperstats.covar_sigy2.value, 0.0)
        apertable.add_row([aperstats.fwhm.value, aperstats.semimajor_sigma.value,
            aperstats.semiminor_sigma.value, aperstats.orientation.value,
            aperstats.eccentricity, aperstats.sum/aperstats.sum_err,
            aperstats.sum, aperstats.sum_err, aperstats.xcentroid,
            aperstats.ycentroid, np.sqrt(covx),
            np.sqrt(covy)])
    
    return(apertable)


def generate_epsf(img_file, x, y, size=11, oversampling=2, maxiters=11,
    log=None):
    """Build an ePSF from cutouts at (x, y) positions in the image.

    Parameters
    ----------
    img_file : str
        Path to FITS with SCI extension.
    x, y : array-like
        Star x,y positions (pixels).
    size : int, optional
        Cutout size in pixels. Default is 11.
    oversampling : int, optional
        EPSFBuilder oversampling. Default is 2.
    maxiters : int, optional
        EPSFBuilder max iterations. Default is 11.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    photutils.psf.EPSFModel
        Fitted ePSF model.
    """
    # Construct stars table from bright
    stars_tbl = Table()
    stars_tbl['x'] = x
    stars_tbl['y'] = y

    with fits.open(img_file) as img_hdu:
        ndimage = NDData(data=img_hdu['SCI'].data)

    stars = extract_stars(ndimage, stars_tbl, size=size)

    if log:
        log.info(f'Extracted {len(stars)} stars.  Building EPSF...')
    else:
        print(f'Extracted {len(stars)} stars.  Building EPSF...')

    epsf_builder = EPSFBuilder(oversampling=oversampling,
        maxiters=maxiters, progress_bar=False, smoothing_kernel='quadratic',
        sigma_clip=SigmaClip(sigma=5, sigma_lower=5, sigma_upper=5, 
            maxiters=20, cenfunc='median', stdfunc='std', grow=False))

    epsf, fitted_stars = epsf_builder(stars)

    return(epsf)

def extract_fwhm_from_epsf(epsf, fwhm_init):
    """Estimate FWHM from ePSF model (Moffat2D fit).

    Parameters
    ----------
    epsf : photutils.psf.EPSFModel
        ePSF model with .data array.
    fwhm_init : float
        Initial FWHM guess for fit (pixels).

    Returns
    -------
    float
        FWHM in pixels from fitted Moffat2D.
    """
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
    """Run PSF photometry and aperture photometry; append result tables to FITS.

    Parameters
    ----------
    img_file : str
        Path to FITS with SCI, MASK, ERROR extensions.
    epsf : photutils.psf.EPSFModel
        ePSF model for PSF fit.
    fwhm : float
        FWHM in pixels (for aperture radius).
    threshold : float
        Detection threshold (unused in this wrapper; for compatibility).
    shape : int
        Fit shape (pixels) for PSFPhotometry.
    stars : astropy.table.Table
        Star table with xcentroid, ycentroid, flux_best for initial guesses.

    Returns
    -------
    tuple
        (result_table, residual_image). result_table has flux_fit, x_fit, y_fit, etc.
    """
    with fits.open(img_file) as img_hdu:
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

    if _HAS_PSF_PHOTOMETRY:
        photometry = PSFPhotometry(psf_model=psf, fit_shape=(shape, shape),
            aperture_radius=int(shape*1.5), progress_bar=False)

        result_tab = photometry(image, mask=mask, error=error,
            init_params=stars_tbl)

        # Also generate a residual image for quality control
        residual_image = photometry.make_residual_image(
            image, psf_shape=(shape, shape), include_localbkg=True
        )
    else:
        # photutils 1.x: BasicPSFPhotometry (no per-pixel error support)
        group_maker = DAOGroup(crit_separation=max(float(fwhm), 2.0))
        bkg_estimator = MMMBackground()
        photometry = BasicPSFPhotometry(
            group_maker=group_maker,
            bkg_estimator=bkg_estimator,
            psf_model=psf,
            fitshape=(shape, shape),
            aperture_radius=int(shape * 1.5),
        )
        result_tab = photometry(image, mask=mask, init_guesses=stars_tbl, progress_bar=False)
        residual_image = photometry.get_residual_image()

        # Normalize output column names to match newer photutils expectations.
        if 'flux_err' not in result_tab.colnames and 'flux_unc' in result_tab.colnames:
            result_tab.rename_column('flux_unc', 'flux_err')
        if 'x_err' not in result_tab.colnames and 'x_unc' in result_tab.colnames:
            result_tab.rename_column('x_unc', 'x_err')
        if 'y_err' not in result_tab.colnames and 'y_unc' in result_tab.colnames:
            result_tab.rename_column('y_unc', 'y_err')

        # If uncertainties are not provided, estimate something finite so downstream
        # filtering doesn't drop all sources.
        if 'flux_err' not in result_tab.colnames:
            flux_fit = np.asarray(result_tab['flux_fit'], dtype=float)
            result_tab['flux_err'] = np.sqrt(np.maximum(np.abs(flux_fit), 1.0))
        if 'x_err' not in result_tab.colnames:
            result_tab['x_err'] = np.full(len(result_tab), 0.1, dtype=float)
        if 'y_err' not in result_tab.colnames:
            result_tab['y_err'] = np.full(len(result_tab), 0.1, dtype=float)

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
    """Detect stars with DAOStarFinder and compute aperture stats; return star table.

    Parameters
    ----------
    img_data : np.ndarray
        Science image data.
    img_mask : np.ndarray
        Boolean mask (True = masked).
    img_error : np.ndarray
        Per-pixel error array.
    fwhm_init : float, optional
        DAOStarFinder FWHM and aperture scale. Default is 5.0.
    threshold : float, optional
        DAOStarFinder detection threshold. Default is 50.0.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    astropy.table.Table
        Star table with centroid, flux, fwhm, flux_best, etc.
    """

    # Construct the finder with input FWHM and threshold.
    # photutils has changed DAOStarFinder kwargs across versions (e.g. older
    # releases do not accept min_separation).
    try:
        daofind = DAOStarFinder(
            fwhm=fwhm_init,
            threshold=threshold,
            exclude_border=True,
            min_separation=fwhm_init,
        )
    except TypeError:  # pragma: no cover (depends on photutils version)
        daofind = DAOStarFinder(
            fwhm=fwhm_init,
            threshold=threshold,
            exclude_border=True,
        )

    # Get initial set of stars in image with iraffind
    if log:
        log.info('Finding stars...')
    else:
        print('Finding stars...')

    # Do the finding...
    stars = daofind(img_data, mask=img_mask)

    if stars is None:
        raise PhotometryError(
            'DAOStarFinder returned no source table (None): no detections passed '
            'DAOStarFinder quality cuts, or the image is shallow / heavily masked. '
            'Try a lower detection threshold or different fwhm_init.'
        )

    stars = _normalize_daofind_catalog(stars)

    # Ignore stars whose peak is below the background-subtracted level
    mask = stars['peak'] > 0.
    stars = stars[mask]

    stars.sort('flux')

    # Limit the number of stars passed into the (expensive) ApertureStats loops.
    # We keep the brightest sources, which are the most useful for PSF building.
    max_sources = 2000
    if len(stars) > max_sources:
        stars = stars[-max_sources:]

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
    """Run full PSF and aperture photometry on a stack; write extensions to FITS.

    Builds ePSF, finds stars, runs PSF and aperture photometry, and appends
    APPPHOT, PSFPHOT, EPSF, and optional residual extensions to the FITS file.

    Parameters
    ----------
    img_file : str
        Path to stack FITS (SCI, MASK, ERROR).
    fwhm_scale_psf : float, optional
        Scale factor for PSF fit shape. Default is 4.5.
    oversampling : int, optional
        ePSF oversampling. Default is 1.
    star_param : dict, optional
        Keys: snthresh_psf, fwhm_init, snthresh_final. Defaults as shown.
    save_psf_img : bool, optional
        If True, save ePSF image extension. Default is False.
    save_residual_hdu : bool, optional
        If True, save residual image extension. Default is False.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
        img_file is updated in place with new extensions.
    """
    img_hdu = fits.open(img_file)
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

    if len(stars) == 0:
        raise PhotometryError(
            'Star detection returned zero sources after DAOStarFinder and peak>0 '
            'cut. Check image depth, MASK coverage, ERROR array, fwhm_init, or '
            'lower the photometry S/N threshold (phot_sn_max / photloop).'
        )

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

    if len(fwhm_stars) == 0:
        raise PhotometryError(
            'No stars left after sharpness/roundness/FWHM sigma-clipping for ePSF '
            'construction. Try a larger fwhm_init or lower S/N thresholds.'
        )
    
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

    if len(bright) == 0:
        raise PhotometryError(
            'No stars passed S/N cuts for ePSF building (bright catalog empty after '
            'PSF-star selection). Lower snthresh_psf / photometry S/N or inspect '
            'image quality and MASK.'
        )

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
    if len(final_stars) == 0:
        raise PhotometryError(
            f"No stars above snthresh_final={star_param['snthresh_final']} for PSF "
            'fitting. Lower the photometry S/N threshold (photloop phot_sn_*).'
        )

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

    final_stars = final_stars[colnames]
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
    flux_err = np.asarray(photometry['flux_err'], dtype=float)
    sn = np.full_like(flux_err, np.nan)
    np.divide(photometry['flux_fit'], flux_err, out=sn, where=flux_err > 0)
    photometry['SN'] = sn
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

    photometry = photometry[colnames]
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
    img_hdu.close()

def photloop(stack, phot_sn_min=3.0, phot_sn_max=40.0, fwhm_init=5.0, log=None):
    """Try PSF photometry with decreasing S/N threshold until do_phot succeeds.

    Parameters
    ----------
    stack : str
        Path to stack FITS (SCI, MASK, ERROR).
    phot_sn_min : float, optional
        Minimum S/N threshold to try. Default is 3.0.
    phot_sn_max : float, optional
        Initial S/N threshold. Default is 40.0.
    fwhm_init : float, optional
        Initial FWHM for star finding. Default is 5.0.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
        Stack FITS is updated in place when do_phot succeeds.

    Raises
    ------
    PhotometryError
        If every S/N attempt fails (no APPPHOT written).
    """
    if phot_sn_max <= phot_sn_min:
        raise PhotometryError(
            f'phot_sn_max ({phot_sn_max}) must be greater than phot_sn_min ({phot_sn_min}).'
        )

    signal_to_noise = phot_sn_max
    last_exc = None
    sn_tried = []

    while signal_to_noise > phot_sn_min:
        sn_tried.append(signal_to_noise)
        if log:
            log.info(
                f'Photometry: trying S/N threshold={signal_to_noise} '
                f'(fwhm_init={fwhm_init}, stack={stack})'
            )
        else:
            print(
                f'Photometry: trying S/N threshold={signal_to_noise} '
                f'(fwhm_init={fwhm_init})'
            )

        star_param = {'snthresh_psf': signal_to_noise * 2.0,
                      'fwhm_init': fwhm_init,
                      'snthresh_final': signal_to_noise}
        try:
            do_phot(stack, star_param=star_param, log=log)
        except Exception as e:
            last_exc = e
            if log:
                log.exception(
                    'Photometry attempt failed for S/N threshold=%s: %s',
                    signal_to_noise, e,
                )
            else:
                traceback.print_exc()
            signal_to_noise = signal_to_noise / 2.0
            continue

        if log:
            log.info(
                f'Photometry finished successfully at S/N threshold={signal_to_noise} '
                f'(APPPHOT written to {stack})'
            )
        else:
            print(f'Photometry succeeded at S/N threshold={signal_to_noise}')
        return

    msg = (
        f'Photometry failed for {stack!r}: no successful run after trying S/N '
        f'thresholds {sn_tried}. The stack has no APPPHOT extension, so downstream '
        f'zeropoint fitting will fail. Last error: {last_exc!r}. '
        f'Try lowering --phot-sn-min (currently {phot_sn_min}), raising '
        f'--phot-sn-max, adjusting --fwhm-init, or inspect SCI/MASK/ERROR and WCS.'
    )
    if log:
        log.error(msg)
    else:
        print(msg, flush=True)
    raise PhotometryError(msg) from last_exc
