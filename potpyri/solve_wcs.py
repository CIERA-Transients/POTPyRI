"Python script for WCS solution."
"Authors: Kerry Paterson, Charlie Kilpatrick"

# Last updated 09/21/2024
__version__ = "2.0"

import numpy as np
import subprocess
import os
import progressbar

import matplotlib.pyplot as plt

from photutils.aperture import ApertureStats
from photutils.aperture import CircularAperture

from astroquery.vizier import Vizier

import astropy.units as u
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table, Column
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points

# Internal dependency
from potpyri import utilities
from potpyri import photometry

#turn Astropy warnings off
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

def run_sextractor(input_file, cat_name, tel, sex_config_dir='./Config', 
    log=None):

    if not os.path.exists(os.path.join(sex_config_dir,'config')):
        e = 'ERROR: Could not find source extractor config file!'
        if log:
            log.error(e)
        raise Exception(e)
    elif not os.path.exists(os.path.join(sex_config_dir,'params')):
        e = 'ERROR: Could not find source extractor param file!'
        if log:
            log.error(e)
        raise Exception(e)

    config = os.path.join(sex_config_dir,'config')
    param_file = os.path.join(sex_config_dir,'params')

    cmd=['sex','-c',config,input_file,'-CATALOG_NAME',cat_name]
    subprocess.call(cmd)

    params = []
    with open(param_file) as f:
        for line in f:
            params.append(line.split()[0].strip())

    if os.path.exists(cat_name):
        table = ascii.read(cat_name, names=params, comment='#')
        if log:
            log.info('SExtractor succesfully run, %d sources extracted.'%len(table))
        return(table)
    else:
        e = 'ERROR: source extractor did not successfully produce {0}'
        if log:
            log.error(e)
        raise Exception(e.format(cat_name))


def get_gaia_catalog(input_file, log=None):

    hdu = fits.open(input_file)
    data = hdu[0].data
    header = hdu[0].header
    hkeys = list(header.keys())

    check_pairs = [('RA','DEC'),('CRVAL1','CRVAL2'),('OBJCTRA','OBJCTDEC')]
    coord = None

    for pair in check_pairs:
        if pair[0] in header.keys() and pair[1] in header.keys():
            ra = header[pair[0]]
            dec = header[pair[1]]
            coord = utilities.parse_coord(ra, dec)
            if coord:
                break

    if not coord:
        if log: log.error(f'Could not parse RA/DEC from header of {file}')
        return(False)

    vizier = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'Plx', 'PSS','PM'])
    vizier.ROW_LIMIT = -1

    cat = vizier.query_region(coord, width=20 * u.arcmin, 
        catalog='I/355/gaiadr3')

    if len(cat)>0:
        cat = cat[0]
    else:
        return(None)

    mask = (cat['PSS'] > 0.99) & (cat['Plx'] < 20) & (cat['PM']<10)
    cat = cat[mask]

    return(cat)

def clean_up_astrometry(directory, file, exten):
    filelist = [file.replace(exten,'.axy'),
                file.replace(exten,'.corr'),
                file.replace(exten,'-indx.xyls'),
                file.replace(exten,'.match'),
                file.replace(exten,'.rdls'),
                file.replace(exten,'.solved'),
                file.replace(exten,'.wcs')]

    filelist = [os.path.join(directory, f) for f in filelist]

    for f in filelist:
        if os.path.exists(f):
            os.remove(f)

def solve_astrometry(file, tel, radius=0.5, replace=True, 
    shift_only=False, log=None):

    # Starting solve, print file and directory for reference
    fullfile = os.path.abspath(file)
    directory = os.path.dirname(file)

    if log: 
        log.info(f'Trying to solve file: {file}')
    else:
        print(f'Trying to solve file: {file}')

    if not os.path.exists(fullfile):
        return(False)

    hdu = fits.open(fullfile)
    data = hdu[0].data
    header = hdu[0].header
    hkeys = list(header.keys())

    exten = '.'+file.split('.')[-1]
    if not replace:
        if os.path.exists(fullfile.replace(exten,'.solved.fits')):
            if log: 
                log.info(f'SUCCESS: solved {fullfile}')
            else:
                print(f'SUCCESS: solved {fullfile}')
            return(True)

    exten = '.'+file.split('.')[-1]

    check_pairs = [('RA','DEC'),('CRVAL1','CRVAL2'),('OBJCTRA','OBJCTDEC')]
    coord = None

    for pair in check_pairs:
        if pair[0] in header.keys() and pair[1] in header.keys():
            ra = header[pair[0]]
            dec = header[pair[1]]
            coord = utilities.parse_coord(ra, dec)
            if coord:
                break

    if not coord:
        if log: log.error(f'Could not parse RA/DEC from header of {file}')
        return(False)

    if 'solved.fits' in fullfile:
        newfile = 'tmp.fits'
    else:
        newfile = fullfile.replace(exten,'.solved.fits')

    # Handle pixel scale guess
    scale = tel.pixscale
    scale_high = float('%.4f'%(scale * 1.2))
    scale_low = float('%.4f'%(scale * 0.8))

    # Get RA and Dec
    ra = float('%.6f'%coord.ra.degree)
    dec = float('%.6f'%coord.dec.degree)

    cmd = 'solve-field'
    args = '--scale-units arcsecperpix '
    args += f'--scale-low {scale_low} --scale-high {scale_high} '
    args += f'--ra {ra} --dec {dec} '
    args += f' --radius {radius} --no-plots -T '
    args += f'--overwrite -N {newfile} --dir {directory} '

    # Test for --use-source-extractor flag
    p = subprocess.run(['solve-field','-h'],capture_output=True)
    data = p.stdout.decode().lower()

    if '--use-source-extractor' in data:
        args += '--use-source-extractor '
    elif '--use-sextractor' in data:
        args += '--use-sextractor '

    extra_opts = '--downsample 2 --no-verify --odds-to-tune-up 1e4 --objs 15'

    tries = 1
    good = False
    while tries < 4 and not good:
        input_args = args + extra_opts

        if log: 
            log.info(f'Try #{tries} with astrometry.net...')
        else:
            print(f'Try #{tries} with astrometry.net...')

        process = [cmd,fullfile]+input_args.split()

        if log: 
            log.info(' '.join(process))
        else:
            print(' '.join(process))

        FNULL = open(os.devnull, 'w')

        p = subprocess.Popen(process, stdout=FNULL, stderr=subprocess.STDOUT)
        try:
            p.wait(90)
        except subprocess.TimeoutExpired:
            p.kill()

        if os.path.exists(newfile):
            good = True
        else:
            tries += 1
            if tries==2:
                extra_opts='--objs 15'
            elif tries==3:
                extra_opts=''


    file_exists=os.path.exists(newfile)

    if log: 
        log.info(f'{newfile} exists: {file_exists}')
    else:
        print(f'{newfile} exists: {file_exists}')

    if os.path.exists(newfile):

        clean_up_astrometry(directory, file, exten)
        if log: 
            log.info(f'SUCCESS: solved {fullfile}')
        else:
            print(f'SUCCESS: solved {fullfile}')

        if replace or newfile=='tmp.fits':
            output_file = fullfile

            # astrometry.net replaces extra extensions, so instead of renaming
            # copy the new header into the old file
            newhdu = fits.open(newfile)
            hdu = fits.open(output_file)

            if shift_only:
                hdu[0].header['CRPIX1'] = newhdu[0].header['CRPIX1']
                hdu[0].header['CRPIX2'] = newhdu[0].header['CRPIX2']
                hdu[0].header['CRVAL1'] = newhdu[0].header['CRVAL1']
                hdu[0].header['CRVAL2'] = newhdu[0].header['CRVAL2']
            else:
                hdu[0].header = newhdu[0].header
                
            hdu.writeto(output_file, overwrite=True, output_verify='silentfix')
            os.remove(newfile)
        else:
            output_file = fullfile.replace(exten,'.solved.fits')

        hdu = fits.open(output_file)
        for i,h in enumerate(hdu):

            if 'COMMENT' in hdu[i].header.keys():
                del hdu[i].header['COMMENT']
            if 'HISTORY' in hdu[i].header.keys():
                del hdu[i].header['HISTORY']

            # Delete other WCS header keys
            wcs_key = ['CSYER1','CSYER2','CRDER1','CRDER2','CD1_1','CD1_2',
                'CD2_1','CD2_2','CRPIX1','CRPIX2','CUNIT1','CUNIT2','EQUINOX',
                'RADESYS','CNAME1','CNAME2','CTYPE1','CTYPE2','WCSNAME',
                'CRVAL1','CRVAL2']

            for key in [w + 'C' for w in wcs_key]:
                if key in list(hdu[i].header.keys()):
                    del hdu[i].header[key]

            for key in [w + 'S' for w in wcs_key]:
                if key in list(hdu[i].header.keys()):
                    del hdu[i].header[key]

            # Delete old WCS keywords starting with '_'
            for key in list(hdu[i].header.keys()):
                if key.startswith('_'):
                    del hdu[i].header[key]

        try:
            hdu.writeto(output_file, overwrite=True, output_verify='silentfix')
        except TypeError:
            if log: 
                log.error(f'FAILURE: could not save file {fullfile}')
            else:
                print(f'FAILURE: could not save file {fullfile}')
            return(False)

        return(True)

    else:
        if log: 
            log.error(f'FAILURE: did not solve {fullfile}')
        else:
            print(f'FAILURE: did not solve {fullfile}')
        clean_up_astrometry(directory, file, exten)
        return(False)

def align_to_gaia(file, tel, radius=0.5, log=None):

    cat = get_gaia_catalog(file, log=log)

    if cat is None or len(cat)==0:
        if log:
            log.error(f'Could not get Gaia catalog or no stars for {file}')
        else:
            print(f'Could not get Gaia catalog or no stars for {file}')
        return(False)

    if log: 
        log.info(f'Got {len(cat)} Gaia DR3 alignment stars')
    else:
        print(f'Got {len(cat)} Gaia DR3 alignment stars')

    # Estimate sources that land within the image using WCS from header
    hdu = fits.open(file)
    w = WCS(hdu[0].header)
    naxis1 = hdu[0].header['NAXIS1']
    naxis2 = hdu[0].header['NAXIS2']

    # Get pixel coordinates of all stars in image
    coords = SkyCoord(cat['RA_ICRS'], cat['DE_ICRS'], unit=(u.deg, u.deg))
    x_pix, y_pix = w.world_to_pixel(coords)

    # Mask to stars that are nominally in image
    mask = (x_pix > 0) & (x_pix < naxis1) & (y_pix > 0) & (y_pix < naxis2)
    x_pix = x_pix[mask] ; y_pix = y_pix[mask]
    # Also mask catalog
    cat = cat[mask]

    if cat is None or len(cat)==0:
        if log: 
            log.error(f'No stars found in {file}')
        else:
            print(f'No stars found in {file}')
        return(False)

    if log: 
        log.info(f'Found {len(cat)} stars in the image')
    else:
        print(f'Found {len(cat)} stars in the image')

    cat = photometry.run_sextractor(file, log=log)

    if log: 
        log.info(f'Centroided on {len(cat)} stars')
    else:
        print(f'Centroided on {len(cat)} stars')

    cat_coords = SkyCoord(cat['ALPHA_J2000'], cat['DELTA_J2000'])
    idx, d2d, d3d = match_coordinates_sky(cat_coords, coords)

    sep_mask = d2d < 1.0 * u.arcsec
    cat_coords_match = cat[sep_mask]
    sky_coords_match = coords[idx[sep_mask]]

    cat_coords_match.add_column(sky_coords_match['RA_ICRS'])
    cat_coords_match.add_column(sky_coords_match['DE_ICRS'])

    if cat_coords_match is None or len(cat_coords_match)==0:
        if log:
            log.error(f'Could not centroid on any stars in {file}')
        else:
            print(f'Could not centroid on any stars in {file}')
        return(False)
    else:
        if log:
            log.info(f'Number of matched coordinates is {len(cat_coords_match)}')
        else:
            print(f'Number of matched coordinates is {len(cat_coords_match)}')

    # Get central coordinate of image based on centroiding
    central_pix = (float(naxis1)/2., float(naxis2)/2.)
    central_coo = w.pixel_to_world(*central_pix)

    # Bootstrap WCS solution and get center pixel
    ra_cent = [] ; de_cent = []
    coords = SkyCoord(cat_coords_match['RA_ICRS'], cat_coords_match['DE_ICRS'], unit=(u.deg, u.deg))

    sip_degree = None
    if len(coords)>20:
        sip_degree = 2
    if len(coords)>60:
        sip_degree = 4
    if len(coords)>200:
        sip_degree = 6

    if log:
        log.info(f'Solving with SIP degree = {sip_degree}')
    else:
        print(f'Solving with SIP degree = {sip_degree}')

    bar = progressbar.ProgressBar(maxval=500)
    bar.start()

    idx = np.arange(len(cat_coords_match))

    for i in np.arange(500):
        bar.update(i+1)
        size = 100

        idxs = np.random.choice(idx, size=size, replace=False)

        xy = (cat_coords_match['X_IMAGE'][idxs], cat_coords_match['Y_IMAGE'][idxs])
        c = coords[idxs]

        new_wcs = fit_wcs_from_points(xy, c, projection=w)

        cent_coord = new_wcs.pixel_to_world(*central_pix)
        ra_cent.append(cent_coord.ra.degree)
        de_cent.append(cent_coord.dec.degree)

    bar.finish()

    ra_cent = np.array(ra_cent)
    de_cent = np.array(de_cent)

    # Estimate systematic precision of new WCS
    mean_ra, med_ra, std_ra = sigma_clipped_stats(ra_cent, maxiters=11,
        sigma=3)
    mean_de, med_de, std_de = sigma_clipped_stats(de_cent, maxiters=11,
        sigma=3)

    mask = (np.abs(ra_cent-med_ra)<3*std_ra) & (np.abs(de_cent-med_de)<3*std_de)

    fig, ax = plt.subplots()
    ax.scatter(ra_cent[mask], de_cent[mask])
    plt.savefig(file.replace('.fits','_centroids.png'))

    # Generate "final" WCS
    xy = (cat_coords_match['X_IMAGE'], cat_coords_match['Y_IMAGE'])
    coords = SkyCoord(cat_coords_match['RA_ICRS'], cat_coords_match['DE_ICRS'], 
        unit=(u.deg, u.deg))
    proj_point = SkyCoord(med_ra, med_de, unit='deg')
    new_wcs = fit_wcs_from_points(xy, coords, projection=w)

    print(new_wcs)

    header = new_wcs.to_header()

    # Delete all old WCS keys
    #delkeys = ['WCSNAME','CUNIT1','CUNIT2','CTYPE1','CTYPE2','CRPIX1','CRPIX2',
    #    'CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2','RADECSYS',
    #    'PC1_1','PC1_2','PC2_1','PC2_2','LONPOLE','LATPOLE','CDELT1','CDELT2']
    #while True:
    #    found_key = False
    #    for key in hdu[0].header.keys():
    #        if any([k in key for k in delkeys]):
    #            found_key=True
    #            del hdu[0].header[key]
    #    if not found_key:
    #        break

    hdu[0].header.update(header)

    ra_disp = std_ra*3600.0*np.cos(np.pi/180.0 * mean_de)
    de_disp = std_de*3600.0

    ra_disp = float('%.6f'%ra_disp)
    de_disp = float('%.6f'%de_disp)

    if log:
        log.info(f'Dispersion in R.A.={ra_disp} arcsec, Dispersion in Decl.={de_disp} arcsec')
    else:
        print(f'Dispersion in R.A.={ra_disp} arcsec, Dispersion in Decl.={de_disp} arcsec')

    # Add header variables for dispersion in WCS solution
    hdu[0].header['RADISP']=(ra_disp, 'Dispersion in R.A. of WCS solution [Arcsec]')
    hdu[0].header['DEDISP']=(de_disp, 'Dispersion in Decl. of WCS solution [Arcsec]')

    hdu.writeto(file, overwrite=True, output_verify='silentfix')

    return(True)

if __name__ == "__main__":
    file='/Users/ckilpatrick/FRB241005/red/workspace/FRB20241005_r.r.ut241007.2.11.1412743557.fits'

    import importlib
    global tel
    tel = importlib.import_module('instruments.BINOSPEC')
    tel = tel.BINOSPEC()

    solve_astrometry(file, tel, radius=0.5)
    align_to_gaia(file, tel, radius=0.5, log=None)
