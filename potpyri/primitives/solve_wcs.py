"""WCS solution using astrometry.net and Gaia DR3 alignment.

Provides coarse solution via solve-field and fine alignment by matching
detected sources to Gaia. Writes RADISP/DEDISP and updated WCS to FITS.
Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from potpyri._version import __version__

import numpy as np
import subprocess
import os
import progressbar
import matplotlib.pyplot as plt

from photutils.aperture import ApertureStats
from photutils.aperture import CircularAperture

import astropy.units as u
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Column
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points

# Internal dependencies
from potpyri.utils import catalogs
from potpyri.utils import utilities
from potpyri.primitives import photometry

#turn Astropy warnings off
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

def _log_fine_align(log, level, message):
    """Small helper to keep fine-alignment logging consistent."""
    if log:
        getattr(log, level)(message)
    else:
        print(message)

def _write_fine_align_fallback_header(hdu, file, reason, log=None, disp_arcsec=1.0):
    """Write fallback astrometric dispersion when fine catalog alignment fails.

    Preserves the coarse astrometry.net WCS while recording the failure in the
    header and logs.
    """
    _log_fine_align(log, 'warning',
        f'Fine alignment fallback for {file}: {reason}. '
        f'Keeping coarse astrometry.net WCS.')
    hdu[0].header['RADISP'] = (disp_arcsec, 'Dispersion in R.A. of WCS [Arcsec]')
    hdu[0].header['DEDISP'] = (disp_arcsec, 'Dispersion in Decl. of WCS [Arcsec]')
    hdu[0].header['ALGNFAIL'] = (reason[:32], 'Fine alignment fallback reason')
    hdu.writeto(file, overwrite=True, output_verify='silentfix')

def _validate_refined_wcs(w, header, tel, log=None):
    """Validate refined WCS against expected pixel scale and anisotropy.

    Returns
    -------
    tuple
        (is_valid, reason_string)
    """
    try:
        m = np.asarray(w.pixel_scale_matrix, dtype=float)
        _, s, _ = np.linalg.svd(m)
        # singular-value pixel scales in arcsec/pixel
        scales_arcsec = s * 3600.0
        max_scale = float(np.max(scales_arcsec))
        min_scale = float(np.min(scales_arcsec))
        aniso = float(max_scale / max(min_scale, 1e-12))
        det = float(np.linalg.det(m))
    except Exception:
        return (False, 'could not compute WCS pixel-scale matrix diagnostics')

    # Expected scale from instrument + binning in header
    try:
        binn = tel.get_binning(header)
        expected = float(tel.get_pixscale() * int(str(binn)[0]))
    except Exception:
        expected = None

    _log_fine_align(log, 'info',
        f'Fine-align WCS diagnostics: scales={scales_arcsec}, anisotropy={aniso:.3f}, det={det:.3e}, expected_scale={expected}')

    # Basic physical sanity
    if min_scale <= 0.0 or max_scale <= 0.0:
        return (False, 'non-positive WCS pixel scale')
    if aniso > 3.0:
        return (False, f'excessive WCS anisotropy ({aniso:.2f} > 3.0)')
    if expected is not None:
        if min_scale < 0.5 * expected:
            return (False, f'WCS min scale too small ({min_scale:.4f} < {0.5*expected:.4f} arcsec/pix)')
        if max_scale > 1.5 * expected:
            return (False, f'WCS max scale too large ({max_scale:.4f} > {1.5*expected:.4f} arcsec/pix)')

    return (True, '')

def validate_existing_wcs_header(header):
    """Check that a FITS header already defines a usable celestial WCS.

    Requires CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, and a linear
    celestial transform via either the CD matrix (CD1_1..CD2_2) or CDELT1
    and CDELT2 (optionally with PC matrix). Confirms :class:`astropy.wcs.WCS`
    can parse the header and evaluate pixel-to-world for one point.

    Parameters
    ----------
    header : astropy.io.fits.Header
        Header to validate (typically primary HDU).

    Returns
    -------
    tuple
        (ok, reason) where *ok* is True if validation passed; *reason* is an
        empty string on success, or a short message on failure.
    """
    required = ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2')
    for k in required:
        if k not in header:
            return (False, f'missing keyword {k}')
    has_cd = all(k in header for k in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'))
    has_cdelt = 'CDELT1' in header and 'CDELT2' in header
    if not has_cd and not has_cdelt:
        return (False, 'need CD1_1..CD2_2 or both CDELT1 and CDELT2')
    try:
        w = WCS(header, naxis=2)
        w.pixel_to_world(1.0, 1.0)
    except Exception as e:
        return (False, f'WCS not usable: {e}')
    return (True, '')

def _field_center_from_fits_primary(input_file, log=None):
    """Parse a field-center SkyCoord from the primary FITS header.

    Returns
    -------
    SkyCoord or False
        ICRS center, or False if RA/Dec could not be parsed.
    """
    with fits.open(input_file) as hdu:
        header = hdu[0].header

    check_pairs = [('CRVAL1', 'CRVAL2'), ('RA', 'DEC'), ('OBJCTRA', 'OBJCTDEC')]
    for pair in check_pairs:
        if pair[0] in header.keys() and pair[1] in header.keys():
            ra = header[pair[0]]
            dec = header[pair[1]]
            coord = utilities.parse_coord(ra, dec)
            if coord:
                return coord

    if log:
        log.error(f'Could not parse RA/DEC from header of {input_file}')
    return False


def get_fine_align_reference_catalog(input_file, catalog='gaia', radius_deg=0.5,
    log=None):
    """Query the configured reference catalog for fine WCS alignment.

    Parameters
    ----------
    input_file : str
        Path to FITS file; primary header used for field center.
    catalog : str, optional
        Catalog key (see :data:`potpyri.utils.catalogs.FINE_ALIGN_CATALOG_CHOICES`).
        Default is ``'gaia'``.
    radius_deg : float, optional
        Half-width of the VizieR query box in degrees. Default 0.5.
    log : ColoredLogger, optional
        Logger.

    Returns
    -------
    astropy.table.Table or None or False
        Table with ``RA_ICRS``, ``DE_ICRS`` in degrees, None if empty or query
        failed, False if field center could not be parsed from the header.
    """
    coord = _field_center_from_fits_primary(input_file, log=log)
    if coord is False:
        return False
    field_width = float(radius_deg) * u.deg
    return catalogs.fetch_astrometry_reference_table(
        coord, catalog, field_width, log=log)

def clean_up_astrometry(directory, file, exten):
    """Remove astrometry.net temporary files (.axy, .corr, .solved, .wcs, etc.).

    Parameters
    ----------
    directory : str
        Directory containing the files.
    file : str
        Base filename (with extension); exten is replaced to form temp names.
    exten : str
        File extension (e.g. '.fits') to replace with .axy, .corr, etc.

    Returns
    -------
    None
    """
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

def _find_astrometry_index_dir():
    """Return a directory that contains astrometry.net ``index-*.fits`` files.

    Prefers common system install locations and ``$ANET_DATA`` when the
    default conda config points at an empty data directory.
    """
    candidates = []
    env = os.environ.get('ANET_DATA') or os.environ.get('ASTROMETRY_DATA')
    if env:
        candidates.append(env)
    candidates.extend([
        '/usr/local/astrometry/data',
        '/usr/share/astrometry',
        '/opt/astrometry/data',
        os.path.expanduser('~/astrometry-data'),
    ])
    # Also probe beside the solve-field that will run
    which = subprocess.run(['which', 'solve-field'], capture_output=True, text=True)
    sf = which.stdout.strip()
    if sf:
        bindir = os.path.dirname(os.path.abspath(sf))
        root = os.path.dirname(bindir)
        candidates.extend([
            os.path.join(root, 'data'),
            os.path.join(root, 'share', 'astrometry'),
        ])
    for path in candidates:
        if not path or not os.path.isdir(path):
            continue
        try:
            if any(name.startswith('index-') and name.endswith('.fits')
                   for name in os.listdir(path)):
                return path
        except OSError:
            continue
    return None

def solve_astrometry(file, tel, binn, paths, radius=0.5, replace=True,
    shift_only=False, index=None, log=None):
    """Run astrometry.net (solve-field) to get coarse WCS; optionally use custom index.

    Parameters
    ----------
    file : str
        Path to FITS file to solve.
    tel : Instrument
        Instrument instance (name, get_pixscale).
    binn : str
        Binning string (e.g. '22') for pixel scale.
    paths : dict
        Paths dict (source_extractor key for solve-field).
    radius : float, optional
        Search radius in degrees. Default is 0.5.
    replace : bool, optional
        If True, write WCS into original file; else write .solved.fits.
    shift_only : bool, optional
        If True, only update CRPIX/CRVAL. Default is False.
    index : str, optional
        Path to custom index file or directory for solve-field.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    bool
        True if solve succeeded (or an existing usable header WCS was kept as
        fallback), False otherwise.
    """
    # Starting solve, print file and directory for reference
    fullfile = os.path.abspath(file)
    directory = os.path.dirname(file)

    if log:
        log.info(f'Trying to solve file: {file}')
    else:
        print(f'Trying to solve file: {file}')

    if not os.path.exists(fullfile):
        return(False)

    with fits.open(fullfile) as hdu:
        header = hdu[0].header

    exten = '.'+file.split('.')[-1]
    if not replace:
        if os.path.exists(fullfile.replace(exten,'.solved.fits')):
            if log:
                log.info(f'SUCCESS: solved {fullfile}')
            else:
                print(f'SUCCESS: solved {fullfile}')
            return(True)

    exten = '.'+file.split('.')[-1]

    # Prefer telescope pointing keywords, but fall back to CRVAL when needed.
    # (BINOSPEC often has a more reliable CRVAL than RA/DEC.)
    if getattr(tel, 'name', '').upper() == 'BINOSPEC':
        check_pairs = [('CRVAL1','CRVAL2'),('RA','DEC'),('OBJCTRA','OBJCTDEC')]
    else:
        check_pairs = [('RA','DEC'),('OBJCTRA','OBJCTDEC'),('TARGRA','TARGDEC'),
            ('CRVAL1','CRVAL2')]

    coord = None
    existing_wcs_ok, existing_wcs_reason = validate_existing_wcs_header(header)

    for pair in check_pairs:
        if pair[0] in header.keys() and pair[1] in header.keys():
            ra = header[pair[0]]
            dec = header[pair[1]]
            coord = utilities.parse_coord(ra, dec)
            if log: log.info(f'Got coord={coord} for RA key {pair[0]}, DEC key {pair[1]}')
            if log: log.info(f'RA value is {ra}, DEC value is {dec}')
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
    scale = tel.get_pixscale() * int(binn[0])
    if scale is not None:
        scale = float(scale)
        scale_high = float('%.4f'%(scale * 1.2))
        scale_low = float('%.4f'%(scale * 0.8))
    else:
        scale_high = 5.0
        scale_low = 0.05

    # Get RA and Dec
    ra = float('%.6f'%coord.ra.degree)
    dec = float('%.6f'%coord.dec.degree)

    # Optional Source Extractor flags (preferred for NIR / crowded residuals).
    # Only pass --source-extractor-path when we have a real executable; an empty
    # value makes solve-field treat the next flag (e.g. --continue) as the path.
    p = subprocess.run(['solve-field','-h'],capture_output=True)
    helptext = p.stdout.decode().lower()
    sex_args = ''
    sex_path = (paths.get('source_extractor') or '').strip()
    if not sex_path or not os.path.exists(sex_path):
        which = subprocess.run(['which', 'sex'], capture_output=True, text=True)
        sex_path = which.stdout.strip()
    if sex_path and os.path.exists(sex_path):
        if '--use-source-extractor' in helptext:
            sex_args = (
                f'--use-source-extractor '
                f'--source-extractor-path {sex_path} '
            )
        elif '--use-sextractor' in helptext:
            sex_args = '--use-sextractor '

    index_args = ''
    if index and os.path.exists(index):
        if os.path.isfile(index):
            index_args = f'--index-file {index} '
        elif os.path.isdir(index):
            index_args = f'--index-dir {index} '
    else:
        # Conda astrometry builds often ship an empty index directory. Prefer an
        # explicit index-dir that actually contains index-*.fits.
        index_dir = _find_astrometry_index_dir()
        if index_dir:
            index_args = f'--index-dir {index_dir} '
            if log:
                log.info(f'Using astrometry index directory: {index_dir}')


    def _base_args(include_radec=True, include_sex=True):
        args = (
            f'--scale-units arcsecperpix '
            f'--scale-low {scale_low} --scale-high {scale_high} '
            f'--no-plots -T --overwrite -N {newfile} --dir {directory} '
        )
        if include_radec:
            args += f'--ra {ra} --dec {dec} --radius {radius} '
        if include_sex:
            args += sex_args
        args += index_args
        return args

    # Retry ladder tuned for instruments with usable header WCS (e.g. MOSFIRE)
    # and for NIR frames where bright residuals make a tiny --objs cut harmful.
    # Empirically, --objs 15 fails on MOSFIRE K-band while --objs 50 / no limit
    # and --continue (existing WCS near the telescope pointing) succeed.
    attempts = []
    continue_ok = False
    if existing_wcs_ok:
        try:
            wcs_coord = utilities.parse_coord(header['CRVAL1'], header['CRVAL2'])
            sep_deg = float(coord.separation(wcs_coord).degree)
            continue_ok = sep_deg <= float(radius)
            if log and not continue_ok:
                log.info(
                    f'Skipping --continue: CRVAL is {sep_deg:.3f} deg from '
                    f'pointing (limit {radius} deg)'
                )
        except Exception:
            continue_ok = False
    if continue_ok and sex_args:
        attempts.append(('continue existing WCS', _base_args(),
                         '--continue --no-verify'))
    elif continue_ok:
        attempts.append(('continue existing WCS', _base_args(include_sex=False),
                         '--continue --no-verify'))
    if sex_args:
        attempts.append(('sextractor objs=50 downsample=2', _base_args(),
                         '--downsample 2 --no-verify --objs 50'))
        attempts.append(('sextractor objs=100', _base_args(),
                         '--no-verify --objs 100'))
        attempts.append(('sextractor no objs limit', _base_args(),
                         '--no-verify'))
    attempts.append(('built-in detection + RA/Dec',
                     _base_args(include_sex=False), '--no-verify'))
    attempts.append(('built-in detection, no RA/Dec',
                     _base_args(include_radec=False, include_sex=False),
                     '--no-verify'))

    for tries, (label, args, extra_opts) in enumerate(attempts, start=1):
        # Avoid reusing a truncated .axy from a prior --objs-limited attempt
        clean_up_astrometry(directory, file, exten)
        if os.path.exists(newfile):
            os.remove(newfile)

        input_args = f'{args} {extra_opts}'
        if log:
            log.info(f'Try #{tries} with astrometry.net ({label})...')
        else:
            print(f'Try #{tries} with astrometry.net ({label})...')

        process = ['solve-field', fullfile] + input_args.split()
        if log:
            log.info(' '.join(process))
        else:
            print(' '.join(process))

        with open(os.devnull, 'w') as FNULL:
            p = subprocess.Popen(process, stdout=FNULL, stderr=subprocess.STDOUT)
            try:
                p.wait(180)
            except subprocess.TimeoutExpired:
                p.kill()

        if os.path.exists(newfile):
            break

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
            with fits.open(newfile) as newhdu:
                with fits.open(output_file) as hdu:
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

        with fits.open(output_file) as hdu:
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
        clean_up_astrometry(directory, file, exten)
        # MOSFIRE (and similar) often ship a usable header WCS with only a small
        # offset; keep it so fine_align_wcs can refine rather than dropping the frame.
        if existing_wcs_ok:
            msg = (
                f'astrometry.net did not solve {fullfile}; keeping existing '
                f'header WCS for fine alignment'
                + (f' ({existing_wcs_reason}).' if existing_wcs_reason else '.')
            )
            if log:
                log.warning(msg)
            else:
                print(msg)
            with fits.open(fullfile, mode='update') as hdu:
                hdu[0].header['ASTNET'] = (
                    False, 'astrometry.net solve failed; kept header WCS')
                hdu.flush()
            return(True)

        if log:
            log.error(f'FAILURE: did not solve {fullfile}')
        else:
            print(f'FAILURE: did not solve {fullfile}')
        return(False)


def fine_align_wcs(file, tel, catalog='gaia', radius=0.5,
    max_search_radius=5.0*u.arcsec, save_centroids=False, min_cat_match=7, log=None):
    """Refine WCS by matching detections to a reference catalog (fine alignment).

    Parameters
    ----------
    file : str
        Path to FITS file with WCS (e.g. after solve_astrometry).
    tel : Instrument
        Instrument instance.
    catalog : str, optional
        Reference catalog: ``gaia`` (default), ``panstarrs``, ``sdss``, ``legacy``,
        ``twomass``, or ``skymapper``. See
        :data:`potpyri.utils.catalogs.FINE_ALIGN_CATALOG_CHOICES`.
    radius : float, optional
        Catalog query half-width in degrees. Default is 0.5.
    max_search_radius : astropy.units.Quantity, optional
        Match radius between detections and catalog positions. Default is 5 arcsec.
    save_centroids : bool, optional
        If True, save centroid table. Default is False.
    min_cat_match : int, optional
        Minimum number of in-frame catalog sources required before matching.
        Default is 7.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    bool
        True if alignment succeeded and WCS was updated, False otherwise.
    """
    catalog = catalogs.normalize_fine_align_catalog(catalog)
    cat = get_fine_align_reference_catalog(
        file, catalog=catalog, radius_deg=radius, log=log)

    with fits.open(file) as hdu:
        header = hdu[0].header

        if cat is False:
            _write_fine_align_fallback_header(hdu, file,
                'could not parse coordinate keywords', log=log)
            return(True)
        if cat is None or len(cat) == 0:
            _write_fine_align_fallback_header(hdu, file,
                f'{catalog} query returned no alignment stars', log=log)
            return(True)

        _log_fine_align(log, 'info', f'Got {len(cat)} {catalog} alignment stars')

        # Estimate sources that land within the image using WCS from header
        w = WCS(header, relax=True)
        naxis1 = header['NAXIS1']
        naxis2 = header['NAXIS2']

        ref_coords = SkyCoord(cat['RA_ICRS'], cat['DE_ICRS'], unit=(u.deg, u.deg))
        x_pix, y_pix = w.world_to_pixel(ref_coords)

        mask = (x_pix > 0) & (x_pix < naxis1) & (y_pix > 0) & (y_pix < naxis2)
        x_pix = x_pix[mask]
        y_pix = y_pix[mask]
        cat = cat[mask]
        ref_coords = ref_coords[mask]

        if cat is None or len(cat) == 0:
            _write_fine_align_fallback_header(hdu, file,
                f'no {catalog} stars project into image bounds', log=log)
            return(True)

        _log_fine_align(log, 'info',
            f'{catalog} stars projected into image bounds: {len(cat)}')

        if len(cat) < min_cat_match:
            _write_fine_align_fallback_header(hdu, file,
                f'too few in-frame catalog stars ({len(cat)} < {min_cat_match})',
                log=log)
            return(True)

        sky_cat = photometry.run_sextractor(file, log=log)
        if sky_cat is None or len(sky_cat) == 0:
            _write_fine_align_fallback_header(hdu, file,
                'Source Extractor returned no sources for fine alignment', log=log)
            return(True)
        _log_fine_align(log, 'info', f'SExtractor sources available for matching: {len(sky_cat)}')

        central_pix = (float(naxis1)/2., float(naxis2)/2.)
        central_coo = w.pixel_to_world(*central_pix)

        _log_fine_align(log, 'info', 'Converging on fine WCS solution')

        bar = progressbar.ProgressBar(maxval=10)
        bar.start()
        succeeded = False
        cat_coords_match = None
        best_n_match = 0
        best_iter = -1
        n_iter_no_match = 0
        for i in np.arange(10):
            bar.update(i+1)

            cat_coords = w.pixel_to_world(sky_cat['X_IMAGE'], sky_cat['Y_IMAGE'])

            idx, d2d, d3d = match_coordinates_sky(ref_coords, cat_coords)

            sep_mask = d2d < max_search_radius
            n_match = int(np.sum(sep_mask))
            if n_match > best_n_match:
                best_n_match = n_match
                best_iter = int(i + 1)
            if n_match == 0:
                n_iter_no_match += 1
                continue

            sky_coords_match = ref_coords[sep_mask]
            cat_coords_match = sky_cat[idx[sep_mask]]

            sip_degree = 0
            if len(cat_coords_match)>100:
                sip_degree = 2
            if len(cat_coords_match)>500:
                sip_degree = 4
            if len(cat_coords_match)>2000:
                sip_degree = 6

            # Make a fresh table each iteration; avoid duplicate-column failures
            cat_coords_match = Table(cat_coords_match)
            cat_coords_match['RA_ICRS'] = sky_coords_match.ra.degree
            cat_coords_match['DE_ICRS'] = sky_coords_match.dec.degree
            fit_coords = SkyCoord(cat_coords_match['RA_ICRS'],
                cat_coords_match['DE_ICRS'], unit=(u.deg, u.deg))

            xy = (cat_coords_match['X_IMAGE'], cat_coords_match['Y_IMAGE'])
            central_coo = w.pixel_to_world(*central_pix)
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore", category=RuntimeWarning,
                        message=".*cdelt will be ignored since cd is present.*")
                    w = fit_wcs_from_points(xy, fit_coords, proj_point=central_coo,
                        sip_degree=sip_degree)
                succeeded = True
            except ValueError:
                # Try with next iteration
                continue

        bar.finish()
        _log_fine_align(log, 'info',
            f'Fine-align match summary for {file} ({catalog}): '
            f'best matches={best_n_match} at iter={best_iter}, '
            f'iters with zero matches={n_iter_no_match}/10, '
            f'min required={min_cat_match}, match radius={max_search_radius}.')

        if not succeeded:
            _write_fine_align_fallback_header(hdu, file,
                'fine WCS fit failed during iterative catalog matching', log=log)
            return(True)

        _log_fine_align(log, 'info', f'Solved with SIP degree = {sip_degree}')

        if cat_coords_match is None or len(cat_coords_match)==0:
            _write_fine_align_fallback_header(hdu, file,
                'zero matched catalog/SExtractor sources after iterations', log=log)
            return(True)
        else:
            _log_fine_align(log, 'info', f'Matched to {len(cat_coords_match)} coordinates')

        if len(cat_coords_match)<min_cat_match:
            _write_fine_align_fallback_header(hdu, file,
                f'too few matched stars ({len(cat_coords_match)} < {min_cat_match})',
                log=log)
            return(True)

        if n_iter_no_match >= 8 and len(cat_coords_match) <= (min_cat_match + 1):
            _write_fine_align_fallback_header(hdu, file,
                f'unstable catalog matching (zero-match iters={n_iter_no_match}/10, matched={len(cat_coords_match)})',
                log=log)
            return(True)

        # Bootstrap WCS solution and get center pixel
        ra_cent = [] ; de_cent = []
        coords = SkyCoord(cat_coords_match['RA_ICRS'], cat_coords_match['DE_ICRS'],
            unit=(u.deg, u.deg))

        bar = progressbar.ProgressBar(maxval=1000)
        bar.start()

        idx = np.arange(len(cat_coords_match))
        for i in np.arange(1000):
            bar.update(i+1)

            if len(idx)>200:
                size = 100
            else:
                size = int(len(idx)/2)

            idxs = np.random.choice(idx, size=size, replace=False)

            xy = (cat_coords_match['X_IMAGE'][idxs],
                  cat_coords_match['Y_IMAGE'][idxs])
            c = coords[idxs]

            central_coo = w.pixel_to_world(*central_pix)
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", category=RuntimeWarning,
                    message=".*cdelt will be ignored since cd is present.*")
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

        # Create a scatter plot of the centroid values to validate the astrometry
        if save_centroids:
            fig, ax = plt.subplots()
            ax.scatter(ra_cent[mask], de_cent[mask])
            plt.savefig(file.replace('.fits','_centroids.png'))

        header = w.to_header()
        header['CRPIX1']=central_pix[0]
        header['CRPIX2']=central_pix[1]
        header['CRVAL1']=med_ra
        header['CRVAL2']=med_de
        header['CUNIT1']='deg'
        header['CUNIT2']='deg'
        if sip_degree>0:
            header['CTYPE1']='RA---TAN-SIP'
            header['CTYPE2']='DEC--TAN-SIP'
            for kw, val in w._write_sip_kw().items():
                header[kw] = val

        # Validate refined WCS; if pathological, fall back to coarse solution.
        is_valid, reason = _validate_refined_wcs(w, hdu[0].header, tel, log=log)
        if not is_valid:
            _write_fine_align_fallback_header(hdu, file, f'pathological refined WCS: {reason}', log=log)
            return(True)

        # Delete all old WCS keys
        delkeys = ['WCSNAME','CUNIT1','CUNIT2','CTYPE1','CTYPE2','CRPIX1','CRPIX2',
            'CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2','RADECSYS',
            'PC1_1','PC1_2','PC2_1','PC2_2','LONPOLE','LATPOLE','CDELT1','CDELT2']
        while True:
            found_key = False
            for key in hdu[0].header.keys():
                if any([k in key for k in delkeys]):
                    found_key=True
                    del hdu[0].header[key]
            if not found_key:
                break

        hdu[0].header.update(header)

        ra_disp = std_ra*3600.0*np.cos(np.pi/180.0 * mean_de)
        de_disp = std_de*3600.0

        ra_disp = float('%.6f'%ra_disp)
        de_disp = float('%.6f'%de_disp)

        _log_fine_align(log, 'info', f'R.A. dispersion ={ra_disp}\", Decl. dispersion={de_disp}\"')

        hdu[0].header['RADISP']=(ra_disp, 'Dispersion in R.A. of WCS [Arcsec]')
        hdu[0].header['DEDISP']=(de_disp, 'Dispersion in Decl. of WCS [Arcsec]')
        if 'ALGNFAIL' in hdu[0].header:
            del hdu[0].header['ALGNFAIL']

        hdu.writeto(file, overwrite=True, output_verify='silentfix')
    return(True)
