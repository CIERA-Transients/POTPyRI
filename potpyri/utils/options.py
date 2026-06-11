"""Pipeline option parsing and path setup.

Provides argument parsing for the main pipeline and builds directory layouts
(raw, red, log, cals, workspace, etc.) for a given data path and instrument.
Authors: Charlie Kilpatrick.
"""
from potpyri._version import __version__

import os
import shutil
import sys
import subprocess

from potpyri import instruments

def init_options():
    """Build and return the argument parser for the main pipeline.

    Returns
    -------
    argparse.ArgumentParser
        Parser with instrument, data_path, target, proc, include-bad,
        file-list-name, photometry, and processing flags.
    """
    import argparse
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('instrument', 
        default=None, 
        type=str.upper,
        choices=instruments.__all__,
        help='''Name of instrument (must be in instruments directory) of data to 
        reduce. Required to run pipeline.''')
    params.add_argument('data_path', 
        default=None, 
        help='''Path of data to reduce. See manual for specific details. 
        Required to run pipeline.''')
    params.add_argument('--target', 
        type=str, 
        default=None, 
        help='''Option to only reduce a specific target. String used here must 
        be contained within the target name in file headers. Optional 
        parameter.''')
    params.add_argument('--proc',
        type=str,
        default=True,
        help='''File-ingestion mode for raw discovery (instrument-specific glob).
        Use ``fits`` (or ``.fits``) on any instrument for uncompressed ``*.fits``.
        Other values depend on the instrument (e.g. GMOS ``dragons``, LRIS
        ``archive``/``raw``, BINOSPEC ``proc``).''')
    params.add_argument('--include-bad','--incl-bad', 
        default=False,
        action='store_true', 
        help='''Include files flagged as bad in the file list.''')
    params.add_argument('--no-redo-sort',
        default=False,
        action='store_true', 
        help='''Do not resort files and generate new file list.''')
    params.add_argument('--file-list-name', 
        type=str, 
        default='file_list.txt', 
        help='''Change the name of the archive file list.''')
    params.add_argument('--phot-sn-min',
        type=float,
        default=3.0,
        help='''Minimum signal-to-noise to try in photometry loop.''')
    params.add_argument('--phot-sn-max',
        type=float,
        default=20.0,
        help='''Maximum signal-to-noise to try in photometry loop.''')
    params.add_argument('--fwhm-init',
        type=float,
        default=5.0,
        help='''Initial FWHM (in pixels) for photometry loop.''')
    params.add_argument('--skip-skysub',
        default=False,
        action='store_true', 
        help='''Do not perform sky subtraction during image processing. Equivalent
        to ``--bkg-sub none``.''')
    params.add_argument('--bkg-sub',
        default='local',
        choices=['local', 'constant', 'none'],
        dest='bkg_sub',
        help='''Per-frame background subtraction for optical data: ``local`` uses a
        spatially varying 2D mesh (default); ``constant`` subtracts one
        sigma-clipped median per frame; ``none`` skips subtraction.''')
    params.add_argument('--fieldcenter',
        default=None,
        nargs=2,
        help='''Align the output images to this center coordinate.  This is 
        useful for difference imaging where the frames need to be a common
        size, pixel scale, and set of coordinates.''')
    params.add_argument('--out-size',
        default=None,
        type=int,
        help='''Output image size (image will be SIZE x SIZE pixels).''')
    params.add_argument('--skip-flatten',
        default=False,
        action='store_true',
        help='Tell the pipeline to skip flattening.')
    params.add_argument('--skip-cr',
        default=False,
        action='store_true',
        help='Tell the pipeline to skip cosmic ray detection.')
    params.add_argument('--skip-fine-align',
        default=False,
        action='store_true',
        dest='skip_fine_align',
        help='Skip fine WCS alignment (catalog matching after astrometry.net).')
    params.add_argument('--skip-gaia',
        default=False,
        action='store_true',
        dest='skip_fine_align',
        help='Deprecated alias for --skip-fine-align.')
    params.add_argument('--fine-align-catalog',
        type=str,
        default='gaia',
        choices=['gaia', 'panstarrs', 'sdss', 'legacy', '2mass', 'skymapper'],
        dest='fine_align_catalog',
        help='Reference catalog for fine WCS alignment after astrometry.net: gaia (default), '
             'panstarrs, sdss (SDSS V/147), legacy, 2mass, or skymapper.')
    params.add_argument('--skip-external-astrometry',
        default=False,
        action='store_true',
        dest='skip_external_astrometry',
        help='Do not run astrometry.net or fine catalog alignment; require a valid WCS '
             'already in the FITS header (CTYPE/CRVAL/CRPIX and CD or CDELT keywords). '
             'Implies --skip-fine-align.')
    params.add_argument('--keep-all-astro',
        default=False,
        action='store_true',
        help='Tell the pipeline to keep all images regardless of astrometric dispersion.')
    params.add_argument('--relative-calibration',
        default=False,
        action='store_true',
        help='Before stacking, calibrate frames to each other using SExtractor catalogs '
             'and RA/Dec cross-matching; use relative source fluxes to set combine scales.')

    return(params)


def add_options():
    """Parse command-line options and return normalized args.

    Normalizes instrument aliases via :func:`potpyri.instruments.resolve_instrument_name`.

    Returns
    -------
    argparse.Namespace
        Parsed and normalized command-line arguments.
    """
    params = init_options()
    args = params.parse_args()

    if args.instrument is not None:
        args.instrument = instruments.resolve_instrument_name(args.instrument)

    if args.skip_skysub:
        args.bkg_sub = 'none'

    return(args)

def test_for_dependencies(skip_external_astrometry=False):
    """Check that astrometry.net and Source Extractor are on the system path.

    Parameters
    ----------
    skip_external_astrometry : bool, optional
        If True, do not require astrometry.net (solve-field); useful when using
        only header WCS with ``--skip-external-astrometry``.

    Raises
    ------
    Exception
        If solve-field (astrometry.net) or sex (sextractor) is not found
        or does not report expected help text.
    """
    if not skip_external_astrometry:
        p = subprocess.run(['solve-field','-h'], capture_output=True)
        data = p.stdout.decode().lower()

        # Check for astrometry.net in output
        if 'astrometry.net' not in data:
            raise Exception(f'''Astrometry.net is a dependency of POTPyRI.  Download
                and install the binaries and required index files from:
                https://astrometry.net/use.html''')

    p = subprocess.run(['sex','-h'], capture_output=True)
    data = p.stderr.decode().lower()
    if 'syntax: sex' not in data:
        raise Exception(f'''source extractor is a dependency of POTPyRI.  
            Install via Homebrew (https://formulae.brew.sh/formula/sextractor), 
            apt-get, or directly from the source code:
            https://github.com/astromatic/sextractor.''')

def add_paths(data_path, file_list_name, tel):
    """Build the directory and path dictionary for a reduction run.

    Creates raw, bad, red, log, cal, work dirs under data_path and resolves
    the path to the Source Extractor binary.

    Parameters
    ----------
    data_path : str
        Top-level data directory.
    file_list_name : str
        Filename of the file list (e.g. files.txt).
    tel : Instrument
        Instrument instance (used for name and caldb).

    Returns
    -------
    dict
        Keys include 'data', 'raw', 'bad', 'red', 'log', 'cal', 'work',
        'filelist', 'caldb', 'source_extractor'.

    Raises
    ------
    Exception
        If data_path does not exist.
    """
    if not os.path.exists(data_path):
        raise Exception(f'Data path does not exist: {data_path}')

    # Get path to code directory
    paths = {'data': os.path.abspath(data_path)}
    paths['abspath']=os.path.abspath(__file__)
    paths['code']=os.path.split(paths['abspath'])[0]
    paths['caldb']=os.path.join(paths['code'], '..', 'data', 'cal', tel.name.upper())
    paths['raw']=os.path.join(paths['data'], 'raw') #path containing the raw data
    paths['bad']=os.path.join(paths['data'], 'bad') #path containing the unused data
    paths['red']=os.path.join(paths['data'], 'red') #path to write the reduced files
    paths['log']=os.path.join(paths['data'], 'red', 'log')
    paths['cal']=os.path.join(paths['data'], 'red', 'cals')
    paths['work']=os.path.join(paths['data'], 'red', 'workspace')
    paths['filelist']=os.path.join(paths['data'], file_list_name)
    
    for key in paths.keys():
        if key in ['caldb','filelist']: continue
        if not os.path.exists(paths[key]): os.makedirs(paths[key])

    p=subprocess.run(['which','sex'], capture_output=True)
    paths['source_extractor']=p.stdout.decode().strip()

    return(paths)

def initialize_telescope(instrument, data_path):
    """Load the instrument class and build paths for that instrument.

    Parameters
    ----------
    instrument : str
        Instrument name (e.g. 'GMOS', 'LRIS').
    data_path : str
        Top-level data directory.

    Returns
    -------
    tuple
        (paths, tel) where paths is the dict from add_paths and tel is the
        instrument instance. Uses default file list name.
    """
    tel = instruments.instrument_getter(instrument)

    # Generate code and data paths based on input path (default file list name)
    paths = add_paths(data_path, 'files.txt', tel)

    return(paths, tel)

