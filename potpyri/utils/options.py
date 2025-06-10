"Functions for initializing pipeline options and data paths."
"Authors: Charlie Kilpatrick"

# Initial version tracking on 09/21/2024
__version__ = "1.0"

import os
import importlib
import shutil
import sys
import subprocess

from potpyri import instruments

def init_options():
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
        help='''Option to specify file processing for data ingestion.''')
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
        default=40.0,
        help='''Maximum signal-to-noise to try in photometry loop.''')
    params.add_argument('--fwhm-init',
        type=float,
        default=5.0,
        help='''Initial FWHM (in pixels) for photometry loop.''')
    params.add_argument('--skip-skysub',
        default=False,
        action='store_true', 
        help='''Do not perform sky subtraction during image processing.''')
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
    params.add_argument('--stages',
        default=['CREATECALS','IMAGEPROC','WCS','DOPHOT','ABSPHOT'],
        nargs='+',
        help='''Stages to execute if running the pipeline in a modular way.''')
    params.add_argument('--skip-flatten',
        default=False,
        action='store_true',
        help='Tell the pipeline to skip flattening.')
    params.add_argument('--skip-cr',
        default=False,
        action='store_true',
        help='Tell the pipeline to skip cosmic ray detection.')
    params.add_argument('--skip-gaia',
        default=False,
        action='store_true',
        help='Tell the pipeline to skip Gaia alignment during WCS.')

    return(params)


def add_options():
    params = init_options()
    args = params.parse_args()

    # Handle/parse options
    if 'BINO' in args.instrument:
        args.instrument = 'BINOSPEC'
    if 'MMIR' in args.instrument:
        args.instrument = 'MMIRS'

    return(args)

def test_for_dependencies():

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

def add_paths(data_path, tel):

    if not os.path.exists(data_path):
        raise Exception(f'Data path does not exist: {data_path}')

    # Get path to code directory
    paths = {'data': os.path.abspath(data_path)}
    paths['abspath']=os.path.abspath(__file__)
    paths['code']=os.path.split(paths['abspath'])[0]
    paths['caldb']=os.path.join(paths['code'], 'data', 'cal', tel.name.upper())
    paths['raw']=os.path.join(paths['data'], 'raw') #path containing the raw data
    paths['bad']=os.path.join(paths['data'], 'bad') #path containing the unused data
    paths['red']=os.path.join(paths['data'], 'red') #path to write the reduced files
    paths['log']=os.path.join(paths['data'], 'red', 'log')
    paths['cal']=os.path.join(paths['data'], 'red', 'cals')
    paths['work']=os.path.join(paths['data'], 'red', 'workspace')
    for key in paths.keys():
        if key in ['caldb']: continue
        if not os.path.exists(paths[key]): os.makedirs(paths[key])

    p=subprocess.run(['which','sex'], capture_output=True)
    paths['source_extractor']=p.stdout.decode().strip()

    return(paths)

def initialize_telescope(instrument, data_path):

    module = importlib.import_module(f'potpyri.instruments.{instrument.upper()}')
    tel = getattr(module, instrument.upper())()

    # Generate code and data paths based on input path
    paths = add_paths(data_path, tel)

    return(paths, tel)

