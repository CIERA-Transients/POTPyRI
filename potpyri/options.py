import os
import shutil
import sys

def add_options():
    import argparse
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('instrument', 
        default=None, 
        choices=['Binospec','DEIMOS','GMOS','LRIS','MMIRS','MOSFIRE','TEST'],
        help='''Name of instrument (must be in params folder) of data to 
        reduce. Required to run pipeline. Use TEST for pipeline test.''')
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
        help='''Option to use the _proc data from MMT. Optional parameter. 
        Default is False.''')
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

    args = params.parse_args()

    # Handle/parse options
    if 'bino' in args.instrument:
        args.instrument = 'binospec'
    if 'mmir' in args.instrument:
        args.instrument = 'mmirs'

    return(args)

def add_paths(data_path):

    if not os.path.exists(data_path):
        raise Exception(f'Data path does not exist: {data_path}')

    # Get path to code directory
    paths = {'data': os.path.abspath(data_path)}
    paths['abspath']=os.path.abspath(sys.argv[0])
    paths['code']=os.path.split(paths['abspath'])[0]
    paths['config']=os.path.join(paths['code'], 'config')
    paths['raw']=os.path.join(data_path, 'raw') #path containing the raw data
    paths['bad']=os.path.join(data_path, 'bad') #path containing the unused data
    paths['red']=os.path.join(data_path, 'red') #path to write the reduced files
    paths['log']=os.path.join(data_path, 'red', 'log')
    paths['cal']=os.path.join(data_path, 'red', 'cals')
    paths['work']=os.path.join(data_path, 'red', 'workspace')
    for key in paths.keys():
        if not os.path.exists(paths[key]): os.makedirs(paths[key])

    # Copy config directory to data path
    if not os.path.exists(os.path.join(paths['data'], 'config')):
        shutil.copytree(paths['config'], os.path.join(paths['data'], 'config'))

    return(paths)