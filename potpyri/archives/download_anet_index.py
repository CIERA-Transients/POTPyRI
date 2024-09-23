#!/usr/bin/env python

import os
import sys
import shutil
import argparse
from astropy.utils.data import download_file

def parse_astrometry_config(config_file):

    data_path = None

    with open(config_file, 'r') as f:
        for line in f:
            if line.strip().startswith('add_path'):
                data = line.strip().split()
                data_path = data[1]
                break

    if data_path is None:
        raise Exception('ERROR: could not parse data path from astrometry.cfg file')

    return(data_path)

def download_index_files(outdir=None):

    # Define outdir based on where solve-field is located in path
    if outdir is None:
        path = shutil.which('solve-field')

        if not path or not os.path.exists(path):
            raise Exception('ERROR: cannot find astrometry.net in path.  Install first.')

        path = os.path.dirname(path)
        config_file = os.path.join(path, '..', 'etc', 'astrometry.cfg')

        if not os.path.exists(config_file):
            raise Exception('ERROR: cannot find astrometry.cfg file.  Specify data path instead.')

        outdir = parse_astrometry_config(config_file)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    from bs4 import BeautifulSoup

    import requests

    baseurl = "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/"
    r  = requests.get(baseurl)
    data = r.text
    soup = BeautifulSoup(data, features="html5lib")

    for link in soup.find_all('a'):
        text = link.get('href')
        if 'index' in text and '.fits' in text:
            fileurl = os.path.join(baseurl, text)
            outfile = os.path.join(outdir, text)
            if os.path.exists(outfile):
                continue
            print(f'Downloading outfile {outfile}')
            data = download_file(fileurl, show_progress=True)
            shutil.move(data, outfile)

def add_options():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('--data_path', default=None, 
        help='''Path for downloading data.''')

    args = params.parse_args()
    return(args)

def main():
    args = add_options()
    download_index_files(outdir=args.data_path)

if __name__=="__main__":
    main()
