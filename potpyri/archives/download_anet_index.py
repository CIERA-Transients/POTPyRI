import os
import sys
import shutil
import argparse
from astropy.utils.data import download_file

def download_index_files(outdir='/opt/homebrew/Cellar/astrometry-net/0.95_2/data'):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    from bs4 import BeautifulSoup

    import requests

    baseurl = "http://broiler.astrometry.net/~dstn/4200/"
    r  = requests.get(baseurl)
    data = r.text
    soup = BeautifulSoup(data, features="html5lib")

    for link in soup.find_all('a'):
        text = link.get('href')
        if 'index' in text:
            fileurl = os.path.join(baseurl, text)
            outfile = os.path.join(outdir, text)
            if os.path.exists(outfile):
                continue
            print(f'Downloading outfile {outfile}')
            data = download_file(fileurl, show_progress=True)
            shutil.move(data, outfile)

def add_options():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('data_path', default='.', 
        help='''Path for downloading data.''')

    args = params.parse_args()
    return(args)

if __name__=="__main__":
    args = add_options()
    download_index_files(outdir=args.data_path)
