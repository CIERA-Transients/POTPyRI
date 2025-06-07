import numpy as np
import sys
import io
import os
from pykoa.koa import Koa 
from astropy.time import Time
from astropy.table import Table
from astropy.table import Column

def login(cookie=None):
    if cookie and os.path.isfile(cookie):
        Koa.login(cookie)

def date_query(date, instr, form='ipac', cookie=None):
    t = Time(date)
    outdir = t.datetime.strftime('%Y%m%d')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile = os.path.join(outdir, f'{instr}.tbl')
    Koa.query_date(instr, date, outfile, overwrite=True, 
                   format=form, cookiepath=cookie)

    if os.path.exists(outfile):
        return(outfile, outdir)
    else:
        return(None, None)

def download_data(recfile, outdir, form='ipac', cookie=None):
    Koa.download(recfile, form, outdir, cookiepath=cookie)

def add_options():
    import argparse
    params = argparse.ArgumentParser()
    params.add_argument('date', 
        help='''Date for which to download data.''')
    params.add_argument('instrument',
        help='''Keck instrument to query in archive.''')
    params.add_argument('--cookie-file', default=None,
        help='''Provide a cookie file if you are downloading proprietary data.''')

    args = params.parse_args()
    return(args)

def main():
    args = add_options()

    if args.cookie_file is not None:
        login(args.cookie_file)

    outfile, outdir = date_query(args.date, args.instrument, cookie=args.cookie_file)
    if outfile is not None:
        download_data(outfile, outdir, cookie=args.cookie_file)

if __name__=="__main__":
    main()
