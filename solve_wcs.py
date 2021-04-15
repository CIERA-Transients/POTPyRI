#!/usr/bin/env python

"Python script for WCS solution."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "3.0" #last updated 14/04/2021

import sys
import numpy as np
import time
import astropy
import argparse
import subprocess
import os
import astropy.units as u
import astropy.wcs as wcs
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from astroquery.gaia import Gaia
from skimage import transform as tf
from scipy import stats
import itertools
import importlib
import tel_params

#turn Astropy warnings off
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

#function to calculate rms on astrometric solution
#David's rms function with Kerry's improvements.
def dvrms(x):
    x = np.asarray(x)
    if len(x) == 0:
        rms = 0
    else:
        rms = np.sqrt(np.sum(x**2)/len(x))
    return rms

#plate solution
def apply_wcs_transformation(header,tform):
    crpix1, crpix2, cd11, cd12, cd21, cd22 = wcs_keyword = tel.WCS_keywords()
    header[crpix1] = header[crpix1]-tform.translation[0]
    header[crpix2] = header[crpix2]-tform.translation[1]
    cd = np.array([header[cd11],header[cd12],header[cd21],header[cd22]]).reshape(2,2)
    cd_matrix = tf.EuclideanTransform(rotation=tform.rotation)
    cd_transformed = tf.warp(cd,cd_matrix)
    header[cd11] = cd_transformed[0][0]
    header[cd12] = cd_transformed[0][1]
    header[cd21] = cd_transformed[1][0]
    header[cd22] = cd_transformed[1][1]
    return header

#function to calculate rms on astrometric solution
#replaced np.median() with dvrms() -David V
def calculate_error(coord1, coord2, hd):
    sep = coord1.separation(coord2)
    d_ra = [s.hms[2] for s in sep]
    d_dec = [s.dms[2] for s in sep]
    rms_ra = dvrms(d_ra)
    rms_dec = dvrms(d_dec)
    hd['RA_RMS'] = (rms_ra, 'RMS of RA fit (arsec).')
    hd['DEC_RMS'] = (rms_dec, 'RMS of Dec fit (arcsec).')
    rms_total = dvrms(sep.arcsec)
    return rms_total

def run_sextractor(input_file, cat_name, tel, sex_config_dir='./Config'):

    if not os.path.exists(sex_config_dir+'/config'):
        e = 'ERROR: Could not find source extractor config file!'
        raise Exception(e)
    elif not os.path.exists(sex_config_dir+'/params'):
        e = 'ERROR: Could not find source extractor param file!'
        raise Exception(e)

    config = sex_config_dir+'/config'
    param_file = sex_config_dir+'/params'

    cmd=['sex','-c',config,input_file,'-CATALOG_NAME',cat_name]
    subprocess.call(cmd)

    params = []
    with open(param_file) as f:
        for line in f:
            params.append(line.split()[0].strip())

    if os.path.exists(cat_name):
        table = ascii.read(cat_name, names=params, comment='#')
        return(table)
    else:
        e = 'ERROR: source extractor did not successfully produce {0}'
        raise Exception(e.format(cat_name))

# Make quads for either x,y values or ra,dec under the assumption that x=ra
def make_quads(x, y, use=None, sky_coords=False):
    if use and use < len(x) and use < len(y):
        x = x[:use]
        y = y[:use]

    ind = list(itertools.combinations(np.arange(4), 2))
    if sky_coords:
        stars = []
        coords = [SkyCoord(x[i],y[i], unit='deg') for i in range(len(x))]
        dis_all = [[coords[j].separation(coords[i]).arcsec for i in range(len(coords))] for j in range(len(coords))]
        for j in range(len(coords)):
            i_n = [x for y,x in sorted(zip(dis_all[j],np.arange(0,len(coords))))][1:5]
            stars.append([(x[i],y[i]) for i in i_n])
            i_n = [x for y,x in sorted(zip(dis_all[j],np.arange(0,len(coords))))][5:9]
            stars.append([(x[i],y[i]) for i in i_n])
        coords = [[SkyCoord(stars[j][i][0], stars[j][i][1], unit='deg')
                 for i in range(4)] for j in range(len(stars))]
        d = [[coords[j][ind[i][0]].separation(coords[j][ind[i][1]]).arcsec for i in range(6)]
                for j in range(len(stars))]
    else:
        stars = list(itertools.combinations(zip(x, y), 4))
        dx = [[stars[j][ind[i][0]][0]-stars[j][ind[i][1]][0] for i in range(6)]
                for j in range(len(stars))]
        dy = [[stars[j][ind[i][0]][1]-stars[j][ind[i][1]][1] for i in range(6)]
            for j in range(len(stars))]
        d = [np.hypot(dx[j], dy[j]) for j in range(len(stars))]
    ds = [sorted(d[j]) for j in range(len(stars))]
    ratios = np.array([[ds[j][i]/ds[j][5] for i in range(5)]
            for j in range(len(stars))])

    return(stars, d, ds, ratios)

# Create 2D list with len(l1) x len(l2) with pairwise quotient of all elements
def broadcast_quotient(l1, l2):
    l1 = list(l1) ; l2 = list(l2)
    out = np.array([[a/b for a in l1] for b in l2])
    return(out)

def mask_catalog_for_wcs(gaiaorig, w, o, data_x, data_y):
    gaiacoord = SkyCoord(gaiaorig['ra'], gaiaorig['dec'], unit='deg')
    gaiax, gaiay = wcs.utils.skycoord_to_pixel(gaiacoord, w, o)
    mask = (gaiax > 50) & (gaiax < data_x-50) & (gaiay > 50) & (gaiay < data_y-50)

    return(mask)

def match_quads(stars,gaiastars,d,gaiad,ds,gaiads,ratios,gaiaratios,sky_coords=True):
    l1 = broadcast_quotient(ratios[:,0], gaiaratios[:,0])
    l2 = broadcast_quotient(ratios[:,4], gaiaratios[:,4])
    l3 = np.hypot(np.abs(l1-1), np.abs(l2-1))

    ind = list(itertools.combinations(np.arange(4), 2))
    gaia_ind = list(itertools.permutations(np.arange(0,4), 4))
    indm = []
    for i in range(len(l3)):
        min_all = np.where(l3[i]==np.min(l3[i]))[0]
        for j in min_all:
            indm.append([i,j])
    dr = np.array([[(gaiads[i[0]][j]/ds[i[1]][j]) for j in range(6)]
        for i in indm])
    if sky_coords:
        indmm = np.array([indm[i] for i in range(len(indm))
            if all(np.abs((dr[i]/tel.pixscale())-1)<0.05)])
    else:
        indmm = np.array([indm[i] for i in range(len(indm))
            if all(np.abs(dr[i]-1)<0.05)])
    starsm = [stars[i[1]] for i in indmm]
    gaiam = [gaiastars[i[0]] for i in indmm]
    dm = [d[i[1]] for i in indmm]
    gaiadm = [gaiad[i[0]] for i in indmm]
    starsx = []
    starsy = []
    gaiastarsra = []
    gaiastarsdec = []
    for k in range(len(starsm)):
        star_d = np.array(dm[k])
        ratio_quads = []
        for j in range(len(gaia_ind)):
            if sky_coords:
                gaia_coords = [SkyCoord(gaiam[k][i][0], gaiam[k][i][1], unit='deg') for i in gaia_ind[j]]
                gaia_d = np.array([[gaia_coords[ind[i][0]].separation(gaia_coords[ind[i][1]]).arcsec for i in range(6)]])
                ratio_quads.append(np.sum(np.abs(((gaia_d/star_d)/tel.pixscale())-1)))
            else:
                gaia_d = [np.hypot(gaiam[k][gaia_ind[i][0]][0]-starsm[k][gaia_ind[i][1]][0], gaiam[k][gaia_ind[i][0]][1]-starsm[k][gaia_ind[i][1]][1]) for i in range(6)]
                ratio_quads.append(np.sum(np.abs((gaia_d/star_d)-1)))
        starsx.append([starsm[k][i][0] for i in range(4)])
        starsy.append([starsm[k][i][1] for i in range(4)])
        gaiastarsra.append([gaiam[k][gaia_ind[np.argmin(ratio_quads)][i]][0] for i in range(4)])
        gaiastarsdec.append([gaiam[k][gaia_ind[np.argmin(ratio_quads)][i]][1] for i in range(4)])
    starsx_unq = list(dict.fromkeys(np.concatenate(starsx)))
    starsy_unq = list(dict.fromkeys(np.concatenate(starsy)))
    gaiastarsra_unq = list(dict.fromkeys(np.concatenate(gaiastarsra)))
    gaiastarsdec_unq = list(dict.fromkeys(np.concatenate(gaiastarsdec)))

    return np.array(starsx_unq[0:len(gaiastarsra_unq)]), np.array(starsy_unq[0:len(gaiastarsdec_unq)]), np.array(gaiastarsra_unq[0:len(starsx_unq)]), np.array(gaiastarsdec_unq[0:len(starsy_unq)])

#main function to calculate astrometric solution
def solve_wcs(input_file, telescope, sex_config_dir='./Config', static_mask=None, proc=None):
    #start time
    t_start = time.time()
    #import telescope parameter file
    global tel
    try:
        tel = importlib.import_module('tel_params.'+telescope)
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)

    #get data and header info
    hdr = fits.open(input_file)
    data = hdr[tel.wcs_extension()].data
    data_y, data_x = np.shape(data)
    header = hdr[tel.wcs_extension()].header
    ra = header['CRVAL1']
    dec = header['CRVAL2']

    #run sextractor
    cat_name = input_file.replace('.fits','.cat')
    table = run_sextractor(input_file, cat_name, tel, sex_config_dir=sex_config_dir)

    #mask and sort table
    table = table[(table['FLAGS']==0)&(table['EXT_NUMBER']==tel.wcs_extension()+1)]
    if static_mask:
        with fits.open(static_mask) as mask_hdu:
            stat_mask = mask_hdu[0].data
        table = table[(stat_mask[table['YWIN_IMAGE'].astype(int)-1,table['XWIN_IMAGE'].astype(int)-1]!=0)]
    table.sort('MAG_BEST')
    table = table[table['MAG_BEST']<(np.median(table['MAG_BEST'])-np.std(table['MAG_BEST']))]

    #make quads using brightest stars
    ind = list(itertools.combinations(np.arange(4), 2))
    stars, d, ds, ratios = make_quads(table['XWIN_IMAGE'],
        table['YWIN_IMAGE'], use=50)

    #query gaia
    Gaia.ROW_LIMIT = -1
    job = Gaia.cone_search_async(SkyCoord(ra*u.deg, dec*u.deg,frame='icrs'),10*u.arcmin,table_name='gaiaedr3.gaia_source')
    gaiaorig = job.get_results()
    gaiaorig.sort('phot_g_mean_mag')

    # Mask catalog so all data are inside nominal image
    mask = mask_catalog_for_wcs(gaiaorig, wcs.WCS(header), 1, data_x, data_y)
    gaia = gaiaorig[mask]

    #make quads using stars brigher than 19 mag
    gaia = gaia[gaia['phot_g_mean_mag']<20]
    gaiastars, gaiad, gaiads, gaiaratios = make_quads(gaia['ra'], gaia['dec'],
        use=50, sky_coords=True)

    #match quads
    starsx, starsy, gaiastarsra, gaiastarsdec = match_quads(stars,gaiastars,d,gaiad,
            ds,gaiads,ratios,gaiaratios,sky_coords=True)

    #calculate transformation
    gaiax, gaiay = wcs.utils.skycoord_to_pixel(SkyCoord(gaiastarsra,gaiastarsdec, unit='deg'), wcs.WCS(header), 1)
    tform = tf.estimate_transform('euclidean', np.c_[starsx, starsy], np.c_[gaiax, gaiay])

    #apply transformation to header
    header_new = apply_wcs_transformation(header,tform)
    header_new['WCS_REF'] = ('GAIA-DR2', 'Reference catalog for astrometric solution.')
    header_new['WCS_NUM'] = (len(starsx), 'Number of stars used for astrometric solution.')

    #calculate error. This error now uses dvrms() instead of np.median.
    stars_ra, stars_dec = (wcs.WCS(header_new)).all_pix2world(starsx,starsy,1)
    stars_radec = SkyCoord(stars_ra*u.deg,stars_dec*u.deg)
    gaia_radec = SkyCoord(gaiastarsra,gaiastarsdec, unit='deg')
    error = calculate_error(stars_radec, gaia_radec, header_new)
    print('%7.4f rms error in WCS solution'%error)

    #write out new fits
    fits.writeto(input_file.replace('.fits','_wcs.fits'),data,header_new,overwrite=True)

    #end and run time
    t_end = time.time()
    print('Code ran in '+str(t_end-t_start)+' sec')

    return 'Median error on astrometry is %.3f arcsec.'%error

def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('input_file', default=None, help='Name of file.')
    params.add_argument('--telescope', default=None, help='')
    params.add_argument('--proc', default=None, help='')
    params.add_argument('--static_mask', default=None, help='')
    args = params.parse_args()

    solve_wcs(args.input_file,telescope=args.telescope,proc=args.proc,static_mask=args.static_mask)

if __name__ == "__main__":
    main()
