#!/usr/bin/env python

"Python script for WCS solution."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "2.7" #last updated 18/02/2021

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
from scipy.optimize import curve_fit
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
def plate_sol(data,*p):
    xi,eta=data
    a,b,c = p
    return a*xi+b*eta+c

#function to remove old WCS keywords from header
def strip_header(hd):
    wcs_obj = wcs.WCS(hd)
    wcs_keywords_old = wcs_obj.to_header()
    for wcs_keyword in wcs_keywords_old:
        try:
            del hd[wcs_keyword]
        except:
            pass
    return hd

#function to add new astrometric solution to header
def add_to_header(hd,crpix1,crpix2,c1,c2,num):
    hd['WCSAXES'] = 2
    hd['CRPIX1'] = crpix1
    hd['CRPIX2'] = crpix2
    hd['CRVAL1'] = c1[2]
    hd['CRVAL2'] = c2[2]
    hd['CUNIT1'] = 'deg'
    hd['CUNIT2'] = 'deg'
    hd['CTYPE1'] = 'RA---TNX'
    hd['CTYPE2'] = 'DEC--TNX'
    hd['CD1_1'] = c1[0]
    hd['CD1_2'] = c1[1]
    hd['CD2_1'] = c2[0]
    hd['CD2_2'] = c2[1]
    WAT0_001, WAT1_001, WAT1_002, WAT1_003, WAT1_004, WAT1_005, \
        WAT2_001, WAT2_002, WAT2_003, WAT2_004, WAT2_005 = tel.WCS_keywords()
    hd['WAT0_001'] = WAT0_001
    hd['WAT1_001'] = WAT1_001
    hd['WAT1_002'] = WAT1_002
    hd['WAT1_003'] = WAT1_003
    hd['WAT1_004'] = WAT1_004
    hd['WAT1_005'] = WAT1_005
    hd['WAT2_001'] = WAT2_001
    hd['WAT2_002'] = WAT2_002
    hd['WAT2_003'] = WAT2_003
    hd['WAT2_004'] = WAT2_004
    hd['WAT2_005'] = WAT2_005
    hd['RADESYSa'] = 'ICRS'
    hd['WCS_REF'] = ('GAIA-DR2', 'Reference catalog for astrometric solution.')
    hd['WCS_QUAD'] = (num, 'Number of stars used for astrometric solution.')
    return hd

#function to calculate rms on astrometric solution
#replaced np.median() with dvrms() -David V
def calculate_error(d, coord1, coord2, hd):
    rms_total = dvrms(d)
    sep = [coord1[i].separation(coord2[i]) for i in range(len(coord1))]
    d_ra = [s.hms[2] for s in sep]
    d_dec = [s.dms[2] for s in sep]
    rms_ra = dvrms(d_ra)
    rms_dec = dvrms(d_dec)
    hd['RA_RMS'] = (rms_ra, 'RMS of RA fit (arsec).')
    hd['DEC_RMS'] = (rms_dec, 'RMS of Dec fit (arcsec).')
    return rms_total

def run_sextractor(input_file, cat_name, tel, sex_config_dir='./Config'):

    if not os.path.exists(sex_config_dir+'/config'):
        e = 'ERROR: Could not find source extractor config file!'
        raise Exception(e)
    elif not os.path.exists(sex_config_dir+'/params'):
        e = 'ERROR: Could not find source extractor param file!'
        raise Exception(e)

    params = []
    with open(sex_config_dir+'/params') as f:
        for line in f:
            params.append(line.split()[0].strip())

    config = sex_config_dir+'/config'

    cmd=['sex','-c',config,input_file,'-CATALOG_NAME',cat_name,'-FLAG_IMAGE',tel.static_mask(),'-FLAG_TYPE','AND']
    subprocess.call(cmd)

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
                gaia_d = [np.hypot(gaiam[k][m_ind[i][0]][0]-stars[k][m_ind[i][1]][0], gaiam[k][m_ind[i][0]][1]-stars[k][m_ind[i][1]][1]) for i in range(6)]
                ratio_quads.append(np.sum(np.abs((gaia_d/star_d)-1)))
        starsx.append([starsm[k][i][0] for i in range(4)])
        starsy.append([starsm[k][i][1] for i in range(4)])
        gaiastarsra.append([gaiam[k][gaia_ind[np.argmin(ratio_quads)][i]][0] for i in range(4)])
        gaiastarsdec.append([gaiam[k][gaia_ind[np.argmin(ratio_quads)][i]][1] for i in range(4)])
    starsx_unq = list(dict.fromkeys(np.concatenate(starsx)))
    starsy_unq = list(dict.fromkeys(np.concatenate(starsy)))
    gaiastarsra_unq = list(dict.fromkeys(np.concatenate(gaiastarsra)))
    gaiastarsdec_unq = list(dict.fromkeys(np.concatenate(gaiastarsdec)))

    return np.array(starsx_unq), np.array(starsy_unq), np.array(gaiastarsra_unq), np.array(gaiastarsdec_unq)

#main function to calculate astrometric solution
def solve_wcs(input_file, telescope, sex_config_dir='./Config'):
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
    table = table[(table['FLAGS']==0)&(table['IMAFLAGS_ISO']==0)
            &(table['EXT_NUMBER']==tel.wcs_extension()+1)]
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

    #load ref pixels
    crpix1, crpix2 = tel.ref_pix()

    #solve plate sol
    c1 = curve_fit(plate_sol,(starsx-crpix1,starsy-crpix2),gaiastarsra,p0=(0,0,ra))[0]
    c2 = curve_fit(plate_sol,(starsx-crpix1,starsy-crpix2),gaiastarsdec,p0=(0,0,dec))[0]

    #strip header of WCS
    hd = strip_header(header)

    #add new WCS header keywords
    header_new = add_to_header(hd,crpix1,crpix2,c1,c2,len(starsx))

    #match stars
    stars_ra, stars_dec = (wcs.WCS(header_new)).all_pix2world(table[:50]['XWIN_IMAGE'],
        table[:50]['YWIN_IMAGE'],1)
    stars = SkyCoord(stars_ra*u.deg,stars_dec*u.deg)
    mask = mask_catalog_for_wcs(gaiaorig, wcs.WCS(header_new), 1, data_x, data_y)
    gaia = gaiaorig[mask]
    gaia_stars = SkyCoord(gaia['ra'],gaia['dec'], unit='deg')

    idx, d2, d3 = gaia_stars.match_to_catalog_sky(stars)
    idx_m = {}
    for i in range(len(idx)):
        if len(np.where(idx==idx[i])[0]) > 1:
            idx_db = np.where(idx==idx[i])[0]
            d2_m = idx_db[d2.arcsec[idx_db]==np.min(d2.arcsec[idx_db])][0]
        else:
            d2_m = i
        if d2[d2_m].arcsec < 15:
            idx_m.update({idx[i]:d2_m})
    d2_m = [d2[idx_m[i]].arcsec for i in idx_m]
    stars_idx = [stars[i] for i in idx_m]
    gaia_idx = [gaia_stars[idx_m[i]] for i in idx_m]

    #calculate error. This error now uses dvrms() instead of np.median. 
    error = calculate_error(d2_m, stars_idx, gaia_stars, header_new)
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
    args = params.parse_args()

    solve_wcs(args.input_file,telescope=args.telescope)

if __name__ == "__main__":
    main()
