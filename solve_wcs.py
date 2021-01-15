#!/usr/bin/env python

"Python script for WCS solution."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "2.0" #last updated 15/01/2021

import sys
import numpy as np
import time
import astropy
from astropy.io import fits
import argparse
import subprocess
import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.wcs as wcs
from astroquery.irsa import Irsa
import warnings
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from scipy.optimize import curve_fit
from scipy import stats
import itertools
import importlib
import tel_params

#turn Astropy warnings off
warnings.simplefilter('ignore', category=AstropyWarning)

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
    WAT0_001, WAT1_001, WAT1_002, WAT1_003, WAT1_004, WAT1_005, WAT2_001, WAT2_002, WAT2_003, WAT2_004, WAT2_005 = tel.WCS_keywords()
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
    hd['WCS_QUAD'] = (num, 'Number of quads used for astrometric solution.')
    return hd

#function to calculate rms on astrometric solution
def calculate_error(d, coord1, coord2, hd):
    rms_total = np.median(d)
    sep = [coord1[i].separation(coord2[i]) for i in range(len(coord1))]
    d_ra = [s.hms[2] for s in sep]
    d_dec = [s.dms[2] for s in sep]
    rms_ra = np.median(d_ra)
    rms_dec = np.median(d_dec)
    hd['RA_RMS'] = (rms_ra, 'RMS of RA fit.')
    hd['DEC_RMS'] = (rms_dec, 'RMS of Dec fit.')
    return rms_total

#main function to calculate astrometric solution
def solve_wcs(input_file,telescope):
    #start time
    t_start = time.time()
    #import telescope parameter file
    global tel
    try:
        tel = importlib.import_module('tel_params.'+telescope) 
    except ImportError:
        print('No such telescope file, please check that you have entered the correct name or this telescope is available.')
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
    subprocess.call(['sex','-c','./Config/config',input_file,'-CATALOG_NAME',cat_name])

    #read in stars and sort
    mag, xwin, ywin, flags = np.loadtxt(cat_name,usecols=(1,3,4,7),unpack=True)
    mags = sorted(mag[(flags==0)])
    xwins = [x for y, x in sorted(zip(mag[(flags==0)],xwin[(flags==0)]))]
    ywins = [x for y, x in sorted(zip(mag[(flags==0)],ywin[(flags==0)]))]
    
    #make quads using 10 brightest
    ind = list(itertools.combinations(np.arange(4), 2))
    stars = list(itertools.combinations(zip(xwins[0:10],ywins[0:10]), 4))
    dx = [[stars[j][ind[i][0]][0]-stars[j][ind[i][1]][0] for i in range(6)] for j in range(len(stars))]
    dy = [[stars[j][ind[i][0]][1]-stars[j][ind[i][1]][1] for i in range(6)] for j in range(len(stars))]
    d = [np.hypot(dx[j],dy[j]) for j in range(len(stars))]
    ds = [sorted(d[j]) for j in range(len(stars))]
    ratios = [[ds[j][i]/ds[j][5] for i in range(5)] for j in range(len(stars))]
    
    #query gaia
    gaia = Irsa.query_region(SkyCoord(ra*u.deg, dec*u.deg,frame='icrs'), catalog="gaia_dr2_source", spatial="Cone",radius=10*u.arcmin)
    gaiara = np.array([x for y, x in sorted(zip(gaia['phot_g_mean_mag'],gaia['ra']))])
    gaiadec = np.array([x for y, x in sorted(zip(gaia['phot_g_mean_mag'],gaia['dec']))])
    gaiax, gaiay = (wcs.WCS(header)).all_world2pix(gaiara,gaiadec,1)
    gaiarain = gaiara[(gaiax>0)&(gaiax<(data_x-30))&(gaiay>0)&(gaiay<(data_y-30))]
    gaiadecin = gaiadec[(gaiax>0)&(gaiax<(data_x-30))&(gaiay>0)&(gaiay<(data_y-30))]
    
    #make quads using 10 brightest
    gaiastars = list(itertools.combinations(zip(gaiarain[0:10],gaiadecin[0:10]), 4))
    gaiadra = [[gaiastars[j][ind[i][0]][0]-gaiastars[j][ind[i][1]][0] for i in range(6)] for j in range(len(gaiastars))]
    gaiaddec = [[gaiastars[j][ind[i][0]][1]-gaiastars[j][ind[i][1]][1] for i in range(6)] for j in range(len(gaiastars))]
    gaiad = [np.hypot(gaiadra[j],gaiaddec[j]) for j in range(len(gaiastars))]
    gaiads = [sorted(gaiad[j]) for j in range(len(gaiastars))]
    gaiaratios = [[gaiads[j][i]/gaiads[j][5] for i in range(5)] for j in range(len(gaiastars))]
    
    #match quads
    l1 = [[ratios[k][0]/gaiaratios[j][0] for j in range(len(gaiastars))] for k in range(len(stars))]
    l2 = [[ratios[k][4]/gaiaratios[j][4] for j in range(len(gaiastars))] for k in range(len(stars))]
    l3 = [np.hypot(np.abs(np.array(l1[k])-1),np.abs(np.array(l2[k])-1)) for k in range(len(stars))]
    indm = np.array([(k,np.where(l3[k]==np.min(l3[k]))[0][0]) for k in range(len(stars))])
    dr = np.array([[(gaiads[i[1]][j]/ds[i[0]][j])*3600 for j in range(6)] for i in indm])
    indmm = np.array([indm[i] for i in range(len(indm)) if all(np.abs((dr[i]/tel.pixscale())-1)<0.05)])
    starsm = [stars[i[0]] for i in indmm]
    gaiam = [gaiastars[i[1]] for i in indmm]
    dm = [d[i[0]] for i in indmm]
    gaiadm = [gaiad[i[1]] for i in indmm]
    starsx = []
    starsy = []
    gaiastarsra = []
    gaiastarsdec = []
    for k in range(len(starsm)):
        indxyrdm_good = True
        starsi = [x for y,x in sorted(zip(dm[k],ind))]
        gaiai = [x for y,x in sorted(zip(gaiadm[k],ind))]
        inds = {0:[],1:[],2:[],3:[]}
        for ij in range(4):
            li = []
            for i in range(6):
                if ij in starsi[i]:
                    [li.append(gaiai[i][j]) for j in range(2)]
            for idxm in inds:
                try:
                    li = list(filter((inds[idxm][0]).__ne__, li))
                except IndexError:
                    pass
            indxyrdm, cnts = stats.mode(li)
            if len(li) > 1 and cnts[0] == 1:
                indxyrdm_good = False
                continue
            inds[ij].append(indxyrdm[0])
        if indxyrdm_good:
            starsx.append([starsm[k][i][0] for i in range(4)])
            starsy.append([starsm[k][i][1] for i in range(4)])
            gaiastarsra.append([gaiam[k][stats.mode(inds[i])[0][0]][0] for i in inds])
            gaiastarsdec.append([gaiam[k][stats.mode(inds[i])[0][0]][1] for i in inds])
    
    #load ref pixels
    crpix1, crpix2 = tel.ref_pix()

    #solve plate sol
    c1 = curve_fit(plate_sol,(np.concatenate(starsx)-crpix1,np.concatenate(starsy)-crpix2),np.concatenate(gaiastarsra),p0=(0,0,0))[0]
    c2 = curve_fit(plate_sol,(np.concatenate(starsx)-crpix1,np.concatenate(starsy)-crpix2),np.concatenate(gaiastarsdec),p0=(0,0,0))[0]
    
    #strip header of WCS
    hd = strip_header(header)

    #add new WCS header keywords
    header_new = add_to_header(hd,crpix1,crpix2,c1,c2,len(starsx))

    #match stars
    stars_ra, stars_dec = (wcs.WCS(header_new)).all_pix2world(xwins[0:50],ywins[0:50],1)
    stars = SkyCoord(stars_ra*u.deg,stars_dec*u.deg)
    gaiaxn, gaiayn = (wcs.WCS(header_new)).all_world2pix(gaiara,gaiadec,1)
    gaiarainn = gaiara[(gaiaxn>0)&(gaiaxn<(data_x-30))&(gaiayn>0)&(gaiayn<(data_y-30))]
    gaiadecinn = gaiadec[(gaiaxn>0)&(gaiaxn<(data_x-30))&(gaiayn>0)&(gaiayn<(data_y-30))]
    gaia_stars = SkyCoord(gaiarainn*u.deg,gaiadecinn*u.deg)
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
    
    #calculate error
    error = calculate_error(d2_m, stars_idx, gaia_stars, header_new)

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