#!/usr/bin/env python

"Python script for WCS solution."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "1.0" #last updated 18/02/2020

import re
import sys
import numpy as np
import os
import datetime
import time
import astropy
import astropy.units.astrophys as u
from astropy.io import fits
import ccdproc
import glob
import argparse
import subprocess
import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, ZScaleInterval
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.io import ascii
import astropy.wcs as wcs
import shutil
from astropy.nddata import CCDData
from astropy.modeling import models
from photutils import make_source_mask, MeanBackground, StdBackgroundRMS, CircularAperture, CircularAnnulus, aperture_photometry, Background2D, MedianBackground, DAOStarFinder
from ccdproc import wcs_project
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.utils import calc_total_error
from astroquery.vizier import Vizier
from astroquery.irsa import Irsa
import warnings
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from scipy.optimize import curve_fit

def plate_sol(data,*p):
    xi,eta=data
    a,b,c = p
    return a*xi+b*eta+c

def poly_fit(data,*p):
    xi,eta=data
    d = [1,xi,eta,np.sqrt(xi**2+eta**2),xi**2,xi*eta,eta**2,xi**3,xi**2*eta,xi*eta**2,eta**3,np.sqrt(xi**2+eta**2)**3]
    if len(d) < len(p):
        return np.sum([p[i]*d[i] for i in range(len(d))])
    if len(d) > len(p):
        return np.sum([p[i]*d[i] for i in range(len(p))])

def ro(x1,y1,theta,x2,y2): 
    x3 = x2*np.cos(theta)-y2*np.sin(theta)-x1*np.cos(theta)+y1*np.sin(theta)+x1 
    y3 = x2*np.sin(theta)+y2*np.cos(theta)-x1*np.sin(theta)-y1*np.cos(theta)+y1 
    return x3,y3

def shift_stars(sources,gaia,wcs_obj,ms):
    offset_median, offset_std, star_idx, gaia_idx = [], [], [], []
    for j in range(len(gaia)):
        for m in range(36):
            x, y = wcs_obj.all_world2pix(gaia['ra'][j],gaia['dec'][j],1)
            dx, dy = sources['xcentroid'][0]-x, sources['ycentroid'][0]-y
            x_all, y_all = wcs_obj.all_world2pix(gaia['ra'],gaia['dec'],1)
            gaia_x_shift, gaia_y_shift = x_all+dx, y_all+dy
            gaia_x, gaia_y = ro(sources['xcentroid'][0], sources['ycentroid'][0], m*2*np.pi/36, gaia_x_shift, gaia_y_shift)
            gaia_x_sub, gaia_y_sub, idxs, mag = [], [], [], []                                                                                             
            for k,n in enumerate(gaia_x): 
                if gaia_x[k]>3 and gaia_x[k]<(np.shape(ms.mask)[1]-3) and gaia_y[k]>3 and gaia_y[k]<(np.shape(ms.mask)[0]-3) and ms.mask[int(gaia_y[k]),int(gaia_x[k])] == False: 
                    gaia_x_sub.append(gaia_x[k]) 
                    gaia_y_sub.append(gaia_y[k])
                    idxs.append(k)
                    mag.append(gaia['phot_g_mean_mag'][k])
            gaia_x_sub = [x for _,x in sorted(zip(mag,gaia_x_sub))]
            gaia_y_sub = [x for _,x in sorted(zip(mag,gaia_y_sub))]
            idxs = np.array([x for _,x in sorted(zip(mag,idxs))])
            if len(idxs) > len(sources):
                gaia_x_sub = gaia_x_sub[0:len(sources)]
                gaia_y_sub = gaia_y_sub[0:len(sources)]
                idxs = idxs[0:len(sources)]
            gaia_ra, gaia_dec = wcs_obj.all_pix2world(gaia_x_sub,gaia_y_sub,1)*u.deg 
            gaia_stars = SkyCoord(gaia_ra,gaia_dec)
            ra, dec = wcs_obj.all_pix2world(sources['xcentroid'],sources['ycentroid'],1)*u.deg 
            stars = SkyCoord(ra,dec)
            if len(gaia_stars) > 3:
                idx, d2, d3 = gaia_stars.match_to_catalog_sky(stars)
                _, median, std = sigma_clipped_stats(d2.deg*3600, sigma=1.5)
                if len(idx[(d2.deg*3600<median+3*std)&(d2.deg*3600>median-3*std)]) > 3:              
                    offset_median.append(median)
                    offset_std.append(std)
                    star_idx.append(idx[(d2.deg*3600<median+3*std)&(d2.deg*3600>median-3*std)])
                    gaia_idx.append(idxs[(d2.deg*3600<median+3*std)&(d2.deg*3600>median-3*std)])
    return offset_median, offset_std, star_idx, gaia_idx

def add_to_header(header,cd1_1,cd1_2,cd2_1,cd2_2,CRPIX1,CRPIX2,crval1,crval2,CRVAL2,p1,p2,rms_ra,rms_dec):
    header['LONPOLE'] = 180
    header['LATPOLE'] = CRVAL2
    header['CRVAL1'] = crval1
    header['CRVAL2'] = crval2
    header['CRPIX1'] = CRPIX1
    header['CRPIX2'] = CRPIX2
    header['CD1_1'] = cd1_1
    header['CD1_2'] = cd1_2
    header['CD2_1'] = cd2_1
    header['CD2_2'] = cd2_2
    # for i in range(len(p1)):
    #     try:
    #         header['PV1_'+str(i)] = p1[i]
    #         header['PV2_'+str(i)] = p2[i]
    #     except:
    #         pass
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['RADESYSa'] = 'ICRS'
    header['WCS_REF'] = ('GAIA-DR2', 'Reference catalog for astrometric solution.')
    header['WCS_STAR'] = (len(p1), 'Number of stars used for astrometric solution.')
    header['RA_RMS'] = (rms_ra, 'RMS of RA fit.')
    header['DEC_RMS'] = (rms_dec, 'RMS of Dec fit.')
    return header

def calculate_error(sources,ms,crval1,crval2,wcs_obj):
    gaia = Irsa.query_region(SkyCoord(crval1*u.deg, crval2*u.deg,frame='icrs'), catalog="gaia_dr2_source", spatial="Cone",radius=5*u.arcmin)
    gaia_x, gaia_y = wcs_obj.all_world2pix(gaia['ra'],gaia['dec'],1)
    gaia_ra, gaia_dec, mag = [], [], []
    for k in range(len(gaia['ra'])):
        if gaia_x[k]>3 and gaia_x[k]<(np.shape(ms.mask)[1]-3) and gaia_y[k]>3 and gaia_y[k]<(np.shape(ms.mask)[0]-3) and ms.mask[int(gaia_y[k]),int(gaia_x[k])] == False: 
            gaia_ra.append(gaia['ra'][k])
            gaia_dec.append(gaia['dec'][k])
            mag.append(gaia['phot_g_mean_mag'][k])
    gaia_ra = [x for _,x in sorted(zip(mag,gaia_ra))]
    gaia_dec = [x for _,x in sorted(zip(mag,gaia_dec))]
    if len(mag) > len(sources):
        gaia_ra = gaia_ra[0:len(sources)]
        gaia_dec = gaia_dec[0:len(sources)]
    gaia_stars = SkyCoord(gaia_ra*u.deg ,gaia_dec*u.deg )
    ra, dec = wcs_obj.all_pix2world(sources['xcentroid'],sources['ycentroid'],1)*u.deg 
    stars = SkyCoord(ra,dec)
    idx, d2, d3 = gaia_stars.match_to_catalog_sky(stars)
    _, median, std = sigma_clipped_stats(d2.deg*3600, sigma=1.5)
    return median

def make_wcs_object(CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip):
        header = fits.Header()
        header['CRVAL1'] = CRVAL1
        header['CRVAL2'] = CRVAL2
        header['CRPIX1'] = CRPIX1
        header['CRPIX2'] = CRPIX2
        header['CD1_1'] = -1*xflip*PIXSCALE*np.cos(np.deg2rad(rot))
        header['CD1_2'] = yflip*PIXSCALE*np.sin(np.deg2rad(rot))
        header['CD2_1'] = xflip*PIXSCALE*np.sin(np.deg2rad(rot))
        header['CD2_2'] = -1*yflip*PIXSCALE*np.cos(np.deg2rad(rot))
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['LONPOLE'] = 180
        header['LATPOLE'] = CRVAL2
        return header, wcs.WCS(header)

def user_inputs():
    CRVAL1 = np.float(input('Please provide RA in degrees: '))
    CRVAL2 = np.float(input('Please provide Dec in degrees: '))
    PIXSCALE = np.float(input('Please provide pixel scale in arcsec: '))/3600
    rot = np.deg2rad(np.float(input('Please provide rotation in degrees (assuming N up and E left): ')))
    xflip = np.float(input('Please provide if there is a flip in x (assuming N up and E left after rotation, -1 for yes, 1 for no): '))
    yflip = np.float(input('Please provide if there is a flip in y (assuming N up and E left after rotation, -1 for yes, 1 for no): '))
    return CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip

def solve_wcs(input_file,ra=None,dec=None,pixscale=None,rot=None,xflip=None,yflip=None):
    user = False
    with fits.open(input_file) as hdr:
        l = 0
        sci_med = hdr[l].data
        try:
            CRVAL1 = hdr[l].header['CRVAL1']
            CRVAL2 = hdr[l].header['CRVAL2']
        except:
            l = 1
            if sci_med is None:
                sci_med = hdr[l].data
            try:
                CRVAL1 = hdr[l].header['CRVAL1']
                CRVAL2 = hdr[l].header['CRVAL2']
            except:
                l = -1
        CRPIX2, CRPIX1 = np.array(np.shape(sci_med))/2
        if ra is not None:
            header, wcs_obj = make_wcs_object(ra,dec,CRPIX1,CRPIX2,pixscale,rot,xflip,yflip)
            user = True
        elif l >= 0:
            header = CCDData.read(input_file,ext=l).header
            wcs_obj = wcs.WCS(hdr[l].header)
        else:
            print('Could not find any coordinates in the header.')
            CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip = user_inputs()
            header, wcs_obj = make_wcs_object(CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip)
            user = True

    sigma_clip = SigmaClip(sigma=3)
    bkg_estimator = MedianBackground()
    _, median, std = sigma_clipped_stats(sci_med, sigma=3.0) 
    m = np.ma.masked_outside(sci_med,median-3*std,median+3*std)
    bkg = Background2D(sci_med, (60, 60), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=m.mask, exclude_percentile=80) 
    ms = np.ma.masked_outside(bkg.background_rms,bkg.background_rms_median-1*np.std(bkg.background_rms),bkg.background_rms_median+1*np.std(bkg.background_rms)) 
    daofind = DAOStarFinder(fwhm=7.0, threshold=5.*std,exclude_border=True,brightest=50)
    sources = daofind(sci_med,mask=ms.mask)

    gaia = []
    while len(gaia) == 0:
        gaia = Irsa.query_region(SkyCoord(CRVAL1*u.deg, CRVAL2*u.deg,frame='icrs'), catalog="gaia_dr2_source", spatial="Cone",radius=5*u.arcmin)
        if len(gaia) == 0:
            if user:
                return 'No GAIA stars found within 5 arcmin of the given RA and Dec for WCS.'
            else:
                print('No GAIA stars found within 5 arcmin of the header RA and Dec for WCS.')
                print('Trying with user values:')
                CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip = user_inputs()
                header, wcs_obj = make_wcs_object(CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip)
                user = True

    offset_median = [5]
    while np.min(offset_median) > 3:
        offset_median, offset_std, star_idx, gaia_idx = shift_stars(sources,gaia,wcs_obj,ms)
        if np.min(offset_median) > 3:
            print('Could not star match within 3 arcsec error with these coordinates.')
            if user:
                print('Trying with best match.')
                break
            else:
                print('Trying with user values:')
                CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip = user_inputs()
                header, wcs_obj = make_wcs_object(CRVAL1,CRVAL2,CRPIX1,CRPIX2,PIXSCALE,rot,xflip,yflip)
                user = True            
                gaia = Irsa.query_region(SkyCoord(CRVAL1*u.deg, CRVAL2*u.deg,frame='icrs'), catalog="gaia_dr2_source", spatial="Cone",radius=5*u.arcmin)

    x = np.array(sources['xcentroid'][star_idx[np.where(offset_median==np.min(offset_median))[0][0]]])
    y = np.array(sources['ycentroid'][star_idx[np.where(offset_median==np.min(offset_median))[0][0]]])
    ra = np.deg2rad(np.array(gaia['ra'][gaia_idx[np.where(offset_median==np.min(offset_median))[0][0]]]))
    dec = np.deg2rad(np.array(gaia['dec'][gaia_idx[np.where(offset_median==np.min(offset_median))[0][0]]]))
    dec_ref = np.deg2rad(CRVAL2)
    xi = np.sin(ra-np.pi)/(np.sin(dec_ref)*np.tan(dec)+np.cos(dec_ref)*np.cos(ra-np.pi))
    eta = (np.tan(dec)-np.tan(dec_ref)*np.cos(ra-np.pi))/(np.tan(dec_ref)*np.tan(dec)+np.cos(ra-np.pi))

    total_error = 0.3
    while np.max(total_error) >= 0.15:
        c1 = curve_fit(plate_sol,(x-CRPIX1,y-CRPIX2),xi-x,p0=(0,0,0))[0]
        c2 = curve_fit(plate_sol,(x-CRPIX1,y-CRPIX2),eta-y,p0=(0,0,0))[0]
        xi_fit = plate_sol((x-CRPIX1,y-CRPIX2),c1[0],c1[1],c1[2])+x
        eta_fit = plate_sol((x-CRPIX1,y-CRPIX2),c2[0],c2[1],c2[2])+y
        error_ra = (np.rad2deg(ra) - np.rad2deg(np.arctan(xi_fit/(np.cos(dec_ref)-eta_fit*np.sin(dec_ref)))+np.pi))*3600
        error_dec = (np.rad2deg(dec) - np.rad2deg(np.arcsin((np.sin(dec_ref)+eta_fit*np.cos(dec_ref))/(np.sqrt(1+xi_fit**2+eta_fit**2)))))*3600
        total_error = np.sqrt(error_ra**2+error_dec**2)
        if len(x[total_error<0.15]) < 3:
            break
        x = x[total_error<0.15]
        y = y[total_error<0.15]
        xi = xi[total_error<0.15]
        eta = eta[total_error<0.15]
        ra = ra[total_error<0.15]
        dec = dec[total_error<0.15]

    xi_0 = plate_sol((0,0),c1[0],c1[1],c1[2])+CRPIX1
    eta_0 = plate_sol((0,0),c2[0],c2[1],c2[2])+CRPIX2
    crval1 = np.rad2deg(np.arctan(xi_0/(np.cos(dec_ref)-eta_0*np.sin(dec_ref)))+np.pi)
    crval2 = np.rad2deg(np.arcsin((np.sin(dec_ref)+eta_0*np.cos(dec_ref))/(np.sqrt(1+xi_0**2+eta_0**2)))) 
    
    xx_shift = plate_sol((1,0),c1[0],c1[1],c1[2])+CRPIX1+1
    xy_shift = plate_sol((0,1),c1[0],c1[1],c1[2])+CRPIX1
    yx_shift = plate_sol((1,0),c2[0],c2[1],c2[2])+CRPIX2
    yy_shift = plate_sol((0,1),c2[0],c2[1],c2[2])+CRPIX2+1
    ra_x = np.rad2deg(np.arctan(xx_shift/(np.cos(dec_ref)-yx_shift*np.sin(dec_ref)))+np.pi)
    ra_y = np.rad2deg(np.arctan(xy_shift/(np.cos(dec_ref)-yy_shift*np.sin(dec_ref)))+np.pi)
    dec_x = np.rad2deg(np.arcsin((np.sin(dec_ref)+yx_shift*np.cos(dec_ref))/(np.sqrt(1+xx_shift**2+yx_shift**2))))
    dec_y = np.rad2deg(np.arcsin((np.sin(dec_ref)+yy_shift*np.cos(dec_ref))/(np.sqrt(1+xy_shift**2+yy_shift**2))))
    cd1_1 = ra_x - crval1
    cd1_2 = ra_y - crval1
    cd2_1 = dec_x - crval2
    cd2_2 = dec_y - crval2

    res1 = xi-xi_fit
    res2 = xi-eta_fit
    p1 = curve_fit(poly_fit,(x-CRPIX1,y-CRPIX2),res1,p0=[0]*len(res1))[0]
    p2 = curve_fit(poly_fit,(x-CRPIX1,y-CRPIX2),res2,p0=[0]*len(res2))[0]

    header = add_to_header(header,cd1_1,cd1_2,cd2_1,cd2_2,CRPIX1,CRPIX2,crval1,crval2,CRVAL2,p1,p2,np.std(error_ra),np.std(error_dec))

    fits.writeto('/Users/kerry/test_wcs3.fits',sci_med,header,overwrite=True)

    error = calculate_error(sources,ms,crval1,crval2,wcs.WCS(header))

    return 'Median error on astrometry is %.3f arcsec.'%error


def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('input_file', default=None, help='Name of file.')
    params.add_argument('--ra', default=None, help='')
    params.add_argument('--dec', default=None, help='')
    params.add_argument('--pixscale', default=None, help='')
    params.add_argument('--rot', default=None, help='')
    params.add_argument('--xflip', default=None, help='')
    params.add_argument('--yflip', default=None, help='')
    args = params.parse_args()

    solve_wcs(args.input_file,ra=args.ra,dec=args.dec,pixscale=args.pixscale,rot=args.rot,xflip=args.xflip,yflip=args.yflip)

if __name__ == "__main__":
    main()