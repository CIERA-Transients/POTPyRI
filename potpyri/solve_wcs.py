#!/usr/bin/env python

"Python script for WCS solution."
"Author: Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "3.12" #last updated 01/10/2021

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
from astropy.visualization import ImageNormalize, ZScaleInterval
from psf import write_out_catalog as woc
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astroquery.gaia import Gaia
from skimage import transform as tf
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from utilities import util
import itertools
import importlib
import params
from utilities import catreg
from colorama import init, Fore, Back, Style
init()

#turn Astropy warnings off
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

def man_wcs(telescope, stack, cat, cat_stars_ra, cat_stars_dec):

    # Load instrument parameters
    global tel
    tel = importlib.import_module(f'params.{telescope}')

    # Get image data and header
    stack_data, stack_header = fits.getdata(stack, header=True)

    # Plot
    fig, ax = plt.subplots(figsize=(7,7))

    # Display image first
    ax.imshow(stack_data, cmap='gray', norm=ImageNormalize(stack_data, interval=ZScaleInterval()))

    # Catalog objects
    x, y = (wcs.WCS(stack_header)).all_world2pix(cat_stars_ra,cat_stars_dec,1)

    for xi,yi, rai, deci, in zip(x,y,cat_stars_ra,cat_stars_dec):
        ax.add_patch(patches.Circle((xi,yi),
                                    radius=6,
                                    edgecolor='g',
                                    alpha=0.8,
                                    facecolor='none',
                                    linewidth=1,
                                    label='Catalog star: RA = %f, Dec = %f'%(rai,deci),
                                    picker=True)
                                    )
    
    # Now extract objects from the image
    _, _, std = sigma_clipped_stats(stack_data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=5.0, threshold=10.*std)
    sources = daofind(np.asarray(stack_data))
    for entry in sources:
        ax.add_patch(patches.Circle((entry['xcentroid'],entry['ycentroid']),
                                    radius=4,
                                    edgecolor='r',
                                    alpha=0.5,
                                    facecolor='none',
                                    linewidth=1,
                                    label='Source star x = %f, y = %f'%(entry['xcentroid'],entry['ycentroid']),
                                    picker=True)
                                    )
    
    # Begin interactive plot
    source_star = []
    cat_star = []
    fig.canvas.mpl_connect('pick_event', lambda event: util.onpick(event,source_star,cat_star,[]))
    print(Fore.RED+'Displaying interactive plot to select star for WCS solution.'+Style.RESET_ALL)
    print(Fore.RED+'First select star (red) and corresponding catalog match (green).'+Style.RESET_ALL)
    print(Fore.RED+'A message will confirm the selection of the star in the terminal.'+Style.RESET_ALL)
    print(Fore.RED+'Note: star selection is turned off in zoom/pan mode (you need to de-select the mode in order to select stars).'+Style.RESET_ALL)
    print(Fore.RED+'Close figure when finished.'+Style.RESET_ALL)
    plt.show()

    # Get user selected stars and save to file
    starsx = np.array([source_star[i][0] for i in range(len(source_star))])
    starsy = np.array([source_star[i][1] for i in range(len(source_star))])
    catstarsra = np.array([cat_star[i][0] for i in range(len(cat_star))])
    catstarsdec = np.array([cat_star[i][1] for i in range(len(cat_star))])
    cat_radec = SkyCoord(catstarsra, catstarsdec,unit=u.deg)
    np.savetxt(stack.replace('.fits','.xy'),np.c_[starsx,starsy])
    np.savetxt(stack.replace('.fits','.radec'),np.c_[catstarsra,catstarsdec])
    check_len = True

    # Prompt user to check
    while check_len:
        con = input(Back.GREEN+'Please check and edit the .xy and .radec files in the case of errors. Hit any key to continue. '+Style.RESET_ALL)
        starsx,starsy = np.loadtxt(stack.replace('.fits','.xy'),unpack=True)
        catstarsra,catstarsdec = np.loadtxt(stack.replace('.fits','.radec'),unpack=True)
        if len(starsx)==len(catstarsra):
            check_len = False
        else:
            print('Length of the x and y star positions and the RA and Dec position are not equal.')
    catx, caty = (wcs.WCS(stack_header)).all_world2pix(catstarsra,catstarsdec,1)

    # This is the new way of doing things, which is hopefully working
    corrected_wcs = wcs.utils.fit_wcs_from_points((starsx, starsy), cat_radec, sip_degree=3)
    corrected_wcshdr = corrected_wcs.to_header(relax=True)
    header_trans = stack_header.copy()
    for kw, val, comment in corrected_wcshdr.cards:
        header_trans.set(kw, value=val,comment=comment)

    stars_ra, stars_dec = (wcs.WCS(header_trans)).all_pix2world(starsx,starsy,1)
    stars_radec = SkyCoord(stars_ra*u.deg,stars_dec*u.deg)
    if 'wcs' not in stack:
        wcs_file = stack.replace('.fits','_wcs.fits')
    else:
        wcs_file = stack.replace('_wcs','').replace('.fits','_wcs.fits')
    error = calculate_error(wcs_file, stars_radec, cat_radec, header_trans)
    make_reg(wcs_file, starsx, starsy)
    wcs_plots(wcs_file, stack_data, starsx, starsy, catx, caty, 'image')
    fits.writeto(wcs_file,stack_data,header_trans,overwrite=True)
    return 'Median rms on astrometry is %.3f arcsec.'%error

def wcs_plots(filename,data,x1,y1,x2,y2,type):
    if type=='image':
        fig, ax = plt.subplots(figsize=(12,12))
        ax.imshow(data, cmap='gray', norm=ImageNormalize(data, interval=ZScaleInterval()))
        [ax.add_patch(patches.Circle((x1[i],y1[i]),radius=7,edgecolor='r',alpha=0.8,facecolor='none',linewidth=1)) for i in range(len(x1))]
        [ax.add_patch(patches.Circle((x2[i],y2[i]),radius=7,edgecolor='g',alpha=0.8,facecolor='none',linewidth=1)) for i in range(len(x2))]
        fig.savefig(filename.replace('.fits','.png'),dpi=1200)
        plt.clf()
    if type=='hist':
        fig = plt.figure(figsize=(12,12))
        left, width = 0.1, 0.75
        bottom, height = 0.1, 0.75
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height , width, 0.12]
        rect_histy = [left + width, bottom, 0.12, height]
        ax_scatter = plt.axes(rect_scatter)
        ax_Histx = plt.axes(rect_histx)
        ax_Histy = plt.axes(rect_histy)
        ax_scatter.tick_params('both',which='both',labelsize=18)
        ax_Histx.tick_params('both',which='both',labelbottom=False,labelsize=16)
        ax_Histy.tick_params('both',which='both',labelleft=False,labelsize=16)
        ax_scatter.plot(x1,y1,marker='o',linestyle='',markersize=10,color='orchid',alpha=0.5)
        ax_Histx.hist(x1,histtype='step',linewidth=2,color='darkorchid')
        ax_Histy.hist(y1,histtype='step',orientation='horizontal',linewidth=2,color='darkorchid')
        ax_scatter.set_xlabel('Delta RA')
        ax_scatter.set_ylabel('Delta Dec')
        fig.savefig(filename.replace('.fits','.hist.png'))
        plt.clf()

def make_reg(filename,xdata,ydata):
    catreg.wcsreg(filename,xdata,ydata)

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
    header_new = header.copy()
    header_new['CRPIX1'] = header[crpix1]-tform.translation[0]
    header_new['CRPIX2'] = header[crpix2]-tform.translation[1]
    try:
        cd = np.array([header_new[cd11],header_new[cd12],header_new[cd21],header_new[cd22]]).reshape(2,2)
    except:
        cd11, cd12, cd21, cd22 = 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'
        cd = np.array([header_new[cd11],header_new[cd12],header_new[cd21],header_new[cd22]]).reshape(2,2)
    cd_matrix = tf.EuclideanTransform(rotation=tform.rotation)
    cd_transformed = tf.warp(cd,cd_matrix)
    header_new['CD1_1'] = cd_transformed[0][0]
    header_new['CD1_2'] = cd_transformed[0][1]
    header_new['CD2_1'] = cd_transformed[1][0]
    header_new['CD2_2'] = cd_transformed[1][1]
    header_new['CTYPE1'] = 'RA---TAN-SIP'
    header_new['CTYPE2'] = 'DEC--TAN-SIP'
    old_keywords = tel.WCS_keywords_old()
    for old in old_keywords:
        try:
            del header_new[old]
        except:
            pass
    return header_new
def ptv_fit(data,*p):
    xi,eta=data.T
    r = np.sqrt(xi**2+eta**2)
    d = [np.ones_like(xi),xi,eta,r,xi**2,xi*eta,eta**2,xi**3,xi**2*eta,xi*eta**2,eta**3,r**3]
    if len(d) <= len(p):
        return np.sum([p[i]*d[i] for i in range(len(d))])
    else:
        return np.sum([p[i]*d[i] for i in range(len(p))])

def apply_wcs_distortion(header,starsx,starsy,gaiastarsra,gaiastarsdec):
    xi = header['CD1_1']*(starsx-header['CRPIX1'])+header['CD1_2']*(starsy-header['CRPIX2'])
    eta = header['CD2_1']*(starsx-header['CRPIX1'])+header['CD2_2']*(starsy-header['CRPIX2'])
    dra, ddec = SkyCoord(header['CRVAL1'],header['CRVAL2'],unit='deg').spherical_offsets_to(SkyCoord(gaiastarsra,gaiastarsdec,unit='deg'))
    if len(xi) > 12:
        p0 = [0,1]+[0]*10
    else:
        p0 = [0,1]+[0]*(len(xi)-2)
    p1 = curve_fit(ptv_fit,np.array((xi,eta)).T,dra.value,p0=np.array(p0))[0]
    p2 = curve_fit(ptv_fit,np.array((eta,xi)).T,ddec.value,p0=np.array(p0))[0]
    header['CTYPE1'] = 'RA---TPV'
    header['CTYPE2'] = 'DEC--TPV'
    del header['PV*']
    for i in range(len(p1)):
        header['PV1_'+str(int(i))] = p1[i]
        header['PV2_'+str(int(i))] = p2[i]
    return header

#function to calculate rms on astrometric solution
#replaced np.median() with dvrms() -David V
def calculate_error(filename, coord1, coord2, hd):
    sep = coord1.separation(coord2)
    d_ra = [s.hms[2] for s in sep]
    d_dec = [s.dms[2] for s in sep]
    wcs_plots(filename,None,d_ra,d_dec,None,None,'hist')
    rms_ra = dvrms(d_ra)
    rms_dec = dvrms(d_dec)
    hd['RA_RMS'] = (rms_ra, 'RMS of RA fit (arsec).')
    hd['DEC_RMS'] = (rms_dec, 'RMS of Dec fit (arcsec).')
    rms_total = dvrms(sep.arcsec)
    return rms_total

def run_sextractor(input_file, cat_name, tel, sex_config_dir='./Config', log=None):

    if not os.path.exists(sex_config_dir+'/config'):
        e = 'ERROR: Could not find source extractor config file!'
        if log:
            log.error(e)
        raise Exception(e)
    elif not os.path.exists(sex_config_dir+'/params'):
        e = 'ERROR: Could not find source extractor param file!'
        if log:
            log.error(e)
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
        if log:
            log.info('SExtractor succesfully run, %d sources extracted.'%len(table))
        return(table)
    else:
        e = 'ERROR: source extractor did not successfully produce {0}'
        if log:
            log.error(e)
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
def solve_wcs(input_file, telescope, sex_config_dir='./Config', static_mask=None, proc=None, log=None, apply_distortion=True):
    #start time
    t_start = time.time()

    if log:
        log.info('Running solve_wcs version: '+str(__version__))

    #import telescope parameter file
    global tel
    try:
        tel = importlib.import_module(f'params.{telescope}')
    except ImportError:
        print('No such telescope file, please check that you have entered the'+\
            ' correct name or this telescope is available.''')
        sys.exit(-1)

    #get data and header info
    hdr = fits.open(input_file)
    data = hdr[tel.wcs_extension()].data
    data_y, data_x = np.shape(data)
    header = hdr[tel.wcs_extension()].header
    ra, dec = (wcs.WCS(header)).all_pix2world(data_x/2,data_y/2,1)

    #run sextractor
    if log:
        log.info('Running SExtractor.')
    cat_name = input_file.replace('.fits','.cat')
    table = run_sextractor(input_file, cat_name, tel, sex_config_dir=sex_config_dir,log=log)

    #mask and sort table
    if log:
        log.info('Masking sources and applying brightness cuts.')
    # Remove bad sources
    table = table[(table['FLAGS']==0)&(table['EXT_NUMBER']==tel.wcs_extension()+1)&(table['MAGERR_BEST']!=99)]
    # Apply static mask?
    if static_mask:
        with fits.open(static_mask) as mask_hdu:
            stat_mask = mask_hdu[0].data
        table = table[(stat_mask[table['YWIN_IMAGE'].astype(int)-1,table['XWIN_IMAGE'].astype(int)-1]!=0)]
    # Sort by brightness
    table.sort('MAG_BEST')
    # Keep only the bright sources
    table = table[table['MAG_BEST']<(np.median(table['MAG_BEST'])-np.std(table['MAG_BEST']))]
    if log:
        log.info('Total of %d usable stars'%len(table))

    #make quads using brightest stars
    if log:
        log.info('Making quads with the brightest 50 stars.')
    ind = list(itertools.combinations(np.arange(4), 2))
    stars, d, ds, ratios = make_quads(table['XWIN_IMAGE'],
        table['YWIN_IMAGE'], use=50)

    #query gaia
    if log:
        log.info('Searching GAIA DR3 for stars within 10 arcmins of WCS in header.')
    Gaia.ROW_LIMIT = -1
    job = Gaia.cone_search_async(SkyCoord(ra*u.deg, dec*u.deg,frame='icrs'),10*u.arcmin,table_name='gaiaedr3.gaia_source')
    gaiaorig = job.get_results()
    gaiaorig.sort('phot_g_mean_mag')
    if log:
        log.info(str(len(gaiaorig))+' GAIA stars found.')

    # Mask catalog so all data are inside nominal image
    if log:
        log.info('Applying image mask to GAIA sources.')
    mask = mask_catalog_for_wcs(gaiaorig, wcs.WCS(header), 1, data_x, data_y)
    gaia = gaiaorig[mask]
    if log:
        log.info(str(len(gaia))+' GAIA stars found within the image footprint.')

    #writing out gaia catalog
    columns = ['ra','dec','phot_g_mean_mag']
    sigfig=[7,7,4]
    if '_wcs.fits' in input_file:
        outfile = input_file.replace('_wcs.fits','_wcs.gaia')
    else:
        outfile = input_file.replace('.fits','_wcs.gaia')
    woc(gaia, input_file, columns, sigfig, outfile, {})

    #make quads using stars brigher than 20 mag
    if log:
        log.info('Making quads with the brightest stars brighter than 20 mag and fainter than 14 mag.')
    gaia = gaia[gaia['phot_g_mean_mag']<20]
    gaia = gaia[gaia['phot_g_mean_mag']>14]
    gaiastars, gaiad, gaiads, gaiaratios = make_quads(gaia['ra'], gaia['dec'],
        use=50, sky_coords=True)

    #match quads
    if log:
        log.info('Matching quads between the detected and cataloged stars.')
    try:
        starsx, starsy, gaiastarsra, gaiastarsdec = match_quads(stars,gaiastars,d,gaiad,
                ds,gaiads,ratios,gaiaratios,sky_coords=True)
        if log:
            log.info('Found '+str(len(starsx))+' unique star matches.')

        #calculate inital transformation
        if log:
            log.info('Calculating the initial transformation between the matched stars.')
        gaiax, gaiay = wcs.utils.skycoord_to_pixel(SkyCoord(gaiastarsra,gaiastarsdec, unit='deg'), wcs.WCS(header), 1)
        tform = tf.estimate_transform('euclidean', np.c_[starsx, starsy], np.c_[gaiax, gaiay])

        #apply initial transformation to header
        if log:
            log.info('Applying the initial transformation to the existing WCS in the header.')
        header_new = apply_wcs_transformation(header,tform)
    except:
        if log:
            log.error('Unique star matching failed.')
        crpix1, crpix2, cd11, cd12, cd21, cd22 = wcs_keyword = tel.WCS_keywords()
        header['CD1_1'] = header[cd11]
        header['CD1_2'] = header[cd12]
        header['CD2_1'] = header[cd21]
        header['CD2_2'] = header[cd22]
        header['CTYPE1'] = 'RA---TAN-SIP'
        header['CTYPE2'] = 'DEC--TAN-SIP'
        old_keywords = tel.WCS_keywords_old()
        for old in old_keywords:
            try:
                del header[old]
            except:
                pass
        header_new = header

    #matching all stars to catalog
    if log:
        log.info('Matching all stars to the catalog.')
    stars_ra, stars_dec = (wcs.WCS(header_new)).all_pix2world(table['XWIN_IMAGE'],table['YWIN_IMAGE'],1)
    stars_radec = SkyCoord(stars_ra*u.deg,stars_dec*u.deg)
    gaia_radec = SkyCoord(gaia['ra'],gaia['dec'], unit='deg')
    idx, d2, _ = gaia_radec.match_to_catalog_sky(stars_radec)
    match = d2<5.0*u.arcsec
    idx = idx[match]
    starx_match = [table['XWIN_IMAGE'][x] for x in idx]
    stary_match = [table['YWIN_IMAGE'][x] for x in idx]
    gaia_radec = SkyCoord(gaia[match]['ra'],gaia[match]['dec'], unit='deg')
    gaiax_match, gaiay_match = wcs.utils.skycoord_to_pixel(gaia_radec, wcs.WCS(header_new), 1)
    if log:
        log.info('Found '+str(len(starx_match))+' star matches within 5".')
    
    #calculate and apply full transformation
    if log:
        log.info('Calculating the full transformation between the matched stars.')
    
    # This is the new way of doing things, which is hopefully working
    corrected_wcs = wcs.utils.fit_wcs_from_points((np.array(starx_match), np.array(stary_match)),
                                                   gaia_radec, sip_degree=3)
    corrected_wcshdr = corrected_wcs.to_header(relax=True)
    header_dist = header_new.copy()
    for kw, val, comment in corrected_wcshdr.cards:
        header_dist.set(kw, value=val,comment=comment)
    header_dist['WCS_REF'] = ('GAIA-DR3', 'Reference catalog for astrometric solution.')
    header_dist['WCS_NUM'] = (len(starx_match), 'Number of stars used for astrometric solution.')

    #write region file and plots
    if '_wcs.fits' in input_file:
        wcs_file = input_file
    else:
        wcs_file = input_file.replace('.fits','_wcs.fits')
    make_reg(wcs_file, starx_match, stary_match)
    wcs_plots(wcs_file, data, starx_match, stary_match, gaiax_match, gaiay_match, 'image')

    #calculate error. This error now uses dvrms() instead of np.median.
    if log:
        log.info('Calculating the rms on the astrometry.')
    stars_ra, stars_dec = (wcs.WCS(header_dist)).all_pix2world(starx_match,stary_match,1)
    stars_radec = SkyCoord(stars_ra*u.deg,stars_dec*u.deg)
    error = calculate_error(wcs_file, stars_radec, gaia_radec, header_dist)

    #write out new fits
    fits.writeto(wcs_file,data,header_dist,overwrite=True)

    #end and run time
    t_end = time.time()
    if log:
        log.info('Solve_wcs ran in '+str(t_end-t_start)+' sec')
    else:
        print('Solve_wcs ran in '+str(t_end-t_start)+' sec')

    return 'Median rms on astrometry is %.3f arcsec.'%error

def clean_up_astrometry(directory, file, exten):
    filelist = [file.replace(exten,'.axy'),
                file.replace(exten,'.corr'),
                file.replace(exten,'-indx.xyls'),
                file.replace(exten,'.match'),
                file.replace(exten,'.rdls'),
                file.replace(exten,'.solved'),
                file.replace(exten,'.wcs')]

    filelist = [os.path.join(directory, f) for f in filelist]

    for f in filelist:
        if os.path.exists(f):
            os.remove(f)

def solve_astrometry(directory, file, tel, replace=True, log=None):

    # Starting solve, print file and directory for reference
    fullfile = directory + '/' + file
    message = 'Trying to solve file={file}'
    if log: log.info(message.format(file=fullfile))

    if not os.path.exists(fullfile):
        return(False)

    hdu = fits.open(fullfile)
    data = hdu[0].data
    header = hdu[0].header
    hkeys = list(header.keys())

    exten = '.'+file.split('.')[-1]
    if not replace:
        if os.path.exists(fullfile.replace(exten,'.solved.fits')):
            if log: log.info('SUCCESS: solved {0}'.format(fullfile))
            return(True)

    exten = '.'+file.split('.')[-1]

    check_pairs = [('RA','DEC'),('CRVAL1','CRVAL2'),('OBJCTRA','OBJCTDEC')]
    coord = None

    for pair in check_pairs:
        if pair[0] in header.keys() and pair[1] in header.keys():
            ra = header[pair[0]]
            dec = header[pair[1]]
            coord = util.parse_coord(ra, dec)
            if coord:
                break

    if not coord:
        error = 'Could not parse RA/DEC from header of {file}'
        if log: log.error(error.format(file=file))
        return(False)

    if 'solved.fits' in fullfile:
        newfile = 'tmp.fits'
    else:
        newfile = fullfile.replace(exten,'.solved.fits')

    scale = tel.pixscale()
    scale_high = scale * 1.2
    scale_low = scale * 0.8

    cmd = 'solve-field'
    args = '--scale-units arcsecperpix '
    args += '--scale-low {lo} --scale-high {hi} --ra {ra} --dec {dec} '
    args += ' --radius 0.5 --no-plots --overwrite -N {newfile} --dir {outdir} '

    args = args.format(ra=coord.ra.degree, dec=coord.dec.degree,
        outdir=directory, newfile=newfile, lo=scale_low, hi=scale_high)

    extra_opts = '--downsample 2 --no-verify --odds-to-tune-up 1e4 --objs 15'

    tries = 1
    good = False
    while tries < 4 and not good:
        if log: log.info('Try #{n} with astrometry.net...'.format(n=tries))
        input_args = args + extra_opts
        if log: log.info(input_args)

        process = [cmd,fullfile]+input_args.split()

        FNULL = open(os.devnull, 'w')

        p = subprocess.Popen(process, stdout=FNULL, stderr=subprocess.STDOUT)
        try:
            p.wait(30)
        except subprocess.TimeoutExpired:
            p.kill()

        if os.path.exists(newfile):
            good = True
        else:
            tries += 1
            if 'downsample' in extra_opts:
                extra_opts='--objs 15'
            else:
                extra_opts=''

    file_exists=os.path.exists(newfile)
    if log: log.info(f'{newfile} exists: {file_exists}')

    if os.path.exists(newfile):
        if log: log.info('SUCCESS: solved {0}'.format(fullfile))

        if replace or newfile=='tmp.fits':
            output_file = fullfile
            os.rename(newfile, output_file)
        else:
            output_file = fullfile.replace(exten,'.solved.fits')

        hdu = fits.open(output_file)
        for i,h in enumerate(hdu):

            if 'COMMENT' in hdu[i].header.keys():
                del hdu[i].header['COMMENT']
            if 'HISTORY' in hdu[i].header.keys():
                del hdu[i].header['HISTORY']

            # Delete other WCS header keys
            wcs_key = ['CSYER1','CSYER2','CRDER1','CRDER2','CD1_1','CD1_2','CD2_1',
                'CD2_2','CRPIX1','CRPIX2','CUNIT1','CUNIT2','EQUINOX','RADESYS',
                'CNAME1','CNAME2','CTYPE1','CTYPE2','WCSNAME','CRVAL1','CRVAL2']

            for key in [w + 'C' for w in wcs_key]:
                if key in list(hdu[i].header.keys()):
                    del hdu[i][key]

            for key in [w + 'S' for w in wcs_key]:
                if key in list(hdu[i].header.keys()):
                    del hdu[i][key]

        try:
            hdu.writeto(output_file, overwrite=True, output_verify='silentfix')
        except TypeError:
            if log: log.error('FAILURE: could not properly save file {0}'.format(
                fullfile))

            clean_up_astrometry(directory, file, exten)
            return(False)

        clean_up_astrometry(directory, file, exten)

        return(True)

    else:
        #if log: log.info('FAILURE: did not solve {0}'.format(fullfile))
        print('FAILURE: did not solve {0}'.format(fullfile))
        clean_up_astrometry(directory, file, exten)
        return(False)

def main():
    params = argparse.ArgumentParser(description='Path of data.')
    params.add_argument('input_file', default=None, help='Name of file.')
    params.add_argument('--telescope', default=None, help='')
    params.add_argument('--proc', default=None, help='')
    params.add_argument('--static_mask', default=None, help='')
    params.add_argument('--log', default=None, help='')
    args = params.parse_args()

    solve_wcs(args.input_file,telescope=args.telescope,proc=args.proc,static_mask=args.static_mask,log=args.log)

if __name__ == "__main__":
    main()
