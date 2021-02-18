#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 14:39:02 2021

@author: dvelasco
"""
import numpy as np
from astropy.coordinates import SkyCoord
    
#To make reg file from test.cat from sextraxt

def secat():
    alpha,delta,flags,flags_m,spread_m = np.loadtxt('./MOSFIRE/MF.20181018.22663_red.cat',usecols=(5,6,7,8,9),unpack=True)
    sex_ra, sex_dec = [],[]
    for k in range(len(alpha)):
        if (flags[k]==0):
            if (flags_m[k]==0):
                if (spread_m[k]>0.0):
                    sex_ra.append(alpha[k])
                    sex_dec.append(delta[k])
    stellae = SkyCoord(sex_ra, sex_dec, unit='deg')
    stellae_ra = [c.to_string(style='hmsdms', sep=':').split()[0] for c in stellae]
    stellae_dec = [c.to_string(style='hmsdms', sep=':').split()[1] for c in stellae]
    f= open("/Users/dvelasco/Imaging_pipelines/secat.reg","w+")  
    line1 = '# Region file format: DS9 version 4.1'
    line2 = 'global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
    line3= 'fk5'
    #print(len(stellae_ra))
    f.write(line1+"\n")
    f.write(line2+"\n") 
    f.write(line3+"\n") 
    for i in range(len(stellae_ra)):
        coords = 'point('+stellae_ra[i]+','+stellae_dec[i]+') # point=circle'
        f.write(coords+"\n") 
    f.close()
    