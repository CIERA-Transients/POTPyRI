#!/usr/bin/env python

"Python script to convert SExtractor files to ds9 region files."
"Author: David Velasco, Kerry Paterson"
"This project was funded by AST "
"If you use this code for your work, please consider citing ."

__version__ = "1.1" #last updated 08/09/2021

import sys
import numpy as np
from astropy.coordinates import SkyCoord
from utilities.util import *

def secat(filename,cattype):
    header, table = import_catalog(filename)
    if cattype == 'zpt':
        x,y = 'Xpos', 'Ypos'
    else:
        x,y = 'XWIN_IMAGE', 'YWIN_IMAGE'
    ext = filename.split('.')[-1]
    f = open(filename.replace(ext,ext+'.reg'),"w+")  
    f.write('# Region file format: DS9 version 4.1\n')
    f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n') 
    f.write('image\n') 
    for i in range(len(table)):
        f.write('circle('+str(table[x][i])+','+str(table[y][i])+',10)\n')
    f.close()

def wcsreg(filename,xdata,ydata):
    ext = filename.split('.')[-1]
    f = open(filename.replace(ext,ext+'.wcs.reg'),"w+")
    f.write('# Region file format: DS9 version 4.1\n')
    f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n') 
    f.write('image\n') 
    for i in range(len(xdata)):
        f.write('circle('+str(xdata[i])+','+str(ydata[i])+',10)\n')
    f.close()