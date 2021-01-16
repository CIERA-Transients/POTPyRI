#!/usr/bin/env python
import glob
import sys
import pdb
import os
import csv
from astropy.io import fits
import numpy as np
from collections import defaultdict
from tmath.pydux.dms_string_to_dd import dms_string_to_dd
from tmath.pydux.airmass import airmass
from tmath.pydux.jdcnv import jdcnv
from tmath.pydux.dateparser import dateparser
from tmath.pydux.parsehourangle import parsehourangle
from tmath.pydux.ddaytodhms import ddaytodhms
from tmath.pydux.hacalc import hacalc
from tmath.pydux.julday import julday
from tmath.pydux.calday import calday
from tmath.pydux.greensid import greensid
from tmath.pydux.dhtohms import dhtohms
from tmath.pydux.parallactic import parallactic
from tmath.pydux.dayofyear import dayofyear
from tmath.pydux.leapyear import leapyear
import inspect

reduxdir=inspect.getfile(calday)
reduxdir=reduxdir.rsplit('/',1)[0]+'/'

degrad=180.0/np.pi

# if (len(sys.argv) < 2):
#     sys.exit('Please include files in the command line (in single quotes)!')
    
# filefilter=sys.argv[1]
filefilter='d*.fits'

# print(sys.argv[1])

obsfile=reduxdir+'observatories'
csvfile=open(obsfile,'r')
csv_reader=csv.reader(csvfile)
obsdict=defaultdict(list)
obsflag=True
for row in csv_reader:
    obsdict[row[0]]=row[1:]

files=glob.glob(filefilter)

print(files)

reconfile=files[0]

recon=fits.open(reconfile)

try:
    observat=recon[0].header['OBSERVAT'].strip().lower()
except KeyError:
    print('OBSERVAT keyword not present in header')
    observat=input('Please enter observatory: ')
    observat=observat.strip().lower()
    obsflag=False

recon.close()
while (observat not in obsdict):
    print('{} is not in the observatory database'.format(observat))
    print('These are the choices: ')
    print(obsdict.keys())
    observat=input('Please enter observatory (q ends program): ')
    observat=observat.strip().lower()
    obsflag=False
    if (observat == 'q'):
        sys.exit()
    

for f in files:
    fitsfile=fits.open(f)
    fitsfile[0].header.set('OBSERVAT',observat)
    latst=obsdict[observat][2]
    lat=dms_string_to_dd(latst)
    latrad=lat/degrad
    lonst=obsdict[observat][1]
    lon=dms_string_to_dd(lonst)
    lonrad=lon/degrad
    decst=fitsfile[0].header['DEC']
    dec=dms_string_to_dd(decst)
    deltarad=dec/degrad
    rast=fitsfile[0].header['RA']
    ra=dms_string_to_dd(rast)
    rarad=ra*15./degrad
    hainput=fitsfile[0].header['HA']
    ha=parsehourangle(hainput)
    harad=ha*15./degrad
    
    # deal with lris header nonsense
    exptime=float(fitsfile[0].header['EXPTIME'])
    dateobs=fitsfile[0].header['DATE-OBS'].strip().split('T')[0]

    year,month,day=dateparser(dateobs)
    if ('UT' in fitsfile[0].header):
        utst=fitsfile[0].header['UT'].strip()
        if ('T' in utst):
            utst=utst.split('T')[1]
    elif ('UT-TIME' in fitsfile[0].header):
            utst=fitsfile[0].header['UT-TIME'].strip()
        
    else:
        if 'T' in fitsfile[0].header['DATE-OBS']:
            utst=fitsfile[0].header['DATE-OBS'].strip().split('T')[1]
        else:
            utst=fitsfile[0].header['TIME-OBS'].strip()
            

            
    hour,minute,second=utst.split(':')
    dhour=int(hour)+(int(minute)+float(second)/60.)/60.
    doy=dayofyear(year,month,day)
    islyr=leapyear(year)
    epoch=year+doy/(365.+islyr)
    jdn=julday(year,month,day)
    jd=jdn+dhour/24.
    jd=jd+exptime/60./60./24./2.
    if (observat == 'flwo'):
        jd = jd - exptime/60./60./24.
    ayear, amonth, adaydec=calday(jd)
    aday,ahour,aminute,asecond=ddaytodhms(adaydec)
    autst=str(ahour)+':'+str(aminute)+':'+str(float(int(asecond*10000))/10000.)
    ajd=julday(ayear, amonth, aday)
    gsid=greensid(ajd,ahour,aminute,asecond)
    localsid=(gsid*15./degrad-lonrad)*degrad/15.
    if (localsid < 0):
        localsid+=24.
    if (localsid > 24.):
        localsid-=24.
    sidh, sidm, sids=dhtohms(localsid)
    sids=float(int(sids*10000))/10000.
    sidst=str(sidh)+':'+str(sidm)+':'+str(sids)
    harad=hacalc(rarad,lonrad,gsid*15./degrad)
    harad=np.fmod(harad,2*np.pi)
    if (harad < -np.pi):
        harad+=2*np.pi
    hah, ham, has=dhtohms(harad*degrad/15.)
    has=float(int(has*10000))/10000.
    hast=str(hah)+':'+str(ham)+':'+str(has)
    el=np.arcsin(np.sin(latrad)*np.sin(deltarad)+np.cos(latrad)*np.cos(deltarad)*np.cos(harad))*degrad
    airm=airmass(el)
    optpa=parallactic(latrad,deltarad,harad)
    print(jd, ajd)
    print('UTMIDDLE',autst)
    print('HA',hast)
    print('ST',sidst)
    print('AIRMASS', airm)
    print('OPT-PA',optpa)
    print('EPOCH',epoch)
    fitsfile[0].header.set('UTMIDDLE',autst)
    fitsfile[0].header.set('HA',hast)
    fitsfile[0].header.set('ST',sidst)
    fitsfile[0].header.set('AIRMASS',airm)
    fitsfile[0].header.set('OPT_PA',optpa)
    fitsfile[0].header.set('EPOCH',epoch)
    fitsfile.writeto(f,overwrite=True)
    fitsfile.close()
# UTMIDDLE
# flwo DISKFILE -> OBSNUM, DISPAXIS=1, st, epoch, optpa, setairmass
#kast EXPOSURE-> EXPTIME, TIME->UT DISPAXIS=1

#mmto EXPOSURE -> EXPTIME st, epoch, ha, optpa, setairmass

#gemini RA/DEC are in decimal degrees

#gemini/lco N&S

# greensid (jd,hour,minute,sec) nutation, obliquity
# julday (year,month,day) 
