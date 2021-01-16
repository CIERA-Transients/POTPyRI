#!/usr/bin/env python
from __future__ import print_function

import sys
from optparse import OptionParser
import util
import quick_reduc
import time
import glob
import matplotlib
import instruments

from pyraf import iraf
import os

matplotlib.use('TkAgg')

description = "> Fast reduction of spectra "
usage = "%prog    \t [option] \n Recommended syntax: %prog -i -c"

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    option, args = parser.parse_args()

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)
    iraf.disp(inlist='1', reference='1')

    iraf.ccdred.verbose = 'no'
    iraf.specred.verbose = 'no'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''

    iraf.longslit.mode = 'h'
    iraf.specred.mode = 'h'
    iraf.noao.mode = 'h'
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"

    if len(args) > 1:
        files=[]
        sys.argv.append('--help')
        option, args = parser.parse_args()
    elif len(args) == 1:
        files = util.readlist(args[0])
    else:
        listfile = glob.glob('*.fits')
        files_science = []
        files_arc = []
        files_dflat = []
        #print 'checking your files ...'
        for img in listfile:
            _type = ''
            hdr0 = util.readhdr(img)
            _type=util.readkey3(hdr0, 'object')
            if _type.lower().startswith("arc"):
                files_arc.append(img)
            elif _type.lower().startswith("dflat"):
                files_dflat.append(img)
            else:
                files_science.append(img)

    list_arc_b = []
    list_arc_r = []
    for arcs in files_arc:
        hdr = util.readhdr(arcs)
        if util.readkey3(hdr, 'VERSION') == 'kastb':
            list_arc_b.append(arcs)
        elif util.readkey3(hdr, 'VERSION') == 'kastr':
            list_arc_r.append(arcs)
        else:
            print(util.readkey3(hdr, 'VERSION') + 'not in database')
            sys.exit()

    list_flat_b = []
    list_flat_r = []
    for dflats in files_dflat:
        hdr = util.readhdr(dflats)
        if util.readkey3(hdr, 'VERSION') == 'kastb':
            list_flat_b.append(dflats)
        elif util.readkey3(hdr, 'VERSION') == 'kastr':
            list_flat_r.append(dflats)
        else:
            print(util.readkey3(hdr, 'VERSION') + 'not in database')
            sys.exit()

    if not os.path.isdir('pre_reduced/'):
        os.mkdir('pre_reduced/')

    mkarc = raw_input("Make arc? ([y]/n)")
    mkflat = raw_input("Make flat? ([y]/n)")

    tlist_files = []
    for img in listfile:
        hdr = util.readhdr(img)
        if util.readkey3(hdr, 'VERSION') == 'kastb':
            inst = instruments.kast_blue
        elif util.readkey3(hdr, 'VERSION') == 'kastr':
            inst = instruments.kast_red
        else:
            print(util.readkey3(hdr, 'VERSION') + 'not in database')
            sys.exit()

        iraf.specred.dispaxi = inst.get('dispaxis')
        iraf.longslit.dispaxi = inst.get('dispaxis')

        _biassec0 = inst.get('biassec')
        _trimsec0 = inst.get('trimsec')

        #OVERSCAN CORRECT AND TRIM ALL IMAGES (OVERSCAN IS CAUSING PROBLEMS)
        # iraf.ccdproc(img, output='t'+img, overscan='yes', trim='yes', zerocor="no", flatcor="no", readaxi='line',
        #              trimsec=str(_trimsec0),biassec=str(_biassec0), Stdout=1)
        if not img.startswith('t'):
            iraf.ccdproc(img, output='t'+img, overscan='no', trim='yes', zerocor="no", flatcor="no", readaxi='line',
                         trimsec=str(_trimsec0), Stdout=1)

    #CREATE RESPONSE FILES, NEED TO IMPLEMENT FLAT COMBINING
    if mkarc != 'n':
        if len(list_arc_b) == 1:
            arc_blue = list_arc_b[0]
            os.system('cp ' + arc_blue + ' ' + 'ARC_blue.fits')
        else:
            arc_str = ''
            for arc in list_arc_b:
                arc_str = arc_str + 't'+ arc + ','
            iraf.imcombine(arc_str, output='ARC_blue')

        if len(list_arc_r) == 1:
            arc_red = list_arc_r[0]
            os.system('cp ' + arc_red + ' ' + 'ARC_red.fits')
        else:
            arc_str = ''
            for arc in list_arc_r:
                arc_str = arc_str + 't'+ arc + ','
            iraf.imcombine(arc_str, output='ARC_red')

        os.system('mv ' + 'ARC_blue.fits' + ' ' + 'pre_reduced' + '/')
        os.system('mv ' + 'ARC_red.fits' + ' ' + 'pre_reduced' + '/')

    if mkflat != 'n':
        inter = 'yes'
        iraf.specred.dispaxi = 1
        if len(list_flat_b) == 1:
            Flat_blue = list_flat_b[0]
        else:
            flat_str = ''
            for flat in list_flat_b:
                flat_str = flat_str + 't'+ flat + ','
            iraf.flatcombine(flat_str, output='tFlat_blue', ccdtype='',rdnoise=3.7, subsets='yes', process='no')
            Flat_blue = 'Flat_blue.fits'

        iraf.specred.response('t'+Flat_blue, normaliz='t'+Flat_blue, response='RESP_blue', interac=inter, thresho='INDEF',
                                                     sample='*', naverage=2, function='legendre', low_rej=3,
                                                     high_rej=3, order=60, niterat=20, grow=0, graphic='stdgraph')

        iraf.specred.dispaxi = 2
        if len(list_flat_r) == 1:
            Flat_red = list_flat_r[0]
        else:
            flat_str = ''
            for flat in list_flat_r:
                flat_str = flat_str + 't'+ flat + ','
            iraf.flatcombine(flat_str, output='tFlat_red', ccdtype='', rdnoise=3.8, subsets='yes', process='no')
            Flat_red = 'Flat_red.fits'

        iraf.specred.response('t'+Flat_red, normaliz='t'+Flat_red, response='RESP_red', interac=inter, thresho='INDEF',
                                                     sample='*', naverage=2, function='legendre', low_rej=3,
                                                     high_rej=3, order=80, niterat=20, grow=0, graphic='stdgraph')

        os.system('mv ' + 'RESP_blue.fits' + ' ' + 'pre_reduced' + '/')
        os.system('mv ' + 'RESP_red.fits' + ' ' + 'pre_reduced' + '/')
    

    #science files should have 't' in front now
    for obj in files_science:
        hdr = util.readhdr(obj)
        if util.readkey3(hdr, 'VERSION') == 'kastb':
            inst = instruments.kast_blue
            flat_file = 'RESP_blue'
        elif util.readkey3(hdr, 'VERSION') == 'kastr':
            inst = instruments.kast_red
            flat_file = 'RESP_red'
        else:
            print(util.readkey3(hdr, 'VERSION') + 'not in database')
            sys.exit()
        iraf.specred.dispaxi = inst.get('dispaxis')

        # iraf.ccdproc('t'+obj, output='ft'+obj, overscan='no', trim='no', zerocor="no", flatcor="yes", readaxi='line', 
        #               flat=flat_file, Stdout=1)
        if not obj.startswith('t'):
            os.system('mv ' + 't'+ obj + ' ' + 'pre_reduced' + '/')
        else:
            os.system('mv ' +  obj + ' ' + 'pre_reduced' + '/')
