from __future__    import print_function
try:      from astropy.io import fits as pyfits
except:   import pyfits
import numpy as np
import glob

#---------------------------------------------------------------------------
# 
# kastbias -- perform overscan bias subtraction from KAST images
#
# Inputs:
#   ifile -- input file to perform overscan bias subtration
#   ofile -- output file
#
# Returns:
#   0, but writes to ofile
#
# Description:
#   Performs KAST bias subtraction.
#   Adapted from Elinor Gates routine. 
#
# Author:
#   J. Brown, UCSC Astronomy
#   brojonat@ucsc.edu
#   2018 Sep 26
#
# Version 2.1 -- Elinor Gates, 2016 Aug 25 - changed math to BITPIX -32
# Version 2.0 -- Elinor Gates, 2016 Aug 11 - improved file reading
# Version 1.0 -- Elinor Gates, 2015 Nov 24
def kastbias(ifile,ofile):
    
    data, header = pyfits.getdata(ifile,header=True)

    # change data to float
    data=data.astype('float32')

    # read necessary keywords from fits header
    xsize = header['NAXIS1']
    ysize = header['NAXIS2']
    xorig = header['CRVAL1U']
    yorig = header['CRVAL2U']
    cdelt1 = header['CDELT1U']
    cdelt2 = header['CDELT2U']
    rover = header['ROVER']
    cover = header['COVER']
    inxsize = header['DNAXIS1']
    inysize = header['DNAXIS2']
    ampsx = header['AMPSCOL']
    ampsy = header['AMPSROW']

    # determine number and sizes of overscan and data regions
    namps = ampsx*ampsy
    if rover > 0:
        over=rover
        sys.exit('Program does not yet deal with row overscans. Exiting.')
    else:
        over = cover
    if over == 0:
        sys.exit('No overscan region specified in FITS header. Exiting.')

    # single amplifier mode
    if namps == 1:
        biassec = data[:,xsize-cover:xsize]
        datasec = data[0:,0:xsize-cover]

        # median overscan section
        bias=np.median(biassec, axis=1) 
        # subtract overscan
        datanew = datasec
        for i in range(datasec.shape[1]):
            datanew[:,i] = datasec[:,i]-bias

    # two amplifier mode
    if namps == 2:
        biasseca = data[:,xsize-cover*2:xsize-cover]
        biassecb = data[:,xsize-cover:xsize]

        # median overscan sections
        biasa=np.median(biasseca,axis=1)
        biasb=np.median(biassecb,axis=1)

        # extract data regions

        #determine size of binned data region
        hsize=abs(inxsize/cdelt1)

        # calculate x origin of readout in binned units if cdelt1 negative or positive
        if cdelt1 < 0:
            xorig=(xorig-(xsize-2*cover)*abs(cdelt1))/abs(cdelt1)
        else:
            xorig=xorig/cdelt1
        x0=xorig+xsize-1-cover*2 # need to test is need -1 because starting counting at 0

        # determine which columns are on which amplifier and subtract proper overscan region

        if x0 < hsize/2: # all data on left amplifier
            datanew=data[:,0:xsize-cover*2]
            m=datanew.shape[1]
            for i in range(0,m):
                datanew[:,i]=datanew[:,i]-biasa

        if xorig >= hsize/2: # all data on right amplifier
            datanew=data[:,0:xsize-cover*2]
            m=datanew.shape[1]
            for i in range(0,m):
                datanew[:,i]=datanew[:,i]-biasb

        if xorig < hsize/2 and x0 > hsize/2:
            x1=hsize/2-xorig
            dataa=data[:,0:x1]
            datab=data[:,x1:-cover*2]
            ma=dataa.shape[1]
            mb=datab.shape[1]
            for i in range(0,ma):
                dataa[:,i]=dataa[:,i]-biasa
            for i in range(0,mb):
                datab[:,i]=datab[:,i]-biasb
            # merge dataa and datab into single image
            datanew=np.hstack([dataa,datab])
        
    if namps > 2: 
        sys.exit('Program does not yet deal with more than two overscan regions. Exiting.')

    # add info to header
    header['HISTORY'] = 'Overscan subtracted'

    # write new fits file
    pyfits.writeto(ofile,datanew,header,clobber=True)

    return 0

def ReadAscii2(ascifile):
    # print "LOGX:: Entering `ReadAscii2` method/function in %(__file__)s" %
    # globals()
    import string

    f = open(ascifile, 'r')
    ss = f.readlines()
    f.close()
    vec1, vec2 = [], []
    for line in ss:
        if line[0] != '#':
            vec1.append(float(string.split(line)[0]))
            vec2.append(float(string.split(line)[1]))
    return vec1, vec2


# ########################################################################
def readlist(listfile):

    import string
    import sys
    import re
    import glob

    if '*' in listfile:
        imglist = glob.glob(listfile)
    elif ',' in listfile:
        imglist = string.split(listfile, sep=',')
    else:
        try:
            hdulist = pyfits.open(listfile)
        except:
            hdulist = []
        if hdulist:
            imglist = [listfile]
        else:
            try:
                ff = open(listfile, 'r')
                files = ff.readlines()
                ff.close()
                imglist = []
                for ff in files:
                    ff = re.sub(' ', '', ff)
                    if not ff == '\n' and ff[0] != '#':
                        ff = re.sub('\n', '', ff)
                        try:
                            hdulist = pyfits.open(ff)
                            imglist.append(ff)
                        except:
                            try:
                                correctcard(ff)
                                hdulist = pyfits.open(ff)
                                imglist.append(ff)
                            except:
                                pass
            except:
                sys.exit('\n##### Error ###\n file ' +
                         str(listfile) + ' do not  exist\n')
    if len(imglist) == 0:
        sys.exit('\n##### Error ###\nIf "' + str(listfile)
                 + '" is an image, it is corrupted \n or is not a list of image\n')
    return imglist


##############################################################################
def delete(listfile):

    import os
    import string
    import re
    import glob

    if listfile[0] == '@':
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files:
            ff = re.sub(' ', '', ff)
            if not ff == '\n' and ff[0] != '#':
                ff = re.sub('\n', '', ff)
                imglist.append(ff)
    elif ',' in listfile:
        imglist = string.split(listfile, sep=',')
    else:
        imglist = [listfile]
    lista = []
    for _file in imglist:
        lista = lista + glob.glob(_file)
    if lista:
        for _file in lista:
            try:
                os.system('rm ' + _file)
            except:
                pass


###############################################################
def readhdr(img):
    try:
        hdr = pyfits.open(img)[0].header
    except:
        try:
            correctcard(img)
        except:
            import sys

            sys.exit('image ' + str(img) +
                     ' is corrupted, delete it and start again')
        hdr = pyfits.open(img)[0].header
    return hdr

###############################################################
def readhdr2(img):
    try:
        # hdr = pyfits.open(img)[0].header
        hdr = pyfits.open(img)[1].header
    except:
        try:
            correctcard(img)
        except:
            import sys

            sys.exit('image ' + str(img) +
                     ' is corrupted, delete it and start again')
        # hdr = pyfits.open(img)[0].header
        hdr = pyfits.open(img)[0].header
    return hdr


def readkey3(hdr, keyword):

    value = hdr.get(keyword)
    return value



################################################
def correctcard(img):
    import numpy as np
    # print "LOGX:: Entering `correctcard` method/function in %(__file__)s" %
    # globals()
    from numpy import asarray
    import re
    import os
    hdulist = pyfits.open(img)
    a = hdulist[0]._verify('exception')
    _header = hdulist[0].header
    hdulist.close()

    ######   change 20161003 
    #print a
    #for i in range(len(a)):
    #    if not a[i]:
    #        a[i] = ['']
    #ww = asarray([i for i in range(len(a)) if (re.sub(' ', '', a[i][0]) != '')])

    ww = np.where(a)[0]
    try:
        if int(re.sub('\.', '', str(pyfits.__version__))[:2]) <= 30:
            aa = 'HIERARCH '
        else:
            aa = ''
    except:
        aa = ''

    # print(aa)
    # print(ww)
    if len(ww) > 0:
        newheader = []
        headername = []
        for j in list(_header.items()):
            headername.append(j[0])
            newheader.append(j[1])
        for i in ww:
            if len(headername[i]) > 8:
                data, hdr = pyfits.getdata(img, 0, header=True)
                comm = hdr.comments[aa + headername[i]]
                hdr.pop(aa + headername[i])
                cc = re.sub('HIERARCH ', '', headername[i])
                hdr.append(
                    ('HIERARCH ' + cc, newheader[i], comm), False, False)
                os.system('rm ' + img)
                pyfits.writeto(img, data, hdr)
            else:
                try:
                    imm = pyfits.open(img, mode='update')
                    imm.close(output_verify='silentfix', verbose=False)
                except:
                    imm = pyfits.open(img, mode='update')
                    _header = imm[0].header
                    for i in ww:
                        if headername[i]:
                            if len(headername[i]) > 8 and 'HIERARCH' not in headername[i]:
                                bb = 'HIERARCH '
                            else:
                                bb = ''
                            try:
                                _header.update(
                                    bb + headername[i], newheader[i])
                            except:
                                _header.update(bb + headername[i], 'xxxx')
                    imm.flush()
                    imm.close()
    else:
        pass


#    data, hdr = pyfits.getdata(img, 0, header=True)
#       print 'correction not needed'
##########################################################################

def updateheader(image, dimension, headerdict):
    # added to cut long header
    while len(max([str(headerdict[i][0]) for i in headerdict], key=len)) > 68:
        key = [i for i in headerdict]
        valori, commenti = list(zip(*[headerdict[i] for i in headerdict]))
        num = valori.index(max([str(headerdict[i][0])
                                for i in headerdict], key=len))
        headerdict[key[num]] = [valori[num][0:68], commenti[num]]
        print('warning: header to long, ', str(key[num]), str(valori[num][0:68]), str(commenti[num]))
        #   keytochange=hdr.keys()[hdr.values().index(max([str(i) for i in hdr.values()],key=len))]
        #   hdr[keytochange]=[str(hdr[keytochange])[0:68]]

    if 'version' in dir(pyfits):
        try:
            imm = pyfits.open(image, mode='update')
            _header = imm[dimension].header
            for i in list(headerdict.keys()):
                _header.update(i, headerdict[i][0], headerdict[i][1])
            imm.flush()
            imm.close()
        except:
            print('warning: problem to update header, try to correct header format ....')
            correctcard(image)
            try:
                print(headerdict)
                imm = pyfits.open(image, mode='update')
                _header = imm[dimension].header
                for i in list(headerdict.keys()):
                    _header.update(i, headerdict[i][0], headerdict[i][1])
                imm.flush()
                imm.close()
            except:
                print('error: not possible update header')
    else:
        #
        #
        # astropy.io.fits requre a tuple to update header
        #
        #
        imm = pyfits.open(image, mode='update')
        _header = imm[dimension].header
        for i in list(headerdict.keys()):
            _header.update( { i : (headerdict[i][0], headerdict[i][1]) } )
        imm.flush()
        imm.close()
        
##########################################################################


def display_image(img, frame, _z1, _z2, scale, _xcen=0.5, _ycen=0.5, _xsize=1, _ysize=1, _erase='yes'):
    # print "LOGX:: Entering `display_image` method/function in %(__file__)s"
    # % globals()
    goon = 'True'
    import glob
    import subprocess
    import time
    import os

    ds9 = subprocess.Popen("ps -U" + str(os.getuid()) + "|grep -v grep | grep ds9", shell=True,
                           stdout=subprocess.PIPE).stdout.readlines()
    if len(ds9) == 0:
        subproc = subprocess.Popen('ds9', shell=True)
        time.sleep(3)

    if glob.glob(img):
        from pyraf import iraf

        iraf.images(_doprint=0)
        iraf.tv(_doprint=0)
        import string
        import os

        if _z2:
            try:
                sss = iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,
                                   fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2, Stdout=1)
            except:
                print('')
                print('### ERROR: PROBLEM OPENING DS9')
                print('')
                goon = 'False'
        else:
            try:
                sss = iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,
                                   fill='yes', Stdout=1)
            except:
                print('')
                print('### ERROR: PROBLEM OPENING DS9')
                print('')
                goon = False

        if scale and goon:
            answ0 = input('>>> Cuts OK ? [y/n] ? [y] ')
            if not answ0:
                answ0 = 'y'
            elif answ0 == 'no' or answ0 == 'NO':
                answ0 = 'n'

            while answ0 == 'n':
                _z11 = float(string.split(string.split(sss[0])[0], '=')[1])
                _z22 = float(string.split(string.split(sss[0])[1], '=')[1])
                z11 = input('>>> z1 = ? [' + str(_z11) + '] ? ')
                z22 = input('>>> z2 = ? [' + str(_z22) + '] ? ')
                if not z11:
                    z11 = _z11
                else:
                    z11 = float(z11)
                if not z22:
                    z22 = _z22
                else:
                    z22 = float(z22)
                print(z11, z22)
                sss = iraf.display(img, frame, fill='yes', xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize,
                                   erase=_erase, zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
                answ0 = input('>>> Cuts OK ? [y/n] ? [y] ')
                if not answ0:
                    answ0 = 'y'
                elif answ0 == 'no' or answ0 == 'NO':
                    answ0 = 'n'
        if goon:
            _z1, _z2 = string.split(string.split(sss[0])[0], '=')[
                1], string.split(string.split(sss[0])[1], '=')[1]
    else:
        print('Warning: image ' + str(img) + ' not found in the directory ')
    return _z1, _z2, goon


###########################################################################

def readspectrum(img):
    # print "LOGX:: Entering `readspectrum` method/function in %(__file__)s" %
    # globals()
    from numpy import array
    import string

    fl = ''
    lam = ''
    graf = 1
    spec = pyfits.open(img)
    head = spec[0].header
    try:
        if spec[0].data.ndim == 1:
            fl = spec[0].data
        elif spec[0].data.ndim == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.ndim == 3:
            fl = spec[0].data[0, 0, :]
    except:
        if spec[0].data.rank == 1:
            fl = spec[0].data
        elif spec[0].data.rank == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.rank == 3:
            fl = spec[0].data[0, 0, :]
    naxis1 = head['naxis1']
    try:
        crpix1 = head['crpix1']
        crval1 = head['crval1']
        try:
            cdelt1 = head['cdelt1']
        except:
            cdelt1 = head['cd1_1']
        pix = array(list(range(1, naxis1 + 1, 1)))
        pix = array(list(range(1, len(fl) + 1, 1)))
        lam = (pix - crpix1) * cdelt1 + crval1
    except:
        try:
            WAT = head['WAT2_001']
            pix = array(list(range(1, naxis1 + 1, 1)))
            crpix1 = string.split(string.split(WAT, '"')[1])[0]
            crval1 = string.split(string.split(WAT, '"')[1])[3]
            cdelt1 = string.split(string.split(WAT, '"')[1])[4]
            lam = (pix - float(crpix1)) * float(cdelt1) + float(crval1)
        except:
            graf = 0
    return lam, fl


###########################################################################
def pval(_xx, p):
    # print "LOGX:: Entering `pval` method/function in %(__file__)s" %
    # globals()
    _y = +p[0] + p[1] * _xx
    return _y


def residual(p, y, x):
    # print "LOGX:: Entering `residual` method/function in %(__file__)s" %
    # globals()
    for i in range(len(p)):
        err = (y - p[i] * x ** i)
    return err



##########################################################################
#Edit this for extraction parameters
def dvex():
    # print "LOGX:: Entering `dvex` method/function in %(__file__)s" %
    # globals()
    dv = {}
    #dv['line'] = {'Gr16': 300, 'Gr11': 430, 'Gr13': 200, 'GR': 150, 'GB': 430}
    dv['std'] = {'_t_order': 4, '_t_niter': 20, '_t_sample': '*', '_t_nlost': 20, 
                 '_width': 10, '_radius': 10,
                 '_weights': 'variance',
                 '_nsum': 30, 
                 '_t_step': 10, '_t_nsum': 10, 
                 '_lower': -10, '_upper': 10, 
                 '_b_sample': '-40:-20,20:40',
                 '_b_order': 2, '_resize': 'no'}
    dv['obj'] = {'_t_order': 4, '_t_niter': 20, '_t_sample': '*', '_t_nlost': 20, 
                 '_width': 10, '_radius': 10,
                 '_weights': 'variance',
                 '_nsum': 40, '_t_step': 10, '_t_nsum': 10, 
                 '_lower': -5, '_upper': 5, 
                 '_b_sample': '-25:-15,15:25',
                 '_b_order': 2,
                 '_resize': 'yes'}
    return dv



###############################################################################


#################################################################
def repstringinfile(filein, fileout, string1, string2):
    # print "LOGX:: Entering `repstringinfile` method/function in
    # %(__file__)s" % globals()
    import re

    f = open(filein, 'r')
    ss = f.readlines()
    f.close()
    f = open(fileout, 'w')
    for n in range(len(ss)):
        if string1 in ss[n]:
            f.write(re.sub(string1, string2, ss[n]))
        else:
            f.write(ss[n])
    f.close()


###################################################

def StoN2(img, show=False):
    # print "LOGX:: Entering `StoN2` method/function in %(__file__)s" %
    # globals()
    import numpy as np
    data, hdr0 = pyfits.getdata(img, header=True)
    yy1 = data[0][0]
    #      yy3=data[2][0]
    #      yy2=data[1][0]
    yy4 = data[3][0]
    xx = np.arange(len(data[0][0]))
    xxmed = hdr0['CRVAL1'] + xx * hdr0['CD1_1']
    sntot = yy1 / yy4
    if show:
        import pylab as pl

        pl.ion()
        pl.plot(xxmed, yy1, '-r', label='Spectrum')
        pl.plot(xxmed, sntot, '-b', label='StoN ' + str(np.mean(sntot)) + ' ')
        pl.legend(numpoints=1, markerscale=1.5, loc=1, ncol=1)
        pl.xlabel('Wavlength')
        pl.ylabel('StoN')
    return np.median(sntot)

################################################


def extractspectrum(img, dv, inst, _interactive, _type, automaticex=False, host_ex = False, match_aperture = 'n'):
    # print "LOGX:: Entering `extractspectrum` method/function in
    # %(__file__)s" % globals()
    import glob
    import os
    import string
    import sys
    import re
    import datetime
    import numpy as np

    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 0o1, 0o1)).days
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.apall', 'specred.transform']
    for t in toforget:
        iraf.unlearn(t)

    hdr = readhdr(img)
    iraf.specred.dispaxi = inst.get('dispaxis')

    imgex = re.sub('.fits', '_ex.fits', img)
    print (img)
    imgfast = re.sub(string.split(img, '_')[-2] + '_', '', img)
    if not os.path.isfile(imgex) and not os.path.isfile(
            'database/ap' + re.sub('.fits', '', img)) and not os.path.isfile(
            'database/ap' + re.sub('.fits', '', imgfast)):
        _new = 'yes'
        _extract = 'yes'
    else:
        if automaticex:
            if _interactive in ['Yes', 'yes', 'YES', 'y', 'Y']:
                answ = 'x'
                while answ not in ['o', 'n', 's']:
                    answ = input(
                        '\n### New extraction [n], extraction with old parameters [o], skip extraction [s] ? [o]')
                    if not answ:
                        answ = 'o'
                if answ == 'o':
                    _new, _extract = 'no', 'yes'
                elif answ == 'n':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'no', 'yes'
        else:
            if _interactive in ['Yes', 'yes', 'YES', 'y', 'Y']:
                answ = 'x'
                while answ not in ['y', 'n']:
                    answ = input(
                        '\n### do you want to extract again [[y]/n] ? ')
                    if not answ:
                        answ = 'y'
                if answ == 'y':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'yes', 'yes'
    if _extract == 'yes':
        delete(imgex)
        dist = 200
        _reference = ''
        _fittrac = 'yes'
        _trace = 'yes'
        if _new == 'no':
            if not os.path.isfile('database/ap' + re.sub('.fits', '', img)):
                repstringinfile('database/ap' + re.sub('.fits', '', imgfast),
                                         'database/ap' +
                                         re.sub('.fits', '', img), re.sub(
                                             '.fits', '', imgfast),
                                         re.sub('.fits', '', img))
            _find = 'no'
            _recenter = 'no'
            _edit = 'no'
            _trace = 'no'
            _fittrac = 'no'
            _mode = 'h'
            _resize = 'no'
            _review = 'no'
            iraf.specred.mode = 'h'
            _interactive = 'no'
        else:
            iraf.specred.mode = 'q'
            _mode = 'q'
            _find = 'yes'
            _recenter = 'yes'
            _edit = 'yes'
            _review = 'yes'
            _resize = dv[_type]['_resize']

        # _interactive = False
        _recenter = 'no'
        _resize = 'no'
        _edit = 'yes'
        _trace = 'yes'
        _fittrace = 'no'
        _review = 'yes'

        _lower=dv[_type]['_lower']
        _upper=dv[_type]['_upper']
        _b_sample = dv[_type]['_b_sample']

        if match_aperture == 'y':

            aps = glob.glob('database/ap*')
            for ap in aps:
                print (ap.split('/')[-1])
            ap_select = raw_input('Choose image to match apertures: ')
            
            #this is just to make things automatic for now
            if 'blue' in ap_select and 'lris' in inst.get('name'):
                ap_match = ap_select.replace('blue','red')
                ap_binning = 1.
            elif 'red' in ap_select and 'lris' in inst.get('name'):
                ap_match = ap_select.replace('red','blue')
                ap_binning = 2.
            elif 'blue' in ap_select and 'kast' in inst.get('name'):
                ap_match = ap_select.replace('blue','red')
                ap_binning = 1.
            elif 'red' in ap_select and 'kast' in inst.get('name'):
                ap_match = ap_select.replace('red','blue')
                ap_binning = 1.


            # img_binning = hdr.get('BINNING', None).strip()
            # if img_binning != None:
            #     img_binning = float(img_binning.split(',')[0])
            img_binning = inst.get('spatial_binning')

            # ap_binning = raw_input('Enter aperture spatial binning [1]: ') or 1
            # ap_binning = float(ap_binning)

            delete('database/'+ap_match)
            new_apfile = 'database/'+ap_match
            ismainap = False

            aps = glob.glob('../master_files/ap*')
            for ap in aps:
                print (ap.split('/')[-1])
            center_ap = raw_input('Choose master aperture for center: ')

            #get center of main reference aperture (should be a std)
            with open('../master_files/'+center_ap) as c_ap:
                c_ap_data = c_ap.readlines()
                for i, c_line in enumerate(c_ap_data):
                    if 'center' in c_line:
                        center_ap_loc = float(c_line.split()[2])

                with open('database/'+ap_select) as old_file:
                    ap_data = old_file.readlines()
                    lows, highs, b_samps = get_relevant_ap_data(ap_data, ap_binning, img_binning, inst)
                    with open(new_apfile,'w') as new_file:
                        for i in range(len(lows)):
                            low = str(lows[i])
                            high = str(highs[i])
                            for c_line in c_ap_data:
                                if 'begin' in c_line and 'aperture' in c_line:
                                    new_file.write(c_line.replace((c_line.split()[2] + ' 1'), c_line.split()[2] + ' ' + str(i+1)))
                                elif 'begin' not in c_line and 'aperture' in c_line:
                                    new_file.write(c_line.replace(c_line.split()[1], str(i+1)))
                                elif 'low' in c_line and 'reject' not in c_line:
                                    new_file.write(c_line.replace(c_line.split()[2], low))
                                elif 'high' in c_line and 'reject' not in c_line:
                                    new_file.write(c_line.replace(c_line.split()[2], high))
                                elif 'sample' in c_line:
                                    b_sample = str(b_samps[i])
                                    new_file.write(c_line.replace(c_line.split('sample')[1].split('\n')[0], ' '+ b_sample))
                                else:
                                    new_file.write(c_line)



            _find = 'yes'
            _recenter = 'yes'
            _reference=ap_match[2:]
            # _reference = 'BD262606_lris_red_1'


        # use_diff_aperture = raw_input('Trace different aperture? y/[n]: ')
        # if use_diff_aperture == 'y':
        #     aps = glob.glob('../master_files/ap*')
        #     for ap in aps:
        #         print (ap.split('/')[-1])
        #     ap_select = raw_input('Choose ref image from list above: ')
        #     os.system('cp ' + ' ../master_files/' + ap_select + ' database/'+ap_select)
        #     _reference = ap_select[2:]
        #     _fittrac = 'no'
        #     _trace = 'no'


        if host_ex:
            _recenter = 'n'
            _trace = 'n'
            _fittrac = 'n'
            _find = 'n'
        
        iraf.specred.apall(img, output=imgex, referen=_reference, trace=_trace, fittrac=_fittrac, find=_find,
                           recenter=_recenter, edit=_edit,
                           nfind=1, backgro='fit', lsigma=3, usigma=3, b_order=dv[_type]['_b_order'],
                           format='multispec', extras='yes',
                           b_function='chebyshev', b_sample=_b_sample, clean='yes', pfit='fit1d',
                           lower=_lower, upper=_upper, t_niter=dv[_type]['_t_niter'],
                           width=dv[_type]['_width'],
                           radius=dv[_type]['_radius'], 
                           line=inst.get('approx_extract_column','INDEF'), nsum=dv[_type]['_nsum'], 
                           t_step=dv[_type]['_t_step'],
                           t_nsum=dv[_type]['_t_nsum'],
                           t_nlost=dv[_type]['_t_nlost'], t_sample=dv[
                               _type]['_t_sample'], resize=_resize,
                           t_order=dv[_type]['_t_order'],
                           weights=dv[_type]['_weights'], 
                           interactive=_interactive, review=_review, mode=_mode)
    else:
        print('\n### skipping new extraction')
    return imgex

def get_relevant_ap_data(ap_data, ap_binning, img_binning, inst):

    num_aps = 0
    lows = []
    highs = []
    b_samps = []

    ap_centers = {}
    for i, line in enumerate(ap_data):
        if len(line.split()) == 2 and 'aperture' in line:
            current_ap = line.split()[1]
            num_aps+=1
        elif 'center' in line:
            ap_centers[current_ap] = float(line.split()[2])

    if 'kast' in inst.get('name'):
        sign = 1.
    else:
        sign = 1.

    for line in ap_data:

        if len(line.split()) == 2 and 'aperture' in line:
            if line.split()[1] == '1':
                ismainap = True
            else:
                ismainap = False
            current_ap = line.split()[1]

        elif 'center' in line:
            if ismainap:
                diff=0
            else:
                diff = (float(line.split()[2]) - ap_centers['1'])*(ap_binning/img_binning)
            if 'lris' in inst.get('name'):
                diff *= -1.

        elif 'low' in line and 'reject' not in line:
            if 'kast_blue' in inst.get('name'): #current inst is blue so want to get correct red ap
                low = str(sign*float(line.split()[1])*(ap_binning/img_binning) + diff) 
            else:
                low = str(sign*float(line.split()[2])*(ap_binning/img_binning) + diff) 
            lows.append(low)

        elif 'high' in line and 'reject' not in line:
            if 'kast_blue' in inst.get('name'): #current inst is blue so want to get correct red ap
                high = str(sign*float(line.split()[1])*(ap_binning/img_binning) + diff)
            else:
                high = str(sign*float(line.split()[2])*(ap_binning/img_binning) + diff)
            highs.append(high)

        elif 'sample' in line:
            _b_sample = line.split('sample')[1].split('\n')[0]
            if ',' in _b_sample:
                b1 = float(_b_sample.split(',')[0].split(':')[0])
                b2 = float(_b_sample.split(',')[0].split(':')[1])
                b3 = float(_b_sample.split(',')[1].split(':')[0])
                b4 = float(_b_sample.split(',')[1].split(':')[1])
            else:
                b1 = float(_b_sample.split()[0].split(':')[0])
                b2 = float(_b_sample.split()[0].split(':')[1])
                b3 = float(_b_sample.split()[1].split(':')[0])
                b4 = float(_b_sample.split()[1].split(':')[1])
            b1 = b1*(ap_binning/img_binning) + diff
            b2 = b2*(ap_binning/img_binning) + diff
            b3 = b3*(ap_binning/img_binning) + diff
            b4 = b4*(ap_binning/img_binning) + diff
            if 'kast' in inst.get('name'):
                _b_sample_new = str(sign*b4)+':'+str(sign*b3)+' '+str(sign*b2)+':'+str(sign*b1)
            else:
                _b_sample_new = str(b1)+':'+str(b2)+' '+str(b3)+':'+str(b4)
            b_samps.append(_b_sample_new)

    return lows, highs, b_samps



def skyfrom2d(fitsfile, skyfile, interac=True):

    from astropy.io import fits as pyfits
    import numpy as np

    hdulist = pyfits.open(fitsfile)
    scidata = hdulist[0].data[:,0]
    scidata_fits=scidata[3,:]

    #yy1 = pyfits.open(fitsfile)[0].data[:, :].mean(1)

    yy1 = scidata_fits

    #print yy1
    crval2 = pyfits.open(fitsfile)[0].header.get('CRVAL1')
    cd2 = pyfits.open(fitsfile)[0].header.get('CD1_1')

    #hdulist = pyfits.open(fitsfile)
    #scidata = hdulist[0].data[:,0]
    #crval2=hdulist[0].header['CRVAL1']
    #cd2=hdulist[0].header['CD1_1']
    #scidata_fits=scidata[3,:]

    delete('new3.fits')
    hdu = pyfits.PrimaryHDU(yy1)
    hdulist = pyfits.HDUList([hdu])
#    hdulist[0].header.update('CRVAL1', crval2)
#    hdulist[0].header.update('CD1_1', cd2)
    hdulist[0].header['CRVAL1']= crval2
    hdulist[0].header['CD1_1'] = cd2
    hdulist.writeto('new3.fits')
    hdulist.close()

    fitsfile = continumsub('new3.fits', 6, 1)
    yy1 = pyfits.open(fitsfile)[0].data
    xx1 = np.arange(len(yy1))
    aa1 = crval2 + (xx1) * cd2
    delete('new3.fits')

    skyff = pyfits.open(skyfile)[0].data
    #print 'here',skyff[100],yy1[100]
    crval1 = pyfits.open(skyfile)[0].header.get('CRVAL1')
    cd1 = pyfits.open(skyfile)[0].header.get('CDELT1')
    skyxx = np.arange(len(skyff))
    skyaa = crval1 + (skyxx) * cd1
    shift = checkwavelength_arc(
        aa1, yy1, skyaa, skyff, 5500, 6500, interac)
    return shift

def continumsub(imagefile, _order1, _order2):
    # print "LOGX:: Entering `continumsub` method/function in %(__file__)s" %
    # globals()
    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.continuum']
    for t in toforget:
        iraf.unlearn(t)
    delete('tsky.fits')
    iraf.specred.continuum(imagefile, output='tsky.fits', type='difference',
                           interact='no', function='legendre', niterat=300, low_rej=3, high_re=2, sample='*',
                           order=_order1, ask='YES')
    delete(imagefile)
    iraf.continuum('tsky.fits', output=imagefile, type='difference',
                   interact='no', function='spline1', overrid='yes', niterat=10, low_rej=3, high_re=1, sample='*',
                   order=_order2, ask='YES')
    delete('tsky.fits')
    return imagefile

def checkwavelength_arc(xx1, yy1, xx2, yy2, xmin, xmax, inter=True):
    # print "LOGX:: Entering `checkwavelength_arc` method/function in
    # %(__file__)s" % globals()
    import numpy as np

    minimo = max(min(xx1), min(xx2)) + 50
    massimo = min(max(xx1), max(xx2)) - 50
    yy1 = [0 if e < 0 else e for e in np.array(yy1)]
    yy2 = [0 if e < 0 else e for e in np.array(yy2)]
    _shift, integral = [], []
    for shift in range(-500, 500, 1):
        xxnew = xx1 + shift / 10.
        yy2interp = np.interp(xxnew, xx2, yy2)
        yy2timesyy = yy2interp * yy1
        xxcut = np.compress((np.array(xxnew) >= minimo) & (
            np.array(xxnew) <= massimo), np.array(xxnew))
        yycut = np.compress((np.array(xxnew) >= minimo) & (
            np.array(xxnew) <= massimo), np.array(yy2timesyy))
        integrale = np.trapz(yycut, xxcut)
        integral.append(integrale)
        _shift.append(shift / 10.)
    result = _shift[integral.index(max(integral))]
    if inter:
        import matplotlib.pyplot as mpl
        #mpl.use("TKAgg")
        #import pylab as pl
        #pl.ion()
        #pl.clf()
        ratio = np.trapz(yy1, xx1) / np.trapz(yy2, xx2)
        yy3 = np.array(yy2) * float(ratio)
        xx4 = xx1 + result
        fig = mpl.figure()
        mpl.plot(xx1, yy1, label='spectrum')
        mpl.plot(xx2, yy3, label='reference sky')
        mpl.plot(xx4, yy1, label='shifted spectrum')
        mpl.legend(numpoints=1, markerscale=1.5)
        mpl.xlim(xmin, xmax)
        mpl.ylim(0,100)
        #if xmin != '' and xmax != '':
        #    mpl.xlim(xmin, xmax)
        #mpl.show()
        fig.savefig('_shift.png')
        mpl.close(fig)
    return result
##########################################################################
def name_duplicate(img, nome, ext):
    # print "LOGX:: Entering `name_duplicate` method/function in %(__file__)s"
    # % globals()
    import glob

    dimg = readkey3(readhdr(img), 'DATE-OBS')
    listafile = glob.glob(nome + '_?' + ext + '.fits') + \
        glob.glob(nome + '_??' + ext + '.fits')
    if len(listafile) == 0:
        nome = nome + "_1" + ext + '.fits'
    else:
        date = []
        for l in listafile:
            date.append(readkey3(readhdr(l), 'DATE-OBS'))
        if dimg in date:
            nome = listafile[date.index(dimg)]
        #         if overwrite:
        #            delete(nome)
        else:
            n = 1
            while nome + '_' + str(n) + str(ext) + '.fits' in listafile:
                n = n + 1
            nome = nome + '_' + str(n) + str(ext) + '.fits'
    return nome

