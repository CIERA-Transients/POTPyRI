from __future__ import print_function
import os,sys,pdb,argparse,glob,shutil,string,time,datetime
import numpy as np 
from astropy.io import fits
from astropy.convolution import convolve
from scipy import signal

def sigclipmedian(datain,sigmahi=3.5,sigmalo=3.5):
    data = 1.*datain
    maxiter = 6
    endclipfrac = 0.01
    ndata = len(data)
    if ndata == 1:
        return data[0]
    if ndata <= 3:
        return np.median(data)

    goodMask = np.isfinite(data)
    goodInds = np.where(goodMask)[0]
    ctgood = len(goodInds)
    if ctgood == 0:
        return np.nan
    if ctgood == 1: 
        return data[goodInds[0]]

    data = data[goodInds]
    clipdata = np.sort(data) # apparently this is for speed, i dunno

    nremoved = ndata
    niter = 0
    while niter < maxiter and nremoved > 0:
        runningmedian = np.median(clipdata)
        runningstdev = np.std(clipdata,ddof=1) # this is how the IDL routine calcs std
        lolimit = runningmedian - sigmalo * runningstdev
        hilimit = runningmedian + sigmahi * runningstdev
        runningn = len(clipdata)

        goodInds = np.where((clipdata >= lolimit) & (clipdata <= hilimit))[0]
        ctgood = len(goodInds)

        if ctgood > 0:
            clipdata = clipdata[goodInds]
        else:
            break


        nremoved = runningn - ctgood
        niter = niter + 1

    return np.median(clipdata)



def writefits(image,outfile,header='HEADER',CLOBBER=False):

    # create hdu for data
    if header == 'HEADER':
        hdu = fits.PrimaryHDU(image)
    else:
        hdu = fits.PrimaryHDU(image,header)

    # check the output location, clear if necessary
    if os.path.isfile(outfile):
        if CLOBBER:
            os.remove(outfile)
        else:
            print('File exists, set CLOBBER=True or move the file')
            return 1
    #Create HDU list and write in one go  
    hdu.writeto(outfile,output_verify='ignore')  
    return 0



def compare_images(img1='',img2='',imgOut='',method='divide'):

    if len(img1) == 0:
        img1 = 'input.fits'
    if len(img2) == 0:
        img2 = 'final.fits'
    if len(imgOut) == 0:
        imgOut = 'comp.fits'

    img1Data = fits.open(img1)[0].data
    img2Data = fits.open(img2)[0].data

    if method == 'divide':
        relative = img1Data / img2Data
    elif method == 'subtract':
        relative = img1Data - img2Data
    else:
        print('I don\'t recognize that method of comparison!')
        pdb.set_trace()
        return 1

    res = writefits(relative,imgOut,CLOBBER=True)

    return 0


def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'A python verison of pzapspec -- CR removal for LRIS'
    parser = argparse.ArgumentParser(description=descStr)

    # required args
    parser.add_argument('infile',type=str,
                       help='the flat field to inspect/mask')

    # optional
    parser.add_argument('-v','--verbose',
                        help='print diagnostic info',action='store_true')
    parser.add_argument('-c','--clobber',action='store_true',
                        help='Clobber files already in pre_reduced/ but not subdirs')
    parser.add_argument('-d','--debug',action='store_true',
                        help='Debug pyzapspec')

    parser.add_argument('--outfile',
                        help='Destination of cleaned file')
    parser.add_argument('--maskfile',
                        help='Destination of the CR mask file')

    parser.add_argument('--boxsize',type=int,
                        help='Kernel size (should approximate seeing in pixels, must be odd)')
    parser.add_argument('--nsigma',type=float,
                        help='Rejection threshold (recommended >10)')
    parser.add_argument('--subsigma',type=float,
                        help='Rejection threshold for neighboring pixels during percolation')

    # parse
    cmdArgs = parser.parse_args()

    # logic mapping to my args/kwargs
    args = (cmdArgs.infile,)
    kwargs = {}
    kwargs['VERBOSE'] = cmdArgs.verbose
    kwargs['CLOBBER'] = cmdArgs.clobber
    kwargs['DEBUG'] = cmdArgs.debug
    kwargs['outfile'] = cmdArgs.outfile 
    kwargs['maskfile'] = cmdArgs.maskfile
    kwargs['boxsize'] = cmdArgs.boxsize
    kwargs['nsigma'] = cmdArgs.nsigma
    kwargs['subsigma'] = cmdArgs.subsigma

    return (args,kwargs)

# def pyzapspec(infile, 
#               outfile='',
#               maskfile='', 
#               WRITE_OUTFILE=False,
#               DEBUG_DIR='../test_data/',DEBUG=False,
#               boxsize=7,nsigma=16.,subsigma=2.7,sfactor=1.0,
#               nzap=0,mask=0,writemodel=0,verbose=0,skysubtract=0,
#               zero=0,method=0,usamp=0,ybin=0,nan=-999,inmaskname=0,**kwargs):

def pyzapspec(infile, 
              outfile='',
              maskfile='', 
              WRITE_OUTFILE=False,
              DEBUG_DIR='../test_data/',DEBUG=False,
              boxsize=9,nsigma=15.,subsigma=2.8,sfactor=1.0,
              nzap=0,mask=0,writemodel=0,verbose=0,skysubtract=0,
              zero=0,method=0,usamp=0,ybin=0,nan=-999,inmaskname=0,**kwargs):
    # defaults
    if len(outfile) == 0:
        dirpath,infile_base = os.path.split(infile)
        outfile_base = 'cr{}'.format(infile_base)
        outfile = os.path.join(dirpath,outfile_base)

    if len(maskfile) == 0:
        dirpath,infile_base = os.path.split(infile)
        maskfile_base = 'cr{}_mask.fits'.format(infile_base.split('.')[0])
        maskfile = os.path.join(dirpath,maskfile_base)

    tstart = time.time()

    # read in a copy of the input image to construct the output instance
    outimgHDU = fits.open(infile)
    outimg = outimgHDU[0].data
    header = outimgHDU[0].header

    dims = outimg.shape
    nx = dims[1] #spectral
    ny = dims[0] #spatial
    outmask = np.full((ny,nx),0)
    ymedimage = np.zeros((ny,nx))
    xmedimage = np.zeros((ny,nx))
    zapimage = np.zeros((ny,nx))

    nzapCount = 0
    iterCount = 0
    nbadCount = 1

    # first do a crude median subtraction of the sky lines
    # for x in range(nx):
    for x in xrange(nx):
        ymedimage[:,x] = np.median(outimg[:,x])
    ysubimage = outimg - ymedimage

    # now subtract traces so we can do a better sky-line subtraction
    # otherwise traces will be killed/wings oversubtracted
    sourceblocksize = 100
    nxblock = int(np.ceil(1.*nx/sourceblocksize))
    realsourcefiltsize = 1*sourceblocksize
    x0 = (np.arange(0,nxblock)*realsourcefiltsize).astype(int)
    x1 = np.append(x0[1:]-1,nx-1)
    xs0 = np.insert(x0[0:nxblock-1],x0[0],0)
    xs1 = np.append(x1[1:nxblock],x1[nxblock-1])

    # for b in range(nxblock):
    for b in xrange(nxblock):
        # for y in range(ny):
        for y in xrange(ny):
            xmedimage[y, x0[b]:x1[b]+1] = np.median(ysubimage[y, xs0[b]:xs1[b]+1])

    kernel = 1.*np.arange(0,realsourcefiltsize/2)
    kernelRev = kernel[::-1]
    kernel = np.append(kernel,kernelRev)
    kernel = kernel / np.sum(kernel)

    # pzap is convolving the 2D xmedimage with a 1D kernel, 
    # I'm pretty sure this is just a row by row convolution.
    # for xrow in range(ny):
    for xrow in xrange(ny):
        newRow = np.convolve(xmedimage[xrow,:],kernel,mode='same')
        xmedimage[xrow,:] = 1.*newRow
    xsubimage = outimg - xmedimage

    # here I'm assuming the method=0 (pzapspec line 226),
    # since I don't feel like coding up the alternative right now
    skyblocksize = 40
    nyblock = int(np.ceil(1.*ny/skyblocksize))
    realskyblocksize = 1.*skyblocksize
    y0 = (np.arange(0,nyblock)*realskyblocksize).astype(int)
    y1 = np.append(y0[1:]-1,ny-1)
    ys0 = np.insert(y0[0:nyblock-1],y0[0],0)
    ys1 = np.append(y1[1:nyblock],y1[nyblock-1])

    # for b in range(nyblock):
    #     for x in range(nx):
    for b in xrange(nyblock):
        for x in xrange(nx):
            scm = sigclipmedian(xsubimage[ys0[b]:ys1[b]+1,x])
            ymedimage[y0[b]:y1[b]+1,x] = scm

    kernel = 1.*np.arange(0,realskyblocksize/2)
    kernelRev = kernel[::-1]
    kernel = np.append(kernel,kernelRev)
    kernel = kernel / np.sum(kernel)

    # pzap is convolving the 2D xmedimage with a 1D kernel, 
    # I'm pretty sure this is just a column by column convolution.
    # for ycol in range(nx):
    for ycol in xrange(nx):
        newCol = np.convolve(ymedimage[:,ycol],kernel,mode='same')
        ymedimage[:,ycol] = 1.*newCol
    ysubimage = outimg - ymedimage


    skysubimage = outimg - ymedimage - xmedimage  # actually subtracts sky AND sources.

    filterimage = signal.medfilt2d(skysubimage,[boxsize,boxsize])
    residualimage = skysubimage - filterimage

    sigmaimage = np.zeros((ny,nx)) + np.nan
    nyb = round(ny / 200)
    yint = np.ceil(ny*1./nyb)

    # for yb in range(nyb):
    #     for x in range(nx):
    for yb in xrange(int(nyb)):
        for x in xrange(nx):
            # select the pixels in this chunk of rows
            selectRowIndsLo = int(yb*yint)
            selectRowIndsHi = int(np.min([(yb+1)*yint,ny])) # this indexing is potentially hazardous
            s = residualimage[selectRowIndsLo : selectRowIndsHi, x]
            goodInds = np.where(np.isfinite(s))[0] # select the pixels with finite values
            goodIndsCount = len(goodInds)
            # make sure we have enough pixels to process
            if goodIndsCount < 2:
                continue
            s = np.sort(s[goodInds]) #now actually select and sort them
            ns = len(s)
            # exclude the tails of the distribution
            keepIndsLo = int(np.floor(0.03*ns))
            keepIndsHi = int(np.ceil(0.97*ns))
            s = s[keepIndsLo:keepIndsHi]
            sigmaimage[selectRowIndsLo : selectRowIndsHi, x] = np.std(s,ddof=1)

    sigmaimage = np.sqrt(sigmaimage**2 + (np.fabs(filterimage*sfactor) + np.fabs(xmedimage)))  # this scaling assumes gain ~ 1
    residualnsigmaimage = residualimage / sigmaimage


    residualnsigmaimage_ravel = residualnsigmaimage.ravel() # this is a view of residualnsigmaimage
    zapimage_ravel = zapimage.ravel() # this is a view of zapimage
    skysubimage_ravel = skysubimage.ravel() # this is a view of zapimage

    crcores = np.where(residualnsigmaimage_ravel > nsigma)[0]
    newzaps = len(crcores)
    if newzaps > 0:
        zapimage_ravel[crcores] = 1


    # #this is hacky but MIGHT needed for host analysis 
    # zapimage[50:80, 1280:1360] = 0 #halpha keck
    # zapimage_ravel = zapimage.ravel()


    # res = writefits(zapimage,'zapimage1.fits',CLOBBER=True)
    if DEBUG:
        res = writefits(zapimage,'{}/zapimage.fits'.format(DEBUG_DIR),CLOBBER=True)


    outStr = 'Flagged {} initial affected pixels before percolation.'.format(newzaps)
    print(outStr)

    nperczap = 0
    # iterCount = 0
    iterCount = 1 #testing
    d0 = nx
    d1 = ny

    # percolate outward to get all pixels covered by each cosmic ray.
    while iterCount < 32:
        nextperc = np.where(zapimage_ravel == iterCount)[0]
        ct = len(nextperc)
        iterCount += 1
        newzaps = 0

        if ct <= 3:
            break

        nrays = len(nextperc)
        # for c in range(nrays):
        for c in xrange(nrays):
            ci = nextperc[c]

            # avoid the detector edges
            if ci < d0-1 or ci > nx*ny-d0-2:
                continue

                 # here's the structure assumed in the IDL
                 # ci is the CR pixel we're looking at, 
                 # this is trying to check the neighboring pixels

                 # ci-d0-1,    ci-d0,    ci-d0+1, 
                 #    ci-1,      ci      ci+1, 
                 # ci+d0-1,    ci+d0,    ci+d0+1 

                 # the way the IDL where() function works in the original pzap code
                 # makes indexing the arrays this way trivial since it does the
                 # ravel implicitly, but np.where() works differently. The easiest way
                 # to handle this is to retain this neighbor indexing from the
                 # original IDL pzap code and just run the np.where() on ravel'd
                 # numpy arrays.

                 # just to be explicit, if you have an N-D numpy array x, then
                 # x.ravel() will return a VIEW of that array mapped to a
                 # 1-D array, so a 2x2 array becomes a 4 element 1-D array. The
                 # important thing to remember is that if you modify a VIEW of
                 # an array, it will modify your original array too.

            coarseblockpos = np.array([       ci-d0, 
                                        ci-1,        ci+1, 
                                              ci+d0])
            newzap = np.where( (np.fabs(residualnsigmaimage_ravel[coarseblockpos]) > subsigma) &
                               (zapimage_ravel[coarseblockpos] == 0) )[0]
            addzap = len(newzap)

            if addzap > 0:
                zapimage_ravel[coarseblockpos[newzap]] = iterCount
                newzaps += addzap
            nperczap += addzap


    # finally, zap anything hiding in a cosmic ray "corner" (three neighbors are cosmic ray pixels)
    countneighbor = np.zeros((ny,nx))
    countneighbor_ravel = countneighbor.ravel()
    nextperc = np.where(zapimage_ravel > 0)[0]
    ct = len(nextperc)
    if ct > 3:
        nrays = len(nextperc)
        # for c in range(nrays):
        for c in xrange(nrays):
            ci = nextperc[c]

            # avoid the detector edges
            # this wasn't present in the IDL code, but I believe it's necessary
            if ci < d0-1 or ci > nx*ny-d0-2:
                continue
            coarseblockpos = np.array([ci-d0-1,    ci-d0,    ci-d0+1, 
                                          ci-1,              ci+1, 
                                       ci+d0-1,    ci+d0,    ci+d0+1])
            countneighbor_ravel[coarseblockpos] = countneighbor_ravel[coarseblockpos] + 1

        newzap = np.where((countneighbor_ravel > 3) &
                          (zapimage_ravel == 0))[0]
        newzaps = len(newzap)
        if newzaps > 0:
            zapimage_ravel[newzap] = iterCount +1


    # actually do the zapping by replacing with the sky plus local median value.
    ibad = np.where(zapimage_ravel != 0)[0]
    nbad = len(ibad)

    if nbad > 0:
        skysubimage_ravel[ibad] = np.nan # NaNs ignored by median() and derivative image modelers 

        # filterimage is really the replacement image (without sky and sources)
        filterimage = signal.medfilt2d(skysubimage,[boxsize,boxsize])

        witer = 0
        ctnan = -1

        while ctnan != 0 and witer <= 2:
            if ctnan > 0:
                # note not active 0th iteration, since we already did the median above
                # so its ok that filterbad is not yet defined...
                filterimage.ravel()[filterbad] = signal.medfilt2d(filterimage, [boxsize,boxsize]).ravel()[filterbad]

            filterbad = np.where( ~np.isfinite(filterimage.ravel()) & 
                                   np.isfinite(xmedimage.ravel()) &
                                   np.isfinite(ymedimage.ravel()))[0]
            ctnan = len(filterbad)
            witer += 1
            if witer == 2 and ctnan > 0:
                filterimage.ravel()[filterbad] = 0   # give up

        outimg.ravel()[ibad] = filterimage.ravel()[ibad] + ymedimage.ravel()[ibad] + xmedimage.ravel()[ibad]
        outmask.ravel()[ibad] = 1

    nzap += nbad

    outStr = 'Zapped {} pixels in {} seconds'.format(nzap,int(time.time()-tstart))
    print(outStr)

    histStr = 'Processed by pyzapspec UT {}'.format(datetime.datetime.now())
    header.add_history(histStr)

    # res = writefits(zapimage,'zapimage2.fits',CLOBBER=True)
    if DEBUG:
        res = writefits(ymedimage,'{}/ymedimage.fits'.format(DEBUG_DIR),CLOBBER=True)
        res = writefits(xmedimage,'{}/xmedimage.fits'.format(DEBUG_DIR),CLOBBER=True)
        res = writefits(ysubimage,'{}/ysubimage.fits'.format(DEBUG_DIR),CLOBBER=True)
        res = writefits(xsubimage,'{}/xsubimage.fits'.format(DEBUG_DIR),CLOBBER=True)
        res = writefits(skysubimage,'{}/fullsubimage.fits'.format(DEBUG_DIR),CLOBBER=True)
        res = writefits(zapimage,'{}/zapimage.fits'.format(DEBUG_DIR),CLOBBER=True)
        res = writefits(outimg,outfile,header=header,CLOBBER=True)
        res = writefits(outmask,maskfile,CLOBBER=True)
        pdb.set_trace()

    if WRITE_OUTFILE:
        res = writefits(outimg,outfile,header=header,CLOBBER=True)
        # res = writefits(outmask,maskfile,CLOBBER=True)

    return outimg,outmask,header

    return 0


def main(*args,**kwargs):
    ''' Run the main pyzapspec machinery '''

    infile = args[0]

    # get rid of NoneTypes in kwargs
    for key,value in kwargs.items():
        if value is None:
            garbage = kwargs.pop(key)

    if not kwargs.get('outfile'):
        dirpath,infile_base = os.path.split(infile)
        outfile_base = 'cosmic_{}'.format(infile_base)
        outfile = os.path.join(dirpath,outfile_base)
        kwargs['outfile'] = outfile


    if not kwargs.get('maskfile'):
        dirpath,infile_base = os.path.split(infile)
        maskfile_base = 'cosmic_{}_mask.fits'.format(infile_base.split('.')[0])
        maskfile = os.path.join(dirpath,maskfile_base)
        kwargs['maskfile'] = maskfile

    outfile = kwargs.get('outfile')
    maskfile = kwargs.get('maskfile')
    CLOBBER = kwargs.get('CLOBBER',False)
    VERBOSE = kwargs.get('VERBOSE',False)
    DEBUG = kwargs.get('DEBUG',False)

    # run it
    outimg,outmask,header = pyzapspec(infile,**kwargs)

    # attempt to output stuff safely
    if not os.path.isfile(outfile):
        res = writefits(outimg,outfile,header=header,CLOBBER=CLOBBER)
        res = writefits(outmask,maskfile,CLOBBER=CLOBBER)

    elif os.path.isfile(outfile):
        promptStr = 'Output file {} exists...overwrite? y/[n]: '.format(outfile)
        usrAns = raw_input(promptStr)
        if usrAns.strip().upper()[0] == 'Y':
            print('Overwritting...')
            # set clobber to true since we've cleared it with the user
            res = writefits(outimg,outfile,header=header,CLOBBER=True)
            res = writefits(outmask,maskfile,header=header,CLOBBER=True)
        else:
            print('Ok, doing nothing and exiting.')


    return 0

if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)

