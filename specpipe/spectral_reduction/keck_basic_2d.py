#!/usr/bin/env python
from __future__ import print_function
import os,sys,pdb,shutil,glob,argparse
import numpy as np

from astropy.io import fits
from astropy.io import ascii
from scipy import signal, ndimage


def load_sections(string, fmt_iraf=True):
    """
    Modified from pypit.core.parse.load_sections -- Jon Brown
    
    From the input string, return the coordinate sections

    Parameters
    ----------
    string : str
      character string of the form [x1:x2,y1:y2]
      x1 = left pixel
      x2 = right pixel
      y1 = bottom pixel
      y2 = top pixel
    fmt_iraf : bool
      Is the variable string in IRAF format (True) or
      python format (False)

    Returns
    -------
    sections : list (or None)
      the detector sections
    """
    xyrng = string.strip('[]()').split(',')
    if xyrng[0] == ":":
        xyarrx = [0, 0]
    else:
        xyarrx = xyrng[0].split(':')
        # If a lower/upper limit on the array slicing is not given (e.g. [:100] has no lower index specified),
        # set the lower/upper limit to be the first/last index.
        if len(xyarrx[0]) == 0: xyarrx[0] = 0
        if len(xyarrx[1]) == 0: xyarrx[1] = -1
    if xyrng[1] == ":":
        xyarry = [0, 0]
    else:
        xyarry = xyrng[1].split(':')
        # If a lower/upper limit on the array slicing is not given (e.g. [5:] has no upper index specified),
        # set the lower/upper limit to be the first/last index.
        if len(xyarry[0]) == 0: xyarry[0] = 0
        if len(xyarry[1]) == 0: xyarry[1] = -1
    if fmt_iraf:
        xmin = max(0, int(xyarry[0])-1)
        xmax = int(xyarry[1])
        ymin = max(0, int(xyarrx[0])-1)
        ymax = int(xyarrx[1])
    else:
        xmin = max(0, int(xyarrx[0]))
        xmax = int(xyarrx[1])
        ymin = max(0, int(xyarry[0]))
        ymax = int(xyarry[1])
    return [[xmin, xmax], [ymin, ymax]]
    
    
def sec2slice(subarray, one_indexed=False, include_end=False, require_dim=None, transpose=False):
    """
    Modified from pypit.core.parse.sec2slice -- Jon Brown
    
    Convert a string representation of an array subsection (slice) into
    a list of slice objects.

    Args:
        subarray (str):
            The string to convert.  Should have the form of normal slice
            operation, 'start:stop:step'.  The parser ignores whether or
            not the string has the brackets '[]', but the string must
            contain the appropriate ':' and ',' characters.
        one_indexed (:obj:`bool`, optional):
            The string should be interpreted as 1-indexed.  Default
            is to assume python indexing.
        include_end (:obj:`bool`, optional):
            **If** the end is defined, adjust the slice such that
            the last element is included.  Default is to exclude the
            last element as with normal python slicing.
        require_dim (:obj:`int`, optional):
            Test if the string indicates the slice along the proper
            number of dimensions.
        transpose (:obj:`bool`, optional):
            Transpose the order of the returned slices.  The
            following are equivalent::
                
                tslices = parse_sec2slice('[:10,10:]')[::-1]
                tslices = parse_sec2slice('[:10,10:]', transpose=True)

    Returns:
        tuple: A tuple of slice objects, one per dimension of the
        prospective array.

    Raises:
        TypeError:
            Raised if the input `subarray` is not a string.
        ValueError:
            Raised if the string does not match the required
            dimensionality or if the string does not look like a
            slice.
    """
    # Check it's a string
    if not isinstance(subarray, basestring):
        raise TypeError('Can only parse string-based subarray sections.')
    # Remove brackets if they're included
    sections = subarray.strip('[]').split(',')
    # Check the dimensionality
    ndim = len(sections)
    if require_dim is not None and ndim != require_dim:
        raise ValueError('Number of slices ({0}) in {1} does not match '.format(ndim, subarray) + 
                         'required dimensions ({0}).'.format(require_dim))
    # Convert the slice of each dimension from a string to a slice
    # object
    slices = []
    for s in sections:
        # Must be able to find the colon
        if ':' not in s:
            raise ValueError('Unrecognized slice string: {0}'.format(s))
        # Initial conversion
        _s = [ None if x == '' else int(x) for x in s.split(':') ]
        if len(_s) > 3:
            raise ValueError('String as too many sections.  Must have format \'start:stop:step\'.')
        if len(_s) < 3:
            # Include step
            _s += [ None ]
        if one_indexed:
            # Decrement to convert from 1- to 0-indexing
            _s = [ None if x is None else x-1 for x in _s ]
        if include_end and _s[1] is not None:
            # Increment to include last 
            _s[1] += 1
        # Append the new slice
        slices += [slice(*_s)]
    return tuple(slices[::-1] if transpose else slices)


def read_lris(raw_file, det=None, TRIM=False):
    """
    Modified from pypeit.spectrographs.keck_lris.read_lris -- Jon Brown
    
    Read a raw LRIS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename
    det : int, optional
      Detector number; Default = both
    TRIM : bool, optional
      Trim the image?

    Returns
    -------
    array : ndarray
      Combined image 
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        #msgs.error("Found {:d} files matching {:s}".format(len(fil)))
        errStr = "Found {:d} files matching {:s}".format(len(fil),raw_file)
        print(errStr)
        sys.exit()

    # Read
    outStr = "Reading LRIS file: {:s}".format(fil[0])
    print(outStr)
    hdu = fits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']
    
    # get the detector
    # this just checks if its the blue one and assumes red if not
    # note the red fits headers don't even have this keyword???
    if head0['INSTRUME'] == 'LRISBLUE':
        redchip = False
    else:
        redchip = True
        

    # Setup for datasec, oscansec
    dsec = []
    osec = []
    nxdata_sum = 0

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # First read over the header info to determine the size of the output array...
    n_ext = len(hdu)-1  # Number of extensions (usually 4)
    xcol = []
    xmax = 0
    ymax = 0
    xmin = 10000
    ymin = 10000
    for i in np.arange(1, n_ext+1):
        theader = hdu[i].header
        detsec = theader['DETSEC']
        if detsec != '0':
            # parse the DETSEC keyword to determine the size of the array.
            x1, x2, y1, y2 = np.array(load_sections(detsec, fmt_iraf=False)).flatten()

            # find the range of detector space occupied by the data
            # [xmin:xmax,ymin:ymax]
            xt = max(x2, x1)
            xmax = max(xt, xmax)
            yt =  max(y2, y1)
            ymax = max(yt, ymax)

            # find the min size of the array
            xt = min(x1, x2)
            xmin = min(xmin, xt)
            yt = min(y1, y2)
            ymin = min(ymin, yt)
            # Save
            xcol.append(xt)

    # determine the output array size...
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    # change size for binning...
    nx = nx // xbin
    ny = ny // ybin

    # Update PRECOL and POSTPIX
    precol = precol // xbin
    postpix = postpix // xbin

    # Deal with detectors
    if det in [1,2]:
        nx = nx // 2
        n_ext = n_ext // 2
        det_idx = np.arange(n_ext, dtype=np.int) + (det-1)*n_ext
        ndet = 1
    elif det is None:
        ndet = 2
        det_idx = np.arange(n_ext).astype(int)
    else:
        raise ValueError('Bad value for det')

    # change size for pre/postscan...
    if not TRIM:
        nx += n_ext*(precol+postpix)
        ny += preline + postline

    # allocate output array...
    array = np.zeros( (nx, ny) )
    gain_array = np.zeros( (nx, ny) )
    order = np.argsort(np.array(xcol))

    # insert extensions into master image...
    for kk, i in enumerate(order[det_idx]):

        # grab complete extension...
        data, gaindata, predata, postdata, x1, y1 = lris_read_amp(hdu, i+1, redchip=redchip)
                            
        # insert components into output array...
        if not TRIM:
            # insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = kk*precol
            xe = xs + nxpre

            array[xs:xe, :] = predata
            gain_array[xs:xe, :] = predata

            # insert data...
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]
            
            # JB: have to track the number of xpixels
            xs = n_ext*precol + nxdata_sum
            xe = xs + nxdata
            
            # now log how many pixels that was
            nxdata_sum += nxdata
            
            # Data section
            #section = '[{:d}:{:d},{:d}:{:d}]'.format(preline,nydata-postline, xs, xe)  # Eliminate lines
            section = '[{:d}:{:d},{:d}:{:d}]'.format(preline,nydata, xs, xe)  # DONT eliminate lines

            dsec.append(section)
            array[xs:xe, :] = data   # Include postlines
            gain_array[xs:xe, :] = gaindata   # Include postlines

            #; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext*postpix + kk*postpix
            xe = xs + nxpost 
            section = '[:,{:d}:{:d}]'.format(xs, xe)
            osec.append(section)

            array[xs:xe, :] = postdata
            gain_array[xs:xe, :] = postdata

        else:
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]

            xs = (x1-xmin)//xbin
            xe = xs + nxdata 
            ys = (y1-ymin)//ybin
            ye = ys + nydata - postline

            yin1 = preline
            yin2 = nydata - postline 

            array[xs:xe, ys:ye] = data[:, yin1:yin2]
            gain_array[xs:xe, ys:ye] = gaindata[:, yin1:yin2]

    # make sure BZERO is a valid integer for IRAF
    obzero = head0['BZERO']
    head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero
    
    # Return, transposing array back to goofy Python indexing
    return array.T,gain_array.T, head0, (dsec, osec)

def lris_read_amp(inp, ext, redchip=False, applygain=True):
    """
    Modified from pypeit.spectrographs.keck_lris.lris_read_amp -- Jon Brown
    
    Read one amplifier of an LRIS multi-extension FITS image

    Parameters
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns
    -------
    data
    predata
    postdata
    x1
    y1

    ;------------------------------------------------------------------------
    function lris_read_amp, filename, ext, $
      linebias=linebias, nobias=nobias, $
      predata=predata, postdata=postdata, header=header, $
      x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
    ;------------------------------------------------------------------------
    ; Read one amp from LRIS mHDU image
    ;------------------------------------------------------------------------
    """
    # Parse input
    if isinstance(inp, str):
        hdu = fits.open(inp)
    else:
        hdu = inp

    # Get the pre and post pix values
    # for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol = head0['precol']
    postpix = head0['postpix']
    postline = head0['postline']
        
    # Deal with binning
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol//xbin
    postpix = postpix//xbin

    # get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
    tsize = temp.shape
    nxt = tsize[0]    

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(load_sections(datasec, fmt_iraf=False)).flatten()
    
    # grab the components...
    predata = temp[0:precol, :]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.  They also appear inconsistent between
    # LRISr and LRISb!
    #data     = temp[xdata1-1:xdata2-1,*]
    #data = temp[xdata1:xdata2+1, :]
    
    # JB: LRIS-R is windowed differently, so the default pypeit checks fail
    # xshape is calculated from datasec.
    # For blue, its 1024, 
    # For red, the chip dimensions are different AND the observations are windowed
    # In windowed mode each amplifier has differently sized data sections
    if not redchip:
        xshape = 1024 // xbin # blue
    else:
        xshape = (xdata2-xdata1 + 1//xbin) # red

    # do some sanity checks
    if (xdata1-1) != precol:
        #msgs.error("Something wrong in LRIS datasec or precol")
        errStr = 'Something wrong in LRIS datasec or precol'
        print(errStr)
        pdb.set_trace()
        sys.exit()
        
    if (xshape+precol+postpix) != temp.shape[0]:
        #msgs.error("Wrong size for in LRIS detector somewhere.  Funny binning?")
        errStr = 'Wrong size for in LRIS detector somewhere.  Funny binning?'
        print(errStr)
        pdb.set_trace()
        sys.exit()
                
    data = temp[precol:precol+xshape,:]
    postdata = temp[nxt-postpix:nxt, :]

    # flip in X as needed...
    if x1 > x2:
        xt = x2
        x2 = x1
        x1 = xt
        data = np.flipud(data) #reverse(temporary(data),1)

    # flip in Y as needed...
    if y1 > y2:
        yt = y2
        y2 = y1
        y1 = yt
        data = np.fliplr(data)
        predata = np.fliplr(predata)
        postdata = np.fliplr(postdata)
        
    # construct gain image
    if redchip:
        # get the amplifier
        if header['extname'].strip().upper() == 'VIDINP1':
            gain = 1.255 # keck page
            #gain = 1.240 # gain page
        elif header['extname'].strip().upper() == 'VIDINP2':
            gain = 1.180 # keck page
            #gain = 1.180 # gain page
        elif header['extname'].strip().upper() == 'VIDINP3':
            gain = 1.191 # keck page
            #gain = 1.250 # gain page
        elif header['extname'].strip().upper() == 'VIDINP4':
            gain = 1.162 # keck page
            #gain = 1.173 # gain page
        else:
            errStr = 'Unsupported extension??? Exiting after debugging.'
            print(errStr)
            pdb.set_trace()
            sys.exit()
    else:
        # get the amplifier
        if header['extname'].strip().upper() == 'VIDINP1':
            gain = 1.55 # keck page
            #gain = 1.758 # gain page
        elif header['extname'].strip().upper() == 'VIDINP2':
            gain = 1.56 # keck page
            #gain = 1.739 # gain page
        elif header['extname'].strip().upper() == 'VIDINP3':
            gain = 1.63 # keck page
            #gain = 1.636 # gain page
        elif header['extname'].strip().upper() == 'VIDINP4':
            gain = 1.70 # keck page
            #gain = 1.681 # gain page
        else:
            errStr = 'Unsupported extension??? Exiting after debugging.'
            print(errStr)
            pdb.set_trace()
            sys.exit()
    gaindata = 0.*data + gain
            
    return data, gaindata, predata, postdata, x1, y1
    
    
def subtract_overscan(rawframe, numamplifiers, datasec, oscansec, gain_image=None,
                      method='savgol', params=[5, 65]):
    """
    Modified from pypeit.core.procimg.subtract_overscan -- Jon Brown
    
    Subtract overscan

    Args:
        frame (:obj:`numpy.ndarray`):
            Frame from which to subtract overscan
        numamplifiers (int):
            Number of amplifiers for this detector.
        datasec (list):
            List of slices, one per amplifier, that contain the data in
            the raw frame.  The slices and be lists of slice ojects or
            strings.  If they are strings, :func:`parse.sec2slice` is
            used to convert them for use in the function.
        oscansec (list):
            List of slices, one per amplifier, that contain the
            overscane regions in the raw frame.  The slices and be lists
            of slice ojects or strings.  If they are strings,
            :func:`parse.sec2slice` is used to convert them for use in
            the function.
        method (:obj:`str`, optional):
            The method used to fit the overscan region.  Options are
            polynomial, savgol, median.
        params (:obj:`list`, optional):
            Parameters for the overscan subtraction.  For
            method=polynomial, set params = order, number of pixels,
            number of repeats ; for method=savgol, set params = order,
            window size ; for method=median, params are ignored.

    Returns:
        :obj:`numpy.ndarray`: The input frame with the overscan region
        subtracted
    """
    # Check input
    if len(datasec) != numamplifiers or len(oscansec) != numamplifiers:
        pdb.set_trace()
        raise ValueError('Number of amplifiers does not match provided image sections.')

    # If the input image sections are strings, convert them
    if isinstance(datasec[0], str):
        _datasec = [None]*len(datasec.copy())
        for i in range(numamplifiers):
            tmp = sec2slice(datasec[i], require_dim=2)
            _datasec[i] = sec2slice(datasec[i], require_dim=2)
    else:
        _datasec = datasec
            
    if isinstance(oscansec[0], str):
        _oscansec = [None]*len(oscansec.copy())
        for i in range(numamplifiers):
            _oscansec[i] = sec2slice(oscansec[i], require_dim=2)
    else:
        _oscansec = oscansec
    
    # Check that there are no overlapping data sections
    testframe = np.zeros_like(rawframe, dtype=int)
    for i in range(numamplifiers):
        testframe[_datasec[i]] += 1
    if np.any(testframe > 1):
        raise ValueError('Image has overlapping data sections!')

    # Copy the data so that the subtraction is not done in place
    nobias = rawframe.copy()

    # Perform the bias subtraction for each amplifier
    for i in range(numamplifiers):
        # Pull out the overscan data
        overscan = rawframe[_oscansec[i]]

        # Shape along at least one axis must match
        data_shape = rawframe[_datasec[i]].shape
        if not np.any([ dd == do for dd, do in zip(data_shape, overscan.shape)]):
            #msgs.error('Overscan sections do not match amplifier sections for'
            #           'amplifier {0}'.format(i+1))
            errStr = 'Overscan sections do not match amplifier sections for '
            errStr += 'amplifier {0}'.format(i+1)
            print(errStr)
            pdb.set_trace()
            sys.exit()
        compress_axis = 1 if data_shape[0] == overscan.shape[0] else 0
        
        # Fit/Model the overscan region
        osfit = np.median(overscan) if method.lower() == 'median' \
                        else np.median(overscan, axis=compress_axis)

        if method.lower() == 'polynomial':
            # TODO: Use np.polynomial.polynomial.polyfit instead?
            c = np.polyfit(np.arange(osfit.size), osfit, params[0])
            ossub = np.polyval(c, np.arange(osfit.size))
        elif method.lower() == 'savgol':
            ossub = signal.savgol_filter(osfit, params[1], params[0])
        elif method.lower() == 'median':
            # Subtract scalar and continue
            nobias[_datasec[i]] -= osfit
            continue
        else:
            raise ValueError('Unrecognized overscan subtraction method: {0}'.format(method))

        # Subtract along the appropriate axis
        nobias[_datasec[i]] -= (ossub[:,None] if compress_axis == 1 else ossub[None,:])
        
        # apply gain correction
        if gain_image is not None:
            nobias[_datasec[i]] *= gain_image[_datasec[i]]
            
    # return bias subtracted gain corrected frame        
    return nobias


def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Main driver for the basic Keck CCD reduction. '
    #descStr += 'Recommended calling sequence: \n \n'
    #descStr += '$ python keck_basic_2d.py -v -c \n'
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # required args
    #parser.add_argument('requried_arg',type=str,
    #                    help='a required arguement')

    # optional
    parser.add_argument('-v','--verbose',
                        help='print diagnostic info',action='store_true')
    parser.add_argument('-c','--clobber',action='store_true',
                        help='Clobber files already in pre_reduced/ but not subdirs')
    parser.add_argument('-f','--full-clean',action='store_true',
                        help='Do a complete wipe of pre_reduced, including subdirs')

    parser.add_argument('--no-reorient',action='store_true',
                        help='Do not reorient to wavelength increasing rightward')
    parser.add_argument('--no-trim',action='store_true',
                        help='Do not trim 2D image to hardcoded section')
    parser.add_argument('--is-dflat',action='store_true',
                        help='The current file is a dome flat')
    parser.add_argument('--mask-middle-blue',action='store_true',
                        help='Mask the middle section of rows on the blue detector')
    parser.add_argument('--mask-middle-red',action='store_true',
                        help='Mask the middle section of rows on the red detector')

    # parse
    cmdArgs = parser.parse_args()

    # package up
    args = () # no args implemented yet
    kwargs = {}
    kwargs['VERBOSE'] = cmdArgs.verbose
    kwargs['CLOBBER'] = cmdArgs.clobber
    kwargs['FULL_CLEAN'] = cmdArgs.full_clean
    kwargs['REORIENT'] = not cmdArgs.no_reorient
    kwargs['TRIM'] = not cmdArgs.no_trim
    kwargs['ISDFLAT'] = cmdArgs.is_dflat
    kwargs['MASK_MIDDLE_BLUE'] = cmdArgs.mask_middle_blue
    kwargs['MASK_MIDDLE_RED'] = cmdArgs.mask_middle_red

    return (args,kwargs)




def main(rawFiles,*args,**kwargs):
    '''
    Run basic 2D CCD reduction on Keck LRIS data

    Parameters
    ----------
    CLOBBER : bool, optional (default=False)
        Overwrite the individual files in pre_reduced, but
        do not wipe subdirectories
    FULL_CLEAN : bool, optional (default=False)
        Completely wipe pre_reduced and all subdirectories
    BPM : bool, optional (default=False)
        Mask bad pixels (not currently implemented)
    PIXEL_FLOOR : bool, optional (default=False)
        Removes negative pixel values
    REORIENT : bool, optional (default=True)
        Transposes to wavelength increasing rightward. LRIS images
        are by default (spatial, spectral), and in general should
        be transposed to (spectral, spatial).
    TRIM : bool, optional (default=True)
        Trim to some hard coded section of the detector
    MASK_MIDDLE_BLUE : bool (default=False)
        Mask the middle section of the blue images
    MASK_MIDDLE_RED : bool (default=False)
        Mask the middle section of the red images (useful if
        there's wildly disparate values that make the iraf
        windowing tedious)
        
    Returns
    -------
    int : 0, and writes files to disk
    '''

    # unpack supported kwargs
    CLOBBER = kwargs.get('CLOBBER',False) # 
    FULL_CLEAN = kwargs.get('FULL_CLEAN',False) # this will completely wipe pre_reduced.
    BPM = kwargs.get('BPM',False) # no bad pixel mask available (at least for our binning)
    PIXEL_FLOOR = kwargs.get('PIXEL_FLOOR',False) 
    REORIENT = kwargs.get('REORIENT',True) 
    TRIM = kwargs.get('TRIM',False) 
    ISDFLAT = kwargs.get('ISDFLAT',False) 
    RED_AMP_BAD = kwargs.get('RED_AMP_BAD',False) 
    MASK_MIDDLE_BLUE = kwargs.get('MASK_MIDDLE_BLUE',False)
    MASK_MIDDLE_RED = kwargs.get('MASK_MIDDLE_RED',False)

    # pre_reduced does not exist, needs to be made
    if not os.path.isdir('pre_reduced/'):
        os.mkdir('pre_reduced/')
        
    # pre_reduced exists, but we want to clobber/do a clean reduction
    elif FULL_CLEAN:
        
        promptStr = 'Do you really want to wipe pre_reduced? [y/n]: '
        usrRespOrig = raw_input(promptStr)
        if usrRespOrig and usrRespOrig[0].strip().upper() == 'Y':

            # remove all pre_reduced files
            shutil.rmtree('pre_reduced')
            os.mkdir('pre_reduced/')
      
    # pre_reduced exists, need to document what is there
    else:
        
        # get files
        preRedFiles = glob.glob('pre_reduced/*.fits')
    
  
    # get all raw files, sort them
    # rawFiles = sorted(glob.glob('???????_????.fits'))

    # loop over raw files, if the destination exists, do nothing
    # otherwise, do the bias/reorient/trim/output
    for i in xrange(len(rawFiles)):
        
        rawFile = rawFiles[i]
        oScanFile = 'pre_reduced/to{}'.format(rawFile)
        oScanFile_amp1 = oScanFile.split('.')[0]+'_amp1.'+oScanFile.split('.')[1]
        oScanFile_amp2 = oScanFile.split('.')[0]+'_amp2.'+oScanFile.split('.')[1]
        
        if (not os.path.isfile(oScanFile) and not os.path.isfile(oScanFile_amp1) and not os.path.isfile(oScanFile_amp2)) or CLOBBER:
            
            # read file
            img,gain_img,head,secs = read_lris(rawFile)
            # img,gain_img,head,secs = read_lris(rawFile, TRIM=TRIM)
            # print (secs)
            dsec,osec = np.array(secs[0]),np.array(secs[1])
            xbin, ybin = [int(ibin) for ibin in head['BINNING'].split(',')]
                        
            # get number of extensions for nAmps
            tmpHDU = fits.open(rawFile)
            nAmps = len(tmpHDU) - 1
            
            # perform oscan/bias subtraction
            # noBiasImg = subtract_overscan(img,nAmps,
            #                               dsec,osec,
            #                               gain_image=gain_img,
            #                               method='median',    # median should be fine
            #                               params=[5,65])       # default savgol params
            noBiasImg = subtract_overscan(img,nAmps,
                                          dsec,osec,
                                          gain_image=gain_img,
                                          method='polynomial',    # median should be fine
                                          params=[5,65])       # default savgol params
                                          
            # mask bad pixels
            if BPM:
                print('Bad pixel masking not yet implemented...')
                
            # apply floor to pixel values    
            if PIXEL_FLOOR:
                noBiasImg[noBiasImg < 0] = 0.
                
                                              
            # rotate/flip/transpose (wavelength increasing w/ increasing row)
            if REORIENT:
                outImg = noBiasImg.T
            
            # trim
            if TRIM:
                if rawFile[0] == 'b':
                    if ISDFLAT:
                        outImg_amp1 = outImg[2260//xbin:2800//xbin,:]
                        outImg_amp2 = outImg[1800//xbin:2260//xbin,:]
                    else:
                        outImg = outImg[1800//xbin:2800//xbin,:]
                else:
                    if nAmps == 2:
                        if ISDFLAT and not RED_AMP_BAD:
                            outImg_amp1 = outImg[290:575,:-55] # trimming for windowed and removes bottom amplifier (assumes xbin = 2)
                            outImg_amp2 = outImg[0:290,:-55]
                        elif ISDFLAT and RED_AMP_BAD:
                            # outImg_amp1 = outImg[290:575,:-55]
                            outImg_amp1 = outImg[290:545,:-55]
                        elif not ISDFLAT and RED_AMP_BAD:
                            # outImg = outImg[290:575,:-55]
                            outImg = outImg[290:545,:-55]
                        else:
                            # outImg = outImg[600//xbin:1100//xbin,:]
                            outImg = outImg[0:575,:-55]

                    else:
                        outImg = outImg[1600//xbin:2600//xbin,:-55]
                        outImg_amp1 = outImg[250:575,:] #12/7/18 ignoring middle of red
                        outImg_amp2 = outImg[0:250,:]
                        # This is a weird case. Red is being read out in full frame
                        # and if its a mask, the slits aren't on the normal longlist
                        # sections; we probably don't want to trim at all. Just save
                        # the full frame and let the use do whatever they want after.
                        
                        # This is only problematic if the also want the blue data,
                        # which is trimmed to longslit by default. The work around is
                        # to just move all the mask data to its own directory and run
                        # this with TRIM=False.
                        #pass # if red is full frame, we probabaly don't want to bin?

            if MASK_MIDDLE_RED and rawFile[0] == 'r':
                # these rows should be chosen programatically
                # this is much more aggressive than necessary
                outImg[190:250,:] = np.median(outImg)
            if MASK_MIDDLE_BLUE and rawFile[0] == 'b':
                # these rows should be chosen programatically
                # this is much more aggressive than necessary
                outImg[380:480,:] = np.median(outImg)
                        
            # adjust the header (these keywords aren't present by default)
            head['EXPTIME'] = head['ELAPTIME']
            head['DATE-OBS'] = head['DATE_BEG']
            head['BASIC-2D'] = 'DONE'
                
            # write the images  
            if ISDFLAT:  
                oScanFile_amp1 = oScanFile.split('.')[0]+'_amp1.'+oScanFile.split('.')[1]
                hdu = fits.PrimaryHDU(outImg_amp1,head)
                if os.path.isfile(oScanFile_amp1):
                    os.remove(oScanFile_amp1)
                hdu.writeto(oScanFile_amp1,output_verify='ignore')  

                if rawFile[0] == 'b':
                    oScanFile_amp2 = oScanFile.split('.')[0]+'_amp2.'+oScanFile.split('.')[1]
                    hdu = fits.PrimaryHDU(outImg_amp2,head)
                    if os.path.isfile(oScanFile_amp2):
                        os.remove(oScanFile_amp2)
                    hdu.writeto(oScanFile_amp2,output_verify='ignore')    
                elif not RED_AMP_BAD:
                    oScanFile_amp2 = oScanFile.split('.')[0]+'_amp2.'+oScanFile.split('.')[1]
                    hdu = fits.PrimaryHDU(outImg_amp2,head)
                    if os.path.isfile(oScanFile_amp2):
                        os.remove(oScanFile_amp2)
                    hdu.writeto(oScanFile_amp2,output_verify='ignore')
            else:
                hdu = fits.PrimaryHDU(outImg,head)
                if os.path.isfile(oScanFile):
                    os.remove(oScanFile)
                hdu.writeto(oScanFile,output_verify='ignore')
                        
        else:
            outStr = 'File exists: {}'.format(oScanFile)
            print(outStr)
   
    return 0
  
if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)

