#!/usr/bin/env python
from __future__ import print_function
import sys, os, glob, numpy as np, shutil, copy, warnings
from astropy.io import fits
from astropy.table import Table,Column
from astropy.time import Time
import progressbar
import astromatic_wrapper as aw
from scipy.odr import *
warnings.filterwarnings('ignore')

# Wrapper for GSAOI tool disco_stu.  We mainly want this for aligning
# and reprojecting images.  We'll use a separate method for optimal stacking
def run_disco_stu(files, options=''):
    cmd = 'disco {o} --no_stack {f}'.format(f=files, o=options)
    os.system(cmd)

def stack_files(files, obj):

    print('Starting stack routine...')
    print('Need to stack {0} images'.format(len(files)))
    print('Splitting data and variance extensions and masking...')

    # Separate data and weight arrays
    for i,file in enumerate(files):

        print('file {0}/{1}:'.format(i+1,len(files)),file)
        hdu = fits.open(file)

        data = hdu[1].data
        var = hdu[2].data
        mask = hdu[3].data

        invvar = 1./var**2
        data[np.where(mask != 0)]=np.nan
        data[np.where(var == 0.)]=np.nan
        invvar[np.where(mask != 0)]=0.
        invvar[np.where(var == 0.)]=0.

        # Create separate data, inverse-variance, and mask images
        dathdu = fits.PrimaryHDU(data)
        varhdu = fits.PrimaryHDU(invvar)
        mashdu = fits.PrimaryHDU(mask)

        dathdu.header = hdu[1].header
        varhdu.header = hdu[2].header

        dathdu.writeto(file.replace('.fits','.data.fits'),
            overwrite=True,output_verify='silentfix')
        varhdu.writeto(file.replace('.fits','.var.fits'),
            overwrite=True,output_verify='silentfix')

    scales=[list(v) for v in zip([None]*len(files),files)]
    reftable=None
    radius=0.1/3600.
    for i,scale in enumerate(scales):
        if i==0:
            scales[i][0]=1.0
            reftable = Table.read(scale[1],format='fits',hdu=4)
            print('file {0}/{1}: {2}, scale={3}'.format(i+1,len(scales),
                scale[1],scales[i][0]))
            continue
        table = Table.read(scale[1],format='fits',hdu=4)
        # Do cross match
        matchtable = Table([[0.],[0.],[0.],[0.]],
            names=('refflux','referr','flux','err')).copy()[:0]
        # We can take advantage of source number for matching
        for row in table:
            match = (reftable['RA']-row['RA'])**2+\
                (reftable['DEC']-row['DEC'])**2<radius**2
            if len(reftable[match])==1:
                matchtable.add_row([reftable[match]['FLUX_AUTO'],
                    reftable[match]['FLUXERR_AUTO'],row['FLUX_AUTO'],
                    row['FLUXERR_AUTO']])

        # Require a minimum of 7 matches to include image
        if len(matchtable)>7:
            # Calculate relative scaling using x-y fit
            def scale_fit(x, p):
                return(p*x)

            model = Model(scale_fit)
            data = RealData(matchtable['flux'], matchtable['refflux'],
                sx=matchtable['err'], sy=matchtable['referr'])
            odr = ODR(data, model, beta0=[1.])
            out = odr.run()
            scales[i][0]=out.beta[0]

        print('file {0}/{1}: {2}, scale={3} with {4} sources'.format(i+1,
            len(scales), scale[1],scales[i][0],len(matchtable)))

    # Iterate through scales to decide if we want to use image for stack
    datafile = open('data.lis','w')
    varfile = open('var.lis','w')
    for scale in scales:
        if scale[0]:
            datafilename = scale[1].replace('.fits','.data.fits')
            datafile.write(datafilename+'\n')
            varfile.write(scale[1].replace('.fits','.var.fits')+'\n')

            # Add FLXSCALE key to data file (preferred SWarp flux scale method)
            hdu = fits.open(datafilename)
            hdu[0].header['FLXSCALE']=scale[0]
            hdu.writeto(datafilename, overwrite=True, output_verify='silentfix')

    datafile.close()
    varfile.close()

    kwargs = {
        'code': 'SWarp',
        'store_output': True,
        'config_file': None,
        'config': {
            'IMAGEOUT_NAME': obj+'.fits',
            'WEIGHTOUT_NAME': obj+'.weight.fits',
            'WEIGHT_IMAGE': '@var.lis',
            'WEIGHT_TYPE': 'MAP_WEIGHT',
            'COMBINE_TYPE': 'MEDIAN',
            'RESCALE_WEIGHTS': True,
            'RESAMPLING_TYPE': 'LANCZOS4',
            'FSCALASTRO_TYPE': 'FIXED',
            'PIXELSCALE_TYPE': 'MANUAL',
            'PIXEL_SCALE': '0.02,0.02',
            'INTERPOLATE': True,
            'WRITE_XML': False,
            'WRITE_FILEINFO': False,
            'BACK_TYPE': 'AUTO',
            'BACK_SIZE': 256,
            'FSCALE_KEYWORD': 'FLXSCALE',
            'VERBOSE_TYPE': 'NORMAL'
        }
    }

    swarp = aw.api.Astromatic(**kwargs)
    os.system('swarp -d > default.swarp')
    print('\n\nStarting SWarp stack...')
    output = swarp.run('@data.lis')

    os.remove('data.lis')
    os.remove('var.lis')
    os.remove('default.swarp')

def cleanup(dirs=['pyraf','uparm'],suffixes=['.log','.cl','.fits'],
    prefixes=['tmp']):

    for directory in dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)
    for suffix in suffixes:
        for file in glob.glob('*'+suffix):
            if suffix=='.fits' and 'flat' in file or 'dark' in file:
                continue
            os.remove(file)
    for prefix in prefixes:
        for file in glob.glob(prefix+'*'):
            os.remove(file)

    return(0)

def clean_amps(file, percentile=5.):
    hdu = fits.open(file)

    # Get the four science amps of the image
    for i,h in enumerate(hdu):
        if h.name=='SCI':
            data = h.data
            # Iterate through rows
            perlen = int(percentile/100.*len(data[0,:]))
            for row in np.arange(data.shape[0]):
                back = np.median(sorted(data[row,:])[0:perlen])
                data[row,:]=data[row,:]-back

            perlen = int(percentile/100.*len(data[:,0]))
            for col in np.arange(data.shape[1]):
                back = np.median(sorted(data[:,col])[0:perlen])
                data[:,col]=data[:,col]-back

            hdu[i].data = data

    hdu.writeto(file, overwrite=True, output_verify='silentfix')

def copy_raw(rawdir='raw/', destdir='./'):

    for file in glob.glob(rawdir+'/*.fits'):
        path, basefile = os.path.split(file)
        shutil.copyfile(file, destdir + '/' + basefile)

def sort_files(files):

    table = None

    for file in sorted(files):
        hdu = fits.open(file)

        obj = hdu[0].header['OBJECT']
        typ = hdu[0].header['OBSTYPE']
        exp = hdu[0].header['EXPTIME']
        xof = hdu[0].header['XOFFSET']
        yof = hdu[0].header['YOFFSET']
        dat = hdu[0].header['DATE-OBS']
        tim = hdu[0].header['TIME-OBS']

        t = Time(dat+'T'+tim)

        imagetyp = ''
        objname = ''
        if 'dark' in obj.lower() or 'dark' in typ.lower():
            imagetyp = 'DARK'
            objname = 'DARK'
        elif ('flat' in obj.lower() or 'flat' in typ.lower()):
            imagetyp = 'FLAT'
            objname = 'FLAT'
        elif ('flat' not in obj.lower() and 'dark' not in obj.lower() and
            'flat' not in typ.lower() and 'dark' not in typ.lower() and
            'object' in typ.lower()):
            imagetyp = 'OBJECT'
            objname = obj

            # Sanitizer object name
            objname = objname.replace('-offset','')
            if objname.upper().startswith('SN'):
                objname = objname[2:]
            elif objname.upper().startswith('AT'):
                objname = objname[2:]

            objname = objname.lower()
        else:
            # Can't identify object type
            continue

        if not table:
            table = Table([['X'*30],['X'*30],['X'*30],[0.0],[0.0],[0.0],[t]],
                names=('filename','object','imagetyp','exptime',
                    'xoffset','yoffset','time')).copy()[:0]
        row = [file,objname,imagetyp,exp,xof,yof,t]
        table.add_row(row)

    return(table)

# Returns a list of alternating on and off files for a given object name
def group_on_off(table, obj, sigoffset=2, maxoffset=15.):

    # Get the list of objects from the table
    objtable = table[table['object']==obj]

    # Sort table by time
    objtable.sort('time')

    # Check the min and max xoff and yoff values.
    min_xoff = np.min(abs(objtable['xoffset']))
    max_xoff = np.max(abs(objtable['xoffset']))
    min_yoff = np.min(abs(objtable['yoffset']))
    max_yoff = np.max(abs(objtable['yoffset']))

    xstd = sigoffset*np.std(objtable['xoffset'])
    ystd = sigoffset*np.std(objtable['yoffset'])

    # Step through the table and create lists of on/off rows
    groups=[]
    newgroup=[]
    typ = 'on'
    for row in objtable:
        if (((abs(min_xoff-row['xoffset'])<xstd and
            abs(min_yoff-row['yoffset'])<ystd)) or
            (abs(min_xoff-row['xoffset'])<maxoffset and
            abs(min_yoff-row['yoffset'])<maxoffset)):

            if typ=='on':
                newgroup.append(row['filename'])
            else:
                typ='on'
                if newgroup:
                    groups.append(('off',newgroup))
                newgroup = [row['filename']]

        else:
            if typ=='off':
                newgroup.append(row['filename'])
            else:
                typ='off'
                if newgroup:
                    groups.append(('on',newgroup))
                newgroup = [row['filename']]

    # Append the last group
    if newgroup:
        groups.append((typ,newgroup))

    return(groups)

# Converts an input list or table of fits files into input for pyraf
def make_inp(input_list, suffix='', add=''):
    input_list = list(input_list)

    add_str=''
    if add: add_str=add+','
    inp = add_str+','.join([suffix+f.replace('.fits','') for f in input_list])
    return(inp)

def unpack_files():
    files = glob.glob('*bz2')

    for file in files:
        cmd = 'bzip2 -d {0}'.format(file)
        os.system(cmd)
        unpacked = file.replace('.bz2','')
        if os.path.exists(unpacked):
            message = 'Successfully unpacked: {0}'.format(unpacked)
            print(message)

def sanitize_objname(obj):
    obj = obj.replace('-obs','')

    if obj.upper().startswith('AT'):
        obj = obj[2:]
    if obj.upper().startswith('SN'):
        obj = obj[2:]

    obj = obj.lower()
    return(obj)

def main():

    if not os.path.exists('raw/'):
        os.makedirs('raw/')
        for file in glob.glob('S2*.fits'):
            shutil.copyfile(file, 'raw/'+file)

    cleanup()
    copy_raw()

    os.system('mkiraf -f')
    from pyraf import iraf
    iraf.gemini(_doprint=0)
    iraf.gsaoi(_doprint=0)

    # Check if we need to unarchive and grab files
    if len(glob.glob('*.bz2'))>0:
        unpack_files()

        files = glob.glob('S2*.fits')
        if not os.path.exists('raw/'):
            os.makedirs('raw/')

        for file in files:
            shutil.copyfile(file, 'raw/'+file)


    files = glob.glob('S2*.fits')
    if len(files)==0:
        print('WARNING: no fits files to reduce!!!')
        print('Exiting...')
        sys.exit()

    allfiles = sort_files(files)

    if len(glob.glob('*_flat.fits'))==0:
        mask = allfiles['imagetyp']=='FLAT'
        iraf.gaflat(make_inp(allfiles[mask]['filename']),
            fl_vardq=True, fl_dqprop=True, use_off='yes')
    flatimg = glob.glob('*_flat.fits')[0]

    # Try parsing objects into on and off groups
    objtable = allfiles[allfiles['imagetyp']=='OBJECT']
    objs = np.unique(objtable['object'])
    objs = [obj for obj in objs if ('91' not in obj)]

    for obj in objs:

        groups = group_on_off(allfiles, obj, sigoffset=2, maxoffset=15.)

        # Assume we're alternating on/off
        # Check how many pairs we have
        num_group_pairs = len(groups)

        if num_group_pairs==1:
            # Assume not grouped together by on-off pairs and all just on source
            all_files = groups[0][1]
            reduce_input = make_inp(all_files)
            print(reduce_input)
            iraf.gareduce(reduce_input, fl_vardq=True, fl_dark=False,
                flatimg=flatimg, fl_sky=False, fl_dqprop=True, fl_flat=True)

            all_obj = make_inp(all_files, suffix='rg')
            print(all_obj)

        elif num_group_pairs>1:

            all_obj=''
            # Iterate through list by groups of 2
            for i,group in enumerate(groups):
                if group[0]=='off':
                    continue

                on_group = group
                if i==len(groups)-1:
                    off_group = groups[i-1]
                else:
                    off_group = groups[i+1]

                # Sky reduction
                iraf.gaprepare(make_inp(off_group[1]), fl_vardq=True)
                iraf.gasky(make_inp(off_group[1], suffix='g'),
                    outimage='sky{0}.fits'.format(i), fl_vardq=True,
                    fl_dqprop=True, flat=flatimg)

                # Now object reduction
                iraf.gareduce(make_inp(on_group[1]), fl_vardq=True,
                    fl_dark=False, flatimg=flatimg, fl_sky=True, fl_dqprop=True,
                    fl_flat=True, skyimg='sky{0}.fits'.format(i))

                all_obj = make_inp(on_group[1], suffix='rg', add=all_obj)

        # Run disco_stu on mosaic files
        all_files = [f+'.fits' for f in all_obj.split(',')]

        # Clean the amplifiers for all files
        print('Cleaning images amplifier by amplifier')
        for i,file in enumerate(all_files):
            print('file {0}/{1}:'.format(i+1,len(all_files)),file)
            clean_amps(file)

        # Now run disco_stu on everything
        run_disco_stu(' '.join(all_files))

        # Get re-projected images
        proj = glob.glob('*_proj.fits')

        if len(proj)>0:
            objname = sanitize_objname(obj)
            stack_files(proj, objname)
        else:
            print('WARNING: no projected files to stack!')

        # Clean up remaining files
        for file in glob.glob('tmp*'):
            os.remove(file)
        for file in glob.glob('*.log'):
            os.remove(file)
        shutil.rmtree('uparm')
        shutil.rmtree('pyraf')

if __name__ == '__main__':
    main()
