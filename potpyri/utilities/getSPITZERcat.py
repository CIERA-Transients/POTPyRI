#!/usr/bin/env python

import sys,os,re,math,optparse,warnings
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
import numpy as np
import astropy.units as u
from astropy import table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Column
warnings.filterwarnings('ignore')

def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)

    return(0)

def makepath4file(filename,raiseError=1):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

def rmfile(filename,raiseError=1,gzip=False):
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    if gzip and os.path.lexists(filename+'.gz'):
        os.remove(filename+'.gz')
        if os.path.isfile(filename+'.gz'):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename+'.gz')
            else:
                return(2)
    return(0)

def rmfiles(filenames,raiseError=1,gzip=False):
    if not (type(filenames) is list):
        raise RuntimeError('List type expected as input to rmfiles')
    errorflag = 0
    for filename in filenames:
        errorflag |= rmfile(filename,raiseError=raiseError,gzip=gzip)
    return(errorflag)

class getSPITZERcatclass:
    def __init__(self):
        self.coord = None

        self.alldata=None
        self.filename = None
        self.filterlist=[]
        self.colorfilterlist=[]
        self.errorcolsFlag= True
        self.Nonestring = '--'
        self.inst=''

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage,
                conflict_handler='resolve')

        parser.add_option('-v', '--verbose', action='count',
            dest='verbose',default=0)
        parser.add_option('-d','--debug', default=False, action='store_true',
            help='Debugging output')
        parser.add_option('--cathtml'  ,
            default='https://irsa.ipac.caltech.edu/SCS?table=slphotdr4',
            type='string', help='root html for SPITZER catalog')
        parser.add_option('--tmpdir'  , default='/tmp', type='string',
            help='Directory for temporary files (default=%default)')
        parser.add_option('-L',  default=False, action='store_true',
            help='u band')
        parser.add_option('-M',  default=False, action='store_true',
            help='g band')
        parser.add_option('-N',  default=False, action='store_true',
            help='r band')
        parser.add_option('-Q',  default=False, action='store_true',
            help='i band')
        parser.add_option('-s','--show',  default=False, action='store_true',
            help='print catalog to screen')
        parser.add_option('--keepnans',  default=False, action='store_true',
            help='keep entries with NaN values')
        parser.add_option('--requiregrizy',  default=False, action='store_true',
            help='only keep objects that have uvgriz measurements')
        parser.add_option('-o','--outfile'  , default=None , type='string',
            help='file name of output file. If not specified, then automatic '+\
            'filename depending on RA,Dec and filters (default=%default)')
        parser.add_option('-c','--clobber',  default=False, action='store_true',
            help='overwrite file if it already exists')
        parser.add_option('--size'  , default=None , type='string',
            help='sidelength of box in degree. If format AxB, then rectangle '+\
            'with sidelength A,B along Ra, Dec, respectively(default=%default)')
        parser.add_option('--skipsave',  default=False, action='store_true',
            help='Don\'t save the catalog')
        parser.add_option('--Mmax'  , default=28.0 , type='float',
            help='cut all objects with M<Mmin for specified filters '+\
            '(default=%default)')
        parser.add_option('--Mmin'  , default=5.0 , type='float',
            help='cut all objects with M>Mmax for specified '+\
            'filters (default=%default)')
        parser.add_option('--dMmax'  , default=1.0 , type='float',
            help='cut all objects with dM<dMmin for specified filters '+\
            '(default=%default)')
        parser.add_option('--dMmin', default=0.0 , type='float',
            help='cut all objects with dM>dMmax for specified filters '+\
            '(default=%default)')
        parser.add_option('--swope'  , default=False , action='store_true',
            help='convert the PS1 mag to SWOPE(CSP)')
        parser.add_option('--transfdir'  , default=False , type='string',
            help='path to the PS1transformations dir (default=%default)')
        parser.add_option('--skip_sc'  , default=False ,action='store_true',
            help='Skip the supercal correction to PS1 (default=%default)')
        parser.add_option('--spline'  , default=False ,action='store_true',
            help='Use spline fit (default=linear)')
        parser.add_option('--sexagesimal'  , default=False ,action='store_true',
            help='Output coordinates in sexagesimal format.'),

        return(parser)

    def check4TRANSFdir(self):

        if self.options.transfdir: self.path2name=self.options.transfdir
        else:
            if 'PIPE_PS1TRANSFDIR' in os.environ:
                if os.path.isdir(os.getenv('PIPE_PS1TRANSFDIR')): self.path2name=os.getenv('PIPE_PS1TRANSFDIR')+'/PS1trans_linear.txt'
                else:
                    print('WARNING: %s does not exist, looking into %s' %(os.getenv('PIPE_PS1TRANSFDIR'),os.getcwd()+'/PS1transformations'))
                    if os.path.isdir('./PS1transformations'):
                        self.path2name='./PS1transformations/PS1trans_linear.txt'
                    else: sys.exit('ERROR: Unable to find PS1transformations dir using the current path: %s. Try again with --transfdir option specifing the new path for the dir' % os.getcwd())
            else:
                print('WARNING: setenv PIPE_PS1TRANSFDIR does not exist, looking into %s' %os.getcwd()+'/PS1transformations')
                if os.path.isdir('./PS1transformations'): self.path2name='./PS1transformations/PS1trans_linear.txt'
                else: sys.exit('ERROR:  Unable to find PS1transformations dir using the current path: %s. Try again  with --transfdir option specifing the new path for the dir' % os.getcwd())


    def getfilterlist_from_options(self):

        if not self.inst:
            if self.options.L: self.filterlist.append('L')
            if self.options.M: self.filterlist.append('M')
            if self.options.N: self.filterlist.append('N')
            if self.options.Q: self.filterlist.append('Q')

        if self.filterlist == []:
            self.filterlist = ['L','M','N','Q']

    def parse_coord(self, ra, dec):
        if (':' in ra and ':' in dec):
            # Input RA/DEC are sexagesimal
            return(SkyCoord(ra, dec, frame='icrs'))
        elif (self.is_number(ra) and self.is_number(dec)):
            # Assume input coordiantes are decimal degrees
            return(SkyCoord(ra, dec, frame='icrs',unit='deg'))
        else:
            # Throw an error and exit
            error = 'ERROR: Cannot parse coordinates ra={ra}, dec={dec}'
            print(error.format(ra=ra,dec=dec))
            return(None)

    def is_number(self, num):
        try:
            num = float(num)
        except ValueError:
            return(False)
        return(True)

    def fluxcolname(self,filt):
        if (filt == 'L'): return('L')
        if (filt == 'M'): return('M')
        if (filt == 'N'): return('N')
        if (filt == 'Q'): return('Q')
        return('L')

    def errcolname(self,filt):
        if (filt == 'L'): return('L_err')
        if (filt == 'M'): return('M_err')
        if (filt == 'N'): return('N_err')
        if (filt == 'Q'): return('Q_err')
        return('e_L')

    def add_magnitude(self, tbdata, fluxcol, magcol, errcol, magerrcol):

        # Check if fluxcol and errcol in tbdata
        if fluxcol not in tbdata.keys() or errcol not in tbdata.keys():
            return(tbdata)

        # Do mag calculation
        mag = Column(np.array(-2.5 * np.log10(tbdata[fluxcol]*1e-6/3631.)),
            name=magcol)
        magerr = Column(np.array(1./1.086 * tbdata[errcol] / tbdata[fluxcol]),
            name=magerrcol)
        tbdata.add_column(mag)
        tbdata.add_column(magerr)
        return(tbdata)

    def getdataforRADEC(self, ra, dec, alldata=None,
        ramin=None, ramax=None, decmin=None, decmax=None, size=None):

        if alldata==None:
            alldata=self.alldata

        if ra<0.0:ra+=360.0
        if ra>=360.0:ra-=360.0

        if ramin!=None:
            if (ra-ramin)<-180:
                ramin-=360.0
                ramax-=360.0
            elif (ra-ramin)>180:
                ramin+=360.0
                ramax+=360.0

        radius = 0.5
        if size!=None:
            radius = np.sqrt(2)/2.0*size

        if self.options.debug:
            print('DEBUG: skipping download!')
        else:
            print('Downloading data')
            import pyvo
            Mmax=21.0
            if (self.options.Mmax!=None):
                Mmax = self.options.Mmax
            objects = pyvo.conesearch(self.options.cathtml, pos=(ra,dec), radius=size)
            tbdata = objects.table

            tbdata = self.add_magnitude(tbdata,
                'i1_f_ap1', 'L', 'i1_df_ap1', 'L_err')
            tbdata = self.add_magnitude(tbdata,
                'i2_f_ap1', 'M', 'i2_df_ap1', 'M_err')
            tbdata = self.add_magnitude(tbdata,
                'i3_f_ap1', 'N', 'i3_df_ap1', 'N_err')
            tbdata = self.add_magnitude(tbdata,
                'i4_f_ap1', 'Q', 'i4_df_ap1', 'Q_err')

        if not tbdata:
            print('Data download failed.  Exiting...')
            return(1)

        # calculate mags
        if self.options.verbose>1:
            print('calculating mags...')

        if ramin!=None:
            mask = ((tbdata['ra']<ramax) & (tbdata['ra']>ramin) &
                (tbdata['dec']<decmax) & (tbdata['dec']>decmin))
            data2keep = tbdata[mask]
        else:
            data2keep = tbdata

        print(data2keep)
        if self.options.requiregrizy:
            mask = ((data2keep[self.fluxcolname('L')] > 0) &
                    (data2keep[self.fluxcolname('M')] > 0) &
                    (data2keep[self.fluxcolname('N')] > 0) &
                    (data2keep[self.fluxcolname('Q')] > 0))
            data2keep = data2keep[mask]

        # check if there are limits on the mags and uncertainties
        for filt in self.filterlist:
            if self.options.Mmin!=None:
                mask = (data2keep[self.fluxcolname(filt)]>=self.options.Mmin)
                data2keep = data2keep[mask]
            if self.options.Mmax!=None:
                mask = (data2keep[self.fluxcolname(filt)]<=self.options.Mmax)
                data2keep = data2keep[mask]
            if self.options.dMmin!=None:
                mask = (data2keep[self.errcolname(filt)]>=self.options.dMmin)
                data2keep = data2keep[mask]
            if self.options.dMmax!=None:
                mask = (data2keep[self.errcolname(filt)]<=self.options.dMmax)
                data2keep = data2keep[mask]

        if alldata is None:
            self.alldata = table.Table(data2keep, masked=True)
        else:
            alldata = vstack(alldata, data2keep)
            self.alldata = table.unique(alldata, keys=['ra', 'dec'])

        return(0)

    def getcatalog(self, switch):
        ra = self.coord.ra.degree
        dec = self.coord.dec.degree
        if self.options.verbose>1:
            print('RA: %.7f, Dec:%.7f' % (ra, dec))
        RAboxsize = DECboxsize = None
        if self.options.size!=None:
            # get the boxsize, and ra/dec limits
            if re.search('\w+x\w+',self.options.size):
                m=re.search('(^\S+)x(\S+$)',self.options.size)
                if m!=None:
                    RAboxsize = float(m.groups()[0])
                    DECboxsize = float(m.groups()[1])
                    print('box sizes: ',RAboxsize,DECboxsize)
                else:
                    error = 'Could not parse size {size}.'
                    raise RuntimeError(error.format(self.options.size))
            else:
                RAboxsize = DECboxsize = float(self.options.size)

            # get the maximum 1.0/cos(DEC) term: used for RA cut
            minDec = self.coord.dec.degree-0.5*DECboxsize
            if minDec<=-90.0:minDec=-89.9
            maxDec = self.coord.dec.degree+0.5*DECboxsize
            if maxDec>=90.0:maxDec=89.9

            invcosdec = max(1.0/math.cos(self.coord.dec.degree*math.pi/180.0),
                            1.0/math.cos(minDec  *math.pi/180.0),
                            1.0/math.cos(maxDec  *math.pi/180.0))

            # get the (conservative) boxlimits
            ramin = self.coord.ra.degree-0.5*RAboxsize*invcosdec
            ramax = self.coord.ra.degree+0.5*RAboxsize*invcosdec
            decmin = self.coord.dec.degree-0.5*DECboxsize
            decmax = self.coord.dec.degree+0.5*DECboxsize
            # check for 4 corners...
            radeclist =[(self.coord.ra.degree, self.coord.dec.degree)]
            size = max(RAboxsize,DECboxsize)
        else:
            ramin = ramax = decmin = decmax = None
            radeclist = [(self.coord.ra.degree, self.coord.dec.degree)]
            size = 0.6

        # Which filters?
        self.getfilterlist_from_options()
        if self.options.verbose>1:
            message = 'Filters: {filt}.'
            print(message.format(filt=self.filterlist))

        for (ra_cat,dec_cat) in radeclist:
            message = '### checking {ra} {dec}.'
            print(message.format(ra=ra_cat, dec=dec_cat))
            (errorflag) = self.getdataforRADEC(ra_cat, dec_cat,
                ramin=ramin, ramax=ramax, decmin=decmin, decmax=decmax,
                size=size)

        return(0)

    def get_output_format_info(self,list_of_filters):
        # get cols ...
        cols=['ra','dec']
        colsformat = ['%11.7f','%11.7f']
        header = '#%11s %11s' % ('ra','dec')
        for filt in list_of_filters:
            filt = filt.strip('_psf')
            colsformat.append('%7.4f')
            cols.append(self.fluxcolname(filt))
            header += ' %7s' % filt
            if self.errorcolsFlag:
                colsformat.append('%7.4f')
                cols.append(self.errcolname(filt))
                header += ' %7s' % 'd'+filt

        return(cols,colsformat,header)

    def getformattedcatalog(self, cols, colsformat, alldata=None,
        indices = None, Nonestring=None, addreturnFlag=False):
        if alldata==None:
            alldata=self.alldata

        if Nonestring==None:
            Nonestring = self.Nonestring

        if self.options.sexagesimal:
            message = 'Converting to sexagesimal...'
            coords = SkyCoord(alldata['ra'], alldata['dec'], unit='deg')
            ra_hms = [c.to_string(style='hmsdms', sep=':').split()[0]
                 for c in coords]
            dec_dms = [c.to_string(style='hmsdms', sep=':').split()[1]
                 for c in coords]
            alldata['ra'] = ra_hms
            alldata['dec'] = dec_dms


        alldata = alldata.filled()
        lines = [' '.join([str(data[c]) for c in cols])+'\n' for data in alldata]
        return(lines)


    def savecatalog(self, outfilename, keys = None, cols = None,
        errorcolsFlag=True, Nonestring='--',clobber=False):
        makepath4file(outfilename)

        # file already exists?
        if os.path.isfile(outfilename):
            if clobber:
                if self.options.verbose>0:
                    warning = 'file {file} exists, clobbering it...'
                    print(warning.format(file=outfilename))
                rmfile(outfilename)
            else:
                warning = 'file {file} already exists, stopping '
                warning += '(if you want to overwrite file, use --clobber)'
                print(warning.format(file=outfilename))
                return(1)

        if self.options.verbose>1:
            print('Getting formatted catalog...')

        (cols,colsformat,header) = self.get_output_format_info(self.filterlist)

        (lines) = self.getformattedcatalog(cols, colsformat,
            Nonestring=Nonestring, addreturnFlag=True)

        if self.options.verbose:
            print('Saving %d entries into %s' % (len(lines),outfilename))
        f=open(outfilename,'w')
        f.writelines([header+'\n'])
        f.writelines(lines)
        f.close()

        # keep track of filename
        self.filename = outfilename

        if not os.path.isfile(outfilename):
            raise RuntimeError('Could not write to %s' % outfilename)

        return(0)

    def showcatalog(self, keys = None, cols = None, errorcolsFlag=True,
        Nonestring='--'):
        (cols,colsformat,header) = self.get_output_format_info()
        (lines) = self.getformattedcatalog(cols, colsformat,
            Nonestring=Nonestring)
        print(header)
        for line in lines: print(line)
        return(0)

if __name__=='__main__':

    spitzer=getSPITZERcatclass()
    parser = spitzer.add_options(usage='getSPITZERcat.py RA Dec')

    if len(sys.argv)<3:
        options,args = parser.parse_args(args=['-h'])
        sys.exit(0)

    # this is a hack to allow negative numbers as arguments.
    # this means that all options must be after the arguments.
    spitzer.options,dummy = parser.parse_args(args=sys.argv[3:])
    args = sys.argv[1:3]
    if len(dummy)!=0:
        print('ERROR: extra arguments!',dummy)
        sys.exit(0)

    # Keep track of RA,DEC
    (ra,dec)=args[:2]
    spitzer.coord = spitzer.parse_coord(ra, dec)

    spitzer.check4TRANSFdir()

    # outfilename
    if spitzer.options.outfile!=None:
        outfilename = spitzer.options.outfile
    else:
        outfilename = spitzer.autooutfilename()

    # get the catalog from the web
    spitzer.getcatalog(switch=0)

    # show the catalog?
    if spitzer.options.show:
        spitzer.showcatalog()

    # save the catalog?
    if not spitzer.options.skipsave:
        spitzer.savecatalog(outfilename,clobber=spitzer.options.clobber)

    print('SUCCESS %s' % os.path.basename(sys.argv[0]))
