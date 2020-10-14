#!/usr/bin/env python3

import sys,os,re,math,optparse,warnings
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
import numpy as np
import astropy.units as u
from astropy import table
from astropy.io import fits,ascii
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from PS1toINST import PS1toINSTclass
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

class getPS1DR2catclass:
    def __init__(self):
        self.coord = None

        self.alldata=None
        self.filename = None
        self.filterlist=[]
        self.colorfilterlist=[]
        self.errorcolsFlag= True
        self.Nonestring = '--'
        self.inst=''

        self.supercal={'g':0.020,'r':0.033,'i':0.024,'z':0.028,'y':0.011}

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage,
                conflict_handler='resolve')

        parser.add_option('-v', '--verbose', action='count',
            dest='verbose',default=0)
        parser.add_option('-d','--debug', default=False, action='store_true',
            help='Debugging output')
        parser.add_option('--url'  ,
            default='http://mastweb.stsci.edu/ps1casjobs/services/jobs.asmx',
            type='string', help='root html for SKYMAPPER catalog')
        parser.add_option('--tmpdir'  , default='/tmp', type='string',
            help='Directory for temporary files (default=%default)')
        parser.add_option('-u',  default=False, action='store_true',
            help='u band')
        parser.add_option('-g',  default=False, action='store_true',
            help='g band')
        parser.add_option('-r',  default=False, action='store_true',
            help='r band')
        parser.add_option('-i',  default=False, action='store_true',
            help='i band')
        parser.add_option('-z',  default=False, action='store_true',
            help='z band')
        parser.add_option('-y',  default=False, action='store_true',
            help='y band')
        parser.add_option('-U',  default=False, action='store_true',
            help='U band')
        parser.add_option('-B',  default=False, action='store_true',
            help='B band')
        parser.add_option('-V',  default=False, action='store_true',
            help='V band')
        parser.add_option('-R',  default=False, action='store_true',
            help='R band')
        parser.add_option('-I',  default=False, action='store_true',
            help='I band')
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
        parser.add_option('--size'  , default='0.5' , type='string',
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
        parser.add_option('--transform'  , default=False , action='store_true',
            help='convert the SKYMAPPER mag to BVRI')
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
                path = os.getenv('PIPE_PS1TRANSFDIR')
                if os.path.isdir(os.getenv('PIPE_PS1TRANSFDIR')):
                    linear_trans = 'PS1trans_linear.txt'
                    self.path2name=path + '/' + linear_trans
                else:
                    warning = 'WARNING: {path} does not exist, looking in {new}'
                    ps1path = os.getcwd()+'/PS1transformations'
                    warning = warning.format(path=path, new=ps1path)
                    print(warning)
                    if os.path.isdir('./PS1transformations'):
                        self.path2name='./PS1transformations/PS1trans_linear.txt'
                    else:
                        sys.exit('ERROR: Unable to find PS1transformations dir using the current path: %s. Try again with --transfdir option specifing the new path for the dir' % os.getcwd())
            else:
                print('WARNING: setenv PIPE_PS1TRANSFDIR does not exist, looking into %s' %os.getcwd()+'/PS1transformations')
                if os.path.isdir('./PS1transformations'): self.path2name='./PS1transformations/PS1trans_linear.txt'
                else: sys.exit('ERROR:  Unable to find PS1transformations dir using the current path: %s. Try again  with --transfdir option specifing the new path for the dir' % os.getcwd())


    def getfilterlist_from_options(self):
        if self.options.swope:
            self.inst='CSP1'

        if not self.inst:
            if self.options.u: self.filterlist.append('u')
            if self.options.g: self.filterlist.append('g')
            if self.options.r: self.filterlist.append('r')
            if self.options.i: self.filterlist.append('i')
            if self.options.z: self.filterlist.append('z')
            if self.options.y: self.filterlist.append('y')
            if self.options.U: self.filterlist.append('U')
            if self.options.B: self.filterlist.append('B')
            if self.options.V: self.filterlist.append('V')
            if self.options.R: self.filterlist.append('R')
            if self.options.I: self.filterlist.append('I')

        if self.filterlist == []:
            self.filterlist = ['g_psf','r_psf','i_psf','z_psf','y_psf']

    def addnewfilt(self, alldata, filt):
        g=alldata[self.fluxcolname('g')]
        r=alldata[self.fluxcolname('r')]
        i=alldata[self.fluxcolname('i')]
        e_g=alldata[self.errcolname('g')]
        e_r=alldata[self.errcolname('r')]
        e_i=alldata[self.errcolname('i')]

        # Transformations in Jester et al. (2005)
        # TODO: Update these with different color corrections
        if filt == 'B':
            alldata['B']=g+0.39*(g-r)+0.21
            alldata['e_B']=np.sqrt((1.39*e_g)**2+(0.39*e_r)**2)
        if filt == 'V':
            alldata['V']=g-0.59*(g-r)-0.01
            alldata['e_V']=np.sqrt((0.41*e_g)**2+(0.59*e_r)**2)
        if filt == 'R':
            alldata['R']=1.09*i+0.41*g-0.50*r-0.23
            alldata['e_R']=np.sqrt((1.09*e_i)**2+(0.50*e_r)**2+(0.41*e_g)**2)
        if filt == 'I':
            alldata['I']=2.09*i+0.41*g-1.50*r-0.44
            alldata['e_I']=np.sqrt((2.09*e_i)**2+(0.41*e_g)**2+(1.50*e_r)**2)

        return(0)

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
        if (filt == 'u'): return('u')
        if (filt == 'g'): return('g')
        if (filt == 'r'): return('r')
        if (filt == 'i'): return('i')
        if (filt == 'z'): return('z')
        if (filt == 'y'): return('y')
        if (filt == 'U'): return('U')
        if (filt == 'B'): return('B')
        if (filt == 'V'): return('V')
        if (filt == 'R'): return('R')
        if (filt == 'I'): return('I')
        return('r_psf')

    def errcolname(self,filt):
        if (filt == 'u'): return('uerr')
        if (filt == 'g'): return('gerr')
        if (filt == 'r'): return('rerr')
        if (filt == 'i'): return('ierr')
        if (filt == 'z'): return('zerr')
        if (filt == 'y'): return('yerr')
        if (filt == 'U'): return('Uerr')
        if (filt == 'B'): return('Berr')
        if (filt == 'V'): return('Verr')
        if (filt == 'R'): return('Rerr')
        if (filt == 'I'): return('Ierr')
        return('e_r_psf')

    def docasjobsphot(self, ra, dec, meta=['nDetections','raMean','decMean'],
        filters=['g','r','i','z','y'], ramin=None, ramax=None, demin=None,
        demax=None, dMmax = 0.2):

        import casjobs

        query = 'select '

        meta_query = ', '.join(['o.' + val for val in meta])
        metacol = ['nDetections', 'ra', 'dec']
        phot = []
        photcol = []
        for val in filters:
            phot.append('m.{filt}MeanPSFMag'.format(filt=val))
            phot.append('m.{filt}MeanKronMag'.format(filt=val))
            phot.append('m.{filt}MeanPSFMagErr'.format(filt=val))
            photcol.append(val)
            photcol.append(val+'Kron')
            photcol.append(val+'err')
        phot_query = ', '.join(phot)

        RAboxsize = DECboxsize = float(self.options.size)

        # get the maximum 1.0/cos(DEC) term: used for RA cut
        if demin is None:
            demin = dec-0.5*DECboxsize
            if demin<=-90.0:demin=-89.9
            demax = dec+0.5*DECboxsize
            if demax>=90.0:demax=89.9

        invcosdec = max(1.0/math.cos(dec*math.pi/180.0),
                        1.0/math.cos(demin  *math.pi/180.0),
                        1.0/math.cos(demax  *math.pi/180.0))

        # get the (conservative) boxlimits
        ramin = ra-0.5*RAboxsize*invcosdec
        ramax = ra+0.5*RAboxsize*invcosdec
        demin = dec-0.5*DECboxsize
        demax = dec+0.5*DECboxsize

        query += meta_query + ', ' + phot_query
        query += ' from ObjectThin o '
        query += ' inner join MeanObject m on o.objID=m.objID '
        query += ' where '
        query += ' o.raMean between {ramin} and {ramax}'.format(ramin=ramin,
            ramax=ramax)
        query += ' and o.decMean between {demin} and {demax}'.format(demin=demin,
            demax=demax)
        query += ' and o.nDetections>2'
        query += ' and m.iMeanPSFMag - m.iMeanKronMag < 0.05'
        query += ' and m.iMeanPSFMagErr < {dMmax}'.format(dMmax=dMmax)

        if self.options.verbose:
            print('running quick query:')
            print(query)

        jobs = casjobs.CasJobs(userid='892987546', password='BossTent1',
            base_url=self.options.url)

        job_output = jobs.quick(query, context='PanSTARRS_DR2',
            task_name='PS1cat_ra%.7f_dec%.7f'%(ra,dec))

        names = metacol + photcol

        table = ascii.read(job_output, names=names)
        return(table)

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
            tbdata = self.docasjobsphot(ra, dec, ramin=ramin, ramax=ramax,
                demin=decmin, demax=decmax, dMmax = 0.2)

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
            mask = ((data2keep[self.fluxcolname('g')] > 0) &
                    (data2keep[self.fluxcolname('r')] > 0) &
                    (data2keep[self.fluxcolname('i')] > 0) &
                    (data2keep[self.fluxcolname('z')] > 0) &
                    (data2keep[self.fluxcolname('y')] > 0))
            data2keep = data2keep[mask]

        if not self.options.skip_sc:
            print('>>>>>> Adding supercal correction')
            for filt in self.filterlist:
                errname = self.errcolname(filt)
                filtname = self.fluxcolname(filt)
                data2keep[errname] = 1.06*data2keep[errname]/data2keep[filtname]
                data2keep[filtname] = -2.5*np.log10(data2keep[filtname])
                data2keep[filtname] = data2keep[filtname] + self.supercal[filt]

        # Add new rows if we need to
        if self.options.transform:
            for filt in self.filterlist:
                self.addnewfilt(data2keep, filt)

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
            print(message)
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

    def converter(self, outfilename):
        selected_band=[]
        if self.options.u:
            selected_band.append('u')
        elif self.options.g:
            selected_band.append('g')
        elif self.options.r:
            selected_band.append('r')
        elif self.options.i:
            selected_band.append('i')
        elif self.options.z:
            selected_band.append('z')
        elif self.options.y:
            selected_band.append('y')
        elif self.options.U:
            selected_band.append('U')
        elif self.options.B:
            selected_band.append('B')
        elif self.options.V:
            selected_band.append('V')
        elif self.options.R:
            selected_band.append('R')
        elif self.options.I:
            selected_band.append('I')
        else:
            selected_band.append('r')
            selected_band.extend(self.filterlist)

        PS1toINST=PS1toINSTclass()
        PS1toINST.convert(self.path2name,
            self.inst,outfilename,selected_band,
            False,self.options.spline)

if __name__=='__main__':

    ps1=getPS1DR2catclass()
    parser = ps1.add_options(usage='getPS1DR2cat.py RA Dec')

    if len(sys.argv)<3:
        options,args = parser.parse_args(args=['-h'])
        sys.exit(0)

    # this is a hack to allow negative numbers as arguments.
    # this means that all options must be after the arguments.
    ps1.options,dummy = parser.parse_args(args=sys.argv[3:])
    args = sys.argv[1:3]
    if len(dummy)!=0:
        print('ERROR: extra arguments!',dummy)
        sys.exit(0)

    # Keep track of RA,DEC
    (ra,dec)=args[:2]
    ps1.coord = ps1.parse_coord(ra, dec)

    ps1.check4TRANSFdir()

    # outfilename
    if ps1.options.outfile!=None:
        outfilename = ps1.options.outfile
    else:
        outfilename = ps1.autooutfilename()

    # get the catalog from the web
    ps1.getcatalog(switch=0)

    # show the catalog?
    if ps1.options.show:
        ps1.showcatalog()

    # save the catalog?
    if not ps1.options.skipsave:
        ps1.savecatalog(outfilename,clobber=ps1.options.clobber)

    if ps1.options.swope:
        ps1.converter(outfilename)

    print('SUCCESS %s' % os.path.basename(sys.argv[0]))
