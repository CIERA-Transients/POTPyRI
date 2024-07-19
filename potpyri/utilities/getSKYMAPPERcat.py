#!/usr/bin/env python

from astropy.io import fits
import sys,os,re,types,string,math,random
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import optparse
import numpy as np
import warnings
warnings.filterwarnings("ignore")

#from SMtoINST import SMtoINSTclass
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

def deg2sex(degrees, ra=False, outputformatRA='%02d:%02d:%06.3f',outputformatDEC='%1s%02d:%02d:%05.2f'):
    if type(degrees) is  bytes:
        degrees=float(degrees)
    if ra:
        # a.k.a. indeg and outhours
        if degrees < 0.0:
            while degrees<0:degrees+=360.
        if degrees > 360.0:
            while degrees>360.0:degrees-=360.
        degrees /= 15.

    if degrees < 0:
        sign = '-'
    else:
        sign = '+'

    degrees = abs(degrees)

    d1  = (degrees - (degrees % 1))
    rem = (degrees % 1) * 60
    d2  = rem - (rem % 1)
    srem = (rem % 1) * 60
    d3 = srem

    if ra:
      return outputformatRA % (d1, d2, d3)
    else:
      return outputformatDEC % (sign, d1, d2, d3)

### Converts sexigesimal notation to decimal degrees or decimal hours (if option 'ra=True' invoked)
def sex2deg(sexigecimal, ra=False):
    ### 2005/12/02 - AR: make sure it is in sexagesimal format!
    # is it a string? if not check if it is None
    if not (type(sexigecimal) is bytes):
        if type(sexigecimal) == None:
            raise RuntimeError("ERROR: sex2deg cannot handle 'None'!")
        return sexigecimal
    # Does it have ':' or ' '? If not, it must be a float in string format, just convert it and return
    if re.search('\:|\s',sexigecimal) == None:
        return(string.atof(sexigecimal))

    s1, s2, s3 = list(map(string.atof, re.split('[DHMSdhms:\s]', sexigecimal.strip())[0:3]))

    # Get the sign
    if re.search('-', sexigecimal):
        sign = -1
    else:
        sign = 1

    deg = abs(s1) + s2 / 60. + s3 / 3600.

    deg *= sign

    if ra:
        deg *= 15.

    return deg

# Returns the passed in RA in decimal degrees
# input RA can be in 'HH:MM:SS.ss', 'HH MM SS.ss' or in decimal degrees
def RaInDeg(Ra):
    import types
    if type(Ra)==bytes:
        if re.search('[DHMShms: ]',Ra.strip()):
            return(sex2deg(Ra,ra=True))
    return(float(Ra))

# Returns the passed in Dec in decimal degrees
# input Dec can be in 'DD:MM:SS.ss', 'DD MM SS.ss' or in decimal degrees
def DecInDeg(Dec):
    import types
    if type(Dec)==bytes:
        if re.search('[DHMShms: ]',Dec.strip()):
            return(sex2deg(Dec,ra=False))
    return(float(Dec))

def rmfile(filename,raiseError=1,gzip=False):
    " if file exists, remove it "
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
        raise RuntimeError("List type expected as input to rmfiles")
    errorflag = 0
    for filename in filenames:
        errorflag |= rmfile(filename,raiseError=raiseError,gzip=gzip)
    return(errorflag)

class getSKYMAPPERcatclass:
    def __init__(self):
        self.ra = None
        self.dec = None

        self.alldata={}
        self.filename = None
        self.filterlist=[]
        self.colorfilterlist=[]
        self.errorcolsFlag= True
        self.Nonestring = 'nan'
        self.inst=''


    def autooutfilename(self):
        if self.ra==None or self.dec==None:
            raise RuntimeError('No RA/Dec defined to be used for auto filename!')

        filename = '%011.7f_%011.7f.SKYMAPPER.' % (self.ra,self.dec)
        filename += ''.join(self.filterlist)
        filename += '.cat'

        return(filename)

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_option('-d','--debug', default=False, action="store_true",
                          help='Debugging output')
        parser.add_option('--cathtml'  , default='http://api.skymapper.nci.org.au/public/tap/', type='string',
                          help='root html for SKYMAPPER catalog')
        parser.add_option('--tmpdir'  , default='/tmp', type="string",
                          help='Directory for temporary files (default=%default)')
        parser.add_option('-u',  default=False, action="store_true",
                          help='u band')
        parser.add_option('-v',  default=False, action="store_true",
                          help='v band')
        parser.add_option('-g',  default=False, action="store_true",
                          help='g band')
        parser.add_option('-r',  default=False, action="store_true",
                          help='r band')
        parser.add_option('-i',  default=False, action="store_true",
                          help='i band')
        parser.add_option('-z',  default=False, action="store_true",
                          help='z band')
        parser.add_option('-U',  default=False, action="store_true",
                          help='U band')
        parser.add_option('-B',  default=False, action="store_true",
                          help='B band')
        parser.add_option('-V',  default=False, action="store_true",
                          help='V band')
        parser.add_option('-R',  default=False, action="store_true",
                          help='R band')
        parser.add_option('-I',  default=False, action="store_true",
            help='I band')
        parser.add_option('--sexagesimal',  default=False, action="store_true",
            help='output RA/Dec are in sexagesimal format')
        parser.add_option('-s','--show',  default=False, action="store_true",
            help='print catalog to screen')
        parser.add_option('--keepnans',  default=False, action="store_true",
            help='keep entries with NaN values')
        parser.add_option('--requiregriz',  default=False, action="store_true",
            help='only keep objects that have uvgriz measurements')
        parser.add_option('-o','--outfile'  , default=None , type="string",
            help='file name of output file. If not specified, then automatic filename depending on RA,Dec and filters (default=%default)')
        parser.add_option('-c','--clobber',  default=False, action="store_true",
            help='overwrite file if it already exists')
        parser.add_option('-s','--size'  , default=None , type="string",
            help='sidelength of box in degree. If format AxB, then rectangle with sidelength A,B along Ra, Dec, respectively(default=%default)')
        parser.add_option('--skipsave',  default=False, action="store_true",
            help='Don\'t save the catalog')
        parser.add_option('--swope',  default=False, action="store_true",
            help='Add Swope corrections to photometry.')
        parser.add_option('--spline',  default=False, action="store_true",
            help='Use spline fit model for transformations.')
        parser.add_option('--skip_sc',  default=False, action="store_true",
            help='Skip supercal corrections.')
        parser.add_option('--Mmax'  , default=None , type="float",
            help='cut all objects with M<Mmin for specified filters (default=%default)')
        parser.add_option('--Mmin'  , default=None , type="float",
                          help='cut all objects with M>Mmax for specified filters (default=%default)')
        parser.add_option('--dMmax'  , default=None , type="float",
                          help='cut all objects with dM<dMmin for specified filters (default=%default)')
        parser.add_option('--dMmin'  , default=None , type="float",
                          help='cut all objects with dM>dMmax for specified filters (default=%default)')
        parser.add_option('--transform'  , default=False , action="store_true",
                          help='convert the SKYMAPPER mag to BVRI')

        return(parser)

    def getfilterlist_from_options(self):

        if not self.inst:
            if self.options.U: self.filterlist.append('U')
            if self.options.u: self.filterlist.append('u')
            if self.options.v: self.filterlist.append('v')
            if self.options.g: self.filterlist.append('g')
            if self.options.r: self.filterlist.append('r')
            if self.options.i: self.filterlist.append('i')
            if self.options.z: self.filterlist.append('z')
            if self.options.B: self.filterlist.append('B')
            if self.options.V: self.filterlist.append('V')
            if self.options.R: self.filterlist.append('R')
            if self.options.I: self.filterlist.append('I')

        if self.filterlist == []:
            self.filterlist = ['u_psf','v_psf','g_psf','r_psf','i_psf','z_psf']

    def add2list(self, alldata, col1, newdata, col2):
        Nrowsnew = len(newdata[col2])
        if col1 in alldata:
            Nrows1 = len(alldata[col1])
            if Nrowsnew>0:
                tmp=np.zeros((Nrows1+Nrowsnew))
                tmp[:Nrows1]=alldata[col1]
                tmp[Nrows1:]=newdata[col2]
                alldata[col1]=tmp
        else:
            if Nrowsnew>0:
                alldata[col1]=newdata[col2]
        return(0)

    def addnewfilt(self,alldata,filt):
        u=alldata[self.fluxcolname('u')]
        g=alldata[self.fluxcolname('g')]
        r=alldata[self.fluxcolname('r')]
        i=alldata[self.fluxcolname('i')]
        e_u=alldata[self.errcolname('u')]
        e_g=alldata[self.errcolname('g')]
        e_r=alldata[self.errcolname('r')]
        e_i=alldata[self.errcolname('i')]

        # Transformations in Jester et al. (2005)
        if filt == 'U':
            alldata['U']=0.78*u-0.78*g-0.88+g+0.39*(g-r)+0.21
            alldata['e_U']=np.sqrt((0.78*e_u)**2+(0.61*e_g)**2+(0.39*e_r)**2)
        if filt == 'V':
            alldata['V']=g-0.587*(g-r)-0.011
            alldata['e_V']=np.sqrt((0.41*e_g)**2+(0.59*e_r)**2)
        if filt == 'B':
            alldata['B']=g+0.327*(g-r)+0.216
            alldata['e_B']=np.sqrt((1.39*e_g)**2+(0.39*e_r)**2)
        if filt == 'R':
            alldata['R']=r-0.272*(r-i)-0.159
            alldata['e_R']=np.sqrt((1.09*e_i)**2+(0.50*e_r)**2+(0.41*e_g)**2)
        if filt == 'I':
            alldata['I']=i-0.337*(r-i)-0.370
            alldata['e_I']=np.sqrt((2.09*e_i)**2+(0.41*e_g)**2+(1.50*e_r)**2)

        return(0)

    def fluxcolname(self,filt):
        if (filt == 'U'): return('U')
        if (filt == 'u'): return('u_psf')
        if (filt == 'g'): return('g_psf')
        if (filt == 'r'): return('r_psf')
        if (filt == 'i'): return('i_psf')
        if (filt == 'z'): return('z_psf')
        if (filt == 'B'): return('B')
        if (filt == 'V'): return('V')
        if (filt == 'R'): return('R')
        if (filt == 'I'): return('I')
        return('r_psf')

    def errcolname(self,filt):
        if (filt == 'u'): return('e_u_psf')
        if (filt == 'g'): return('e_g_psf')
        if (filt == 'r'): return('e_r_psf')
        if (filt == 'i'): return('e_i_psf')
        if (filt == 'z'): return('e_z_psf')
        if (filt == 'U'): return('e_U')
        if (filt == 'B'): return('e_B')
        if (filt == 'V'): return('e_V')
        if (filt == 'R'): return('e_R')
        if (filt == 'I'): return('e_I')
        return('e_r_psf')

    def getdataforRADEC(self,ra,dec,alldata=None,ramin=None,
        ramax=None,decmin=None,decmax=None):

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

        if self.options.debug:
            print('DEBUG: skipping download!')
        else:
            print('Downloading data')
            if (self.options.Mmax!=None):
                Mmax = self.options.Mmax
            else:
                Mmax = 30.0
            vquery = Vizier(columns=['RAICRS', 'DEICRS',
                             'uPSF', 'e_uPSF',
                             'vPSF', 'e_vPSF',
                             'gPSF', 'e_gPSF',
                             'rPSF', 'e_rPSF',
                             'iPSF', 'e_iPSF',
                             'zPSF', 'e_zPSF'],
                    column_filters={'gmag':
                                    ('<%f' % Mmax)},
                    row_limit=100000)

            coord = SkyCoord(ra, dec, unit='deg', frame='icrs')
            w = ('%fd' % radius)
            cat = 'II/358/smss'
            tbdata = vquery.query_region(coord, width=w, catalog=cat)[0]
            tbdata.rename_column('RAICRS', 'raj2000')
            tbdata.rename_column('DEICRS', 'dej2000')
            tbdata.rename_column('uPSF', 'u_psf')
            tbdata.rename_column('e_uPSF', 'e_u_psf')
            tbdata.rename_column('vPSF', 'v_psf')
            tbdata.rename_column('e_vPSF', 'e_v_psf')
            tbdata.rename_column('gPSF', 'g_psf')
            tbdata.rename_column('e_gPSF', 'e_g_psf')
            tbdata.rename_column('rPSF', 'r_psf')
            tbdata.rename_column('e_rPSF', 'e_r_psf')
            tbdata.rename_column('iPSF', 'i_psf')
            tbdata.rename_column('e_iPSF', 'e_i_psf')
            tbdata.rename_column('zPSF', 'z_psf')
            tbdata.rename_column('e_zPSF', 'e_z_psf')

        if not tbdata:
            print('Data download failed.  Exiting...')
            return(1)

        # calculate mags
        if self.options.verbose>1:
            print('calculating mags...')

        if ramin!=None:
            mask = ((tbdata['raj2000']<ramax) & (tbdata['raj2000']>ramin) &
                (tbdata['dej2000']<decmax) & (tbdata['dej2000']>decmin))
            data2keep = tbdata[mask]
        else:
            data2keep = tbdata
        if self.options.requiregriz:
            mask =  ((data2keep[self.fluxcolname('u')] != '--') &
                    (data2keep[self.fluxcolname('v')] != '--') &
                    (data2keep[self.fluxcolname('g')] != '--') &
                    (data2keep[self.fluxcolname('r')] != '--') &
                    (data2keep[self.fluxcolname('i')] != '--') &
                    (data2keep[self.fluxcolname('z')] != '--'))
            data2keep = data2keep[mask]
            for elno in range(len(data2keep)):
                print((data2keep[elno][0],data2keep[elno][1],
                       data2keep[elno][2],data2keep[elno][3],
                       data2keep[elno][4],data2keep[elno][5]))

        # check if there are limits on the mags and uncertainties
        for filt in self.filterlist:
            if self.options.Mmin!=None:
                mask =  (data2keep[self.fluxcolname(filt)]>=self.options.Mmin)
                data2keep = data2keep[mask]
            if self.options.Mmax!=None:
                mask =  (data2keep[self.fluxcolname(filt)]<=self.options.Mmax)
                data2keep = data2keep[mask]
            if self.options.dMmin!=None:
                mask =  (data2keep[self.errcolname(filt)]>=self.options.dMmin)
                data2keep = data2keep[mask]
            if self.options.dMmax!=None:
                mask =  (data2keep[self.errcolname(filt)]<=self.options.dMmax)
                data2keep = data2keep[mask]

         # now add the data to self.alldata
        self.add2list(alldata,'raj2000',data2keep,'raj2000')
        self.add2list(alldata,'dej2000',data2keep,'dej2000')

        if self.options.transform:
            for filt in self.filterlist:
                self.addnewfilt(data2keep,filt)

        for filt in self.filterlist:
            self.add2list(alldata,filt,data2keep,self.fluxcolname(filt))
            self.add2list(alldata,'d'+filt,data2keep,self.errcolname(filt))

        del tbdata,data2keep#,color2keep

        return(0)

    def getcatalog(self, ra, dec, switch):
        # keep track of RA,Dec
        self.ra  = RaInDeg(ra)
        self.dec = DecInDeg(dec)
        if self.options.verbose>1:
            print('RA: %.7f, Dec:%.7f' % (self.ra,self.dec))
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
                    raise RuntimeError('Could not parse size %s for' % self.options.size)
            else:
                RAboxsize = DECboxsize = float(self.options.size)

            # get the maximum 1.0/cos(DEC) term: used for RA cut
            minDec = self.dec-0.5*DECboxsize
            if minDec<=-90.0:minDec=-89.9
            maxDec = self.dec+0.5*DECboxsize
            if maxDec>=90.0:maxDec=89.9

            invcosdec = max(1.0/math.cos(self.dec*math.pi/180.0),
                            1.0/math.cos(minDec  *math.pi/180.0),
                            1.0/math.cos(maxDec  *math.pi/180.0))

            # get the (conservative) boxlimits
            ramin = self.ra-0.5*RAboxsize*invcosdec
            ramax = self.ra+0.5*RAboxsize*invcosdec
            decmin = self.dec-0.5*DECboxsize
            decmax = self.dec+0.5*DECboxsize
            # check for url for center and all 4 corners...
            radeclist =[(self.ra,self.dec),
                        (ramin,decmin),
                        (ramin,decmax),
                        (ramax,decmin),
                        (ramax,decmax)]


            if ramax-ramin > 1 and decmax-decmin > 1:
                radeclist=[]
                for n in range(int(ramax-ramin)+1):
                    app_ra=ramin+n
                    for l in range(int(decmax-decmin)+1):
                        app_dec=decmin+l
                        radeclist.extend([(app_ra,app_dec)])

            elif ramax-ramin > 1 and decmax-decmin < 1:
                radeclist=[]
                for n in range(int(ramax-ramin)+1):
                    app_ra=ramin+n
                    for dec in [decmin,decmax]:
                        app_dec=dec
                        radeclist.extend([(app_ra,app_dec)])

            elif ramax-ramin < 1 and decmax-decmin > 1:
                radeclist=[]
                for ra in [ramin,ramax]:
                    app_ra=ra
                    for l in range(int(decmax-decmin)+1):
                        app_dec=decmin+l
                        radeclist.extend([(app_ra,app_dec)])


        else:
            ramin = ramax = decmin = decmax = None
            radeclist = [(self.ra,self.dec)]

        # Which filters?
        self.getfilterlist_from_options()
        if self.options.verbose>1:
            print('Filters:',self.filterlist)

        for (ra_cat,dec_cat) in radeclist:
            print('### checking %f %f'% (ra_cat,dec_cat))
            (errorflag) = self.getdataforRADEC(ra_cat,dec_cat,ramin=ramin,
                ramax=ramax,decmin=decmin,decmax=decmax)

        return(0)

    def get_output_format_info(self,list_of_filters):
        # get cols ...
        cols=['raj2000','dej2000']
        colsformat = ['%11.7f','%11.7f']
        header = '#%11s %11s' % ('ra','dec')
        for filt in list_of_filters:
            colsformat.append('%7.4f')
            cols.append(filt)
            header += ' %7s' % (filt)
            if self.errorcolsFlag:
                colsformat.append('%7.4f')
                cols.append('d'+filt)
                header += ' %7s' % ('d'+filt)

        return(cols,colsformat,header)

    def getformattedcatalog(self,cols,colsformat, alldata=None, indices = None,
        Nonestring=None, addreturnFlag=False):
        if alldata==None:
            alldata=self.alldata

        if indices == None:
            indices = range(len(alldata['raj2000']))

        if Nonestring==None:
            Nonestring = self.Nonestring

        cs = range(len(cols))

        lines=[]
        for i in indices:
            line = ''
            nanflag = False
            for c in cs:
                #print len(alldata[cols[c]]),i
                if math.isnan(alldata[cols[c]][i]) or math.isinf(alldata[cols[c]][i]):
                    nanflag=True
                    line += ' %7s' %  Nonestring
                else:
                    if cols[c] in ['raj2000','dej2000'] and self.options.sexagesimal:
                        line += ' '+deg2sex(alldata[cols[c]][i],
                            ra=(cols[c]=='raj2000'))
                    else:
                        line += ' '+colsformat[c] % alldata[cols[c]][i]

            # skip bad lines if not --keepnans
            if nanflag and (not self.options.keepnans):
                continue

            if addreturnFlag:
                line += '\n'

            lines.append(line)

        return(lines)


    def savecatalog(self, outfilename, keys = None, cols = None,
        errorcolsFlag=True, Nonestring='NaN',clobber=False):
        makepath4file(outfilename)

        # file already exists?
        if os.path.isfile(outfilename):
            if clobber:
                if self.options.verbose>0:
                    print('file %s exists, clobbering it...' % outfilename)
                rmfile(outfilename)
            else:
                print('file %s already exists, stopping (if you want to overwrite file, use --clobber)' % outfilename)
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

    def showcatalog(self, keys = None, cols = None, errorcolsFlag=True, Nonestring='NaN'):
        (cols,colsformat,header) = self.get_output_format_info()
        (lines) = self.getformattedcatalog(cols, colsformat, Nonestring=Nonestring)
        print(header)
        for line in lines: print(line)
        return(0)

    def converter(self,outfilename):
        selected_band=[]
        if self.options.B:
            selected_band.append('B')
        elif self.options.U:
            selected_band.append('U')
        elif self.options.V:
            selected_band.append('V')
        elif self.options.R:
            selected_band.append('R')
        elif self.options.I:
            selected_band.append('I')
        elif self.options.u:
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
        else:
            selected_band.append('r')
            selected_band.extend(self.filterlist)

if __name__=='__main__':

    getSKYMAPPERcat=getSKYMAPPERcatclass()
    parser = getSKYMAPPERcat.add_options(usage='getSKYMAPPERcat.py RA Dec')

    if len(sys.argv)<3:
        options,args = parser.parse_args(args=['-h'])
        sys.exit(0)

    # this is a hack to allow negative numbers as arguments.
    # this means that all options must be after the arguments.
    getSKYMAPPERcat.options,dummy = parser.parse_args(args=sys.argv[3:])
    args = sys.argv[1:3]
    if len(dummy)!=0:
        print('ERROR: extra arguments!',dummy)
        sys.exit(0)

    (ra,dec)=args[:2]

    # outfilename
    if getSKYMAPPERcat.options.outfile!=None:
        outfilename = getSKYMAPPERcat.options.outfile
    else:
        outfilename = getSKYMAPPERcat.autooutfilename()

    # get the catalog from the web
    getSKYMAPPERcat.getcatalog(ra,dec,switch=0)

    # show the catalog?
    if getSKYMAPPERcat.options.show:
        getPS1cat.showcatalog()

    # save the catalog?
    if not getSKYMAPPERcat.options.skipsave:
        getSKYMAPPERcat.savecatalog(outfilename,clobber=getSKYMAPPERcat.options.clobber)

    print('SUCCESS %s' % os.path.basename(sys.argv[0]))
