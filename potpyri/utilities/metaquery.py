import numpy as np
import warnings,sys
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io.votable import parse
import os,requests,sys,io,shutil,time,gzip,ftplib,copy
warnings.filterwarnings('ignore')

instrument_defaults = {
    'standard_keys': ['obsid','ra','dec','instrument',
                      'filter','exptime','jpg','image'],
    'standard_units': [None, u.degree, u.degree, None,
                       None, u.second, None, None],
    'hst': {
        'radius': 1 * u.arcsec,
        'mask': {'instrument_name': ['WFPC2/WFC','PC/WFC','ACS/WFC','ACS/HRC',
                                     'ACS/SBC','WFC3/UVIS','WFC3/IR'],
                 't_exptime': 40,
                 'obs_collection': ['HST'],
                 'filters': ['F220W','F250W','F330W','F344N','F435W','F475W',
                      'F550M','F555W','F606W','F625W','F658N','F660N','F660N',
                      'F775W','F814W','F850LP','F892N','F098M','F105W','F110W',
                      'F125W','F126N','F127M','F128N','F130N','F132N','F139M',
                      'F140W','F153M','F160W','F164N','F167N','F200LP','F218W',
                      'F225W','F275W','F280N','F300X','F336W','F343N','F350LP',
                      'F373N','F390M','F390W','F395N','F410M','F438W','F467M',
                      'F469N','F475X','F487N','F502N','F547M','F600LP','F621M',
                      'F625W','F631N','F645N','F656N','F657N','F658N','F665N',
                      'F673N','F680N','F689M','F763M','F845M','F953N','F122M',
                      'F160BW','F185W','F218W','F255W','F300W','F375N','F380W',
                      'F390N','F437N','F439W','F450W','F569W','F588N','F622W',
                      'F631N','F673N','F675W','F702W','F785LP','F791W','F953N',
                      'F1042M','F502N']
                },
        'keymap': [['s_ra','ra'],
                   ['s_dec','dec'],
                   ['instrument_name','instrument'],
                   ['filters','filter'],
                   ['t_exptime','exptime'],
                   ['dataURL','image'],
                   ['jpegURL','jpg'],
                   ['t_min','date']],
    },
    'chandra': {
        'url': 'https://cxcfps.cfa.harvard.edu/cgi-bin/cda/footprint/'+\
                'get_vo_table.pl?pos={ra},{dec}&size={radius}',
        'options': '&inst=ACIS-I&grating=NONE',
        'jpeg_url_base': 'https://cda.harvard.edu/chaser/viewerImage.do?',
        'jpeg_url_options': 'obsid={obsid}&filename=acisf{fobsid}'+\
                            'N004_cntr_img2.jpg&filetype=hiresimg_jpg',
        'radius': 12 * u.arcmin,
        'keymap': [['ObsId','obsid'],
                   ['RA','ra'],
                   ['Dec','dec'],
                   ['Instrument','instrument'],
                   ['Grating','filter'],
                   ['Exposure','exptime'],
                   ['obs_date','date']],
        'mask': {'Instrument': ['ACIS-I','ACIS-S'],
                 'Grating': ['NONE'],
                 'Exposure': 20.0}
    },
    'spitzer': {
        'url': 'http://sha.ipac.caltech.edu/applications/'+\
               'Spitzer/SHA/servlet/DataService?',
        'radius': 500 * u.arcsec,
        'keymap': [['reqkey','obsid'],
                   ['accessUrl','image'],
                   ['begintime','date']],
        'mask': {'instrument': ['IRAC'],
                 'externalname': ['maic.fits']}
    }
}

class metaquery():
    def __init__(self, ra, dec):
        # Assume ':' means formatted hms, dms
        # Otherwise, assume decimal degrees
        if (':' in str(ra) and ':' in str(dec)):
            self.coord = SkyCoord(ra, dec, unit = (u.hour, u.deg))
        else:
            self.coord = SkyCoord(ra, dec, unit = (u.deg, u.deg))

        self.table = None
        self.options = instrument_defaults

    def map_dtypes(self, type_names, field_widths):
        dtypes = []
        for i, name in enumerate(type_names):
            if name == 'int':
                dtypes.append('i8')
            elif name == 'double':
                dtypes.append('f8')
            elif name == 'char':
                dtypes.append('a{0}'.format(field_widths[i]))
            else:
                raise ValueError('Unexpected type name: {0}.'.format(name))
        return dtypes

    def get_table(self, telescope):
        options = self.options[telescope]
        if (telescope.upper() == 'HST'):
            from astroquery.mast import Observations

            # Get table
            table = Observations.query_region(self.coord,
              radius = options['radius'])

            # HST-specific masks
            filmask = [table['filters'] == good
              for good in options['mask']['filters']]
            filmask = [any(l) for l in list(map(list,zip(*filmask)))]
            expmask = table['t_exptime'] > options['mask']['t_exptime']
            obsmask = [table['obs_collection'] == good
              for good in options['mask']['obs_collection']]
            obsmask = [any(l) for l in list(map(list,zip(*obsmask)))]
            detmask = [table['instrument_name'] == good
              for good in options['mask']['instrument_name']]
            detmask = [any(l) for l in list(map(list,zip(*detmask)))]

            # Construct and apply mask
            mask = [all(l) for l in zip(filmask,expmask,obsmask,detmask)]
            table = table[mask]

        if (telescope.upper() == 'SPITZER'):
          import requests
          from dateutil.parser import parse
          params = {'RA': self.coord.ra.degree, 'DEC': self.coord.dec.degree,
                    'DATASET': 'ivo://irsa.ipac.spitzer.level2',
                    'SIZE': options['radius'].to(u.degree).value,
                    'VERB': 3}

          url = options['url']

          r = requests.get(url, params = params)
          if 'No Matches' in r.text:
              return None

          if r.status_code!=200:
              message = 'status message: {message}.'
              print(message.format(message=r.text))
              error = 'ERROR: could not get url {url}, status code {stat}.'
              raise RuntimeError(error.format(url=url, stat=r.status_code))
              return None

          # Parse table from html text
          # This is borrowed from astroquery.sha
          raw_data = [line for line in r.text.split('\n')]
          field_widths = [len(s) + 1 for s in raw_data[1].split('|')][1:-1]
          col_names = [s.strip() for s in raw_data[1].split('|')][1:-1]
          type_names = [s.strip() for s in raw_data[2].split('|')][1:-1]
          cs = [0] + np.cumsum(field_widths).tolist()

          def parse_line(line, cs = cs):
              return [line[a:b] for a, b in zip(cs[:-1], cs[1:])]

          data = [parse_line(row) for row in raw_data[4:-1]]
          # Parse type names
          dtypes = self.map_dtypes(type_names, field_widths)
          # To table
          # transpose data for appropriate table instance handling
          table = Table(list(zip(*data)), names=col_names, dtype=dtypes)

          # Split wavelength into instrument and filter
          instrument = []
          filt = []
          for w in table['wavelength']:
            dat = w.split()
            inst = dat[0].strip()
            if 'IRS' in inst:
              f = dat[2]
            elif 'IRAC' in inst:
              f = dat[1]
            elif 'MIPS' in inst:
              if 'SED' in dat[1]:
                f = dat[2]
              else:
                f = dat[1]
            else:
              print('Unknown instrument: ',inst)
            instrument.append(inst)
            filt.append(f)

          table.add_column(astropy.table.Column(name = 'instrument',
            data = instrument))
          table.add_column(astropy.table.Column(name = 'filter',
            data = filt))

          detmask = [table['instrument'] == good
            for good in options['mask']['instrument']]
          detmask = [any(l) for l in list(map(list,zip(*detmask)))]
          promask = [[good in x for x in table['externalname']]
            for good in options['mask']['externalname']]
          promask = [any(l) for l in list(map(list,zip(*promask)))]

          # Construct and apply mask
          mask = [all(l) for l in zip(detmask,promask)]
          table = table[mask]

          table.add_column(astropy.table.Column(name = 'jpg',
            data = [None] * len(table)))
          # Estimate exposure time.  This is only defined for Level 1 products
          # since all other images are stacks.  But we can guesstimate the
          # exposure time by looking at the beginning and end time of the
          # observation
          exptime = [(parse(end) - parse(start)).seconds
              for start,end in zip(table['begintime'],table['endtime'])]
          table.add_column(astropy.table.Column(name = 'exptime',
            data = exptime))

        if (telescope.upper() == 'CHANDRA'):
            import requests
            from astropy.io.votable import parse
            url = 'https://cxcfps.cfa.harvard.edu/cgi-bin/'+\
              'cda/footprint/get_vo_table.pl'
            ra  = self.coord.ra.degree
            dec = self.coord.dec.degree
            params = {
                'pos': '{0},{1}'.format(ra, dec),
                'size': 0.2,
                'grating': 'NONE',
                'inst': 'ACIS-I,ACIS-S'
            }
            if True:
                # Try to get a response from the URL
                r = requests.get(url, params=params)

                # Create a temporary file and then read it into a votable object.
                # This is terrible, but best way I've found to complete this part
                fname = 'test.xml'
                f = open(fname, 'w')
                f.write(r.text)
                f.close()
                votable = parse(fname)
                os.remove(fname)
            #except:
            #    return(None)

            tbdata = votable.get_first_table().to_table()

            # Although using a large search radius, check to make sure that the
            # input coordinates are < 5 arcmin off axis
            remove = []
            for i,row in enumerate(copy.copy(tbdata)):
                coord = SkyCoord(row['RA'],row['Dec'],
                    unit=(u.deg, u.deg), frame='icrs')
                if self.coord.separation(coord).arcmin > 5:
                    remove.append(i)
            tbdata.remove_rows(remove)

            # Require that all observations have SI Mode = TE (timed exposure)
            # rather than CC = continuous clocking (and thus not imaging)
            remove = []
            for i,row in enumerate(copy.copy(tbdata)):
                if not row['SIMode'].decode('utf-8').startswith('TE'):
                    remove.append(i)
            tbdata.remove_rows(remove)

            table = tbdata

        # Rename column names to standard names
        for keypair in options['keymap']:
            oldname = keypair[0]
            newname = keypair[1]
            table.rename_column(oldname, newname)

        # Reorder table and make sure only standard keys are used
        return table

if __name__=='__main__':
    if (len(sys.argv) < 2):
        print("Usage: python metaquery.py RA Dec")
        sys.exit(1)

    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    query = metaquery(ra,dec)
    #chandra = query.get_table('chandra')
    hst = query.get_table('hst')
    #spitzer = query.get_table('spitzer')
    print(len(hst))
