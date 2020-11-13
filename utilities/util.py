from astropy.coordinates import SkyCoord
from astropy import utils, units as u
from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
from astropy.table import Table, Column
from astropy.io.votable import parse as votable_parse
from astropy.io import ascii, fits
from astroquery.mast import Observations, Catalogs
from astroquery.vizier import Vizier
import os
import pandas
import requests
import numpy as np
import shutil
import casjobs
import math
import copy
import requests
import warnings
from contextlib import contextmanager

warnings.filterwarnings('ignore')

viziercat = {
    'sdssdr12': {'name':'V/147',
        'columns': ['RAJ2000', 'DEJ2000','class','umag', 'e_umag',
            'gmag','e_gmag', 'rmag','e_rmag', 'imag','i_mag', 'zmag',
            'e_zmag', 'zph']
    },
    '2mass': {'name':'II/246',
        'columns':['RAJ2000', 'DEJ2000', 'Jmag','e_Jmag','Hmag','e_Hmag',
            'Kmag', 'e_Kmag']
    },
    'unwise': {'name':'II/363',
        'columns': ['RAJ2000', 'DEJ2000', 'FW1','e_FW1', 'FW2','e_FW2']
    },
    'glade': {'name':'VII/281',
        'columns': ['RAJ2000', 'DEJ2000', 'Dist', 'e_Dist', 'Bmag', 'Jmag',
            'Hmag', 'Kmag', 'z']
    },
    'des': {'name':'II/357',
        'columns': ['RAJ2000', 'DEJ2000', 'S/Gg', 'S/Gr', 'S/Gi', 'S/Gz',
            'gmag','e_gmag', 'rmag','e_rmag', 'imag','e_imag', 'zmag','e_zmag']
    }
}

@contextmanager
def suppress_stdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in ra and ':' in dec):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)


# Get ramin, ramax, demin, demax corners of a box centered at coord with size
def get_ra_dec_box(coord, size):
    ra = coord.ra.degree
    dec = coord.dec.degree

    RAboxsize = DECboxsize = float(size)

    # get the maximum 1.0/cos(DEC) term: used for RA cut
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

    return(ramin, ramax, demin, demax)

def download_url_to_file(url, filename):

    outdir, file = os.path.split(filename)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.path.exists(filename):
        return(0)
    dat = utils.data.download_file(url, cache=False, show_progress=True)
    shutil.move(dat, filename)
    return(0)

def downloadPS1(coord,  filt, outdir='data/', tempdir='data/delme/'):

    # Check for filter
    if not filt: filt='r'

    obs = Observations.query_region(coord, radius=5*u.arcmin)

    masks = []
    masks.append([str(t).upper()=='PS1' for t in obs['obs_collection']])
    masks.append([str(p).upper()=='IMAGE' for p in obs['dataproduct_type']])
    masks.append([str(r).upper()=='PUBLIC' for r in obs['dataRights']])
    masks.append([str(f).lower()==filt for f in obs['filters']])

    masks = list(map(list, zip(*masks)))
    mask = np.array([all(l) for l in masks])

    obsTable = obs[mask]

    if len(obsTable)==0:
        error = 'ERROR: no matching PS1 images for ra={ra}, dec={dec}.'
        print(error.format(ra=ra, dec=dec))
        return(1)
    else:
        prodList = Observations.get_product_list(obsTable)
        prodList = prodList[prodList['productType']=='SCIENCE']
        prodList = prodList[prodList['description']=='stack data image']

        message = 'We need to download {n} PS1 images.'
        print(message.format(n=len(prodList)))

        for product in prodList:
            filename = product['productFilename']
            outfile = filename.replace('.fz','')
            outname = outfile.replace('.fits', '_mod.fits')
            url = 'https://mast.stsci.edu/api/v0/download/file?uri='
            url += product['dataURI']
            if os.path.isfile(tempdir+outname):
                message = '{filename} exists.  Continuing...'
                print(message.format(filename=tempdir+outname))
            else:
                outfile = tempdir+filename+'.fz'
                errflag = download_url_to_file(url, outfile)

                hdu  = fits.open(outfile)

                # Need to reweight the data for PS1 by effective zpt
                exptime = hdu[1].header['EXPTIME']
                boffset = hdu[1].header['BOFFSET']
                bsoften = hdu[1].header['BSOFTEN']
                a = 1.0857362

                # Get the image data from the file
                img  = hdu[1].data
                mask = np.isnan(hdu[1].data)

                # Adjust data values for asinh compression
                data  = boffset + bsoften * (np.exp(img/a) - np.exp(-img/a))
                zpt   = 25.0 + 2.5 * np.log10(exptime)
                data  = data * 10**(0.4*(27.5-zpt))
                hdu[1].data = data

                # Reset the mask values
                hdu[1].data[mask] = np.nan

                # Write out modified file
                outname = outfile.replace('.fits', '_mod.fits')
                outname = outname.replace('.fz','')

                outname = os.path.basename(outname)
                fulloutname = outdir + outname
                hdu.writeto(fulloutname, overwrite=True)
                print('Wrote out: {0}'.format(outname))
                return fulloutname

def docasjobsstrm(coord, size=0.1, mask=True, verbose=False,
    meta=['raMean','decMean','z_phot','z_photErr']):

        ra = coord.ra.degree ; dec = coord.dec.degree
        ramin, ramax, demin, demax = get_ra_dec_box(coord, size)

        query = 'select '
        meta_query = ', '.join(['c.' + val for val in meta])
        query += meta_query
        query += ' from catalogRecordRowStore c '
        query += ' where '
        query += ' c.raMean between {ramin} and {ramax}'.format(ramin=ramin,
            ramax=ramax)
        query += ' and c.decMean between {demin} and {demax}'.format(demin=demin,
            demax=demax)

        if verbose:
            print('running query:')
            print(query)

        jobs = casjobs.CasJobs(userid='892987546', password='BossTent1',
            base_url='http://mastweb.stsci.edu/ps1casjobs/services/jobs.asmx')

        job_output = jobs.quick(query, context='HLSP_PS1_STRM',
            task_name='PS1cat_ra%.7f_dec%.7f'%(ra,dec))

        names = meta

        table = ascii.read(job_output, names=names)

        # Add an offset key to the table
        coords = SkyCoord(table['raMean'].data, table['decMean'].data, unit='deg')
        sep = coord.separation(coords).arcsec
        col = Column(sep, name='separation')
        table.add_column(col)

        if mask and 'z_phot' in table.keys():
            # Mask values with missing photo-z and anomalously negative values
            table = table[table['z_phot']!=-999.0]
            table = table[table['z_phot']>-0.01]

        return(table)

def get_grb_data(grb_name, verbose=True):
    name = copy.copy(grb_name)
    grb_name = grb_name.replace('GRB','').strip().upper()
    r = requests.get('https://swift.gsfc.nasa.gov/archive/grb_table/fullview/')

    if verbose: print('Successfully downloaded Swift GRB data...')

    if grb_name not in r.text:
        print('WARNING: {0} not in Swift GRB data!'.format(name))
        return(None)

    # This is slowest step - maybe a way to speed up?
    df = pandas.read_html(r.content)

    table = None
    header = [key[0] for key in df[0].columns]

    # Sanitize header for astropy
    astrohead = []
    for key in header:
        key = key.split('[')[0]
        key = key.replace('90%','')
        key = key.split('(')[0]
        key = key.strip()
        astrohead.append(key)

    for i in np.arange(df[0].shape[0]):
        if df[0]['GRB'].iloc[i][-1]!=grb_name: continue
        if verbose: print('Got GRB metadata')
        row = []
        for key in header:
            d = df[0][key].iloc[i][-1]
            # Parse nan values
            if str(d)=='nan': d = ''
            row.append(d)
        # Sanitize
        table = Table([[r] for r in row], names=astrohead)
        table.rename_column('TriggerNumber','Trigger Number')
        break

    table[0]['BAT RA']  = table[0]['BAT RA'][-10:]
    table[0]['BAT Dec'] = table[0]['BAT Dec'][-11:]

    if verbose:
        for key in ['Time','Trigger Number','BAT T90','XRT RA','XRT Dec']:
            fmt='{0}: {1}'
            print(fmt.format(key, table[0][key]))

    return(table[0])

# Function for the handlinng the downloaded data text function
def readData(String):
    file = os.path.abspath(String)
    # Open the text file
    infile = open(file, 'r')
    data = infile.read().split(',' and '\n')
    infile.close()
    lst =data[1:-1]

    return(lst)

# Function to sort through our table now that we have it in a list
def sorting(lst):
    s = ''
    no_punct = ""
    string = s.join(lst)
    punctuation = '()"'
    sorted_list =[]
    # Split the big list into a nested list
    for i in range(len(lst)):
        sorted_list.append(lst[i].split(','))

    # Combine everything in the nested list into a list of strings
    string_list = []
    for i in sorted_list:
        for j in i:
            string = s.join(j)
            string_list.append(string)

    # Take away the punctuations in the list so we don't have to worrry
    #   about them later on in the data
    big_list = []
    for words in string_list:
        for char in words:
            if char in punctuation:
                words = words.replace(char, '')
        big_list.append(words)

    # Big list is the list that contains all the values without the punctuations


    # You can decide how many variables you want to obtain from the catalogs,
    # although the current variables should be comprehensive for Prospector use
    # or any sorts of plotting.
    num_of_variable = 11

    # Split those numbers into nested lists and store them in thet 'nest' list
    #   and then we will beready to analyze the actual data!
    count = 1
    small_list = []
    nest = []
    for i in big_list:
        small_list.append(i)
        if count % num_of_variable == 0:
            #small_list.append(distance[(count/num_of_variable)-1])
            nest.append(small_list)
            small_list = []
        count += 1
    return(nest)

def data(nest):
# Find all galaxy-only photometry and store that in a list
    num_of_variable = 11
    galaxy_list = []
    for i in range(len(nest)):
        for j in range(len(nest[i])):
            if nest[i][j] == ' 3':
                galaxy_list.append(nest[i])

    # We want to do a quick analysis on the photometric redshift so we will
    # make a list of just redshifts here
    phot_z_list = []
    for i in range(len(galaxy_list)):
        phot_z_list.append(galaxy_list[i][(num_of_variable - 1)])


    # Getting rid of all the observations that have no z_phot (nan and -9999.0)
    c_list = []
    galaxy_list_org = copy.deepcopy(galaxy_list)

    for i in range(len(galaxy_list_org)):
        if phot_z_list[i] != ' nan':
            c_list.append(galaxy_list[i])

    c_list_org = copy.deepcopy(c_list)
    for i in c_list_org:
        if i[(num_of_variable-1)] == ' -9999.0':
            c_list.remove(i)

    # EVERTHING WE NEED IS NOW IN THIS INFO_LIST IN NUMBERS
    info_list = []
    for i in range(len(c_list)):
        for j in range(len(c_list[i])):
            c_list[i][j].strip()
            info_list.append(float(c_list[i][j].strip()))

    count = 1
    small = []
    final = []
    for i in info_list:
        small.append(i)
        if count % num_of_variable == 0:
            final.append(small)
            small = []
        count += 1
    return final

def ordering(final):
    #Make a r_mag list
    r_mag = []
    for i in range(len(final)):
        for j in range(len(final[i])):
            if j == 7:
                r_mag.append(final[i][7])

    #Make a distance list
    dist = []
    for i in range(len(final)):
        for j in range(len(final[i])):
            if j == 0:
                dist.append(final[i][0])

    # Here, we want to sort through our photometry
    # You can sort the photometry either based on r_mag or distance. Just simply
    # Change dist to r_mag or r_mag to dist and taa daa!

    # BubbleSort
    n = len(dist)

    # Traverse through all array elements
    for i in range(n-1):
    # range(n) also work but outer loop will repeat one time more than needed.

        # Last i elements are already in place
        for j in range(0, n-i-1):
            # traverse the array from 0 to n-i-1
            # Swap if the element found is greater
            # than the next element
            if dist[j] > dist[j+1]:
                dist[j], dist[j+1] = dist[j+1], dist[j]
                final[j], final[j+1] = final[j+1], final[j]
    for i in final:
        print(i)

    #print(len(final))


    # Write the ordered photometry as a text file using pickle so they can be in floats and not integers
    with open(object_name+'.txt', 'wb') as fp:
        pickle.dump(final,fp)

    # Use this to open the pickled text file
    with open(object_name+'.txt', 'rb') as fp:
        b = pickle.load(fp)
    #print(b)

    # You can also write your result using this method since pickle does not allow you to read it in a text editor
    #np.savetxt('200522A_result_np.txt', np.c_[final])

    # Then you can use this to unpack the text file :)
    #final = np.loadtxt('200522A_result_np.txt',unpack=True)

def getSDSS(coord, impix=1024, imsize=12*u.arcmin, outfile=''):
    from urllib.parse import urlencode
    from urllib.request import urlretrieve

    scale = imsize.to(u.arcsec).value/impix
    cutoutbaseurl = 'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx'
    query_string = urlencode(dict(ra=coord.ra.deg, dec=coord.dec.deg,
        width=impix, height=impix, scale=scale))

    url = '?'.join([cutoutbaseurl, query_string])

    outdir, basename = os.path.split(outfile)
    if outdir and not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            warning = 'WARNING: could not make outdir={0}'
            warning += '\nSetting outdir to current directory'
            print(warning.format(outdir))
            outdir = ''

    fulloutfile = outfile

    # this downloads the image to your disk
    urlretrieve(url, fulloutfile)

    if os.path.exists(fulloutfile):
        return(fulloutfile)
    else:
        error = 'ERROR: could not generate file={0}'
        print(error.format(fulloutfile))
        return('')

def onclick(event,source_list):
    if event.dblclick:
        ra_click, dec_click = event.xdata, event.ydata
        print('RA = %.3f, Dec = %.3f' %(ra_click, dec_click))

def onpick(event,source_list,gaia_list,std_list):
    label = event.artist.get_label()
    print(label)

def get_coord_from_table(table1, cat1ra='ra', cat1dec='dec', cat1coord=False):

    coords1=None
    if cat1coord:
        if not cat1coord in table1.keys():
            return(None)
        try:
            coords1 = SkyCoord(table1[cat1coord])
        except:
            print('ERROR: problem parsing coords in table1')
            return(None)
    elif cat1ra not in table1.keys() or cat1dec not in table1.keys():
        print('ERROR: could not get coordinates from table1')
        return(None)
    else:
        try:
            unit=(u.hour, u.deg)
            if is_number(table1[cat1ra][0]) and is_number(table1[cat1dec][0]):
                unit=(u.deg, u.deg)
            coords1 = SkyCoord(table1[cat1ra], table1[cat1dec], unit=unit)
        except:
            print('ERROR: problem parsing coords in table1')
            return(None)

    return(coords1)

def make_blank_row(table):

    row = []
    for key in table.keys():
        if 'str' in table[key].dtype.name:
            row.append('')
        elif 'float' or 'int' in table[key].dtype.name:
            row.append(-999)
        else:
            row.append(None)

    return(row)

def search_catalogs(coord, catnames, search_radius=2.564*u.arcmin,
    match_radius=2.0*u.arcsec, outfile=''):

    if outfile and os.path.exists(outfile):
        m = '{0} exists.  Do you want to use this catalog?'
        m = m.format(outfile) + '([y]/n): '
        y = input(m)
        if not y or y=='y' or y=='yes':
            outtable = ascii.read(outfile)
            return(outtable)
        else:
            print('Redoing search...')

    outtable = None
    for name in catnames:

        table = None
        print('Searching {0} catalog...'.format(name))
        if name in viziercat.keys():
            table = searchVizier(coord, name, radius=search_radius)
            table_match_ra = 'RAJ2000'
            table_match_dec = 'DEJ2000'
        elif name=='ps1dr2':
            table = Catalogs.query_region(coord, catalog='Panstarrs',
                radius=search_radius.to_value('degree'),
                data_release='dr2', table='mean')
            table_match_ra = 'raMean'
            table_match_dec = 'decMean'
        elif name=='strm':
            table = docasjobsstrm(coord, size=search_radius.to_value('degree'))
            table_match_ra = 'raMean'
            table_match_dec = 'decMean'

        if table and len(table)>0:
            print('{0} records in {1} catalog'.format(len(table), name))
            if not outtable:
                outtable = table
                for key in table.keys(): table.rename_column(key, name+'_'+key)
            else:
                out_match_ra = ''
                out_match_dec = ''
                if 'ra' in outtable.keys() and 'dec' in outtable.keys():
                    out_match_ra = 'ra'
                    out_match_dec = 'dec'
                elif any(['ps1dr2' in key for key in outtable.keys()]):
                    out_match_ra = 'ps1dr2_raMean'
                    out_match_dec = 'ps1dr2_decMean'
                elif any(['strm' in key for key in outtable.keys()]):
                    out_match_ra = 'strm_raMean'
                    out_match_dec = 'strm_decMean'
                else:
                    for key in outtable.keys():
                        if 'raj2000' in key.lower() and not out_match_ra:
                            out_match_ra = key
                        if 'dej2000' in key.lower() and not out_match_dec:
                            out_match_dec = key

                outtable = crossmatch_tables(outtable,table,'',name,
                    radius=match_radius,
                    cat1ra=out_match_ra, cat1dec=out_match_dec,
                    cat2ra=table_match_ra, cat2dec=table_match_dec)

            n=len(outtable)
            print('{0} total records in catalog'.format(n))

        else:
            print('0 records in {0} catalog'.format(name))

    if not outfile:
        if not os.path.exists('data/'):
            os.makedirs('data/')
        # Definitely gonna want to write this out, so write out to outfile.cat
        ascii.write(outtable, 'data/outfile.cat')
    else:
        outdir, base = os.path.split(outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        ascii.write(outtable, outfile)

    return(outtable)


def searchVizier(coord, catname, radius=2.564*u.arcmin):

    catid = viziercat[catname]['name']
    columns = viziercat[catname]['columns']

    v = Vizier(columns=columns, catalog=[catname], row_limit = 1000000000)
    result = v.query_region(coord, radius=2.564*u.arcmin)

    if len(result)==1:
        return(result[0])
    else:
        return(None)

def crossmatch_tables(table1, table2, t1name, t2name, radius=2.0 * u.arcsec,
    cat1ra='ra', cat1dec='dec', cat2ra='ra', cat2dec='dec', cat1coord='',
    cat2coord=''):

    coords1 = get_coord_from_table(table1, cat1ra=cat1ra, cat1dec=cat1dec,
        cat1coord=cat1coord)
    coords2 = get_coord_from_table(table2, cat1ra=cat2ra, cat1dec=cat2dec,
        cat1coord=cat2coord)

    if not coords1 or not coords2:
        print('ERROR: could not get coordinates from one of the tables')
        return(None)

    # Make a new table with all columns from table1 and table2
    if t1name: t1name+='_'
    if t2name: t2name+='_'
    header = [t1name+k for k in table1.keys()] +\
             [t2name+k for k in table2.keys()]
    # Make a dummy row with first row from both tables
    row = [[table1[0][k]] for k in table1.keys()] +\
          [[table2[0][k]] for k in table2.keys()]

    # This makes an empty table with the appropriate header and column types
    mtable = Table(row, names=header).copy()[:0]

    # Crossmatch all of the coordinates from coords1 and coords2
    idx, dist2d, dist3d = match_coordinates_sky(coords1, coords2)

    # Deal with duplicate matches
    close = idx[dist2d < radius]
    if len(close) > len(np.unique(close)):
        # TODO - fix this method.  Currently we just throw out all but the
        # closest for duplicate matches.  Ideally we'd like to keep track of
        # these and iteratively rematch if there is a better candidate (but
        # that would also mean keeping track of whether or not the rematch also
        # has duplicates, etc. until there are no duplicates)
        for i, dist in zip(idx, dist2d):
            matches = np.where(idx==i)[0]
            if len(matches)>1:
                distance = dist
                closest = i
                for j in matches:
                    if distance > dist2d[j]:
                        # Throw out this match by resetting the distance
                        dist2d[closest] = radius
                        closest = j
                    else:
                        dist2d[j] = radius

    t2copy = copy.copy(table2)
    # Combine table1 + table2 data if they match closer than radius and keep
    # track of table2 rows we've already checked off
    radat = [] ; decdat = []
    for j, val in enumerate(zip(idx, dist2d)):
        i, dist = val
        row1 = table1[j]
        if dist < radius:
            row2 = table2[i]
            ridx = np.where(row2[cat2ra]==t2copy[cat2ra])[0][0]
            radat.append(0.5*(row2[cat2ra]+row1[cat1ra]))
            decdat.append(0.5*(row2[cat2dec]+row1[cat1dec]))
            t2copy.remove_row(ridx)
            mtable.add_row(list(row1) + list(row2))
        else:
            radat.append(row1[cat1ra])
            decdat.append(row1[cat1dec])
            mtable.add_row(list(row1) + make_blank_row(table2))

    # Now handle remaining table2 rows
    if len(t2copy)>0:
        for i, row in enumerate(t2copy):
            radat.append(row[cat2ra])
            decdat.append(row[cat2dec])
            mtable.add_row(make_blank_row(table1) + list(row))

    racol = Column(radat, name='ra')
    deccol = Column(decdat, name='dec')
    if 'ra' in mtable.keys():
        mtable['ra'] = racol
    else:
        mtable.add_column(racol)
    if 'dec' in mtable.keys():
        mtable['dec'] = deccol
    else:
        mtable.add_column(deccol)

    # Final mtable should contain all information from table1 and table2
    return(mtable)










