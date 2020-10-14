import requests,os,sys,zipfile,shutil
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import utils
from astropy.table import Table
from astropy import units as u

uri = 'http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/DataService'

# Color strings for download messages
green = '\033[1;32;40m'
red = '\033[1;31;40m'
end = '\033[0;0m'

def post_request(coord, size=0.5):

    ra = coord.ra.degree
    dec = coord.dec.degree

    params = {'RA': ra, 'DEC': dec, 'SIZE': size, 'VERB': 3,
        'DATASET':'ivo%3A%2F%2Firsa.ipac%2Fspitzer.level1'}

    req = requests.get(uri, params=params)
    table = ascii.read(req.text)
    return(table)

# Input a table row from post_request, download the corresponding data product
# and unpack to target directory
def get_cbcd(row, outrootdir = './', tempdir = './', download_zip=True):

    url = row['accessWithAnc1Url']
    if 'NONE' in url or not url:
        return(1)
    unpack_file = tempdir + row['externalname'].replace('bcd.fits','cbcd.fits')
    base_file = os.path.basename(unpack_file)
    ref_file = outrootdir + base_file

    # This speeds up process if downloading multiple SNe where you might have
    # overlapping files from one object to next.  Also if you have to redo
    # the download process.
    if download_zip and os.path.exists(unpack_file):
        # Move downloaded file to ref_file
        print(unpack_file)
        shutil.copyfile(unpack_file, ref_file)
        return(0)

    message = 'Downloading file: {url}'

    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    if download_zip:
        download_file = tempdir + 'download.zip'
    else:
        download_file = ref_file

    message = 'Downloading file: {url}'
    sys.stdout.write(message.format(url=url))
    sys.stdout.flush()
    dat = utils.data.download_file(url, cache=False,
        show_progress=False, timeout=120)
    shutil.copyfile(dat, download_file)
    #os.chmod(download_file, 0775)
    message = '\r' + message
    message += green+' [SUCCESS]'+end+'\n'
    sys.stdout.write(message.format(url=url))

    # Unpack the file if zip
    if download_zip:
        # Unzip file
        zip_ref = zipfile.ZipFile(download_file, 'r')
        zip_ref.extractall(tempdir)
        zip_ref.close()

        # Move downloaded file to ref_file
        if os.path.exists(unpack_file):
            shutil.copyfile(unpack_file, ref_file)
        os.remove(download_file)

    return(0)

def download_from_file(file):
    for row in ascii.read(file):
        # Get SN specific info
        name = row['name']
        ra = row['ra']
        dec = row['dec']
        date = row['ut-date']
        mjd = Time(date).mjd
        coord = SkyCoord(ra, dec, unit='deg', frame='icrs')

        # Get the table of all Spitzer data
        table = post_request(coord)

        # Get only IRAC rows of table
        instru_mask = ['IRAC' in row['modedisplayname'] for row in table]
        before_mask = [Time(row['scet']).mjd < mjd for row in table]
        templa_mask = [Time(row['scet']).mjd > (mjd+900) for row in table]

        # First download the pre-explosion data
        mask1=[all(l) for l in zip(instru_mask,before_mask)]
        for row in table[mask1]:
            outdir = outrootdir + name + '/'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            get_cbcd(row, outrootdir=outdir, tempdir='delme/')

        # Now download post-explosion data
        mask2=[all(l) for l in zip(instru_mask,templa_mask)]
        for row in table[mask2]:
            outdir = outrootdir + name + 'tmpl' + '/'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            get_cbcd(row, outrootdir=outdir, tempdir='delme/')

def download_from_coord(coord, before=None, after=None, outdir='./',
    tempdir='delme/'):
    # Get the table of all Spitzer data
    table = post_request(coord)

    # Get only IRAC rows of table
    instru_mask = ['IRAC' in row['modedisplayname'] for row in table]
    mask = None
    if before:
        before_mask = [Time(row['scet']).mjd < before for row in table]
        mask=[all(l) for l in zip(instru_mask, before_mask)]
    elif after:
        after_mask = [Time(row['scet']).mjd > after for row in table]
        mask=[all(l) for l in zip(instru_mask, after_mask)]
    else:
        mask = instru_mask

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for row in table[mask]:
        get_cbcd(row, outrootdir=outdir, tempdir=tempdir)

def get_bool_from_coord(coord, before=None, after=None):
    # Get the table of all Spitzer data
    table = post_request(coord)

    # Get only IRAC rows of table
    instru_mask = ['IRAC' in row['modedisplayname'] for row in table]
    mask = None
    if before:
        before_mask = [Time(row['scet']).mjd < before for row in table]
        mask=[all(l) for l in zip(instru_mask, before_mask)]
    elif after:
        after_mask = [Time(row['scet']).mjd > after for row in table]
        mask=[all(l) for l in zip(instru_mask, after_mask)]
    else:
        mask = instru_mask

    if len(table[mask])>0:
        return(True)
    else:
        return(False)

obj=sys.argv[1]
ra=sys.argv[2]
dec=sys.argv[3]
outdir=sys.argv[4]

outrootdir = outdir+'/'+obj+'/'
coord = SkyCoord(ra,dec,unit=(u.hour,u.deg))
download_from_coord(coord, outdir=outrootdir)
