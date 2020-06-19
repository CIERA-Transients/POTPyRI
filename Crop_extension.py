# UKIRT pipeline function
# Owen Eskandari
# 6/16/20

# Simple function to crop a multi-extension file

# Imports
from astropy.io import fits
from astropy.nddata import CCDData


def crop_extension(f, ext):

    '''

    Simple function to crop a multi-extension file into a single extension file.
    No return value, but new file is written to same location as old file with ``_c`` + ``ext`` + ``.fits``.
    Only problem with this is that information in the primary header of the multi-extension file may be lost.

    Parameters
    ---------

    :param f: str
        Filename and path to filename of multi-extension file

    :param ext: int
        Number of the extension to use when cropping

    '''

    with fits.open(f) as hdr:
        header = hdr[ext].header
        data = hdr[ext].data
    new_file = CCDData(data, meta=header, unit='adu')  # Make a CCDData object with same header and data in adu
    new_file.write(f.replace(('.' + f.split('.')[len(f.split('.')) - 1]), '_c' + str(ext) + '.fits'), overwrite=True)
