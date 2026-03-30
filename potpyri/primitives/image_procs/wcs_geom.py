"""WCS geometry helpers: PV removal, field center, synthetic TAN WCS."""
from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from potpyri.utils import utilities


def remove_pv_distortion(header):
    """Remove PV* distortion keywords from a FITS header (modifies in place).

    Parameters
    ----------
    header : astropy.io.fits.Header
        Header to modify.

    Returns
    -------
    astropy.io.fits.Header
        The same header reference (modified in place).
    """
    done = False
    while not done:
        bad_key = False
        for key in header.keys():
            if key.startswith('PV'):
                if key in header.keys():
                    bad_key = True
                    del header[key]

        if not bad_key: done = True

    return(header)

def get_fieldcenter(images):
    """Compute mean RA/Dec of image centers from a list of FITS files.

    Parameters
    ----------
    images : list of str
        Paths to FITS files with WCS in primary header.

    Returns
    -------
    list
        [mean_ra_deg, mean_dec_deg].
    """
    ras = [] ; des = []
    for file in images:
        with fits.open(file) as hdu:
            w = WCS(hdu[0].header)

            data_shape = hdu[0].data.shape
            center_pix = (data_shape[1]/2., data_shape[0]/2.)

            coord = w.pixel_to_world(*center_pix)

            ras.append(coord.ra.degree)
            des.append(coord.dec.degree)

    mean_ra = np.mean(ras)
    mean_de = np.mean(des)

    return([mean_ra, mean_de])

def generate_wcs(tel, binn, fieldcenter, out_size):
    """Build a TAN WCS for a given field center and output size.

    Parameters
    ----------
    tel : Instrument
        Instrument instance (for pixel scale).
    binn : str
        Binning string (e.g. '22').
    fieldcenter : sequence
        [ra_deg, dec_deg] or parseable coordinate.
    out_size : int
        Output image size (out_size x out_size pixels).

    Returns
    -------
    astropy.wcs.WCS
        WCS instance.
    """
    w = {'NAXES':2, 'NAXIS1': out_size, 'NAXIS2': out_size,
        'EQUINOX': 2000.0, 'LONPOLE': 180.0, 'LATPOLE': 0.0}

    pixscale = tel.get_pixscale()
    cdelt = pixscale * int(str(binn)[0])/3600.0

    w['CDELT1'] = -1.0 * cdelt
    w['CDELT2'] = cdelt

    w['CTYPE1']='RA---TAN'
    w['CTYPE2']='DEC--TAN'
    w['CUNIT1']='deg'
    w['CUNIT2']='deg'

    w['CRPIX1']=float(out_size)/2 + 0.5
    w['CRPIX2']=float(out_size)/2 + 0.5

    coord = utilities.parse_coord(fieldcenter[0], fieldcenter[1])

    w['CRVAL1']=coord.ra.degree
    w['CRVAL2']=coord.dec.degree

    w = WCS(w)

    return(w)
