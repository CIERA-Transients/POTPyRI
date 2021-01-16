# UKIRT pipeline function
# Owen Eskandari
# 6/16/20

# This function finds the extension (tile) to use based off coordinates given (RA, Dec)

# Imports
import numpy as np
from astropy.io import fits
import astropy.wcs as wcs
from astropy.nddata import CCDData


# Finds the correct tile to make stack and find photometry
# Inputs: filename (string), RA & Dec (floats)
def find_tile(file, ra, dec):

    '''

    Function used to find the extension which contains a set of coordinates (RA, Dec).

    Parameters
    ----------

    :param file: string
        Name of a multi-extension file.

    :param ra: float
        RA of the position of the target.

    :param dec: float
        Dec of the position of the target.

    Returns
    -------

    :return: int, tuple
        Returns the number of the extension the coordinates are in
        and a tuple of the pixel coordinates at this position.

    Raises
    ------

    ValueError: `ValueError`
        If the none of the extensions of the file contain the (RA, Dec) coordinates inputted.

    '''

    coords = np.array([[ra, dec]], np.float)  # Array of coordinates: [[RA, Dec]] in [deg, deg]
    with fits.open(file) as hdr:
        for i in range(1, 5):        # Skip HDU header (i == 0)
            header = hdr[i].header
            data = hdr[i].data
            w = wcs.WCS(header)     # Parse the WCS keywords in the primary HDU
            # Find the pixel coordinates in the image
            # The second argument is "origin" -- in this case we're declaring we have 1-based (Fortran-like) coordinates
            pixel_coords = w.wcs_world2pix(coords, 1)
            shape = CCDData(data, unit='adu').shape     # To set parameters when determining

            if (0 <= pixel_coords[0][0] <= shape[0]) and (0 <= pixel_coords[0][1] <= shape[1]):  # If in this extension
                return i, pixel_coords
    raise ValueError("RA, Dec not in any of the extensions in the file")
