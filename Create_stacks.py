# UKIRT pipeline function
# Owen Eskandari
# 6/16/20

# Function to create a stack from a dictionary of files sorted by filter,
# a path to write them to, and an extension to use
# Can create stacks in three different ways: no alignment, wcs alignment, and pixel alignment
# Pixel alignment function is in a separate folder

# Imports
from astropy.io import fits
import astropy.wcs as wcs
import numpy as np
from astropy.nddata import CCDData
from photutils import MeanBackground
from astropy.stats import SigmaClip
import astropy.units as u
import ccdproc
from ccdproc import wcs_project

# Function imports
# from Pixel_alignment import pixel_alignment       # TODO: When this function is completed

from astropy.utils.exceptions import AstropyWarning
import warnings

warnings.simplefilter('ignore', category=AstropyWarning)        # Removes deprecation warnings


# Creates stacks by filter: subtracts the bkg, and re-bins appropriately for each file; combines into a stack
# mult_ext is True by default (False if files only have a single extension)
# new align variable (None, 'wcs', or 'pixel')
def create_stacks(filter_list, old_path, path, ext, log=None, log2=None, mult_ext=True, align=None):

    '''

    Function used to create one or more stacks of a specific extension from a dictionary of files sorted by filter.
    Currently supports stacking with ``wcs``-alignment, soon will also support ``pixel``-alignment.
    Also supports stacking without any alignment.
    Stacks created have name of event, filter, extension, and type of alignment (if used) in their name and are written
    to ``path``.


    Parameters
    ----------

    :param filter_list: python dictionary
        Dictionary of files. Key is the filter of the files, values are the file names themselves.

    :param old_path: str
        Path to the original data.

    :param path: str
        Path to the reduced data; usually ``old_path`` + '/red/ex' + str(``ext``).

    :param ext: int
        The extension of each original file to use.

    :param log: log, optional
        In-depth log.
        If no log is inputted, information is printed out instead of being written to ``log``.
        Default is ``None``.

    :param log2: log, optional
        Overview log.
        If no log is inputted, information is printed out instead of being written to ``log2``.
        Default is ``None``.

    :param mult_ext: boolean, optional
        Boolean as to whether the original files have multiple extensions. Set to ``False`` if only one extension.
        Default is ``True``.

    :param align: str, optional
        Type of alignment to perform. Currently supports ``wcs`` (and ``pixel``).
        Pixel alignment is still in the process of being written. ``TODO`` comments highlight where to put in new code
        to add in pixel alignment functionality. Remove this line and ``TODO`` comments after adding functionality in.
        Default is ``None`` (no alignment performed).

    Returns
    -------

    :return: python dictionary
        Dictionary of stacks. Key is the filter of the stack, values are the stack names themselves.

    '''

    if log2 is not None:
        log2.info("Filename & Filter")
    else:
        print("Filename & Filter")

    stack_dict = {}        # Dictionary for stack files
    delta_t = 10           # Exposure time, in sec

    for fil in filter_list:
        x_shape, y_shape = 1000000000, 1000000000       # For shape (initially set to a very large value)
        red_list = []
        reprojected = []        # For WCS projections

        for i, f in enumerate(filter_list[fil]):

            if log2 is not None:
                log2.info("%s %s" % (f.split("/")[len(f.split("/")) - 1], fil.upper()))
            else:
                print("%s %s" % (f.split("/")[len(f.split("/")) - 1], fil.upper()))
            if log is not None:
                log.info(f + " " + fil)

            with fits.open(f) as hdr:               # Open the file and extract the header and the data
                if i == 0:                          # Get the event name
                    hdr0 = hdr[0].header
                    event_name = hdr0['OBJECT']

                if mult_ext:
                    header = hdr[ext].header
                    data = hdr[ext].data
                else:                           # In the case that the files only have one extension
                    header = hdr[0].header
                    data = hdr[0].data

            if i == 0 and align == "wcs":
                wcs_object = wcs.WCS(header)     # Create object for WCS projection

            # Find background
            data[np.isnan(data)] = np.nanmedian(data)           # Get rid of any nan data
            sci_data = CCDData(data, meta=header, unit='adu')   # Make a CCDData object with same header and data in adu
            bkg = MeanBackground(SigmaClip(sigma=3.))           # Use sigma clipping to find the mean background
            bkg_value = bkg.calc_background(sci_data)           # Find the background of the CCDData object

            # Subtract background
            red = sci_data.subtract(bkg_value * np.ones(np.shape(sci_data)) * u.adu, propagate_uncertainties=True,
                                    handle_meta='first_found')

            # Divide by exposure time to get all files on the same scale
            red = red.divide(delta_t * u.second, propagate_uncertainties=True, handle_meta='first_found')

            # Setting the shape (finding the minimum to rebin later)
            if red.shape[0] < x_shape:
                x_shape = red.shape[0]
            if red.shape[1] < y_shape:
                y_shape = red.shape[1]

            if align == "wcs":      # Reprojecting the image onto a common WCS
                red.wcs = wcs.WCS(header)
                reprojected.append(wcs_project(red, wcs_object))

            if align == "pixel":                   # TODO: Add pixel alignment function here in the future
                print("Not available yet. Images will not be aligned")

            red_list.append(red)        # Add to a list to combine later

            # Write to a new directory
            red.write(f.replace(old_path, path).replace(('.'+f.split('.')[1]),
                                                        '_c' + str(ext) + '_' + str(fil) + '_red.fits'), overwrite=True)

        # Re-bin so all of the files are of the same shape
        final_list = []
        for red in red_list:
            red1 = ccdproc.rebin(red, (x_shape, y_shape))
            final_list.append(red1)

        # Stack the reduced files that have been re-bined so all have the same dimensions
        stack_hbu = ccdproc.combine(final_list, method='median', sigma_clip=True, sigma_clip_func=np.ma.median)
        stack_hbu.write(path + event_name + "_ext_" + str(ext) + "_" + fil + '_stack.fits', overwrite=True)

        if align == "wcs":
            # Re-bin so all of the files are of the same shape
            final_list = []
            for red in reprojected:
                red1 = ccdproc.rebin(red, (x_shape, y_shape))
                final_list.append(red1)
            # Stack, WCS aligned
            stack_hbu = ccdproc.combine(final_list, method='median', sigma_clip=True, sigma_clip_func=np.ma.median)
            stack_hbu.write(path + event_name + "_ext_" + str(ext) + "_" + fil + '_stack_wcs.fits', overwrite=True)

        if align == "pixel":    # TODO: Add pixel alignment function here in the future
            print("Not available yet. No pixel-aligned stacks")

            # # Re-bin so all of the files are of the same shape
            # final_list = []
            # for red in pixel_reprojected:
            #     red1 = ccdproc.rebin(red, (x_shape, y_shape))
            #     final_list.append(red1)
            # # Stack, pixel aligned
            # stack_hbu = ccdproc.combine(final_list, method='median', sigma_clip=True, sigma_clip_func=np.ma.median)
            # stack_hbu.write(path + event_name + "_ext_" + str(ext) + "_" + fil + '_stack_pixel.fits', overwrite=True)

        try:
            stack_dict[fil]
        except KeyError:
            stack_dict.update({fil: []})
        stack_dict[fil].append(path + event_name + "_ext_" + str(ext) + "_" + fil + '_stack.fits')
        if align == "wcs":
            stack_dict[fil].append(path + event_name + "_ext_" + str(ext) + "_" + fil + '_stack_wcs.fits')
        if align == "pixel":     # TODO: Uncomment when pixel-aligning function has been added
            print("Not available yet. No pixel-aligned stacks")
            # stack_dict[fil].append(path + event_name + "_ext_" + str(ext) + "_" + fil + '_stack_pixel.fits')

    if align is not None and log2 is not None:
        if align != "pixel":       # TODO: Take out later when pixel-aligning function has been added
            log2.info("%s-aligned images" % str(align))

    return stack_dict
