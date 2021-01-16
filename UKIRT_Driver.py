# UKIRT pipeline
# Owen Eskandari
# 6/16/20

# Imports
import os
import datetime
import time
import glob
import argparse
import logging
from astropy.io import fits

# Function imports
from Sort_files import sort_files
from Create_stacks import create_stacks
from Find_tile import find_tile
from Find_zero_point import find_zero_point
from Find_target_phot import find_target_phot

from astropy.utils.exceptions import AstropyWarning
import warnings

warnings.simplefilter('ignore', category=AstropyWarning)        # Removes deprecation warnings

# Set of parameters to run the program
params = argparse.ArgumentParser(description='Path of data.')
params.add_argument('data_path', default=None, help='Path of data')     # Path to the UKIRT data
# Set to yes if you don't yet have a stack
params.add_argument('--create_stack', type=str, default=None, help='Option to create stack ("yes").')
# Set to no if you would not like to see plots/cannot see plots
params.add_argument('--show_plots', type=str, default=True, help='Option to show plots (default is "yes").')
# If known, the extension of the target
params.add_argument('--ext', type=int, default=None, help='Extension of target.')
# If known, the RA of the target
params.add_argument('--ra', type=float, default=None, help='RA of target, in decimal format.')
# If known, the Dec of the target
params.add_argument('--dec', type=float, default=None, help='Dec of target, in decimal format.')
params.add_argument('--use_ext', type=int, default=None, help='Set to ext # if using files with only one extension')
params.add_argument('--fil', type=str, default=None, help='Type the filter if not in image headers')
params.add_argument('--align', type=str, default=None, help='"wcs" or "pixel" (default: None)')
params.add_argument('--date', type=str, default=None, help='Date of observation (puts in log)')


args = params.parse_args()

# Define paths (for sci, red, and extensions)
data_path = args.data_path
sci_path = data_path
red_path = sci_path + '/red/'
ex1_path, ex2_path, ex3_path, ex4_path = red_path + '/ex1/', red_path + '/ex2/', red_path + '/ex3/', red_path + '/ex4/'
if not os.path.exists(red_path):    # Create file paths if they don't exist
    os.makedirs(red_path)
if not os.path.exists(ex1_path):
    os.makedirs(ex1_path)
if not os.path.exists(ex2_path):
    os.makedirs(ex2_path)
if not os.path.exists(ex3_path):
    os.makedirs(ex3_path)
if not os.path.exists(ex4_path):
    os.makedirs(ex4_path)

# To initiate the in-depth log
log_file_name = red_path+'UKIRT_log_'+datetime.datetime.utcnow().strftime('%m%d_%H%M%S')+'.log'  # create log filename
log = logging.getLogger(log_file_name)                  # create logger
log.setLevel(logging.INFO)                              # set level of logger
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")  # set format of logger
logging.Formatter.converter = time.gmtime               # convert time in logger to UCT
filehandler = logging.FileHandler(log_file_name, 'w+')  # create log file
filehandler.setFormatter(formatter)                     # add format to log file
log.addHandler(filehandler)                             # link log file to logger
streamhandler = logging.StreamHandler()                 # add log stream to logger
streamhandler.setFormatter(formatter)                   # add format to log stream
log.addHandler(streamhandler)                           # link logger to log stream

log.info("UKIRT/WFCAM Pipeline")
log.info('Reducing data from '+data_path)
log.info("Arguments inputted: %s" % args)

# To initiate the overview log
log2_file_name = red_path+'UKIRT_log_overview_'+datetime.datetime.utcnow().strftime('%m%d_%H%M%S')+'.log'  # create log
log2 = logging.getLogger(log2_file_name)                  # create logger
log2.setLevel(logging.INFO)                              # set level of logger
filehandler2 = logging.FileHandler(log2_file_name, 'w+')  # create log file
log2.addHandler(filehandler2)                             # link log file to logger
log2.addHandler(streamhandler)                           # link logger to log stream

log2.info("UKIRT/WFCAM Pipeline")
log2.info('Reducing data from '+data_path)
log2.info('Pipeline user: ' + input("Who is running this? Type your name: "))

show_plots = True
if args.show_plots == "no":
    show_plots = False

# Finds the files to use
mult_ext = True            # Assumes the files have multiple extensions
if args.use_ext is not None:
    files = glob.glob(args.data_path + '/sw*_sf_st*.fits')  # Originally .fz
    if len(files) == 0:
        files = glob.glob(args.data_path + '/sw*_sf_st*.fit')  # Originally .fz
    args.ext = args.use_ext
    mult_ext = False
else:
    files = glob.glob(args.data_path+'/w*_sf_st.fits')     # Originally .fz
    if len(files) == 0:
        files = glob.glob(args.data_path+'/w*_sf_st.fit')     # Originally .fz

# Stack needs to be performed
if args.create_stack == "True":
    log.info("Number of files to use in stack: %s" % len(files))
    log2.info("Number of files to use in stack: %s" % len(files))
    filter_list = sort_files(files, args.fil, log2=log2, date=args.date)   # Make a dict of files sorted by fil

    # If correct extension given
    if args.ext is not None:
        ex_path = red_path + '/ex' + str(args.ext) + '/'
        log.info(create_stacks(filter_list, sci_path, ex_path, args.ext, log=log, log2=log2, mult_ext=mult_ext,
                               align=args.align))  # Saves stack files in a dictionary sorted by filter

    # If the coordinates are given
    elif (args.ra is not None) & (args.dec is not None):
        log.info("RA inputted: " + str(args.ra))
        log.info("Dec inputted: " + str(args.dec))
        log2.info("RA inputted: " + str(args.ra))
        log2.info("Dec inputted: " + str(args.dec))

        # Pick out the first file in the dictionary to use to find the extension
        i = 0
        orig_file = None
        for fil in filter_list:
            orig_file = filter_list[fil][0]
            break

        if orig_file is None:
            log.info("No file found")

        else:
            ext, pixel_coords = find_tile(orig_file, args.ra, args.dec)     # Find the correct extension

            if ext is not None:
                log.info("Given the RA and Dec inputted, the object is at pixel (" + str(int(pixel_coords[0][0])) + ", "
                         + str(int(pixel_coords[0][1])) + ") in extension " + str(ext) + ".")
                log2.info("Target at pixel (" + str(int(pixel_coords[0][0])) + ", "
                          + str(int(pixel_coords[0][1])) + ") in extension " + str(ext) + ".")
                ex_path = red_path + '/ex' + str(ext) + '/'
                log.info(create_stacks(filter_list, sci_path, ex_path, ext, log=log, log2=log2, mult_ext=mult_ext,
                                       align=args.align))    # Saves stack files in a dictionary sorted by filter
            else:
                log.info("Extension not found")

    else:   # If neither extension nor coordinates are given
        log.info("You have not inputted coordinates of the target or its extension.")
        proceed = input("Would you like to find the stacks of all extensions? Type 'yes' if so. If not, press return. ")

        if proceed == 'yes':
            log.info(create_stacks(filter_list, sci_path, ex1_path, 1, log=log, log2=log2, mult_ext=mult_ext,
                                   align=args.align))  # Saves stack files in a dictionary sorted by filter
            log.info(create_stacks(filter_list, sci_path, ex2_path, 2, log=log, log2=log2, mult_ext=mult_ext,
                                   align=args.align))  # Saves stack files in a dictionary sorted by filter
            log.info(create_stacks(filter_list, sci_path, ex3_path, 3, log=log, log2=log2, mult_ext=mult_ext,
                                   align=args.align))  # Saves stack files in a dictionary sorted by filter
            log.info(create_stacks(filter_list, sci_path, ex4_path, 4, log=log, log2=log2, mult_ext=mult_ext,
                                   align=args.align))  # Saves stack files in a dictionary sorted by filter
            log.info("Stacks of all extensions have been made")

        else:
            ext = input("Which extension would you like to use? Type the number you would like to use (1, 2, 3, 4): ")
            if ext == "1" or ext == "2" or ext == "3" or ext == "4":
                ex_path = red_path + '/ex' + ext + '/'
                log.info(create_stacks(filter_list, sci_path, ex_path, int(ext), log=log, log2=log2, mult_ext=mult_ext,
                                       align=args.align))    # Saves stack files in a dictionary sorted by filter
            else:
                log.info("Invalid input.")

else:   # Log the event name
    for i, f in enumerate(files):
        if i == 0:
            with fits.open(f) as hdr:
                header = hdr[0].header
            try:
                event_name = header['OBJECT']
                log2.info("Event: " + event_name)
            except KeyError:
                log2.info("Event not in file header")
            if args.date is not None:
                log2.info("Date observed: %s" % args.date)
            break

find_zp = input("Stack has been made. Are you ready to find the zero point of the stack? Type 'yes' to continue: ")
if find_zp != "yes":
    find_zp = input("You did not type yes. Type 'yes' to continue on to find the zero point or 'no' to quit: ")

# Now find the zero point
if find_zp == "yes":
    log.info("Stack has already been made. Time to find the zero point of the field. This will take about 20 seconds")
    stacks = glob.glob(args.data_path+'/red/ex*/*stack*.fits')       # Get the stack
    if len(stacks) == 1:
        log.info("There is 1 stack")
    else:
        log.info("There are %s stacks" % len(stacks))
    show_fwhm, show_bkg, show_contours = False, False, False        # Default is not to show plots
    if show_plots:
        show = input("Would you like to see the fwhm plots? Type 'yes' or 'no': ")
        if show == "yes":
            show_fwhm = True
        show = input("Would you like to see the background plots? Type 'yes' or 'no': ")
        if show == "yes":
            show_bkg = True
        show = input("Would you like to see the contour plots? Type 'yes' or 'no': ")
        if show == "yes":
            show_contours = True

    zp_list = []
    zp_err_list = []
    fwhm_list = []
    for stack in stacks:
        ext, fil = stack.split("_")[2], stack.split("_")[3]            # Extension & Filter
        log.info("Using extension %s and %s filter" % (ext, fil))
        # Note: Zero point and zero point error are in Vega, not AB...converted in photometry
        zp, zp_err, fwhm = find_zero_point(stack, fil, ext, '2MASS', data_path,
                                           show_fwhm=show_fwhm, show_bkg=show_bkg, show_contours=show_contours, log=log)
        zp_list.append(zp)
        zp_err_list.append(zp_err)
        fwhm_list.append(fwhm)
    log.info("Phew we made it! Time to do the interesting stuff")

    # Time to perform photometry on the target itself
    perform_phot = input("Ready to perform the photometry on the target? Type 'yes' to continue: ")
    if perform_phot != "yes":
        perform_phot = input("You did not type yes. Type 'yes' to continue or 'no' to quit: ")

    if perform_phot == "yes":
        log.info("Performing photometry")
        ra, dec = args.ra, args.dec
        stacks = glob.glob(args.data_path+'/red/ex*/*stack*.fits')
        if len(stacks) == 1:
            log.info("There is 1 stack")
        else:
            log.info("There are %s stacks" % len(stacks))

        mag_list = []
        mag_err_list = []
        for i, stack in enumerate(stacks):
            ext, fil = stack.split("_")[2], stack.split("_")[3]             # Extension & Filter
            log2.info("\nUsing stack %s, extension %s, filter %s" %
                      (stack.split("/")[len(stack.split("/")) - 1], ext, fil))
            log.info("Using stack %s" % stack)
            log.info("Using extension %s and %s filter with RA = %s and Dec = %s" % (ext, fil, ra, dec))

            phot = True             # TODO: NEW
            while phot:
                mag, mag_err = find_target_phot(stack, fil, fwhm=fwhm_list[i], zp=zp_list[i], zp_err=zp_err_list[i],
                                                show_phot=show_plots, log=log, log2=log2, ra=ra, dec=dec)
                if mag is not None:
                    mag_list.append(mag)         # In AB magnitudes
                    mag_err_list.append(mag_err)

                # TODO: NEW
                if input("Would you like to perform photometry again on this stack? Type 'yes' or 'no': ") == 'no':
                    phot = False
