# UKIRT pipeline function
# Owen Eskandari
# 6/16/20

# This function takes in a list of files (and a filter if the filter is not in the header)
# and returns a dictionary of the files sorted by filter

# Imports
from astropy.io import fits
import shutil

# Sort the calibration files:
def sort_files(files, manual_filter=None, log2=None, date=None):

    '''

    Function used to sort a list of files into a dictionary of files sorted by filter.

    Parameters
    ----------

    :param files: list (string)
        List of strings of files (path should be included).

    :param manual_filter: string, optional
        Filter name if filter is not given in header of files.
        Default is ``None``.

    :param log2: log, optional
        Overview log used to write the object and date observed (if ``date`` parameter is not ``None``).
        If no log is inputted, information is printed out instead of being written to ``log2``.
        Default is ``None``.

    :param date: string, optional
        String of the date the observations were taken (to be recorded in ``log2`` if it is not ``None``).
        Default is ``None``.

    Returns
    -------

    :return: python dictionary
        Dictionary of files. Key is the filter of the file, values are the file names themselves.

    '''

    filter_list = {}        # Create new dictionaries for the files filtered by filter

    for i, f in enumerate(files):
        with fits.open(f) as file_open:
            hdr = file_open[0].header

        if i == 0:
            event_name = hdr['OBJECT']
            if log2 is not None:
                log2.info("Event: " + event_name)
                if date is not None:
                    log2.info("Date observed: %s" % date)
            else:
                print("Event: " + event_name)
                if date is not None:
                    print("Date observed: %s" % date)

        try:
            fil = hdr['FILTER']         # add filter
        except KeyError:
            fil = manual_filter.upper()
        try:
            filter_list[fil]
        except KeyError:
            filter_list.update({fil: []})

        filter_list[fil].append(f)

    return filter_list

def sort_files_bino(files,raw_path): #sort the calibration files: 
    cal_list = []
    sci_list = []
    for f in files:
        with fits.open(f) as file_open:
            hdr = file_open[1].header
            imtype = hdr['MASK']
            if imtype == 'imaging':
                fil = hdr['FILTER'] #sort by filter
                if hdr['SCRN']=='deployed': #only collect dome flats for imaging    
                    cal_list.append(f)
                elif hdr['SCRN']=='stowed':
                    target = hdr['OBJECT']
                    sci_list.append(f)
                else:
                    shutil.move(f,raw_path+'bad/')
            else:
                shutil.move(f,raw_path)
    return cal_list, sci_list