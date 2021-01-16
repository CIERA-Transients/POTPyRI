from __future__ import print_function

import sys,os,pdb,shutil,argparse
import numpy as np 

from astropy.io import fits
from astropy.io import ascii
from scipy import signal, ndimage

from matplotlib import pyplot as plt


def guess_lamps(file_list):
    return None


def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Generate a helpful plot for line IDs '
    #descStr += 'Recommended calling sequence: \n \n'
    #descStr += '$ python keck_basic_2d.py -v -c \n'
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # required args
    #parser.add_argument('requried_arg',type=str,
    #                    help='a required arguement')

    # optional
    parser.add_argument('-v','--verbose',
                        help='print diagnostic info',action='store_true')
    parser.add_argument('-c','--clobber',action='store_true',
                        help='Clobber files already in pre_reduced/ but not subdirs')

    # parse
    cmdArgs = parser.parse_args()

    # logic mapping to my args/kwargs
    VERBOSE = cmdArgs.verbose
    CLOBBER = cmdArgs.clobber

    # package up
    args = () # no args implemented yet
    kwargs = {}
    kwargs['VERBOSE'] = VERBOSE
    kwargs['CLOBBER'] = CLOBBER

    return (args,kwargs)

def generate_plot(file_list,*args,**kwargs):
    ''' Manage the plotting window '''

    # basic initializations
    master_coadd = np.array([])
    img_secs_dict = {}
    ax_dict = {}
    color_list = ['#2577fa','#18a300','#a30000']

    # configure the plotting window
    fig=plt.figure(figsize=(8,5))

    # set up plotting grid
    nGrid = 32
    nLamps = len(file_list)
    nRowsIntensityPlot = nGrid // 2
    nRowsPerLamp = (nGrid - nRowsIntensityPlot) // nLamps
    lampAspectFactor = nRowsPerLamp / nRowsIntensityPlot
    intensityAspect = 30
    ax_dict['main'] = plt.subplot2grid((nGrid,nGrid), (0,0), 
                            rowspan=nGrid, colspan=nGrid)
    ax_dict['coadd'] = plt.subplot2grid((nGrid,nGrid), (0,0), 
                            rowspan=nRowsIntensityPlot, colspan=nGrid)
    ax_dict['coadd'].set_ylabel('Intensity')
    ax_dict['coadd'].set_xticklabels([])
    ax_dict['coadd'].set_yticklabels([])
    ax_dict['coadd'].format_coord = lambda x, y: 'x={:.2f}, y={:.2f}'.format(x,y)



    # determine lamps/plot labels
    lamp_list = kwargs.get('lamp_list',None)
    if lamp_list is None:
        lamp_list = guess_lamps(file_list)

    # read in the data
    for i,file in enumerate(file_list):
        # open
        hdu = fits.open(file)
        data = hdu[0].data
        header = hdu[0].header

        # allocate space in the plot and set labels
        plotGridRowStart = nRowsIntensityPlot + i*nRowsPerLamp
        ax_dict[lamp_list[i]] = plt.subplot2grid((nGrid,nGrid), (plotGridRowStart,0), 
                                    rowspan=nRowsPerLamp, colspan=nGrid,
                                    sharex=ax_dict['coadd'])
        ax_dict[lamp_list[i]].set_ylabel('row')
        ax_dict[lamp_list[i]].set_xticklabels([])
        ax_dict[lamp_list[i]].set_yticklabels([])
        xticklabels = [1000*x for x in np.arange(12) if x < data.shape[1]//1000]
        ax_dict[lamp_list[i]].set_xticks(xticklabels)
        ax_dict[lamp_list[i]].format_coord = lambda x, y: 'x={:.2f}, y={:.2f}'.format(x,y)


        if i == len(file_list)-1:
            ax_dict[lamp_list[i]].set_xlabel('column')
            ax_dict[lamp_list[i]].set_xticklabels(xticklabels)

        # if we want
        #ax_dict[lamp_list[i]].set_xticks([0,100,200,300,400])
        #ax_dict[lamp_list[i]].set_xticklabels([0,100,200,300,400])


        # cut out img sections (grab a chunk in the upper half)
        nrows = data.shape[0]
        ncols = data.shape[1]
        row_start = nrows // 2 - 200
        row_end = row_start + 30
        img_sec = data[row_start:row_end,:]

        # store the 1D profile
        img_secs_dict[lamp_list[i]] = np.median(img_sec,axis=0)

        # throw up the data
        ax_dict[lamp_list[i]].imshow(img_sec,
                                #aspect=intensityAspect*lampAspectFactor,
                                aspect='auto',
                                origin='lower',
                                vmin=10.,
                                vmax=1000.)
        # annotate here
        xTextCoord = 10.
        yTextCoord = 10.
        ax_dict[lamp_list[i]].annotate(lamp_list[i],
                                       (xTextCoord,yTextCoord),
                                       xycoords='axes points',
                                       color=color_list[i],ha='left',
                                       weight = 'bold',fontsize=15.
        )

        # throw up the profile on the coadd plot
        ax_dict['coadd'].plot(np.arange(img_sec.shape[1]),
                              np.median(img_sec,axis=0),
                              color=color_list[i],
                              lw = 0.5,
                              ls = '-'
        )

        # coadd sections for intensity plot
        if len(master_coadd) > 0:
            master_coadd += img_sec
        else:
            master_coadd = np.copy(img_sec)

    # throw up the master_coadd
    # ax_dict['coadd'].imshow(master_coadd,
    #                     aspect=intensityAspect,
    #                     origin='lower',
    #                     vmin=10.,
    #                     vmax=1000.)

    ax_dict['coadd'].plot(np.arange(master_coadd.shape[1]),
                          np.median(master_coadd,axis=0),
                          color='k',
                          lw = 0.5,
                          ls = '-'
    )
    ax_dict['coadd'].set_xlim([0,master_coadd.shape[1]])


    # finally, show
    plt.show()


    # return axis? no, lets have this manage the plot
    return 0



def main(*args,**kwargs):
    '''
    Generate a useful set of plots to help the user ID arc lines



    '''
    file_list = [
        'test_data/tor190630_1054.fits', 
        'test_data/tor190630_1057.fits'
    ]

    lamp_list = ['CdNe','HgAr']

    plot_args = (file_list,)
    plot_kwargs = {'lamp_list':lamp_list,}
    res = generate_plot(*plot_args,**plot_kwargs)

    return 0


if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)