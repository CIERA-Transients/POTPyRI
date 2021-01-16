#!/usr/bin/env python
import getpass
import logging
import copy
# need copy.deepcopy() for Spec class
import os
import datetime
from astropy.io import fits
from tmath.wombat.specclass import Spec
import tmath.wombat as wombat
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.widgets import Cursor
import matplotlib.pyplot as plt

# needs extinction, astropy.io, scipy

matplotlib.rcParams["savefig.directory"] = "."
matplotlib.rcParams["xtick.minor.visible"] = True
matplotlib.rcParams["ytick.minor.visible"] = True
matplotlib.rcParams["lines.linewidth"] = 0.5



#global nsplinepoints, tmpsplptsx, tmpsplptsy, pflag
WOMVERSION = 0.1

def main():
    # logging straight from docs.python.org cookbook page
    #  INFO level to screen and wombat.log
    #  DEBUG level only to wombat.log
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename='wombat.log',
                        filemode='a')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)


    user = getpass.getuser()
    print('Hello, {} \n'.format(user))
    hop = [Spec() for i in range(21)]

    print('Welcome to Python Wombat, V{}.  ? or help lists commands'.format(WOMVERSION))

    plt.ion()
    screen_width, screen_height = wombat.get_screen_size()
    screenpos = '+{}+{}'.format(int(screen_width*0.2), int(screen_height*0.05))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cursor = Cursor(ax, useblit=True, color='k', linewidth=1)
    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('Wombat')
    fig.set_size_inches(8, 5)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    # this should bring the window to the top, but it doesn't
    wm = plt.get_current_fig_manager()
    #fig.canvas.manager.window.attributes('-topmost', 1)
    blah = wm.window.attributes('-topmost', 1)

    #plt.pause(0.2)
    #fig.canvas.manager.window.attributes('-topmost', 0)
    blah = wm.window.attributes('-topmost', 0)
    #ax.plot([0],[1])

    logging.debug('Wombat starts')
    logging.debug('{}'.format(str(datetime.datetime.now())))
    logging.debug('{} running things'.format(user))

    command_dict = {'atmdisp':wombat.wommkatmdisp,
                    'hop':wombat.womhop,
                    'rh':wombat.womrdhop,
                    'rhop':wombat.womrdhop,
                    'buf':wombat.wombuf,
                    'rpl':wombat.womrpl,
                    'rtxt':wombat.womrpl,
                    'wtxt':wombat.womwpl,
                    'wpl':wombat.womwpl,
                    'p':wombat.womplot,
                    'plot':wombat.womplot,
                    'apl':wombat.womapl,
                    'oplot':wombat.womapl,
                    'bb':wombat.wommkbb,
                    'planck':wombat.wommkbb,
                    'rfits':wombat.womrdfits,
                    'wfits':wombat.womwfits,
                    'cho':wombat.womcho,
                    'choose':wombat.womcho,
                    'bin':wombat.womscipyrebinselector,
                    'ashrebin':wombat.womashrebinselector,
                    'stat':wombat.womstat,
                    'e':wombat.womexam,
                    'examine':wombat.womexam,
                    'smo':wombat.womsmo,
                    'smooth':wombat.womsmo,
                    'plotlog':wombat.womplotlog,
                    'plotl':wombat.womplotlog,
                    'pl':wombat.womplotlog,
                    'sca':wombat.womscale,
                    'scale':wombat.womscale,
                    'avt':wombat.womscale2match,
                    'scale2':wombat.womscale2match,
                    'scalemany':wombat.womscalemany,
                    'avmany':wombat.womscalemany,
                    'w':wombat.womredshift,
                    'redshift':wombat.womredshift,
                    'blueshift':wombat.womblueshift,
                    'ms':wombat.womms,
                    'arith':wombat.womms,
                    'fnuflam':wombat.womfnuflam,
                    'flux':wombat.womfnuflam,
                    'lup':wombat.womlup,
                    'bluen':wombat.wombluen,
                    'relvel':wombat.womrelvel,
                    'xshift':wombat.womxshift,
                    'yshift':wombat.womyshift,
                    'zap':wombat.womzap,
                    'com':wombat.womcomselector,
                    'combine':wombat.womcomselector,
                    'cmw':wombat.womcmwselector,
                    'com3':wombat.womcom3selector,
                    'commany':wombat.womcommanyselector,
                    'commanyw':wombat.womcommanywselector,
                    'join':wombat.womjoin,
                    'red':wombat.womreddening,
                    'head':wombat.womheader,
                    'header':wombat.womheader,
                    'rmsfits':wombat.womrmsfits,
                    'filter':wombat.womfilters,
                    'irfilter':wombat.womirfilters,
                    'hertz':wombat.womhertz,
                    'sqr':wombat.womsqr,
                    'sqrt':wombat.womsqrt,
                    'velcho':wombat.womvelcho,
                    'int':wombat.womint,
                    'gau':wombat.womgau,
                    'wavescale':wombat.womwavescale,
                    'cat':wombat.womcat,
                    'help':wombat.womhelp,
                    '?':wombat.womhelp
    }
    commandarr = []
    done = False
    while (not done):
        print(' ')
        command = input('Enter command: ')
        command = command.strip()
        commandarr.append(command)
        if len(command) > 0:
            if command[0] == '$':
                os.system(command[1:])
                command = ''
        if len(commandarr) > 2:
            if (commandarr[-1] == commandarr[-2] == commandarr[-3]):
                print("\nMonotonous, isn't it?")
        if (command == 'q') or (command == 'quit'):
            answer = ''
            while (answer != 'y') and (answer != 'n'):
                reply = input('Really? (y/n, default y) ')
                reply = reply.strip()
                if len(reply) == 0:
                    reply = 'y'
                answer = reply.lower()[0]
            if (answer == 'n'):
                command = ''
            else:
                print('So long\n')
                done = True
        if (command == 'h') or (command == 'hist'):
            for i, _ in enumerate(commandarr):
                print(commandarr[i])
            command = ''

# spline and blotch (which uses splines) have to have 'fig'
# passed for the connection to the mouse click to work

        if (command == 'spl') or (command == 'spline'):
            hop = wombat.womspl(hop, fig)
        if (command == 'b') or (command == 'blotch'):
            hop = wombat.womblo(hop, fig)
        if (command in command_dict):
            hop = command_dict[command](hop)

    logging.info('moooo')
    plt.close()
    plt.ioff()
    #fig.delaxes(ax)

if __name__ == '__main__':
    main()
