def womcatfinal(blue_data, red_data):
    """concatenate data with overlapping wavelength regions"""
    import numpy as np
    import math
    import logging
    import matplotlib.pyplot as plt
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.yesno import yesno
    from tmath.wombat.get_screen_size import get_screen_size
    from tmath.wombat.womwaverange2 import womwaverange2
    from tmath.wombat import HOPSIZE
    from matplotlib.widgets import Cursor
    from tmath.wombat.womspectres import spectres

    plt.ion()
    screen_width, screen_height = get_screen_size()
    screenpos = '+{}+{}'.format(int(screen_width*0.2), int(screen_height*0.05))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cursor = Cursor(ax, useblit=True, color='k', linewidth=1)
    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('Cat')
    fig.set_size_inches(9, 6)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
    # print("\nThis will combine blue and red pieces from two hoppers\n")
    # hopnum1=0
    # hopnum2=0
    # while (hopnum1 < 1) or (hopnum1 > HOPSIZE):
    #     hopnum1=inputter('Enter first hopper: ','int',False)
    # while (hopnum2 < 1) or (hopnum2 > HOPSIZE):
    #     hopnum2=inputter('Enter second hopper: ','int',False)
    # if (hop[hopnum1].wave[0] > hop[hopnum2].wave[0]):
    #     hopnum1,hopnum2=hopnum2,hopnum1
    # wdelt1 = hop[hopnum1].wave[1]-hop[hopnum1].wave[0]
    # wdelt2 = hop[hopnum2].wave[1]-hop[hopnum2].wave[0]
    # # check if wavelength dispersion same
    # if (abs(wdelt1 - wdelt2) > 0.00001):
    #     print('Spectra do not have same Angstrom/pixel')
    #     print('Blue side: {}'.format(wdelt1))
    #     print('Red side: {}'.format(wdelt2))
    #     return hop
    # if hop[hopnum1].wave[-1] < hop[hopnum2].wave[0]:
    #     print('Spectra do not overlap\n')
    #     return hop
    # print("\nOverlap range is {} to {}".format(hop[hopnum2].wave[0],
    #                                            hop[hopnum1].wave[-1]))
    # print("\nPlotting blue side as blue, red side as red\n")

    waveblue=np.asarray(blue_data[0])
    fluxblue=np.asarray(blue_data[1])
    varblue=np.asarray(blue_data[2])
    wavered=np.asarray(red_data[0])
    fluxred=np.asarray(red_data[1])
    varred=np.asarray(red_data[2])

    blue_dw = waveblue[1] - waveblue[0]
    red_dw = wavered[1] - wavered[0]
    if red_dw > blue_dw:
        print ("Interpolating to {} A/pix".format(red_dw))
        interp_wave = np.arange(math.ceil(waveblue[0])+1.*red_dw, math.floor(wavered[-1])-1.*red_dw, dtype=float, step=red_dw)
        blue_range = np.where((interp_wave >= waveblue[0]+1.*red_dw) & (interp_wave <= waveblue[-1]-1.*red_dw))
        red_range = np.where((interp_wave >= wavered[0]+1.*red_dw) & (interp_wave <= wavered[-1]-1.*red_dw))
        blue_interp = spectres(interp_wave[blue_range], waveblue, fluxblue, spec_errs=np.sqrt(varblue))
        red_interp = spectres(interp_wave[red_range], wavered, fluxred, spec_errs=np.sqrt(varred))

        waveblue = blue_interp[0]
        fluxblue = blue_interp[1]
        varblue = blue_interp[2]**2.
        wavered = red_interp[0]
        fluxred = red_interp[1]
        varred = red_interp[2]**2.
    else:
        print ("Interpolating to {} A/pix".format(blue_dw))
        interp_wave = np.arange(math.ceil(waveblue[0])+1.*blue_dw, math.floor(wavered[-1])-1.*blue_dw, dtype=float, step=blue_dw)
        blue_range = np.where((interp_wave >= waveblue[0]+1*blue_dw) & (interp_wave <= waveblue[-1]-1*blue_dw))
        red_range = np.where((interp_wave >= wavered[0]+1*blue_dw) & (interp_wave <= wavered[-1]-1*blue_dw))
        blue_interp = spectres(interp_wave[blue_range], waveblue, fluxblue, spec_errs=np.sqrt(varblue))
        red_interp = spectres(interp_wave[red_range], wavered, fluxred, spec_errs=np.sqrt(varred))

        waveblue = blue_interp[0]
        fluxblue = blue_interp[1]
        varblue = blue_interp[2]**2.
        wavered = red_interp[0]
        fluxred = red_interp[1]
        varred = red_interp[2]**2.

    indexblue=womget_element(waveblue,wavered[0])
    indexred=womget_element(wavered,waveblue[-1])
    fluxcor=1.0
    blue_mean=np.mean(fluxblue[indexblue:])
    red_mean=np.mean(fluxred[0:indexred+1])
    if (blue_mean/red_mean < 0.8) or (blue_mean/red_mean > 1.2):
        fluxcor=blue_mean/red_mean
        print("Averages very different, scaling red to blue for plot")
        print("Red multiplied by {}".format(fluxcor))
    plt.cla()
    plt.plot(waveblue[indexblue:],fluxblue[indexblue:],drawstyle='steps-mid',color='b')
    plt.plot(wavered[0:indexred],fluxred[0:indexred]*fluxcor,drawstyle='steps-mid',color='r')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.pause(0.01)
    print('Change scale?')
    answer=yesno('n')
    if (answer == 'y'):
        xmin_old,xmax_old=plt.xlim()
        ymin_old,ymax_old=plt.ylim()
        done=False
        while (not done):
            plt.xlim([xmin_old,xmax_old])
            plt.ylim([ymin_old,ymax_old])
            print('Click corners of box to change plot scale')
            newlims=plt.ginput(2,timeout=-1)
            xmin=newlims[0][0]
            ymin=newlims[0][1]
            xmax=newlims[1][0]
            ymax=newlims[1][1]
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            print('Is this OK?')
            loopanswer=yesno('y')
            if (loopanswer == 'y'):
                done=True
    print('\nEnter method to select wavelength ranges\n')
    mode=inputter_single('Enter (w)avelengths or mark with the (m)ouse? (w/m) ','wm')
    print('\nChoose end points of region to compute average\n')
    
    waveb,waver,mode=womwaverange2(waveblue[indexblue:],fluxblue[indexblue:],
                                   wavered[0:indexred],fluxred[0:indexred]*fluxcor
                                   ,mode)
    indexblueb=womget_element(waveblue,waveb)
    indexbluer=womget_element(waveblue,waver)
    indexredb=womget_element(wavered,waveb)
    indexredr=womget_element(wavered,waver)
    mean_blue=np.mean(fluxblue[indexblueb:indexbluer+1])
    mean_red=np.mean(fluxred[indexredb:indexredr+1])
    print("\nAverage for {}:{}".format(waveb,waver))
    print("Blue side: {}".format(mean_blue))
    print("Red side:  {}\n".format(mean_red))
    brscale=inputter_single('Scale to blue or red (b/r)? ','br')
    if (brscale == 'b'):
        brscalefac=mean_blue/mean_red
        logging.info('Cat scaling to blue by {}'.format(brscalefac))
        fluxred=fluxred*brscalefac
        varred=varred*brscalefac**2
    else:
        brscalefac=mean_red/mean_blue
        logging.info('Cat scaling to red by {}'.format(brscalefac))
        fluxblue=fluxblue*brscalefac
        varblue=varblue*brscalefac**2

    # print("\nPlotting blue side as blue, red side as red\n") 
    # plt.cla()
    # plt.plot(waveblue[indexblueb:indexbluer+1],fluxblue[indexblueb:indexbluer+1],drawstyle='steps-mid',color='b')
    # plt.plot(wavered[indexredb:indexredr+1],fluxred[indexredb:indexredr+1],drawstyle='steps-mid',color='r')
    # plt.xlabel('Wavelength')
    # plt.ylabel('Flux')
    # plt.pause(0.01)
    # print('Change scale?')
    # answer=yesno('n')
    # if (answer == 'y'):
    #     xmin_old,xmax_old=plt.xlim()
    #     ymin_old,ymax_old=plt.ylim()
    #     done=False
    #     while (not done):
    #         plt.xlim([xmin_old,xmax_old])
    #         plt.ylim([ymin_old,ymax_old])
    #         print('Click corners of box to change plot scale')
    #         newlims=plt.ginput(2,timeout=-1)
    #         xmin=newlims[0][0]
    #         ymin=newlims[0][1]
    #         xmax=newlims[1][0]
    #         ymax=newlims[1][1]
    #         plt.xlim([xmin,xmax])
    #         plt.ylim([ymin,ymax])
    #         print('Is this OK?')
    #         loopanswer=yesno('y')
    #         if (loopanswer == 'y'):
    #             done=True

    # print('\nChoose end points of region to compute average\n')
    
    # waveb,waver,mode=womwaverange2(waveblue[indexblueb:indexbluer+1],
    #                                fluxblue[indexblueb:indexbluer+1],
    #                                wavered[indexredb:indexredr+1],
    #                                fluxred[indexredb:indexredr+1],mode)
    # indexblueb=womget_element(waveblue,waveb)
    # indexbluer=womget_element(waveblue,waver)
    # indexredb=womget_element(wavered,waveb)
    # indexredr=womget_element(wavered,waver)
    # ewadd=inputter_single('Add overlap region (e)qually or with (w)eights (e/w)?','ew')
    # if (ewadd == 'e'):
        # overflux=(fluxblue[indexblueb:indexbluer+1]+fluxred[indexredb:indexredr+1])/2.
        # overvar=(varblue[indexblueb:indexbluer+1]+varred[indexredb:indexredr+1])

        #replacing with inverse variance weighted average
    overflux = np.average([fluxblue[indexblueb:indexbluer+1], fluxred[indexredb:indexredr+1]], 
                            weights = [1./varblue[indexblueb:indexbluer+1], 1./varred[indexredb:indexredr+1]], axis = 0)
    overvar = np.sum([varblue[indexblueb:indexbluer+1], varred[indexredb:indexredr+1]], axis = 0)
    # logging.info('Cat combined sides with inverse variance weighted average')
    # else:
    #     wei_done = False
    #     while (not wei_done):
    #         weiblue=inputter('Enter fractional weight for blue side: ','float',False)
    #         weired=inputter('Enter fractional weight for red side: ','float',False)
    #         if (np.abs((weiblue+weired)-1.0) > 0.000001):
    #             print('Weights do not add to 1.0')
    #         else:
    #             wei_done = True
    #     overflux=(fluxblue[indexblueb:indexbluer+1]*weiblue+
    #               fluxred[indexredb:indexredr+1]*weired)
    #     overvar=(varblue[indexblueb:indexbluer+1]*weiblue**2+
    #               varred[indexredb:indexredr+1]*weired**2)
        # logging.info('Cat adds blue side with weight {} and red side with weight {}'.format(weiblue, weired))
    newbluewave=waveblue[:indexblueb]
    newblueflux=fluxblue[:indexblueb]
    newbluevar=varblue[:indexblueb]
    overwave=waveblue[indexblueb:indexbluer+1]
    newredwave=wavered[indexredr+1:]
    newredflux=fluxred[indexredr+1:]
    newredvar=varred[indexredr+1:]
    newwave=np.concatenate([newbluewave, overwave, newredwave])
    newflux=np.concatenate([newblueflux, overflux, newredflux])
    newvar=np.concatenate([newbluevar, overvar, newredvar])
    # logging.info('File {} and'.format(hop[hopnum1].obname))
    # logging.info('File {} concatenated'.format(hop[hopnum2].obname))
    # logging.info('over wavelength range {} to {}'.format(waveblue[indexblueb],
    #                                                      waveblue[indexbluer]))
    plt.clf()
    axarr=fig.subplots(2)
    
    axarr[0].plot(overwave,overflux,drawstyle='steps-mid',color='k')
    axarr[0].plot(waveblue[indexblueb:indexbluer+1], fluxblue[indexblueb:indexbluer+1],
             drawstyle='steps-mid',color='b')
    axarr[0].plot(wavered[indexredb:indexredr+1], fluxred[indexredb:indexredr+1],
             drawstyle='steps-mid',color='r')
    axarr[0].set_title('Overlap region, with inputs and combination')
    axarr[1].set_ylabel('Flux')
    axarr[1].plot(newwave,newflux,drawstyle='steps-mid',color='k')
    axarr[1].plot(newwave[indexblueb:indexbluer+1],
                  newflux[indexblueb:indexbluer+1],
                  drawstyle='steps-mid',color='r')
    plt.pause(0.01)
    check=inputter('Check plot [enter when done]: ','string',False)
    # hopout=0
    # while (hopout < 1) or (hopout > HOPSIZE):
    #     hopout=inputter('Enter hopper to store combined spectrum: ','int',False)
    # hop[hopout].wave=newwave
    # hop[hopout].flux=newflux
    # hop[hopout].var=newvar
    # hop[hopout].obname=hop[hopnum1].obname
    # hop[hopout].header=hop[hopnum1].header
    # plt.close()
    #fig.clf()
    # #plt.cla()
    # return hop
    return newwave, newflux, newvar

