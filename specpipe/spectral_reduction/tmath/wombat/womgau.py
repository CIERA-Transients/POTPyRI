def womgau(hop):
    """fit gaussian to line"""
    import numpy as np
    import logging
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.gauss import gauss
    from tmath.wombat.gauss_cont import gauss_cont
    from tmath.wombat.yesno import yesno
    print(' ')
    logging.info('Object is {}'.format(hop[0].obname))
    print(' ')
    print('Spectrum runs from {} to {}'.format(hop[0].wave[0],hop[0].wave[-1]))
    print(' ')
    print('This routine expects the spectrum to be in flambda units.')
    print('It also expects a linear wavelength scale.')
    print(' ')
    print('Choose general region of spectrum\n')
    nwave,nflux,mode=womwaverange(hop[0].wave,hop[0].flux,'none')
    print('\nNow pick the exact range for the fit')
    waveint,fluxint,mode=womwaverange(nwave,nflux,mode)
    indexblue=womget_element(nwave, waveint[0])
    indexred=womget_element(nwave,waveint[-1])
    if (mode == 'w'):
        done = False
        while (not done):
            print(' ')
            wavecenter=inputter('Enter approximate center of Gaussian : ','float',False)
            indexcenter=womget_element(waveint,wavecenter)
            if (indexcenter <= 0) or (wavecenter > waveint[-1]):
                print('Bad central wavelength, try again')
            else:
                done = True
    else:
        done=False
        while (not done):
            print('Mark the approximate center of the Gaussian')
            pickcent=plt.ginput(1,timeout=-1)
            indexcenter=womget_element(waveint,pickcent[0][0])
            print('\nApproximate center at {}'.format(waveint[indexcenter]))
            print('\nIs this OK?')
            answer=yesno('y')
            if (answer == 'y'):
                done=True
    weights=np.sqrt(hop[0].var[indexblue:indexred+1])
    print(' ')
    continuum=inputter_single('Do you want to fit gaussian with (c)ontinuum, or (n)o continuum? ','cn')
    if (continuum == 'c'):
        p=[fluxint[indexcenter], waveint[indexcenter],3.0,1.0,waveint[0]]
        result=curve_fit(gauss_cont,waveint,fluxint,sigma=weights,p0=p,absolute_sigma=True,full_output=True)
    else:
        p=[fluxint[indexcenter], waveint[indexcenter],3.0]
        result=curve_fit(gauss,waveint,fluxint,sigma=weights,p0=p,absolute_sigma=True,full_output=True)
    coefferr=np.sqrt(np.diag(result[1]))
    coeff=result[0]
        # make 'finer-grained' version of fit, 0.2A/pix for calculations
    wavecalc=np.arange(2*5*50*abs(coeff[2]))*0.2+coeff[1]-0.2*5*50*abs(coeff[2])
    calccenter=womget_element(wavecalc,coeff[1])
    if (continuum == 'c'):
        fluxcalc=gauss_cont(wavecalc,*coeff)
        fluxcont=wavecalc*coeff[3]+coeff[4]
        fluxgaussian=fluxcalc-fluxcont
        linecont=fluxcont[calccenter]
    else:
        fluxcalc=gauss(wavecalc,*coeff)
    
    
    deltafit=wavecalc[1]-wavecalc[0]
    calcindexblue=womget_element(wavecalc,waveint[0])
    calcindexred=womget_element(wavecalc,waveint[-1])
    sumfluxcalc=np.sum(fluxcalc[calcindexblue:calcindexred+1]*deltafit)
    sumallfluxcalc=np.sum(fluxcalc*deltafit)
    chi=(result[2]['fvec']**2).sum()
    redchi=chi/(len(waveint)-len(coeff))
    if (continuum == 'c'):
        sumfluxgaussian=np.sum(fluxgaussian[calcindexblue:calcindexred+1]*deltafit)
        sumallfluxgaussian=np.sum(fluxgaussian*deltafit)
        sumfluxcont=np.sum(fluxcont[calcindexblue:calcindexred+1]*deltafit)
        sumallfluxcont=np.sum(fluxcont*deltafit)
        # propagate uncertainty (from old version) not sure this is correct
        height_pct=coefferr[0]/coeff[0]
        sigma_pct=coefferr[2]/coeff[2]
        flux_pct=np.sqrt(height_pct**2+sigma_pct**2)
        sumfluxgaussiansig=sumfluxgaussian*flux_pct
        sumallfluxgaussiansig=sumallfluxgaussian*flux_pct
    plt.cla()
    plt.plot(nwave,nflux,drawstyle='steps-mid',color='k')
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    xmin,xmax=plt.xlim()
    ymin,ymax=plt.ylim()
    plt.plot(wavecalc,fluxcalc,drawstyle='steps-mid',color='b')
    if (continuum == 'c'):
        plt.plot(wavecalc,fluxgaussian,drawstyle='steps-mid',color='r')
        plt.plot(wavecalc,fluxcont,drawstyle='steps-mid',color='g')
    plt.plot([waveint[0],waveint[0]],[ymin,ymax],color='k',linestyle='--')
    plt.plot([waveint[-1],waveint[-1]],[ymin,ymax],color='k',linestyle='--')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    logging.info('For object {} Gaussian fit'.format(hop[0].obname))
    if (continuum == 'c'):
        print('\nData = Black, Fit = Blue, Continuum = Green, Fit-Continuum = Red\n')
    else:
        print('\nData = Black, Fit = Blue\n')
    logging.info('Height      {:16.8f}+/-{:16.8f}'.format(coeff[0],coefferr[0]))
    logging.info('Center      {:16.8f}+/-{:16.8f}'.format(coeff[1],coefferr[1]))
    logging.info('Sigma       {:16.8f}+/-{:16.8f}'.format(coeff[2],coefferr[2]))
    if (continuum == 'c'):
        logging.info('Slope       {:16.8f}+/-{:16.8f}'.format(coeff[3],coefferr[3]))
        logging.info('Y-intercept {:16.8f}+/-{:16.8f}'.format(coeff[4],coefferr[4]))
        logging.info('FWHM        {:16.8f}+/-{:16.8f}'.format(2.35482*np.abs(coeff[2]),2.35482*coefferr[2]))
        logging.info('Flux between dotted lines (Gaussian): {:16.8f}+/-{:16.8f}'.format(sumfluxgaussian, sumfluxgaussiansig))
        logging.info('EW between dotted lines (Gaussian): {:16.8f}'.format(sumfluxgaussian/linecont))
        logging.info('Flux for full (Gaussian): {:16.8f}+/-{:16.8f}'.format(sumallfluxgaussian, sumallfluxgaussiansig))
        logging.info('EW for full (Gaussian): {:16.8f}'.format(sumallfluxgaussian/linecont))
        logging.info('Continuum flux at line center: {:16.8f}'.format(linecont))
    logging.info('Chi^2: {}'.format(chi))
    logging.info('Reduced chi^2: {}'.format(redchi))
    logging.info('All fluxes might need to be scaled by 1e-15')
    print(' ')
    return hop

