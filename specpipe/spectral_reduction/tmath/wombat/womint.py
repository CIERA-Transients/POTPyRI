def womint(hop):
    """calculate intensity in given wavelength range"""
    import logging
    import numpy as np
    from tmath.wombat.womwaverange import womwaverange
    from tmath.wombat.womget_element import womget_element
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
    print('\nNow pick the exact range for the intensity calculation')
    wavint,fluxint,mode=womwaverange(nwave,nflux,mode)
    indexblue=womget_element(nwave, wavint[0])
    indexred=womget_element(nwave,wavint[-1])
    wdelt=nwave[1]-nwave[0]
    lineflux=np.sum(nflux[indexblue:indexred+1])*wdelt
    linefluxin=np.sum(nflux[indexblue+1:indexred])*wdelt
    linefluxout=np.sum(nflux[indexblue-1:indexred+2])*wdelt
    print(' ')
    logging.info('FWZI (approximate): {}'.format(nwave[indexred]-nwave[indexblue]))
    logging.info('Line flux (ergs/sec/cm^2): {}'.format(lineflux))
    logging.info('Line flux one pixel in: {}'.format(linefluxin))
    logging.info('Line flux one pixel out: {}'.format(linefluxout))
    logging.info('Note that flux may need to be scaled by 1e-15')
    logging.info('Average difference (between line flux and one pixel')
    avgdiff=(np.abs(linefluxin - lineflux) + np.abs(linefluxout-lineflux))/2.0
    logging.info('in or out): {}'.format(avgdiff))
    logging.info('As a percentage of line flux: {}'.format(100.0*avgdiff/lineflux))
    return hop


