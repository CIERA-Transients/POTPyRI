def womirfilters(hop):
    """calculate photometric values from spectra"""
    import numpy as np
    import logging
    from tmath.wombat.filtermag import filtermag
    from tmath.wombat.yesno import yesno
    from tmath.wombat.inputter import inputter
    from tmath.wombat.inputter_single import inputter_single
    print('NOTE:  The routine expects an f_lambda spectrum')
    print('       I will try to guess if the spectrum')
    print('       has been scaled by 1E15')
    print(' ')
    print('       Check this before believing fluxes')
    print(' ')
    print('NOTE Also:  These are the 2MASS filter curves')
    print(' ')
    flux=hop[0].flux.copy()
    if (np.mean(flux) > 0.00001):
        flux = flux *1.e-15

    filtwave=np.zeros((109,3))
    filttran=np.zeros((109,3))

    filtwave[:,0]=[1.050, 1.051, 1.062, 1.066, 1.070, 1.075, 1.078, 1.082, \
        1.084, 1.087, 1.089, 1.093, 1.096, 1.102, 1.105, 1.107, 1.109, 1.112, \
        1.116, 1.117, 1.120, 1.123, 1.128, 1.129, 1.132, 1.134, 1.138, 1.140, \
        1.143, 1.147, 1.154, 1.159, 1.164, 1.167, 1.170, 1.173, 1.175, 1.179, \
        1.182, 1.186, 1.188, 1.192, 1.195, 1.199, 1.202, 1.209, 1.216, 1.221, \
        1.227, 1.231, 1.236, 1.240, 1.244, 1.247, 1.253, 1.255, 1.258, 1.260, \
        1.265, 1.270, 1.275, 1.279, 1.286, 1.292, 1.297, 1.302, 1.305, 1.307, \
        1.310, 1.313, 1.316, 1.319, 1.323, 1.326, 1.330, 1.333, 1.334, 1.336, \
        1.339, 1.343, 1.346, 1.349, 1.353, 1.355, 1.360, 1.363, 1.370, 1.373, \
        1.377, 1.383, 1.388, 1.392, 1.395, 1.396, 1.397, 1.398, 1.400, 1.401, \
        1.402, 1.404, 1.406, 1.407, 1.410, 1.412, 1.416, 1.421, 1.426, 1.442, \
        1.450]


    filttran[:,0]=[0.0000, 0.0000, 0.0000, 0.0023, 0.0087, 0.0150, 0.0309, 0.0690, \
        0.1136, 0.1709, 0.2282, 0.2886, 0.3491, 0.4255, 0.4668, 0.5209, \
        0.5687, 0.6228, 0.6546, 0.6864, 0.7150, 0.7437, 0.7595, 0.7595, \
        0.7435, 0.7276, 0.6861, 0.6575, 0.6224, 0.5873, 0.5649, 0.5840, \
        0.6157, 0.6571, 0.6857, 0.7271, 0.7685, 0.8162, 0.8416, 0.8511, \
        0.8447, 0.8256, 0.7937, 0.7554, 0.7172, 0.6757, 0.6629, 0.6883, \
        0.7391, 0.7869, 0.8505, 0.8823, 0.8950, 0.8854, 0.8471, 0.8184, \
        0.7802, 0.7324, 0.6845, 0.6239, 0.5889, 0.5729, 0.5728, 0.5918, \
        0.6172, 0.6681, 0.6968, 0.7286, 0.7667, 0.7954, 0.8431, 0.8813, \
        0.9194, 0.9353, 0.9257, 0.9225, 0.9129, 0.8906, 0.8524, 0.8141, \
        0.7854, 0.7599, 0.7439, 0.7375, 0.7247, 0.7183, 0.7087, 0.7023, \
        0.7022, 0.7181, 0.7339, 0.7147, 0.6829, 0.6446, 0.6160, 0.5873, \
        0.5172, 0.4662, 0.3770, 0.2305, 0.1350, 0.1126, 0.0712, 0.0362, \
        0.0170, 0.0042, 0.0009, 0.0007, 0.0000]


    filtwave[0:57,1]=[1.315, 1.341, 1.368, 1.397, 1.418, 1.440, 1.462, 1.478, \
        1.486, 1.493, 1.504, 1.515, 1.528, 1.539, 1.546, 1.551, 1.556, 1.565, \
        1.572, 1.577, 1.583, 1.592, 1.597, 1.602, 1.613, 1.619, 1.628, 1.633, \
        1.642, 1.648, 1.657, 1.659, 1.671, 1.684, 1.701, 1.715, 1.727, 1.739, \
        1.746, 1.751, 1.753, 1.756, 1.764, 1.775, 1.785, 1.790, 1.796, 1.803, \
        1.810, 1.813, 1.818, 1.828, 1.835, 1.850, 1.871, 1.893, 1.914]


    filttran[0:57,1]=[0.0014, 0.0014, 0.0000, 0.0000, 0.0014, 0.0028, 0.0070, \
        0.0252, 0.0700, 0.1807, 0.3529, 0.4972, 0.6527, 0.7591, 0.8109, \
        0.8319, 0.8403, 0.8389, 0.8305, 0.8235, 0.8193, 0.8277, 0.8347, \
        0.8375, 0.8319, 0.8193, 0.8081, 0.8053, 0.8095, 0.8165, 0.8263, \
        0.8305, 0.8375, 0.8431, 0.8501, 0.8529, 0.8543, 0.8529, 0.8445, \
        0.8305, 0.8151, 0.7927, 0.7255, 0.6275, 0.5084, 0.4258, 0.3291, \
        0.2101, 0.1275, 0.0882, 0.0560, 0.0294, 0.0154, 0.0070, 0.0028, \
        0.0014, 0.0000]


    filtwave[0:76,2]=[1.900, 1.915, 1.927, 1.934, 1.939, 1.948, 1.957, 1.962, \
        1.969, 1.976, 1.981, 1.989, 1.990, 1.998, 2.008, 2.014, 2.019, 2.028, \
        2.037, 2.045, 2.061, 2.072, 2.075, 2.082, 2.089, 2.099, 2.106, 2.113, \
        2.120, 2.124, 2.138, 2.145, 2.155, 2.169, 2.176, 2.185, 2.197, 2.208, \
        2.213, 2.218, 2.232, 2.237, 2.248, 2.256, 2.260, 2.263, 2.265, 2.270, \
        2.272, 2.276, 2.277, 2.281, 2.284, 2.286, 2.291, 2.293, 2.295, 2.297, \
        2.299, 2.306, 2.311, 2.316, 2.320, 2.325, 2.328, 2.335, 2.339, 2.344, \
        2.346, 2.352, 2.361, 2.363, 2.370, 2.375, 2.384, 2.399]

    filttran[0:76,2]=[0.0000, 0.0013, 0.0027, 0.0040, 0.0082, 0.0153, 0.0293, \
        0.0462, 0.0743, 0.1222, 0.1714, 0.2672, 0.3517, 0.4263, 0.6262, \
        0.6797, 0.7487, 0.7853, 0.8120, 0.8303, 0.8485, 0.8513, 0.8583, \
        0.8597, 0.8667, 0.8751, 0.8765, 0.8835, 0.8891, 0.8863, 0.8848, \
        0.8819, 0.8805, 0.8748, 0.8804, 0.8818, 0.8902, 0.8986, 0.9014, \
        0.8999, 0.8999, 0.8956, 0.8913, 0.8969, 0.8997, 0.8997, 0.9053, \
        0.9109, 0.9166, 0.9109, 0.9025, 0.8870, 0.8686, 0.8433, 0.7714, \
        0.7292, 0.6650, 0.5950, 0.5333, 0.4094, 0.3108, 0.2234, 0.1544, \
        0.1234, 0.0896, 0.0599, 0.0416, 0.0320, 0.0300, 0.0162, 0.0063, \
        0.0007, 0.0034, 0.0020, 0.0006, 0.0000]

    filtwave=filtwave*10000.0
    
    filtsize = [109, 57, 76]
    #		Holds the filter zero-points as determined from
#		Vega model by Dreiling & Bell (ApJ, 241,736, 1980)
#
#		B	6.268e-9   erg cm-2 s-1 A-1
#		V	3.604e-9
#		R	2.161e-9
#		I	1.126e-9
#
#		The following zero-points are from Lamla
#		(Landolt-Boernstein Vol. 2b, eds. K. Schaifer & 
#		H.H. Voigt, Berlin: Springer, p. 73, 1982 QC61.L332)
#
#		U	4.22e-9   erg cm-2 s-1 A-1
#
#		J	3.1e-10
#		H	1.2e-10
#		K	3.9e-11
#
#               U        B          V        R         I

    zeropoint = [3.1e-10, 1.2e-10,3.9e-11]

    mag=np.zeros(3)
    filtflux=mag.copy()
    coverage=mag.copy()
    efflambda=mag.copy()
    totflux=mag.copy()
    filtername = ['J', 'H', 'K']
    for i,_ in enumerate(filtername):
        filtw=filtwave[0:filtsize[i],i]
        filtt=filttran[0:filtsize[i],i]
        mag[i], filtflux[i], coverage[i], efflambda[i], totflux[i]= \
              filtermag(hop[0].wave,flux, filtw, filtt, \
              zeropoint[i])                                                            
    logging.info('For object {}'.format(hop[0].obname))
    logging.info('Filter magnitude  Flux(erg/s/cm^2/A) Flux(erg/s/cm^2)  Coverage(%)  Eff. Lambda')
    for i in range(0,3):
        if (mag[i] > 99):
            logging.info('  {:1s}        FILTER AND SPECTRUM DO NOT OVERLAP'.format(filtername[i]))
        else:
            logging.info('  {:1s}     {:6.3f}      {:10.4e}        {:10.4e}         {:5.1f}         {:7.1f}'.format(filtername[i],mag[i],filtflux[i],totflux[i],coverage[i]*100.,efflambda[i]))


            
    print(' ')
    logging.info('Colors')
    colortab=[[0,1],[1,2]]
    for i in range(0,2):
        if (mag[colortab[i][0]] > 99) or (mag[colortab[i][1]] > 99):
            logging.info('{}-{}    ONE OR BOTH FILTERS DO NOT OVERLAP SPECTRUM'.format(filtername[colortab[i][0]],filtername[colortab[i][1]]))
        else:
            logging.info('{:1s}-{:1s}    {:12.4f}'.format(filtername[colortab[i][0]],filtername[colortab[i][1]],mag[colortab[i][0]]-mag[colortab[i][1]]))


    print('\nWould you like to scale the spectrum to match photometry?\n')
    answer=yesno('n')
    if (answer == 'y'):
        print('\nWhich filter do you have?')
        scalefilt=inputter_single_mix('J/H/K: ','JHK')
        filtindex=filtername.index(scalefilt)
        scalemag=inputter('Enter your value for filter {}: '.format(filtername[filtindex]),'float',False)
        print(' ')
        logging.info('Scaling {} from {}={:6.3f} to {}={}'.format(hop[0].obname,filtername[filtindex],mag[filtindex],filtername[filtindex],scalemag))
        logging.info('Multiplying by {:.3f}'.format(10**(0.4*(mag[filtindex]-scalemag))))
        hop[0].flux=hop[0].flux*10**(0.4*(mag[filtindex]-scalemag))
    

    return hop

