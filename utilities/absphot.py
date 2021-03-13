import numpy as np
import os
import sys

class absphot(object):
    def __init__(self, iterations=5, sigma=5):

        self.iterations = iterations
        self.sigma = sigma

        self.odr_iter = 10

    def get_zeropoint(self, flux, fluxerr, mag, magerr):

        from scipy.odr import ODR, Model, Data, RealData
        def magnitude(zpt, flux):
            return(zpt - 2.5*np.log10(flux))

        zpt_guess = np.median(mag + 2.5*np.log10(flux))

        data = RealData(flux, mag, fluxerr, magerr)
        model = Model(magnitude)
        odr = ODR(data, model, [zpt_guess], maxit=self.odr_iter)

        output = odr.run()

        zpt = output.beta[0]
        zpterr = output.sd_beta[0]

        return(zpt, zpterr)

    def zpt_iteration(self, flux, fluxerr, mag, magerr):

        flux = np.array(flux)
        fluxerr = np.array(fluxerr)
        mag = np.array(mag)
        magerr = np.array(magerr)

        for i in np.arange(self.iterations):

            nobj = len(flux)
            zpt, zpterr = self.get_zeropoint(flux, fluxerr, mag, magerr)

            mag_deriv = -2.5*np.log10(flux)+zpt
            magerr_deriv = np.array(2.5/np.log(10) * fluxerr/flux)

            total_err = np.sqrt(magerr**2+magerr_deriv**2+zpterr**2)

            mask = np.abs(mag - mag_deriv) < self.sigma * total_err

            flux=flux[mask] ; fluxerr=fluxerr[mask]
            mag=mag[mask] ; magerr=magerr[mask]

            m='Iteraction {i}: {N} obj, zpt = {zpt}+/-{zpterr}, {s}-sigma '+\
                'clip to {M} obj'
            print(m.format(i=i, N=nobj, zpt='%2.4f'%zpt, zpterr='%2.4f'%zpterr,
                s=self.sigma, M=len(flux)))

        return(zpt, zpterr)




