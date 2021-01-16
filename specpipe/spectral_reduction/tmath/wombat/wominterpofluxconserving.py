import numpy as np
import math
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import scipy.optimize as opt
def interpo_flux_conserving(wave, flux, ivar, waveb, waver, dw=2, testing=False):
    """This is a an interpolation algorithm that does trapezoidal integration
        to conserve flux. The variance is then propagated correctly. Since 
        interpolation always introduces correlation in the variance spectrum,
        we ignore  this correlation byscale the original variance spectrum 
        to the new variance spectrum.
    """
    var = 1./ivar
    pixel_scale = np.median(np.diff(wave))
    # print pixel_scale

    wave_min = 999
    wave_max = 12001
    wavelength_mids = np.arange(math.ceil(wave_min), math.floor(wave_max),
                           dtype=int, step=dw)
    lower = wave[0]
    upper = wave[-1]

    good_data = np.where((wave >= waveb) & (wave <= waver))
    influx = inter.splrep(wave[good_data], flux[good_data])
    invar = inter.splrep(wave[good_data], var[good_data])

    low_int = math.ceil(lower)
    up_int = math.ceil(upper)
    wave_final = []
    flux_final = []
    var_final = []
    for mid_point in wavelength_mids:
        inter_flux_mid_left = float(inter.splev(mid_point, influx, ext = 3))
        inter_var_mid_left = float(inter.splev(mid_point, invar, ext = 3))
        inter_flux_mid_right = float(inter.splev(mid_point+dw, influx, ext = 3))
        inter_var_mid_right = float(inter.splev(mid_point+dw, invar, ext = 3))
        waves_to_sum = np.where((wave > mid_point) & (wave < mid_point+dw))[0]
        wave_final.append((mid_point + mid_point +dw)/2.)
        if (mid_point >= lower and (mid_point+dw) <= upper):
            waves = [mid_point]
            fluxes = [inter_flux_mid_left]
            variances = [inter_var_mid_left]
            for i, w in enumerate(waves_to_sum):
                waves.append(wave[w])
                fluxes.append(flux[w])
                variances.append(var[w])
            waves.append(mid_point+dw)
            fluxes.append(inter_flux_mid_right)
            variances.append(inter_var_mid_right)
            new_point = np.trapz(fluxes, x=waves)
            flux_final.append(new_point)
            diffs = np.diff(waves)
            var_tot = 0.
            for i in range(len(variances)-1):
                v1 = variances[i]
                v2 = variances[i+1]
                var_tot = var_tot + (diffs[i]**2.)*(v1+v2)
            var_tot = var_tot*.25
            var_final.append(var_tot)
        else:
            flux_final.append(0.)
            var_final.append(0.)

    inter_var = inter.splev(wave_final, invar, ext = 3)
    s = scale_together(var_final, inter_var)[0]
    scale = 1./s
    var_uncorrelated = scale*inter_var
    ivar_final = 1./var_uncorrelated

    wave_final = np.asarray(wave_final)
    flux_final = np.asarray(flux_final)
    ivar_final = np.asarray(ivar_final)
    if testing:
        var_final = np.asarray(var_final)

    data = np.where((wave_final > lower+dw) & (wave_final < upper-dw))# padding to avoid interpolation errors near edges
    # flux_final[missing_data] = float('NaN')
    # ivar_final[missing_data] = float('NaN')
    # if testing:
    #     var_final[missing_data] = float('NaN')

    output = np.array([wave_final[data], flux_final[data], ivar_final[data]])

    if testing:
        interp_wave = output[0,:]
        interp_flux = output[1,:]
        interp_ivar = output[2,:]
        # print scale
        plt.plot(wave,var)
        plt.plot(interp_wave,var_final)
        plt.plot(interp_wave,1./interp_ivar)
        plt.xlim([7000,7100])
        plt.ylim([-.05e-32,2.e-32])
        plt.show()

    return output

def scale_together(data, comp):
    """Finds the scale factor the minimizes the difference between two 
        spectra (data and comp)
    """
    scales = []
    guess = 1.
    s = opt.minimize(sq_residuals_in_range, guess, args = (data, comp), 
                 method = 'Nelder-Mead').x
    return s

def sq_residuals_in_range(s, data, comp):
    """Calculates the sum of the square residuals between arrays data and comp
    """
    data = s*data
    res = data - comp
    sq_res = res*res
    return np.sum(sq_res)