# UKIRT pipeline function
# Owen Eskandari
# 6/16/20

# Different Gaussian approximation functions

# Imports
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt


# Fits a Gaussian to a radial profile
def gaus(x, a, x0, sigma, y0):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+y0


# Gaussian but fix x0 = 0
def fix_x0(x, a, sigma, y0):
    return a*np.exp(-(x-np.float64(0))**2/(2*sigma**2))+y0


# Two dimensional Gaussian
# Used for determining if the stack was correctly made and for adding in a fake source
def twoD_Gaussian(data, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = data
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2)))
    return g.ravel()


# Shows a 2D cutout of a source with Gaussian contours overlaid
def cutouts_2d_gaussian(stack, coords, d_x_y, fwhm, show_cutouts=False, log=None):

    '''

    Function to show a 2D cutout of a source with Gaussian contours overlaid. Shows the cutouts if
    ``show_cutouts`` is ``True``.
    Returns the ratios of y to x sigma values for the Gaussian of each source (a ratio of ``1`` would be circular).

    :param stack: str
        Name of the stack to use (include path to the stack).

    :param coords: python list of tuples
        List of pixel coordinates of sources to use.

    :param d_x_y: int
        Width of cutout, in pixels.

    :param fwhm: float
        FWHM of the image.

    :param show_cutouts: boolean, optional
        Option to visually see the cutouts (set to ``True`` to see cutouts).
        Default is ``False``.

    :param log: log, optional
        In-depth log.
        If no log is inputted, information is printed out instead of being written to ``log``.
        Default is ``None``.

    :return: list of floats
        Returns the list of the y to x sigma ratios for each source. A ratio of 1 would be perfectly round.

    '''
    with fits.open(stack) as hdr:
        data = hdr[0].data

    field_ratios = []
    for i, j in enumerate(coords):
        d = data[int(j[1])-d_x_y:int(j[1])+d_x_y, int(j[0])-d_x_y:int(j[0])+d_x_y]      # Cutout of source
        x, y = np.meshgrid(np.linspace(0, np.shape(d)[1] - 1, np.shape(d)[1]),
                           np.linspace(0, np.shape(d)[0] - 1, np.shape(d)[0]))

        initial_guess = (100, d_x_y, d_x_y, fwhm/2.35482, fwhm/2.35482, 0, 0)     # Initial guess is VERY IMPORTANT
        popt, pcov = curve_fit(twoD_Gaussian, (x, y), d.ravel(), p0=initial_guess)      # Fit of 2D Gaussian
        field_ratios.append(popt[4] / popt[3])

        if i == 0 and show_cutouts:     # Show an example (not random because I want to show the same example each time)
            if log is not None:
                log.info('Sample: x sigma = %.3f y sigma = %.3f ratio = %.3f' % (popt[3], popt[4], popt[4] / popt[3]))
            else:
                print('Sample: x sigma = %.3f y sigma = %.3f ratio = %.3f' % (popt[3], popt[4], popt[4] / popt[3]))
            fitted_data = twoD_Gaussian((x, y), *popt)
            contours = np.arange(np.min(d), np.max(d), 30.)
            norm = simple_norm(data, 'sqrt', percent=99)
            plt.imshow(d, norm=norm)
            plt.colorbar()
            plt.contour(x, y, fitted_data.reshape(np.shape(d)), contours, colors='w')
            plt.show()

    return field_ratios
