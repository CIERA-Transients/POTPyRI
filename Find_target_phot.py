# UKIRT pipeline function
# Owen Eskandari
# 6/16/20

# Main function performs photometry on the target in a stack given its position, zp, and zp_err

# Imports
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.wcs as wcs
from photutils import make_source_mask, CircularAperture, CircularAnnulus, aperture_photometry, \
    Background2D, MedianBackground
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.utils import calc_total_error
from scipy.optimize import curve_fit
from astropy.visualization import simple_norm
import random

# Function imports
from Gaussians import fix_x0, twoD_Gaussian
from utilities import util

# Radial profile function; currently only used in find_target_phot; compact since it's outlined earlier in zp function
def radial_profile(data, x, y, step_size, fwhm, rad):

    '''

    Function used to create a plot of the radial profile of a source at the coordinates (``x``, ``y``).
    Performs aperture photometry on concentric circles of radii ``step_size``.
    This function does not show the plot but does everything up to that point.

    :param data: CCDData
        Data of the image.

    :param x: float
        Pixel coordinate.

    :param y: float
        Pixel coordinate.

    :param step_size: float
        Radial size of each ring (concentric circles). Anything less than ``0.7`` is too pixelated and does not
        represent the data accurately

    :param fwhm: float
        FWHM of the image.

    :param rad: float
        Radius to extent to (actually is ``rad`` * ``fwhm``).

    :return: nothing
        Makes the plot but does not return anything.
    '''

    radii = np.arange(step_size, rad*fwhm, step_size)  # Radii of apertures
    x_data = np.arange(0, rad*fwhm, step_size)  # Steps (starting at zero, to 2.5 * fwhm, by step_size)
    apers_area = [np.pi * (step_size ** 2)]  # Circle area = pi*(r**2)
    for r in radii:
        apers_area.append(np.pi * ((r + step_size) ** 2 - r ** 2))  # Annulus area = pi*(r_out**2 - r_in**2)
    apertures = [CircularAperture((x, y), r=step_size)]  # For circle aperture around center
    for r in radii:
        apertures.append(CircularAnnulus((x, y), r_in=r, r_out=r + step_size))  # Annuli apertures
    rad_prof_large = aperture_photometry(data, apertures)  # Radial profile photometry of new fake source
    sums = [rad_prof_large[0][k] / apers_area[k - 3] for k in range(3, len(rad_prof_large[0]))]  # For radial profile
    w_fit = x_data
    f_fit = sums
    g, _ = curve_fit(fix_x0, w_fit, f_fit)
    plt.figure(1)
    # plt.plot(x_data, sums, 'o-', alpha=0.5, lw=2, markersize=4)
    plt.plot(x_data, sums, 'o-', alpha=0.5, lw=2, markersize=4)  # Plot the radial profiles
    plt.plot(x_data, fix_x0(x_data, *g), '^-', alpha=0.5, lw=2, markersize=4)
    plt.title("Radial profile of target")
    plt.xlabel("Radial distance from centroid of source (pixels)")
    plt.ylabel("Count (adu per second per pixel)")


# This function finds the target's photometry and thus its magnitude, but in the case that the target is not found,
# the function finds the limiting magnitude using a false source at the target's location
def find_target_phot(stack, fil, fwhm, zp, zp_err, pixscale, show_phot=False, log=None, log2=None,
                     ra=None, dec=None, x=None, y=None):

    '''

    Function to perform photometry on a target at a specific RA and Dec. Usually used after `find_zero_point`
    (this finds the ``zp``, ``zp_err``, and ``fwhm``).
    During this process, the user has the options to change the target's location and radii used in
    aperture photometry (if ``show_plots`` is ``True``, user can see cutouts of target and its radial profile).
    Attempts to find the magnitude of the target using aperture photometry.


    If magnitude cannot be found at a 3 sigma confidence interval, limiting magnitude is found at the same position.
    Limiting magnitude is found by adding fake sources of increasing magnitude at the target's position until it is
    more than 3 sigma above the background.

    :param stack: str
        Name of the stack to use (include path to the stack).

    :param fil: str
        Name of the filter of the stack.

    :param fwhm: float
        FWHM of image.

    :param zp: float
        Zero point of image in AB mag.

    :param zp_err: float
        Error on the zero point of image.

    :param show_phot: boolean, optional
        Option to see the cutout of the source with apertures along with its radial profile.
        Default is ``False``.

    :param log: log, optional
        In-depth log.
        If no log is inputted, information is printed out instead of being written to ``log``.
        Default is ``None``.

    :param log2: log, optional
        Overview log.
        If no log is inputted, information is printed out instead of being written to ``log2``.
        Default is ``None``.

    :param ra: float, optional
        RA coordinate of target.
        Default is ``None``.

    :param dec: float, optional
        Dec coordinate of target.
        Default is ``None``.

    :param x: float, optional
        X pixel coordinate (not RA).
        Default is ``None``.

    :param y: float, optional
        Y pixel coordiante (not Dec).
        Default is ``None``.

    Returns
    -------

    :return: float, float
        Three options:
            * Returns the magnitude of the target (in AB) and its error.
            * Returns tne limiting magnitude at the target's location (in AB) and its error.
            * Returns ``None``, ``None`` if a fake source of amplitude 5 (adu / sec / pixel) is not three
                sigma above the background.

    '''

    with fits.open(stack) as hdr:
        header, data = hdr[0].header, hdr[0].data
        w = wcs.WCS(header)  # Parse the WCS keywords in the primary HDU

    if ra is not None:
        coords = np.array([[ra, dec]], np.float)  # Array of coordinates: [[RA, Dec]] in [deg, deg]
        pixel_coords = w.wcs_world2pix(coords, 1)[0]     # Find the pixel coordinates in the image
        x = pixel_coords[0]
        y = pixel_coords[1]
    elif x is not None:
        x-=1
        y-=1
        coords = np.array([[x, y]], np.float)  # Array of coordinates: [[RA, Dec]] in [deg, deg]
        pixel_coords = w.wcs_pix2world(coords, 1)[0]     # Find the pixel coordinates in the image
        ra = pixel_coords[0]
        dec = pixel_coords[1]
    else:
        print("No coordinates inputted")

    d_x_y = 125     # Pixels above and below the target when showing image...125 corresponds to a 50" x 50" cutout
    if log is not None:
        log.info("Pixel coordinates: (%.3f, %.3f)" % (x, y))
    else:
        print("Pixel coordinates: (%.3f, %.3f)" % (x, y))
    if input('Would you like to choose the cutout size? Default is 50"x50". ') == "yes":
        try:
            d_x_y = int(input("Choose the radius, in arcsec: "))*2.5
        except TypeError:
            pass

    sigma_clip = SigmaClip(sigma=3)  # This section is to find the error on the background
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (120, 120), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    error = calc_total_error(data, bkg.background_rms, 1)  # Used later as well for limiting magnitude

    if log is not None:
        log.info("You will now get to pick the position you would like. "
                 "Decide what pixel coordinates you would like to use")
        log.info("Blue circle in image is the distance the radial profile extends to")
    else:
        print("You will now get to pick the position you would like. "
              "Decide what pixel coordinates you would like to use")
        print("Blue circle in image is the distance the radial profile extends to")
    correct_position = "no"
    while correct_position == "no":
        if show_phot:
            radial_profile(data, x, y, 0.2, fwhm, rad=6)        # Radial profile
            plt.axvline(6*fwhm)

            print("Double click to set the new center. Do nothing if you are ok with current coordinates.")
            fig = plt.figure(2)
            new_coords = []
            fig.canvas.mpl_connect('button_press_event', lambda event: util.onclick(event, new_coords=new_coords))
            ap_in = CircularAperture((x, y), 1)  # Aperture to perform photometry
            ap_out = CircularAperture((x, y), 6*fwhm)  # Aperture to perform photometry
            ap_in.plot(color='r', lw=1)      # Plot of target with apertures and annulus
            ap_out.plot(color='b', lw=2)
            norm = simple_norm(data, 'sqrt', percent=99)
            plt.imshow(data, norm=norm)
            plt.xlim(int(x) - d_x_y, int(x) + d_x_y)
            plt.ylim(int(y) - d_x_y, int(y) + d_x_y)
            plt.colorbar()
            plt.show()

            if len(new_coords) != 0:
                x, y = new_coords[0]
                real_coords = w.wcs_pix2world(np.array([[x, y]], np.float), 1)[0]  # Find the pixel coords in the image
                ra, dec = real_coords[0], real_coords[1]

                radial_profile(data, x, y, 0.2, fwhm, rad=6)        # Radial profile
                plt.axvline(6*fwhm)

                print("Showing new position.")
                fig = plt.figure(2)
                ap_in = CircularAperture((x, y), 1)  # Aperture to perform photometry
                ap_out = CircularAperture((x, y), 6*fwhm)  # Aperture to perform photometry
                ap_in.plot(color='r', lw=1)      # Plot of target with apertures and annulus
                ap_out.plot(color='b', lw=2)
                norm = simple_norm(data, 'sqrt', percent=99)
                plt.imshow(data, norm=norm)
                plt.xlim(int(x) - d_x_y, int(x) + d_x_y)
                plt.ylim(int(y) - d_x_y, int(y) + d_x_y)
                plt.colorbar()
                plt.show()
        
        if input("Would you like to use a centroid? Type 'yes' or 'no': ") == "yes":
            d_x_y = 50
            d = data[int(y)-d_x_y:int(y)+d_x_y, int(x)-d_x_y:int(x)+d_x_y]      # Cutout of source; upside down from normal
            x_mesh, y_mesh = np.meshgrid(np.linspace(0, np.shape(d)[1] - 1, np.shape(d)[1]),
                                np.linspace(0, np.shape(d)[0] - 1, np.shape(d)[0]))
            popt, pcov = curve_fit(twoD_Gaussian, (x_mesh, y_mesh), d.ravel(),
                                p0=(np.max(d), 50, 50, fwhm/2.35482, fwhm/2.35482, 0, 0))
            x+=popt[1]-d_x_y
            y+=popt[2]-d_x_y
            if log:
                log.info('Centroid calculated position: (%.3f, %.3f)'%(x, y))
            else:
                print('Centroid calculated position: (%.3f, %.3f)'%(x, y))
            
            real_coords = w.wcs_pix2world(np.array([[x, y]], np.float), 1)[0]  # Find the pixel coords in the image
            ra, dec = real_coords[0], real_coords[1]

            radial_profile(data, x, y, 0.2, fwhm, rad=6)        # Radial profile
            plt.axvline(6*fwhm)

            print("Showing centroid position.")
            fig = plt.figure(2)
            ap_in = CircularAperture((x, y), 1)  # Aperture to perform photometry
            ap_out = CircularAperture((x, y), 6*fwhm)  # Aperture to perform photometry
            ap_in.plot(color='r', lw=1)      # Plot of target with apertures and annulus
            ap_out.plot(color='b', lw=2)
            norm = simple_norm(data, 'sqrt', percent=99)
            plt.imshow(data, norm=norm)
            plt.xlim(int(x) - d_x_y, int(x) + d_x_y)
            plt.ylim(int(y) - d_x_y, int(y) + d_x_y)
            plt.colorbar()
            plt.show()

        if input("Are you ok with this position? Type 'yes' or 'no': ") != "yes":
            pass
        else:
            correct_position = "yes"
            if log is not None:
                log.info("Coordinates chosen: (%.3f, %.3f) at RA = %.5f and Dec = %.5f" % (x, y, ra, dec))
            else:
                print("Coordinates chosen: (%.3f, %.3f) at RA = %.5f and Dec = %.5f" % (x, y, ra, dec))
            if log2 is not None:
                log2.info("Final coordinates: (%.3f, %.3f) at RA = %.5f and Dec = %.5f" % (x, y, ra, dec))

    print("You will now get to choose the radii for the circular aperture and the r_in and r_out of the annulus")
    correct_radii = "no"
    while correct_radii == "no":
        if log is not None:
            log.info("Automatic radii picked by comparing to FWHM of field: rad = %.3f, r_in = %.3f, r_out = %.3f" %
                     (2.5 * fwhm, 2.5 * fwhm, 4.5 * fwhm))
        else:
            print("Automatic radii picked by comparing to FWHM of field: rad = %.3f, r_in = %.3f, r_out = %.3f" %
                  (2.5 * fwhm, 2.5 * fwhm, 4.5 * fwhm))
        if input("Would you like to use these radii? Type 'yes or 'no': ") == "yes":
            rad, r_in, r_out = 2.5*fwhm, 2.5*fwhm, 4.5*fwhm
        else:
            if log:
                log.info("FWHM = %.3f pixels" % fwhm)
            rad = float(input("Pick a radius (in pixels) for the circular aperture: "))
            r_in = float(input("Pick an inner radius (in pixels) for the background annulus: "))
            r_out = float(input("Pick an outer radius (in pixels) for the background annulus: "))
            if r_in >= r_out:
                r_out = r_in + 1
                if log is not None:
                    log.info("You selected an invalid r_out value. Automatically set to %s" % r_out)
                else:
                    print("You selected an invalid r_out value. Automatically set to %s" % r_out)

        aperture = CircularAperture((x, y), rad)  # Aperture to perform photometry
        annulus_aperture = CircularAnnulus((x, y), r_in=r_in, r_out=r_out)       # Annulus of target
        annulus_mask = annulus_aperture.to_mask(method='center')  # Create masks to highlight pixels in annuli

        annulus_data = annulus_mask.multiply(data)
        annulus_data_1d = annulus_data[annulus_mask.data > 0]
        _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)

        phot_table = aperture_photometry(data, aperture, error=error)  # Photometry, error accounted for
        bkg_aper = median_sigclip * aperture.area

        if show_phot:
            # Radial profile out to 6 * fwhm (large just to be safe)
            radial_profile(data, x, y, 0.2, fwhm, rad=6)
            plt.axvline(rad, c='r')
            plt.axvline(r_in)
            plt.axvline(r_out)
            plt.figure(2)
            aperture.plot(color='white', lw=2)      # Plot of target with apertures and annulus
            annulus_aperture.plot(color='b', lw=2)
            norm = simple_norm(data, 'sqrt', percent=99)
            plt.imshow(data, norm=norm)
            plt.xlim(int(x) - d_x_y, int(x) + d_x_y)
            plt.ylim(int(y) - d_x_y, int(y) + d_x_y)
            plt.colorbar()
            plt.show()

        correct_radii = input("Are you ok with the previously selected radii? Type 'yes' or 'no': ")
    if log is not None:
        log.info("Radii chosen in pixels: %.3f, %.3f, %.3f" % (rad, r_in, r_out))
    else:
        print("Radii chosen in pixels: %.3f, %.3f, %.3f" % (rad, r_in, r_out))
    if log2 is not None:
        log2.info("Final radii: %.3f, %.3f, %.3f" % (rad, r_in, r_out))

    three_sigma = 3  # Change according to what sigma you'd like to use
    if (phot_table['aperture_sum']) - bkg_aper > 0:
        mag = -2.5 * np.log10((phot_table['aperture_sum']) - bkg_aper)  # Instrumental magnitude
        mag_err = 2.5 / np.log(10.) * (phot_table['aperture_sum_err'])/(phot_table['aperture_sum'])
        mag += zp  # Apparent magnitude (Vega)
        mag_err = np.sqrt(mag_err ** 2 + zp_err ** 2)
        if log is not None:
            log.info("Magnitude of target = %.3f +/- %.3f AB mag" % (mag, mag_err))  # Final mag
        if log2 is not None:
            log2.info("Magnitude of target = %.3f +/- %.3f AB mag" % (mag, mag_err))  # Final mag
        else:
            print("Magnitude of target = %.3f +/- %.3f AB mag" % (mag, mag_err))  # Final magnitude

        if np.abs(mag_err) < 1/three_sigma:     # Good! No limiting magnitude needed
            return float(mag), float(mag_err)
        if log is not None:
            log.info("Target was not %s sigma above the background." % three_sigma)
        else:
            print("Target was not %s sigma above the background." % three_sigma)
        if log2 is not None:
            log2.info("Target not found")

    else:
        if log is not None:
            log.info("Target's aperture sum was less than background aperture sum. Target was not found")
        else:
            print("Target's aperture sum was less than background aperture sum. Target was not found")
        if log2 is not None:
            log2.info("Target not found")

    # Finding the limiting magnitude in the area 30" around where the target is supposed to be located
    if log is not None:
        log.info("Now finding the limiting magnitude to %s sigma." % three_sigma)
    else:
        print("Now finding the limiting magnitude to %s sigma." % three_sigma)
    d_x_y = int(np.round(30/pixscale,-1)/2)
    d = data[int(y)-d_x_y:int(y)+d_x_y, int(x)-d_x_y:int(x)+d_x_y]       # Cutout of source; upside down from normal
    d_error = error[int(y)-d_x_y:int(y)+d_x_y, int(x)-d_x_y:int(x)+d_x_y]
    x_mesh, y_mesh = np.meshgrid(np.linspace(0, np.shape(d)[1] - 1, np.shape(d)[1]),
                                 np.linspace(0, np.shape(d)[0] - 1, np.shape(d)[0]))
    
    source_mask = make_source_mask(d, nsigma=3, npixels=5)
    ind = np.ma.nonzero(~source_mask)
    inds = [(x,y) for x, y in zip(ind[0],ind[1]) if x > 2*fwhm and y > 2*fwhm and x < 2*d_x_y-2*fwhm and y < 2*d_x_y-2*fwhm]
    mag_list = []
    mag_err_list = []
    for i in range(50):
        xy_fake = random.choice(inds)
        gauss_data = twoD_Gaussian((x_mesh, y_mesh), 3*error[xy_fake[1],xy_fake[0]], xy_fake[0],xy_fake[1], fwhm/2.35482, fwhm/2.35482, 0, 0)
        d += gauss_data.reshape(np.shape(d))
        aperture = CircularAperture((xy_fake[0],xy_fake[1]), 2.5*fwhm)
        annulus_aperture = CircularAnnulus((xy_fake[0],xy_fake[1]), r_in=2.5*fwhm, r_out=3.5*fwhm)
        annulus_mask = annulus_aperture.to_mask(method='center')
        annulus_data = annulus_mask.multiply(d)
        annulus_data_1d = annulus_data[annulus_mask.data > 0]
        _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)
        phot_table = aperture_photometry(d, aperture, error=d_error)
        bkg_aper = median_sigclip * aperture.area
        if (phot_table['aperture_sum']) - bkg_aper > 0:
            mag = -2.5 * np.log10((phot_table['aperture_sum']) - bkg_aper)
            mag_err = 2.5 / np.log(10.) * (phot_table['aperture_sum_err'])/(phot_table['aperture_sum'])
            mag_list.append(mag[0])
            mag_err_list.append(mag_err[0])
        d -= gauss_data.reshape(np.shape(d))
    mag_list = [mi for i,mi in enumerate(mag_list) if mag_err_list[i] > 0.2 and mag_err_list[i] < 0.5]
    mag_err_list = [mi for i,mi in enumerate(mag_err_list) if mag_err_list[i] > 0.2 and mag_err_list[i] < 0.5]
    mag_list = [x for _, x in sorted(zip(mag_err_list, mag_list))]
    mag_err_list = sorted(mag_err_list)
    limit3 = np.interp(0.333,mag_err_list,mag_list)+zp
    if log is not None:
        log.info("Limiting magnitude = %.3f +/- %.3f AB mag" % (limit3, 0.333))
    else:
        print("Limiting magnitude = %.3f +/- %.3f AB mag" % (limit3, 0.333))
    return float(limit3), float(0.333)