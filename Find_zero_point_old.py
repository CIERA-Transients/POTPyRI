# Two functions:
# Zero point: Function that finds the zero point of a stack; returns the zero point, the zero point error, and the FWHM
# Revisit stack: Function that revisits images that make up a stack; determines if stack has been combined correctly
# It uses isophotes of "standard" stars to check to see how the individual images look and in relation to the stack

# TODO: ***NOTE*** This uses Irsa, not Vizier, and also does not have functionality for Y-band

# Imports (a lot)
import numpy as np
import astropy.units.astrophys as u
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import astropy.wcs as wcs
from photutils import CircularAperture, CircularAnnulus, aperture_photometry, \
    Background2D, MedianBackground, DAOStarFinder
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.utils import calc_total_error
from astroquery.irsa import Irsa
import ccdproc
from scipy.stats import sigmaclip
from scipy.optimize import curve_fit
from astropy.visualization import simple_norm
import glob

# Function imports
from Gaussians import fix_x0, twoD_Gaussian, cutouts_2d_gaussian

from astropy.utils.exceptions import AstropyWarning
import warnings

warnings.simplefilter('ignore', category=AstropyWarning)        # Removes deprecation warnings


# Finds the zero point of a singular stack (need to run for each filter)
# Inputs: stack (filename as a string), filter as a string, extension as an int, path to the data
# 'show' parameters show various plots, r_annulus_in/out are x*fwhm value, log is log
def find_zero_point(stack, fil, ext, data_path, show_fwhm=False, show_bkg=False, show_contours=False,
                    r_annulus_in=3.5, r_annulus_out=4.5, log=None):
    if log is not None:
        log.info("Computing FWHM of sources in the stack %s" % stack)
    else:
        print("Computing FWHM of sources in the stack %s" % stack)

    with fits.open(stack) as hdr:
        header, data = hdr[0].header, hdr[0].data
    w = wcs.WCS(header)

    # Center in pixel coordinates, change to real coordinates
    real_coords = w.wcs_pix2world(np.array([[header['NAXIS1']/2, header['NAXIS2']/2]], np.float), 1)
    coords = SkyCoord(real_coords[0][0]*u.deg, real_coords[0][1]*u.deg, frame='fk5')

    # Information about the sources from 2MASS found in the extension (radially, not full area)
    two_mass = Irsa.query_region(coords, catalog='fp_psc', spatial='Cone', radius=7*u.arcmin)

    # Find the sources above 10 std in the image
    _, median, std = sigma_clipped_stats(data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=7, threshold=10. * std)  # Looking for best sources, so raise threshold way up
    dao_sources = daofind(np.asarray(data))     # DAOStarFinder sources == dao_sources
    if log is not None:
        log.info("%s objects found in 2MASS and %s sources above 10 std found in the stack" %
                 (len(two_mass), len(dao_sources)))
        log.info("DAOStarFinder found %s sources that are 10 std above the background" % len(dao_sources))
    else:
        print("%s objects found in 2MASS and %s sources above 10 std found in the stack" %
              (len(two_mass), len(dao_sources)))
        print("DAOStarFinder found %s sources that are 10 std above the background" % len(dao_sources))

    ra_dao = []         # Change pixel coordinates to real coordinates for DAOStarFinder sources
    dec_dao = []
    for i, j in enumerate(dao_sources):
        pixel_coords = np.array([[j['xcentroid'], j['ycentroid']]], np.float)
        dao_coords = w.wcs_pix2world(pixel_coords, 1)
        ra_dao.append(dao_coords[0][0])
        dec_dao.append(dao_coords[0][1])

    # Coordinates need to be in real, not pixel
    coords_dao = SkyCoord(ra_dao*u.deg, dec_dao*u.deg, frame='fk5')      # Coordinates from sources
    coords_two_mass = SkyCoord(two_mass['ra'], two_mass['dec'], frame='fk5')    # Coordinates from 2MASS catalog

    # Match each 2MASS sources to the closest dao_source
    idx, d2, d3 = coords_two_mass.match_to_catalog_sky(coords_dao)      # Match smaller to larger
    if log is not None:
        log.info("%s objects have been matched" % len(idx))
    else:
        print("%s objects have been matched" % len(idx))

    # Next step: find the radial profiles from intermediate-brightness sources and then find the fwhm of the field
    step_size = 0.7
    fwhm_orig = 5

    coords = [w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0] for j in idx]    # Positions of DAO

    radii = np.arange(step_size, 2.5*fwhm_orig, step_size)       # Radii of apertures
    x_data = np.arange(0, 2.5*fwhm_orig, step_size)              # Steps (starting at zero, to 2.5 * fwhm, by step_size)

    apers_area = [np.pi*(step_size**2)]                     # Circle area = pi*(r**2)
    for r in radii:
        apers_area.append(np.pi*((r+step_size)**2 - r**2))  # Annulus area = pi*(r_out**2 - r_in**2)

    apertures = [CircularAperture(coords, r=step_size)]     # For circle aperture around center
    for r in radii:
        apertures.append(CircularAnnulus(coords, r_in=r, r_out=r+step_size))  # Annuli apertures

    phot_table = aperture_photometry(data, apertures)   # Each row is a source, with len(x_data) photometry measurements

    if show_fwhm:
        fig, ax = plt.subplots()
    sigmas = []
    new_idx = []
    for i, source in enumerate(phot_table):               # Find the sums per area (adu/pix/sec) for each of the sources
        sums = [source[i] / apers_area[i - 3] for i in range(3, len(source))]
        # Only consider moderately bright sources; too bright are saturated and too dim can have larger error values
        if 100 <= sums[0] <= 300 and sums[len(sums) - 1] < 15:
            try:    # Try unless it takes too much time...if that happens, it's probably a bad fit and we don't want it
                w_fit = x_data
                f_fit = sums
                g, _ = curve_fit(fix_x0, w_fit, f_fit)
                if show_fwhm:       # Add two plots: radial profile and Gaussian fit
                    ax.plot(x_data, sums, 'o-', alpha=0.5, lw=2, markersize=4)  # Plot the radial profiles
                    ax.plot(x_data, fix_x0(x_data, *g), '^-', alpha=0.5, lw=2, markersize=4)
                sigmas.append(np.abs(g[1]))     # Sigma can be either positive or negative; fix_x0 allows both
                new_idx.append(idx[i])  # Append to a list that will be used later for future sources
            except RuntimeError:
                pass

    fwhm_values = [i*2.35482 for i in sigmas]
    sigmas_post_clipping, _, _ = sigmaclip(sigmas, low=3, high=3)       # Clip outlier sigma values
    med = np.median(sigmas_post_clipping)       # Median sigma value
    fwhm = med * 2.35482            # FWHM to be used to find the zero point and make cuts
    if log is not None:
        log.info("Median sigma: %.3f Median FWHM: %.3f Number of sources counted/total found in field: %s/%s" %
                 (med, fwhm, len(sigmas_post_clipping), len(dao_sources)))
    else:
        print("Median sigma: %.3f Median FWHM: %.3f Number of sources counted/total found in field: %s/%s" %
              (med, fwhm, len(sigmas_post_clipping), len(dao_sources)))

    clipped_fwhm_values = [i*2.35482 for i in sigmas_post_clipping]

    if show_fwhm:
        ax.set_title("Radial profiles of sources 10 std above background")  # Plot of radial profiles and Gaussian fits
        ax.set_xlabel("Radial distance from centroid of source (pixels)")
        ax.set_ylabel("Count (adu per second per pixel)")
        plt.show()

        fig, ax = plt.subplots()       # Histogram of fwhm_values after sigma clipping
        plt.hist(clipped_fwhm_values, bins=30)
        ax.set_title("FWHM values for Gaussian-fit radial profiles")
        ax.set_xlabel("FWHM value for best-fit Gaussian (pixels)")
        plt.show()

    # Calculate the zero point
    if log is not None:
        log.info("Computing the zero point of the stack %s with filter %s and FWHM of %.3f" % (stack, fil, fwhm))
    else:
        print("Computing the zero point of the stack %s with filter %s and FWHM of %.3f" % (stack, fil, fwhm))

    orig_coords = [w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0] for j in new_idx]  # DAO pos

    # First cut: removing sources near the edge
    edge_error = 100
    first_cut_coords, first_cut_idx, first_fwhm_values = [], [], []
    for i, j in enumerate(orig_coords):
        if not (j[0] < edge_error or data.shape[0] - j[0] < edge_error or
                j[1] < edge_error or data.shape[0] - j[1] < edge_error):
            first_cut_coords.append(orig_coords[i])
            first_cut_idx.append(new_idx[i])
            first_fwhm_values.append(fwhm_values[i])

    # Second cut: removing sources that are not isolated
    # 1: Make sure no two sources are within (r_annulus_out + r_annulus_in) * fwhm of each other
    # 2: If a second source is present in the annulus, the sigma clipping for the mask will cut out the additional
    # counts from these pixels so that the second source will be cut out and not interfere with the zero point
    # If the sources are close enough together, they will either be cut out here or will look oblong and be cut out next
    overlap = []
    for i, j in enumerate(first_cut_coords):
        for k in range(i + 1, len(first_cut_coords)):
            # (r_ann_out + r_ann_in)*fwhm so that the annulus of one doesn't overlap the circular aperture of the other
            if np.sqrt((j[0] - first_cut_coords[k][0])**2 + (j[1] - first_cut_coords[k][1])**2) < \
                    (r_annulus_out + r_annulus_in) * fwhm:
                overlap.append(i)
                overlap.append(k)

    second_cut_coords, second_cut_idx, second_fwhm_values = [], [], []
    for i, j in enumerate(first_cut_coords):
        if i not in overlap:
            second_cut_coords.append(j)
            second_cut_idx.append(first_cut_idx[i])
            second_fwhm_values.append(first_fwhm_values[i])

    # Third cut: removing sources that are not circular
    # Part 1: remove sources with a fwhm 1.2 * median fwhm of the field (1.2 determined by considering 5 stacks)
    third_cut_coords, third_cut_idx, third_fwhm_values = [], [], []
    for i, j in enumerate(second_fwhm_values):
        if not (j > 1.2 * fwhm):
            third_cut_coords.append(second_cut_coords[i])
            third_cut_idx.append(second_cut_idx[i])
            third_fwhm_values.append(j)

    if log is not None:
        log.info("Checking if the stack was performed correctly with %s sources" % len(third_cut_coords))  # Final check
    else:
        print("Checking if the stack was performed correctly with %s sources" % len(third_cut_coords))  # Final check

    d_x_y = 40
    field_sigmas = cutouts_2d_gaussian(stack, third_cut_coords, d_x_y, fwhm, show_cutouts=show_contours)
    field_ratio = np.median(field_sigmas)
    if log is not None:
        log.info("Median field ratio (y_sigma/x_sigma): %.3f" % field_ratio)
    else:
        print("Median field ratio (y_sigma/x_sigma): %.3f" % field_ratio)

    if not((1 / 1.2) > field_ratio or 1.2 < field_ratio):
        if log is not None:
            log.info("Good to go! Now calculating the zero point with %s sources" % len(third_cut_coords))
        else:
            log.info("Good to go! Now calculating the zero point with %s sources" % len(third_cut_coords))

    else:      # Stack is slightly off. Check to see what the issue is
        if log is not None:
            log.info("Stack seems to be off. Sample cutout shown")
        else:
            print("Stack seems to be off. Sample cutout shown")

        if show_contours:
            for i, j in enumerate(third_cut_coords):
                d = data[int(j[1]) - d_x_y:int(j[1]) + d_x_y, int(j[0]) - d_x_y:int(j[0]) + d_x_y]  # Cutout of source
                x, y = np.meshgrid(np.linspace(0, np.shape(d)[1] - 1, np.shape(d)[1]),
                                   np.linspace(0, np.shape(d)[0] - 1, np.shape(d)[0]))

                initial_guess = (100, d_x_y, d_x_y, fwhm / 2.35482, fwhm / 2.35482, 0, 0)  # Initial guess is IMPORTANT
                popt, pcov = curve_fit(twoD_Gaussian, (x, y), d.ravel(), p0=initial_guess)  # Fit of 2D Gaussian

                if np.abs((popt[4] / popt[3]) - field_ratio) < 0.01:
                    if log is not None:
                        log.info('x sigma = %.3f y sigma = %.3f ratio = %.3f' % (popt[3], popt[4], popt[4] / popt[3]))
                    else:
                        print('x sigma = %.3f y sigma = %.3f ratio = %.3f' % (popt[3], popt[4], popt[4] / popt[3]))

                    fitted_data = twoD_Gaussian((x, y), *popt)
                    contours = np.arange(np.min(d), np.max(d), 30.)
                    norm = simple_norm(data, 'sqrt', percent=99)
                    plt.imshow(d, norm=norm)
                    plt.colorbar()
                    plt.contour(x, y, fitted_data.reshape(np.shape(d)), contours, colors='w')
                    plt.show()
                    break       # Breaks out of loop

        if input("Would you like to revisit the stack and look at isophotes within it? Type 'yes' or 'no': ") == "yes":
            try:
                zp, zp_err, fwhm = revisit_stack(stack, fil, ext, third_cut_coords, fwhm,
                                                 show_fwhm, show_bkg, show_contours, data_path, log)
                return zp, zp_err, fwhm
            except TypeError:
                pass
        elif log is not None:
            log.info("User decided not to remake the stack even though isophotes are not circular")

    final_coords, final_idx = third_cut_coords, third_cut_idx        # Final coordinates and idx to use

    # Annuli of DAO sources
    annulus_apertures = CircularAnnulus(final_coords, r_in=fwhm*r_annulus_in, r_out=fwhm*r_annulus_out)
    annulus_masks = annulus_apertures.to_mask(method='center')              # Create masks to highlight pixels in annuli

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    bkg_median = np.array(bkg_median)

    if show_bkg:            # Finds the sources for each upper median background limit and plots them over the image
        groupings = [[], [], [], [], [], [], [], []]
        for i, j in enumerate(bkg_median):
            if j < 0:
                groupings[0].append(i)
            elif j < 0.1:
                groupings[1].append(i)
            elif j < 0.25:
                groupings[2].append(i)
            elif j < 0.5:
                groupings[3].append(i)
            elif j < 0.75:
                groupings[4].append(i)
            elif j < 1:
                groupings[5].append(i)
            elif j < 2.5:
                groupings[6].append(i)
            else:
                groupings[7].append(i)

        positions = [[], [], [], [], [], [], [], []]
        for i, j in enumerate(final_idx):
            if i in groupings[0]:  # Positions of DAO
                positions[0].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[1]:
                positions[1].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[2]:
                positions[2].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[3]:
                positions[3].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[4]:
                positions[4].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[5]:
                positions[5].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[6]:
                positions[6].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
            elif i in groupings[7]:
                positions[7].append(w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0])
        colors = ['w', 'r', 'orange', 'lime', 'deepskyblue', 'b', 'purple', 'k']
        bkg_med_list = []
        for i, j in enumerate(positions):
            if len(positions[i]) != 0:
                apertures = CircularAperture(positions[i], 2.5 * fwhm)      # Aperture to perform photometry
                # Annuli of DAO sources
                annulus_apertures = CircularAnnulus(positions[i], r_in=fwhm*r_annulus_in, r_out=fwhm*r_annulus_out)
                annulus_masks = annulus_apertures.to_mask(method='center')  # Create masks to highlight pixels in annuli

                bkg_median_1 = []
                for mask in annulus_masks:
                    annulus_data = mask.multiply(data)
                    annulus_data_1d = annulus_data[mask.data > 0]
                    _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)
                    bkg_median_1.append(median_sigclip)
                bkg_median_1 = np.array(bkg_median_1)
                bkg_med_list.append(bkg_median_1)
                apertures.plot(color='white', lw=2)
                annulus_apertures.plot(color=colors[i], lw=2)

        if log is not None:
            log.info("Plotting the sources by median background value")
            log.info("Number of sources per list (<0, 0.1, 0.25, 0.5, 0.75, 1, 2.5, >2.5)")
        else:
            print("Plotting the sources by median background value")
            print("Number of sources per list (<0, 0.1, 0.25, 0.5, 0.75, 1, 2.5, >2.5)")
        for i in range(len(bkg_med_list)):
            if log is not None:
                log.info(len(bkg_med_list[i]))
            else:
                print(len(bkg_med_list[i]))

        norm = simple_norm(data, 'sqrt', percent=99)
        plt.suptitle('w < 0 | r < 0.1 | o < 0.25 | g < 0.5 | lb < 0.75 | db < 1 | p < 2.5 | b > 2.5')
        plt.imshow(data, norm=norm)
        plt.colorbar()
        plt.show()  # Plot of matched sources with annuli, color coded by median background value in annulus

    mag_two_mass = []       # Moving past show_bkg
    for i, j in enumerate(final_idx):
        mag_two_mass.append(two_mass[i][fil.lower() + '_m'])  # Reference magnitudes of 2MASS objects

    final_apertures = CircularAperture(final_coords, 2.5 * fwhm)      # Aperture to perform photometry
    # Annuli of DAO sources
    final_ann_aps = CircularAnnulus(final_coords, r_in=fwhm * r_annulus_in, r_out=fwhm * r_annulus_out)
    final_ann_masks = final_ann_aps.to_mask(method='center')  # Create masks to highlight pixels in annuli

    final_bkg_median = []
    for mask in final_ann_masks:     # For each masked annulus, find the median background value and add to list
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, stddev = sigma_clipped_stats(annulus_data_1d)
        final_bkg_median.append(median_sigclip)
    final_bkg_median = np.array(final_bkg_median)

    sigma_clip = SigmaClip(sigma=3)         # This section is to find the error on the background
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (120, 120), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    error = calc_total_error(data, bkg.background_rms, 1)

    phot_table = aperture_photometry(data, final_apertures, error=error)  # Photometry, error accounted for
    bkg_aper = final_bkg_median * final_apertures.area

    mag_inst = -2.5*np.log10((np.asarray(phot_table['aperture_sum']) - bkg_aper))  # Instrumental magnitudes (and error)

    # Find zero point
    if log is not None:
        log.info('Calculating zero point from %i stars with sigma clipping and a maximum std error of 0.1.'
                 % len(mag_inst))
    else:
        print('Calculating zero point from %i stars with sigma clipping and a maximum std error of 0.1.'
              % len(mag_inst))

    for i in np.flip(np.linspace(1, 3, num=10)):
        _, zp, zp_err = sigma_clipped_stats(mag_two_mass - mag_inst, sigma=i)
        if zp_err < 0.1:
            break

    if log is not None:
        log.info('zp = %.3f +/- %.3f.' % (zp, zp_err))
    else:
        print('zp = %.3f +/- %.3f.' % (zp, zp_err))

    # Online zp (vega): https://www.ukirt.hawaii.edu/instruments/wfcam/user_guide/performance.html
    return zp, zp_err, fwhm


# Check to make sure the stack has been aligned correctly; if not, remakes stack from good images, returns zp function
# If show_contours is False, then plots aren't shown
def revisit_stack(stack, fil, ext, coords, fwhm, show_fwhm, show_bkg, show_contours, data_path, log=None):

    with fits.open(stack) as hdr:
        header, data = hdr[0].header, hdr[0].data
    w_stack = wcs.WCS(header)
    if log is not None:
        log.info("Revisiting the stack")
    else:
        print("Revisiting the stack")
    images = glob.glob(data_path + '/red/ex' + str(ext) + '/*' + fil + '_red.fits')
    d_x_y = 40

    coords_ra = []
    coords_dec = []
    for i, j in enumerate(coords):
        pixel_coords = np.array([[j[0], j[1]]], np.float)
        pix_coords = w_stack.wcs_pix2world(pixel_coords, 1)       # w_stack
        coords_ra.append(pix_coords[0][0])
        coords_dec.append(pix_coords[0][1])
    coords_stack = SkyCoord(coords_ra * u.deg, coords_dec * u.deg, frame='fk5')  # Coordinates from sources

    if log is not None:
        log.info("There are %s sources that will be used in this process" % len(coords))
    else:
        print("There are %s sources that will be used in this process" % len(coords))

    field_ratios = cutouts_2d_gaussian(stack, coords, d_x_y, fwhm, show_cutouts=show_contours)
    stack_med_rat = np.median(field_ratios)
    if log is not None:
        log.info("The median value for the ratio (sigma y / sigma x) for the stack is %.3f" % stack_med_rat)
    else:
        print("The median value for the ratio (sigma y / sigma x) for the stack is %.3f" % stack_med_rat)

    if show_contours:
        plt.hist(field_ratios, bins=np.max([int(len(coords)/3), 10]))
        plt.xlabel("Ratio of sigma y / sigma x")
        plt.show()

    if log is not None:
        log.info("This process will take about 10 seconds per image. There are %s images in the stack" % len(images))
    else:
        print("This process will take about 10 seconds per image. There are %s images in the stack" % len(images))

    median_ratios = []
    for num, image in enumerate(images):        # Same process but for the individual images

        with fits.open(image) as hdr:     # Most of this code is copied (with small adjustments) from find_zero_point()
            image_header = hdr[0].header  # Comments for this section are in the later section
            image_data = hdr[0].data
        w = wcs.WCS(image_header)

        _, median, std = sigma_clipped_stats(data, sigma=3.0)
        daofind = DAOStarFinder(fwhm=7, threshold=10. * std)  # Looking for best sources, so raise threshold way up
        dao_sources = daofind(np.asarray(image_data))  # DAOStarFinder

        # Change pixel coordinates to real coordinates for DAOStarFinder sources
        ra_dao = []
        dec_dao = []
        for i, j in enumerate(dao_sources):
            pixel_coords = np.array([[j['xcentroid'], j['ycentroid']]], np.float)
            dao_coords = w.wcs_pix2world(pixel_coords, 1)
            ra_dao.append(dao_coords[0][0])
            dec_dao.append(dao_coords[0][1])

        # Coordinates need to be in real, not pixel
        coords_dao = SkyCoord(ra_dao * u.deg, dec_dao * u.deg, frame='fk5')  # Coordinates from sources
        idx, d2, d3 = coords_stack.match_to_catalog_sky(coords_dao)  # Match smaller to larger
        image_coords = [w.wcs_world2pix(np.array([[ra_dao[j], dec_dao[j]]], np.float), 1)[0] for j in idx]

        image_field_ratios = cutouts_2d_gaussian(image, image_coords, d_x_y, fwhm, show_cutouts=show_contours)
        med_rat = np.median(image_field_ratios)
        if log is not None:
            log.info("The median value for the ratio (sigma y / sigma x) for the field %s is %.3f" % (num, med_rat))
        else:
            print("The median value for the ratio (sigma y / sigma x) for the field %s is %.3f" % (num, med_rat))

        if show_contours:
            plt.hist(field_ratios, bins=np.max([int(len(image_coords)/3), 10]))
            plt.xlabel("Ratio of sigma y / sigma x")
            plt.show()
        median_ratios.append(med_rat)

    im_remove = []      # Make a list of the bad images (ratios more than 20% off of 1)
    for i, j in enumerate(median_ratios):
        if j > 1.2 or j < 1/1.2:
            im_remove.append(i)

    remove = input("The following images are recommended to be removed: %s. Would you like to remove these? "
                   "Type 'yes or 'no': " % im_remove)
    if remove != 'yes':
        if input("Would you like to continue without removing any images? Type 'yes' or 'no': ") == 'yes':
            if log is not None:
                log.info("User decided not to remake the stack even though isophotes are not circular")
            else:
                print("User decided not to remake the stack even though isophotes are not circular")
            return
        if log is not None:
            log.info("Select the images you would like to remove. The stack median ratio is %.3f" % stack_med_rat)
            log.info("The median ratios for the images are %s" % median_ratios)
        else:
            print("Select the images you would like to remove. The stack median ratio is %.3f" % stack_med_rat)
            print("The median ratios for the images are %s" % median_ratios)
        im_remove = [int(i) for i in input("Input indices by list position, separated by a comma: ").split(",")]
    if log is not None:
        log.info("Images removed: %s" % im_remove)
    else:
        print("Images removed: %s" % im_remove)

    new_images = []
    for i in range(len(images)):
        if i not in im_remove:
            new_images.append(images[i])
    stack_hbu = ccdproc.combine(new_images, method='median', sigma_clip=True, sigma_clip_func=np.ma.median)
    path = data_path + '/red/ex' + str(ext) + '/'
    stack = path + "Extension_" + str(ext) + "_" + fil + '_stack.fits'
    stack_hbu.write(stack, overwrite=True)

    return find_zero_point(stack, fil, ext, data_path,
                           show_fwhm=show_fwhm, show_bkg=show_bkg, show_contours=show_contours, log=log)
