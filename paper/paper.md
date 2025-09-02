---
title: 'Pipeline for Optical/infrared Telescopes in Python for Reducing Images (POTPyRI)'

tags:
  - Python
  - astronomy
  - imaging
  - data reduction

authors:
  - name: K. Paterson^[corresponding author]
    orcid: 0000-0001-8340-3486
    affiliation: "1, 2" #
  - name: C. D. Kilpatrick
    affiliation: 1
  - name: O. Eskandar
    affiliation: 1
  - name: D. Velsaco
    affiliation: 1

affiliations:
 - name: Center for Interdisciplinary Exploration and Research in Astrophysics and Department of Physics and Astronomy, Northwestern University, 1800 Sherman Ave, Evanston, IL 60201, USA
   index: 1
 - name: Max-Planck-Institut für Astronomie, Königstuhl 17, 69117 Heidelberg, Germany
   index: 2

date: 2 Sep 2025

bibliography: paper.bib

---

# Summary

This pipeline was developed for the reduction and stacking of imaging data for a number of telescopes and instruments in the optical and near-infrared (NIR) bands. The purpose of the pipeline is to provide an automated way to reduce imaging data, create a median stack with reliable astrometry and photometric calibration, and perform aperture and PSF photometry on sources within the stacked image.

# Statement of need

The pipeline is written in Python (currently deployed and tested on Python 3.12) and uses packages from Astropy [@astropy:2013, @astropy:2018], ccdproc, astroquery [@Ginsburg2019] and photutils [@photutils:2024]; with the additional use of external dependencies including SExtractor [@Bertin1996] and Astrometry.net [@Lang2010]. The code is available on GitHub at https://github.com/CIERA-Transients/POTPyRI, where instructions on the installation and detailed use can be found. The code is also installable with pip and is release at PyPI at https://pypi.org/p/potpyri.

Currently available instruments (as of June 2025 - see the GitHub for the most recent list) include:
 - MMIRS (MMT) [@McLeod2012]
 - Binospec (MMT) [@Fabricant2019]
 - MOSFIRE (Keck) [@McLean2008]
 - DEIMOS (Keck) [@Cowley1997]
 - LRIS (Keck) [@Oke1995]
 - GMOS (Gemini-N and Gemini-S) [@Davies1997]
 - Flamingos2 (Gemini) [@Eikenberry2012]
 - FourStar (Magellan) [@Persson2013]
 - IMACS (Magellan) [@Bigelow2003]

Since the pipeline is meant to provide the community with a way to reduce imaging data from these instruments; science applications include the rapid reduction and identification of transients, such as Gamma-Ray Burst (GRB) afterglows, as well as the reduction and stacking of follow-up observations of transients to identify potential host galaxies for associated and in depth study (see [@Paterson2020,@Rastinejad2021,@Fong2021,@Rastinejad2025,@Caleb2025]).

# Method

The pipeline requires two parameters to run. The first is the name of the instrument in order to load the correct settings file. This settings file allows new and user instruments to be added as needed. The second is the full data path to the data to be reduced. The data are expected to be in a particular format, including the file name and the number of extensions, as specified in the setting file. For example, for Keck, the pipeline expects the data in the format in which it is download from the Keck Observatory Archive.  In this case, the user can specify that the input data conform to the archive format using *--proc archive*. The details on data formats, including scripts to download and sort data, can be found on GitHub. 

![Flow diagram showing the basic steps taken by the pipeline. The dashed box shows all steps taken under the process science function, while all steps after the calibration creation is done per target.\label{fig:flow}](../images/Pipeline_flow_diagram.pdf)

Figure \autoref{fig:flow} shows the basic pipeline operations. The pipeline will first sort the files in the given path, and create a file list which contains information about each file such as the type, filter, configuration, and target. Using the file list and the settings file, the pipeline will then create calibration files (e.g. bias, dark, flat) and performs the necessary basic reductions to science files, grouped by target, including gain correction, bias subtraction, dark subtraction, flat correction, and trimming. During this step, a static mask, if available, is used to identify bad pixels. This mask is later extended to include image-specific issues such as satellites trails, cosmic rays (cleaned from images) and saturated stars. A 2D sky background is created for each image and subtracted. Astrometry.net is used on each image to provide a rough WCS solution, after which a fine astrometric solution is computed using Gaia-DR3 catalogs [@Gaia2021]. Images are then aligned using the new WCS. After basic reductions, additional corrections where applicable are applied. For the NIR, the sky background can change during the night. Thus, for NIR observations, an additional NIR sky map is created for each science image, based on observations for each target. In the optical, redder filters such as $i$ and $z$-bands, can also suffer from fringing, which introduces low-level structure in the background. Thus, for filters that require this correction (based on each instrument and thus set in the settings file for the instrument), the pipeline creates a fringe map based on all science images for a target for that filter. An uncertainty map for each image is created to include noise added by the reduction such as the background subtraction.

After solving for all the astrometry, the pipeline will remove any images which high astrometric dispersion (i.e. images that have dispersion $>5\sigma$ either in RA or DEC), except if the user specified the *--keep\_all\_astro* option, before stacking images. Images are median stacked, using the uncertainty masks and exposure times as weights, and a final stack mask is created (see Figure \autoref{fig:stacks} for examples of median stacks from MOSFIRE, MMIRS, LRIS (blue channel), BINOSPEC (left side) and DEIMOS).

The pipeline will then perform both aperture and Point Spread Function (PSF) photometry on sources in the stacked image. For the aperture photometry, the pipeline will detect sources within the image and determine statistics within a fiducial aperture radius, which is written as an extension to the stacked image. For the PSF photometry, the PSF is determined using a list of bright unsaturated stars, determined based on cuts on roundness, Full Width Half Maximum (FWHM), and signal-to-noise.
Next, the pipeline uses the PSF to calculate PSF photometry for all the originally extracted sources. After the PSF photometry has been calculated, the pipeline will then calculate a zero point based on aperture photometry by downloading a catalog of standard stars. Currently supported catalogs are the Sloan Digital Sky Survey DR12 [@Eisenstein2011], Pan-STARRS [@Tonry2012], SkyMapper [@Onken2024],  and 2MASS [@Skrutskie2006] covering most optical and infrared bands in the north and south. The zero point is then transformed to the AB system and written to the header.

In addition, an approximate limiting magnitude is calculated for each stack using the FWHM and the approximate standard deviation in the sky background per pixel.  These values are propagated into the image header as 3-, 5-, and 10-$\sigma$ limiting magnitudes as *M3SIGMA*, *M5SIGMA*, and *M10SIGMA*.

For more detailed instructions on running the pipeline and a description of pipeline outputs, please see the documentation available on the GitHub. For user support or errors running the pipeline, please open a Github issue.

While the pipeline runs, the pipeline will also created a detailed log of operations. Along with the log, the pipeline will also write information such as the readnoise, effective gain, and zero point to the stack. The pipeline also produces a number of other output files, including source catalogs, star lists, error plots etc. For more details on these outputs, refer to the documentation of the GitHub.

![Examples of stacks produced by the pipeline from MOSFIRE, MMIRS, LRIS (blue side), BINOSPEC (left side) and DEIMOS data. The position of Gaia DR3 stars are shown by the green circles.\label{fig:stacks}](../images/Stack_examples.pdf)

# Acknowledgements

KP acknowledges the help of the Fong group at Center for Interdisciplinary Exploration and Research in Astrophysics, including W. Fong, J. Rastinejad, A. Rouco Escorial, G. Schroeder, A. E. Nugent, A. Gordon for helping to test the pipeline and providing feedback. KP acknowledges support by the National Science Foundation under grant Nos. AST-1814782, AST-1909358 and CAREER grant No. AST-2047919.

# References
