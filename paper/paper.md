---
title: 'Pipeline for Optical/infrared Telescopes in Python for Reducing Images (POTPyRI)'

tags:
  - Python
  - astronomy
  - imaging
  - data reduction

authors:
  - name: K. Paterson[^corresponding author]
    orcid: 0000-0001-8340-3486
    affiliation: 1
  - name: C. D. Kilpatrick
    affiliation: 1
    orcid: 0000-0002-5740-7747
  - name: O. Eskandar
    affiliation: 1
  - name: D. Velsaco
    affiliation: 1

affiliations:
 - name: Center for Interdisciplinary Exploration and Research in Astrophysics and Department of Physics and Astronomy, Northwestern University, 1800 Sherman Ave, Evanston, IL 60201, USA
   index: 1

date: 28 Feb 2022

bibliography: paper.bib

---

# Summary

The Pipeline for Optical/infrared Telescopes in Python for Reducing Images (POTPyRI) was developed for flexible, open-source processing of imaging data from several optical/infrared (OIR) on large-aperture telescopes.  The goals of this project were to provide: 1) an open-source code base with generic methods that can be adapted to data from any OIR telescope and instrument, 2) a fully automated pipeline with minimal user input that can be used as a "black box" and handle the most common failure modes in reduction of each instrument, and 3) consistent outputs that can be easily interpreted by the user and applied to a wide range of astronomical science cases requiring deep limiting magnitudes (23--27 AB mag), especially in time-domain astronomy.  

# Statement of need

The pipeline is written in Python (currently deployed and tested on Python 3.11) and uses packages from Astropy [@astropy:2013], [@astropy:2018], ccdproc, astroquery [@Ginsburg2019] and photutils; with the additional use of SExtractor [@Bertin1996] and Astrometry.net [@Lang2010]. The code is available on GitHub at https://github.com/CIERA-Transients/POTPyRI, where instructions on the installation and detailed use can be found.

Currently available instruments (as of Feb 2022 - see the GitHub for the most recent list) include:

 - MMIRS (MMT) [@McLeod2012]
 - Binospec (MMT) [@Fabricant2019]
 - MOSFIRE (Keck) [@McLean2008]
 - DEIMOS (Keck) [@Cowley1997]
 - LRIS (Keck) [@Oke1995]
 - GMOS (Gemini-N and Gemini-S) [@Davies1997]

Since the pipeline is meant to provide the community with a way to reduce imaging data from these instruments, one science application has been the rapid reduction and identification of gamma-ray burst (GRB) afterglows, as well as the reduction and stacking of follow-up observations to identify potential host galaxies for associated and in depth study (see [@Paterson2020], [@Rastinejad2021], [@Fong2021]).

# Method

The pipeline requires two parameters to run. The first is the name of the instrument in order to load the correct setting file. This setting file allows new and user instruments to be added as needed. The second is the full data path to the data to be reduced. The data are expected to be in a particular format, including the file name and the number of extensions, as specified in the setting file. For example, for Keck, the pipeline expects the data in the format in which it is download from the Keck Observatory Archive). The details on data formats, including scripts to download and sort data, can be found on GitHub. 

![Flow diagram of the steps taken by the pipeline. Points where user input is required is shown in green. A number of different options are shown depending on the additional optional parameters selected when running the pipeline.\label{fig:flow}](../images/Pipeline_flow_diagram.pdf)

Figure~\autoref{fig:flow} shows the basic pipeline operations. The pipeline will first sort the files in the given path, and create a file list which contains information about each file such as the type, filter, configuration, and target. Using the file list and the setting file, the pipeline will then create calibration files (e.g. bias, dark, flat) and then performs the necessary basic reductions to science files, grouped by target, including gain correction, bias subtraction, dark subtraction, flat correction, and trimming. A source mask, which masks out sources on the image, is created for internal use in downstream steps such as the background, NIR and fringe map creation. A mask is created for each image, which shows the position of bad pixels, cosmic rays, and satellite trails. In order to prepare images for median stacking, the 2D background of each science image is determined and subtracted. Science images are also divided by their exposure time, to allow the median stacking of images with different exposure times, and results in a reduced image with pixel values in e$^{-}$/s. For instruments that do not have astrometry i.e. a World Coordinate System (WCS) in the header, a WCS is added to the header based on the position and orientation of the telescope and instrument.

After basic reductions, additional corrections are applied. For the NIR, the sky background changes on the order of a few minutes or less. Thus, for NIR observations, it is necessary to create an NIR sky map for each science image, based on observations occurring close in time. In the optical, redder filters such as $i$ and $z$-bands, can also suffer from fringing, which introduces low-level structure in the background. For filters that require this correction (which is set in the setting file for the instrument), the pipeline creates a fringe map based on all science images for a target for that filter.

Next, the pipeline prepares the reduced science images (again grouped by target) for aligning and stacking. Images are aligned by calculating the relative shift and transformation between the first image in the observation sequence and each subsequently image using unique star quads matches. Spatial flips to the image are determined automatically based on the WCS header keywords. After the science images are aligned, the pipeline performs some basic quality checks on the images before stacking. These quality checks include removing images of poor quality based on the number of stars detected (using SExtractor) in the image (e.g. no stars found due to poor conditions), the Full-Width Half Maximum (FWHM) and elongation of sources in the image (based on 3$\sigma$ cuts from the median value), and poor alignment. There are no limits on the number of images that can be removed due to poor conditions or alignment, while images removed due to a FWHM or elongation cut are restriction to 10\% of the total number of images for that target. After the quality checks, the pipeline then stacks the remaining images using a median combine (see Figure \ref{fig:stacks} for examples of median stacks from MOSFIRE, MMIRS, LRIS (blue channel), BINOSPEC (left side) and DEIMOS).

Next, the pipeline will attempt to solve the WCS of the stack to improve the astrometry of the image. The pipeline uses unique star quads to match against stars from the Gaia-DR3 catalog [@Gaia2021] to calculate the initial shift and pixel scale between the native WCS and the true astrometry of the image. After the initial solution has been applied, a star to star match between the extracted sources from the image and the Gaia catalog is used to calculate the full WCS of the image in the TPV (TAN, or tangential projection plus distortion polynomials) system.

After the astrometry, the pipeline performs Point Spread Function (PSF) photometry. First, the PSF is determined using a sample of unsaturated star-like sources, determined based on 3$\sigma$ cuts from the median values and flux (to remove saturated and faint sources). Next, the pipeline uses the PSF to calculate PSF photometry for the full source catalog. The pipeline then matches the extracted sources against a photometric catalog (set in the setting file for the instrument) and uses the PSF photometry to calculated the zero point through iterative sigma clipping for the image. Zero points can be calculated using the Sloan Digital Sky Survey DR12 [@Eisenstein2011], Pan-STARRS [@Tonry2012] and 2MASS [@Skrutskie2006] catalogs in the $u,g,r,i,z,Y,J,H,K$ filters. The zero point is then transformed to the AB system and applied to the photometric catalog.

This concludes the automatic reductions that the pipeline performs on input data. However, several times during the course of the pipeline operation, the pipeline will ask the user to confirm changes or the quality of the outputs before continuing. A basic flow diagram of the pipeline steps, including some different ways the pipeline can process and highlighting the points where user input is required, is shown in Figure \autoref{fig:stacks}. These prompts from the pipeline allow the user to make changes to the file list, perform manual astrometry, or abort the reduction for a particular target if the quality of the outputs are not satisfactory. There are also several additional options the user can specify when running the pipeline to perform additional operations such as skipping the creation of the calibration files, only selecting a particular target for reduction, skipping the creation of the stack, and performing manual aperture photometry on a single target. For more details on these additional options, and detailed instructions on running the pipeline and pipeline outputs, please see the documentation available on the GitHub. For user support or errors running the pipeline, please open a Issue on GitHub.

While the pipeline runs, the pipeline will also created a detailed log of operations. Along with the log, the pipeline will write information such as the readnoise, effective gain, zero point, and limiting mag to the stack. The pipeline also produced a number of other output files, including source catalogs, star lists, error plots etc. For more details on these outputs, refer to the documentation of the GitHub.

![Examples of stacks produced by the pipeline from MOSFIRE, MMIRS, LRIS (blue side), BINOSPEC (left side) and DEIMOS data. The position of Gaia DR3 stars are shown by the green circles.\label{fig:stacks}](../images/Stack_examples.pdf)

# Acknowledgements

KP acknowledges the help of the Fong group at Center for Interdisciplinary Exploration and Research in Astrophysics, including W. Fong, J. Rastinejad, A. Rouco Escorial, G. Schroeder, A. E. Nugent, A. Gordon for helping to test the pipeline and providing feedback. KP acknowledges support by the National Science Foundation under grant Nos. AST-1814782, AST-1909358 and CAREER grant No. AST-2047919.

# References
