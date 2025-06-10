# POTPyRI: Pipeline for Optical/infrared Telescopes in Python for Reducing Images

[![GitHub Actions](https://github.com/CIERA-Transients/POTPyRI/actions/workflows/build-test.yml/badge.svg)](https://github.com/CIERA-Transients/POTPyRI/actions/workflows/build-test.yml)
[![Docs](https://img.shields.io/badge/docs-latest-blue)](https://ciera-transients.github.io/POTPyRI/)


Data reduction pipeline for imaging from large aperture telescopes

## Usage

If you use this code, please reference Paterson et al. (in prep.)

## Installation

The recommended installation is via conda and pip.  First create a conda environment:

```conda create -n potpyri python=3.12 pip-tools```

Then to install POTPyRI, run:

```
conda activate potpyri
pip install git+https://github.com/CIERA-Transients/POTPyRI.git
```

which will install the remaining dependencies from the requirements.txt file.

### Non-Python Dependencies

POTPyRI has two non-python dependencies: `astrometry.net` and `source extractor`.  The package file will install these dependencies on `linux-64` and `osx-64` environments for which there are working conda repositories.  On newer Mac OS systems (`osx-arm64`), conda will not be able to install these packages and they must be built manually.  The recommended installation method is with Homebrew (https://brew.sh/) via:

```
brew install sextractor
brew install astrometry-net
```

Finally, note that for `astrometry.net`, index files are required.  There is a utility script in POTPyRI that will attempt to install the latest index files.  Once you have successfully installed `astrometry.net` in your path, run:

```
download_anet_index
```

## Supported Instruments

- BINOSPEC (MMT)
- DEIMOS (Keck)
- Flamingos2 (Gemini)
- FourStar (Magellan)
- GMOS (Gemini)
- IMACS (Magellan)
- LRIS (Keck)
- MMIRS (MMT)
- MOSFIRE (Keck)

## Running

Once POTPyRI and external dependencies are installed, the basic syntax to run the pipeline is:

```main_pipeline instrument data_path```

where `instrument` is the name of the instrument you wish to reduce data from (see **Supported Instruments**), and `data_path` is the full path of the data you wish you reduce. Both the `instrument` and `data_path` are required parameters. Optional parameters are discussed in **Pipeline parameters**.

The pipeline is designed to run with no user input, including file ingestion, pixel calibration with provided calibration files, astrometry, image stacking, photometry, and flux calibration. The pipeline will display the log to the terminal (also saved as a file under `red/log`) to provide feedback and indicate when an error is encountered. For more details on each of these sections, including how the pipeline performs each task, read the sections below. 

All processed data will be written out in the `data_path` under the `red` subdirectory.  Raw data are stored in `raw` and unused files are stored in `bad`.  Processed calibration files will be stored in `red/cals`, logs in `red/log`, and interstitial image files before stacking are stored in `log/workspace`. A description of each type of output file can be found in **Outputs**.

## Pipeline parameters

In addition to the required parameters `instrument` and `data_path`, the following optional parameters can be passed to the pipeline (also viewable via `main_pipeline --help`):

```
options:
  -h, --help            show this help message and exit
  --target TARGET       Option to only reduce a specific target. String used here must be contained within the target name in file headers. Optional
                        parameter.
  --proc PROC           Option to specify file processing for data ingestion.
  --include-bad, --incl-bad
                        Include files flagged as bad in the file list.
  --no-redo-sort        Do not resort files and generate new file list.
  --file-list-name FILE_LIST_NAME
                        Change the name of the archive file list.
  --phot-sn-min PHOT_SN_MIN
                        Minimum signal-to-noise to try in photometry loop.
  --phot-sn-max PHOT_SN_MAX
                        Maximum signal-to-noise to try in photometry loop.
  --fwhm-init FWHM_INIT
                        Initial FWHM (in pixels) for photometry loop.
  --skip-skysub         Do not perform sky subtraction during image processing.
  --fieldcenter FIELDCENTER FIELDCENTER
                        Align the output images to this center coordinate. This is useful for difference imaging where the frames need to be a common size,
                        pixel scale, and set of coordinates.
  --out-size OUT_SIZE   Output image size (image will be SIZE x SIZE pixels).
  --stages STAGES [STAGES ...]
                        Stages to execute if running the pipeline in a modular way.
  --skip-flatten        Tell the pipeline to skip flattening.
  --skip-cr             Tell the pipeline to skip cosmic ray detection.
  --skip-gaia           Tell the pipeline to skip Gaia alignment during WCS.
```

## Outputs

All outputs from the pipeline are written out to the `red` folder and various subdirectories. Calibration files such as the master bias, flat and dark are saved using `mbias`, `mflat` and `mdark` prefixes, along with additional information such as the filter, amplifiers and binning in the file names.  These calibration files are stored under `red/cals`.

All processed science files are renamed with following the basic format: `{object}.{filter}.{ut_date}.{amplifier}.{binning}.{unique_number}*.fits`.  All data are processed together within a common `TargType`, defined as image frames with the same object, filter, amplifier setup, and binning.  This variable is defined in the `file_list.txt` located in `data_path` and that is generated from `potpyri/sort_files.py` during the initial file processing performed by POTPyRI.  A more detailed description of `file_list.txt` is provided below.

The following table provides a basic description of the naming format, location, and brief description for each science output file that POTPyRI will produce:

```
Filename         Location       Description
*.fits           red/workspace  Initial processed science files, containing SCI (image data), MASK (mask), and UNCERT (error) extensions
*_bkg.fits       red/workspace  A background image calculated during image processing, containing only a PRIMARY (background data) extension
*_reproj.fits    red/workspace  Reprojection of the science image to the final data frame common to that `TargType`, containing only a SCI (image data) extension
*_data.fits      red/workspace  Generated directly from the corresponding `*_reproj.fits` file, containing the SCI (image data), MASK (mask), and UNCERT (error) extensions
*stk.fits        red            The stacked image data for the corresponding `TargType`, containing SCI (image data), MASK (mask), and ERR (error) extensions
```

Once photometry is performed, additional extensions are added to each `*stk.fits` file containing fits-formatted tables with aperture phohotometry, PSF stars, and PSF photometry of identified sources in the image.  Currently, only the aperture photometry table is used by `potpyri/absphot.py` for flux calibration.

## File list

The pipeline will sort through all the fits with the correct format for a given instrument (given by the **raw_format** function in the setting file) and create a file list with the file type, target name, exposure time, observation time, and instrument setup such as number of amps and binning.

This process is designed to be automatic and account for common observation or archiving errors by checking various header keywords to classify each file.  If you find your files are being misclassified for a POTPyRI-supported instrument, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

## Image Calibration

Bias, dark, and flat-field images are defined by the corresponding keywords in each `potpyri/instruments/*.py` file.  The methods for generating and applying each calibration frame are generic and defined entirely within the general instrument class in `potpyri/instruments/instrument.py`.  Combined with instrument-specific differences in overscan, trimming, and static masking, we find that these methods provide good quality initial processed images with few bad pixels and generally Gaussian noise for the typical calibration frames that are available by instrument.

If you find that your initial processed images (see **Outputs**) are noisy, contain a large number of bad pixels, or contain other artifacts due to pixel-level processing, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

## Aligning and Astrometry

Images are aligned in a two-step process using `astrometry.net` followed by centroiding to Gaia DR3 astrometric standard stars.  This process has been tested across all POTPyRI-supported instruments in a variety of instrument setups, read out modes, filters, and fields with varying densities of stars.  We find that it generally results in 0.1-1 pixel dispersion in the alignment solution per image, resulting in good image alignment for stacking and preservation of the PSF shape.

If you find that you have good quality images that consistently fail image alignment and astrometry, have extremely large dispersions, or have poor WCS solution, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

## Image Masking and Stacking

POTPyRI has been significantly updated to implement optimal masking for satellite trails, cosmic rays, bad pixels that can be statically masked or are introduced by the pixel-level calibration, and saturation/non-linear effects.  These pixels are tracked throughout the pipeline from calibration through image stacking and are stored within the `*.fits`, `*_data.fits`, and `*stk.fits` using a bitwise image mask.  In general, the following schema is used by `potpyri/image_procs.py` to flag bad pixels:

```
Bit     Value
1       Bad pixel (e.g., static mask, NaN during image processing)
2       Cosmic ray flagged pixel
4       Saturated pixel
8       Neighbor of a bad pixel set within a distance determined by "grow" parameter
```

These pixels will be ignored during image stacking, and only pixels that are not flagged by these parameters will be used to generate the final stack.

In addition, error/uncertainty images are generated for each frame, accounting for the expected read noise, Poisson noise, and an additional empirical noise term defined by the standard deviation in the sky background.  These error images set the `weight` term for each input frame when performing stacking.

Final image stacking is performed by `ccdproc.combine` with the individual image data, mask, and error frames for each science image.  The stacking method is generally `median`, but can be changed via the **stack_method** value in each `potpyri/instruments/*.py` parameter file.  Images are scaled by the exposure time within each header to account for variable depth between frames.

## Automatic Aperture and PSF photometry

The pipeline will perform both automatic aperture and PSF photometry of sources in the stacked image using `potpyri/photometry.py`. Firstly, the pipeline will detect sources within the image and determine statistics within a fiducial aperture radius using `photutils.aperture.ApertureStats`.  This initial table of aperture photometry is saved within the `*stk.fits` image as the extension `APPPHOT`.

Next, based on cuts on roundness, FWHM, and signal-to-noise, the pipeline will define a list of bright stars with which it calculates an empirical PSF model.  The final list of PSF stars is saved in the `*stk.fits` file as the `PSFSTARS` extension.  The PSF itself will be generated from the extracted data around these PSF stars using `photutils.psf.EPSFBuilder`.  A stamp of the final effective PSF is saved in the `*stk.fits` file as the `PSF` extension.  A final FWHM is empirically calculated from the effective PSF model using a `astropy.modeling.functional_models.Moffat2D` defined profile and saved to the `*stk.fits` header.

Once the pipeline has determined the PSF, it will then calculate PSF photometry for all the originally extracted sources with `photutils.psf.PSFPhotometry` and write them to the `PSFPHOT` extension of the `*stk.fits` file.  In addition, as a data quality check a residual image is saved to the `*stk.fits` file showing the original science data with each of these sources subtracted using the PSF model as the `RESIDUAL` extension.

## Flux Calibration

After the PSF photometry has been calculated, the pipeline will then calculate a zero point based aperture photometry by downloading a catalog of standard stars.  Currently supported catalogs are Pan-STARRS, SDSS, SkyMapper, and 2MASS, covering most optical and infrared bands in the north and south.

The zero point and associated uncertainty are saved in the `*stk.fits` files as `ZPTMAG` and `ZPTMUCER`, respectively.  In addition, an approximate limiting magnitude is calculated for each `*stk.fits` file using the `FWHM` and `SKYSIG` (the approximate standard deviation in sky background per pixel) header keywords and saved to the header as `M3SIGMA`, `M5SIGMA`, and `M10SIGMA` (3-, 5-, and 10-sigma limiting magnitudes).

## Quality checks

The pipeline will automatically create a log file in the `red/log` folder named after the instrument, followed by the date and time at the start of the run (UTC). A lot of information is written to the log file about the various steps step in the pipeline and the reduction process. Since this log file can also contain the reduction of multiple targets, it is suggested to make a personal log for each target with key information that can be used for the data analysis. Before, during and after the pipeline the has run, there are some key aspect you should check in terms of quality control that should be noted in the target log. These steps are list and explained below:

### Data Processing Checks

- Check the file list: The first thing the pipeline will do is to create a file list containing information about the files. Files can be misclassified, which can be corrected via direct editing of the file list and reprocessing with `--no-redo-sort`, direct editing of your file headers to correct the misclassification, or submitting an issue to debug POTPyRI-related classification issues.

- Check the stacks: After the pipeline has reduced and stacked the files for a target, it is recommended that you check your final stacked images to assess data quality. Check the fits image named after the target name. If the image is of poor quality with very few stars, WCS and flux calibration will be affected.  If you find your stack has a large number of image artifacts, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

- Check WCS solution: In additional to individual image WCS solutions, the pipeline will calculate an automatic WCS solution on the final stack. The rms on the astromtery is typically $<$1 pixel and the stars aligned well near your target as well as across the image.  Check the image alignment quality via the `RADISP` and `DEDISP` keywords in each image header.  If the WCS is off by a large number of pixels, the zero point will not be calculated correctly and the reduction should be aborted. In this case, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

- Check the PSF: The PSF calculated by the pipeline will determine if the PSF photometry, as well as the zero point can be trusted. The PSF is saved in the `*.epsf.png` file, and will usually look like a 2D Gaussian. The pipeline will also report the x and y sigma of a 2D Gaussian fitted to the PSF. If the PSF looks good and the x and y sigma values reeported are similar, the PSF photometry can be trusted. If the PSF doesn't look Gaussian-like or the reported x and y sigma values are very different or negative, then there was an issue determining the PSF. If there are a lot of good stars in the final stack that could have been used to determine the PSF, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

- Check the zero point calculation: If the PSF photometry is reliable, the zero point calculation should be good if a large number of stars is used to calculate it. In general, zero points calculated with $>$10 stars will provide relatively small (10$--$20 mmag level or smaller) zero point uncertainty. If there are more stars in the field that the pipeline should have used, or the pipeline reported that no stars were found but there is coverage in the corresponding catalog, please contact us or submit a Github Issue with the instrument and specific error/problem indicated.

### Data over multiple nights

POTPyRI can stack data from multiple nights, although it will not use unique calibration frames for each stacked image.  To stack multiple nights of night, move all files to a single `raw` directory and run the pipeline.

## Issues and error

Should you encounter an error in the pipeline, have a special data setup that doesn't run through the pipeline, or wish to add an additional instrument, please contact the developers at `ckilpatrick@northwestern.edu` or open an Issue on GitHub.
