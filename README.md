# POTPyRI: Pipeline for Optical/infrared Telescopes in Python for Reducing Images

[![GitHub Actions](https://github.com/CIERA-Transients/POTPyRI/actions/workflows/build-test.yml/badge.svg)](https://github.com/CIERA-Transients/POTPyRI/actions/workflows/build-test.yml)
[![Docs](https://img.shields.io/badge/docs-latest-blue)](https://ciera-transients.github.io/POTPyRI/)

Data reduction pipeline for imaging from large aperture telescopes. Supports automatic file ingestion, pixel calibration, astrometry (astrometry.net + Gaia DR3), stacking, aperture/PSF photometry, and flux calibration.

## Installation

**Requirements:** Python 3.11 or later.

### From PyPI (recommended)

```bash
pip install potpyri
```

### With conda (recommended for managing Python and non-Python deps)

```bash
conda create -n potpyri python=3.12 pip
conda activate potpyri
pip install potpyri
```

### From source (development)

```bash
git clone https://github.com/CIERA-Transients/POTPyRI
cd POTPyRI
pip install -e .                    # editable install
pip install -e ".[test]"             # add pytest, pytest-cov
pip install -e ".[docs]"             # optional: Sphinx docs
```

### Non-Python dependencies

The pipeline needs **astrometry.net** and **Source Extractor** (SExtractor):

| Platform | Install method |
|----------|----------------|
| **Linux (x86_64)** | `conda install conda-forge::astromatic-source-extractor conda-forge::astrometry` |
| **macOS (Intel, osx-64)** | Same as Linux: both packages are available from conda-forge. |
| **macOS (Apple Silicon, arm64)** | **astromatic-source-extractor** is available from conda-forge. **astrometry** has no conda build for osx-arm64; install it (and optionally both) via Homebrew: `brew install sextractor astrometry-net`. |

So: use the `environment.yml` (or the conda commands above) on Linux and Intel Mac. On Apple Silicon, either create the conda env for Python + source-extractor and add astrometry via Homebrew, or install both with `brew install sextractor astrometry-net` and use a separate conda env with only `python=3.12` and `pip` for POTPyRI.

- **Astrometry.net index files** (~59 GB): after `astrometry.net` is on your PATH, run `download_anet_index` (or use the script in the package to fetch indices).

## Usage

Run the full pipeline from the command line:

```bash
main_pipeline INSTRUMENT DATA_PATH
```

- **INSTRUMENT:** One of the supported instrument names (e.g. `GMOS`, `LRIS`, `BINOSPEC`). See **Supported instruments** below.
- **DATA_PATH:** Top-level directory containing your data (e.g. a `raw` folder with FITS files).

Example:

```bash
main_pipeline GMOS /path/to/my/run
```

All options: `main_pipeline --help`. Processed data are written under `DATA_PATH/red/`; logs go to `red/log/`.

**Citation:** If you use this code, please cite Paterson et al. (in prep.).

## Testing

From the project root (with test deps installed: `pip install -e ".[test]"`):

```bash
pytest tests
```

- **Unit tests only (no network):**  
  `pytest tests -m "not integration"`  
  This skips the five integration tests that download fixtures from external storage. CI runs this by default.
- **All tests (including integration):**  
  `pytest tests`  
  Requires network access for fixture downloads.
- **Verbose:**  
  `pytest tests -v`

## Required Disk Space

The POTPyRI repository itself has a size of: ![Size](https://img.shields.io/github/repo-size/CIERA-Transients/POTPyRI).

In addition, the required index files installed by `download_anet_index` require a large footprint, currently about 59 GB of disk space.

## Supported instruments

- BINOSPEC (MMT)
- DEIMOS (Keck)
- Flamingos2 (Gemini)
- FourStar (Magellan)
- GMOS (Gemini)
- IMACS (Magellan)
- LRIS (Keck)
- MMIRS (MMT)
- MOSFIRE (Keck)

## Running and outputs

The pipeline runs non-interactively: file ingestion, pixel calibration, astrometry, stacking, photometry, and flux calibration. The log is printed to the terminal and saved under `data_path/red/log`.

Directory layout under `data_path`:

- **raw** — raw data
- **bad** — files excluded from reduction
- **red** — all processed outputs
- **red/cals** — master bias, flat, dark
- **red/log** — log files
- **red/workspace** — intermediate per-frame images before stacking

Details of output file names and contents are in **Outputs**.

## Pipeline parameters

Required: `instrument`, `data_path`. All other options are optional. Run `main_pipeline --help` for the full list.

| Option | Default | Description |
|--------|---------|-------------|
| `--target TARGET` | — | Only reduce targets whose name contains this string. |
| `--proc PROC` | — | File processing mode for data ingestion. |
| `--include-bad` / `--incl-bad` | off | Include files flagged as bad in the file list. |
| `--no-redo-sort` | off | Do not re-sort files or regenerate the file list. |
| `--file-list-name NAME` | `file_list.txt` | Name of the generated file list. |
| `--phot-sn-min N` | 3.0 | Minimum S/N to try in the photometry loop. |
| `--phot-sn-max N` | 20.0 | Maximum S/N to try in the photometry loop. |
| `--fwhm-init N` | 5.0 | Initial FWHM (pixels) for the photometry loop. |
| `--skip-skysub` | off | Skip sky subtraction during image processing. |
| `--fieldcenter RA DEC` | — | Align outputs to this center (e.g. for difference imaging). |
| `--out-size N` | — | Output image size (N×N pixels). |
| `--skip-flatten` | off | Skip flat-field correction. |
| `--skip-cr` | off | Skip cosmic-ray detection. |
| `--skip-gaia` | off | Skip Gaia alignment in the WCS step. |
| `--keep-all-astro` | off | Keep all images regardless of astrometric dispersion (default: drop high-dispersion frames). |

## Outputs

All outputs from the pipeline are written out to the `red` folder and various subdirectories. Calibration files such as the master bias, flat and dark are saved using `mbias`, `mflat` and `mdark` prefixes, along with additional information such as the filter, amplifiers and binning in the file names.  These calibration files are stored under `red/cals`.

All processed science files are renamed following the basic format: `{object}.{filter}.{ut_date}.{amplifier}.{binning}.{unique_number}*.fits`.  All data are processed together within a common `TargType`, defined as image frames with the same object, filter, amplifier setup, and binning.  This variable is defined in the `file_list.txt` located in `data_path` and that is generated from `potpyri/primitives/sort_files.py` during the initial file processing performed by POTPyRI.  A more detailed description of `file_list.txt` is provided below.

The following table provides a basic description of the naming format, location, and brief description for each science output file that POTPyRI will produce:

```
Filename         Location       Description
*.fits           red/workspace  Initial processed science files, containing SCI (image data), MASK (mask), and UNCERT (error) extensions
*_bkg.fits       red/workspace  A background image calculated during image processing, containing only a PRIMARY (background data) extension
*_reproj.fits    red/workspace  Reprojection of the science image to the final data frame common to that `TargType`, containing only a SCI (image data) extension
*_data.fits      red/workspace  Generated directly from the corresponding `*_reproj.fits` file, containing the SCI (image data), MASK (mask), and UNCERT (error) extensions
*stk.fits        red            The stacked image data for the corresponding `TargType`, containing SCI (image data), MASK (mask), and ERR (error) extensions
```

Once photometry is performed, additional extensions are added to each `*stk.fits` file containing FITS-formatted tables with aperture photometry, PSF stars, and PSF photometry of identified sources in the image. Currently, only the aperture photometry table is used by `potpyri/primitives/absphot.py` for flux calibration.

## File list

The pipeline will sort through all FITS files with the correct format for a given instrument (given by the **raw_format** function in the settings file) and create a file list with the file type, target name, exposure time, observation time, and instrument setup such as number of amps and binning.

This process is designed to be automatic and account for common observation or archiving errors by checking various header keywords to classify each file. Classification logic lives in `potpyri/primitives/sort_files.py`. If you find your files are being misclassified for a POTPyRI-supported instrument, please contact us or open a GitHub issue with the instrument and specific error indicated.

## Image Calibration

Bias, dark, and flat-field images are defined by the corresponding keywords in each instrument module under `potpyri/instruments/`. The methods for generating and applying each calibration frame are generic and defined in the base class in `potpyri/instruments/instrument.py`; calibration orchestration is in `potpyri/primitives/calibration.py`.  Combined with instrument-specific differences in overscan, trimming, and static masking, we find that these methods provide good quality initial processed images with few bad pixels and generally Gaussian noise for the typical calibration frames that are available by instrument.

If you find that your initial processed images (see **Outputs**) are noisy, contain a large number of bad pixels, or contain other artifacts due to pixel-level processing, please contact us or open a GitHub issue with the instrument and specific error indicated.

## Aligning and Astrometry

Images are aligned in a two-step process using `astrometry.net` followed by centroiding to Gaia DR3 astrometric standard stars.  This process has been tested across all POTPyRI-supported instruments in a variety of instrument setups, read out modes, filters, and fields with varying densities of stars.  We find that it generally results in 0.1-1 pixel dispersion in the alignment solution per image, resulting in good image alignment for stacking and preservation of the PSF shape.

If you find that you have good quality images that consistently fail image alignment and astrometry, have extremely large dispersions, or have a poor WCS solution, please contact us or open a GitHub issue with the instrument and specific error indicated.

## Image Masking and Stacking

POTPyRI has been significantly updated to implement optimal masking for satellite trails, cosmic rays, bad pixels that can be statically masked or are introduced by the pixel-level calibration, and saturation/non-linear effects.  These pixels are tracked throughout the pipeline from calibration through image stacking and are stored within the `*.fits`, `*_data.fits`, and `*stk.fits` using a bitwise image mask.  In general, the following schema is used by `potpyri/primitives/image_procs.py` to flag bad pixels:

```
Bit     Value
1       Bad pixel (e.g., static mask, NaN during image processing)
2       Cosmic ray flagged pixel
4       Saturated pixel
8       Neighbor of a bad pixel set within a distance determined by "grow" parameter
```

These pixels are ignored during image stacking; only pixels not flagged are used to generate the final stack.

Error/uncertainty images are generated for each frame (read noise, Poisson noise, and empirical sky noise). These error images set the `weight` term for each input frame when stacking.

Stacking is performed by `ccdproc.combine` in `potpyri/primitives/image_procs.py` with the individual image data, mask, and error frames for each science image.  The stacking method is generally `median`, but can be changed via the **stack_method** value in each `potpyri/instruments/*.py` parameter file.  Images are scaled by the exposure time within each header to account for variable depth between frames.

## Automatic Aperture and PSF photometry

The pipeline performs both automatic aperture and PSF photometry of sources in the stacked image using `potpyri/primitives/photometry.py`. Firstly, the pipeline will detect sources within the image and determine statistics within a fiducial aperture radius using `photutils.aperture.ApertureStats`.  This initial table of aperture photometry is saved within the `*stk.fits` image as the extension `APPPHOT`.

Next, based on cuts on roundness, FWHM, and signal-to-noise, the pipeline will define a list of bright stars with which it calculates an empirical PSF model.  The final list of PSF stars is saved in the `*stk.fits` file as the `PSFSTARS` extension.  The PSF itself will be generated from the extracted data around these PSF stars using `photutils.psf.EPSFBuilder`.  A stamp of the final effective PSF is saved in the `*stk.fits` file as the `PSF` extension.  A final FWHM is empirically calculated from the effective PSF model using an `astropy.modeling.functional_models.Moffat2D` profile and saved to the `*stk.fits` header.

Once the pipeline has determined the PSF, it will then calculate PSF photometry for all the originally extracted sources with `photutils.psf.PSFPhotometry` and write them to the `PSFPHOT` extension of the `*stk.fits` file.  In addition, as a data quality check a residual image is saved to the `*stk.fits` file showing the original science data with each of these sources subtracted using the PSF model as the `RESIDUAL` extension.

## Flux Calibration

After the PSF photometry has been calculated, the pipeline will then calculate a zero point based on aperture photometry by downloading a catalog of standard stars.  Currently supported catalogs are Pan-STARRS, SDSS, SkyMapper, and 2MASS, covering most optical and infrared bands in the north and south.

The zero point and associated uncertainty are saved in the `*stk.fits` files as `ZPTMAG` and `ZPTMUCER`, respectively.  In addition, an approximate limiting magnitude is calculated for each `*stk.fits` file using the `FWHM` and `SKYSIG` (the approximate standard deviation in sky background per pixel) header keywords and saved to the header as `M3SIGMA`, `M5SIGMA`, and `M10SIGMA` (3-, 5-, and 10-sigma limiting magnitudes).

## Quality checks

The pipeline will automatically create a log file in the `red/log` folder named after the instrument, followed by the date and time at the start of the run (UTC). The log file records each step in the pipeline and the reduction process. Because it can cover multiple targets, we recommend keeping a separate log per target with key information for your analysis. Before, during, and after the pipeline has run, perform the following quality-control checks and note them in the target log:

### Data Processing Checks

- Check the file list: The first thing the pipeline will do is to create a file list containing information about the files. Files can be misclassified, which can be corrected via direct editing of the file list and reprocessing with `--no-redo-sort`, direct editing of your file headers to correct the misclassification, or submitting an issue to debug POTPyRI-related classification issues.

- Check the stacks: After the pipeline has reduced and stacked the files for a target, it is recommended that you check your final stacked images to assess data quality. Check the FITS image named after the target name. If the image is of poor quality with very few stars, WCS and flux calibration will be affected.  If you find your stack has a large number of image artifacts, please contact us or open a GitHub issue with the instrument and specific error indicated.

- Check WCS solution: In addition to individual image WCS solutions, the pipeline calculates an automatic WCS solution on the final stack. The rms on the astrometry is typically $<$1 pixel, and the stars are aligned well near your target and across the image.  Check the image alignment quality via the `RADISP` and `DEDISP` keywords in each image header.  If the WCS is off by a large number of pixels, the zero point will not be calculated correctly and the reduction should be aborted. In that case, please contact us or open a GitHub issue with the instrument and specific error indicated.

- Check the PSF: The PSF calculated by the pipeline will determine whether the PSF photometry and the zero point can be trusted. The PSF is saved in the `*.epsf.png` file, and will usually look like a 2D Gaussian. The pipeline will also report the x and y sigma of a 2D Gaussian fitted to the PSF. If the PSF looks good and the x and y sigma values reported are similar, the PSF photometry can be trusted. If the PSF doesn't look Gaussian-like or the reported x and y sigma values are very different or negative, there was an issue determining the PSF. If there are many good stars in the final stack that could have been used to determine the PSF, please contact us or open a GitHub issue with the instrument and specific error indicated.

- Check the zero point calculation: If the PSF photometry is reliable, the zero point calculation should be good if a large number of stars is used to calculate it. In general, zero points calculated with $>$10 stars will provide relatively small (10$--$20 mmag level or smaller) zero point uncertainty. If there are more stars in the field that the pipeline should have used, or the pipeline reported that no stars were found but there is coverage in the corresponding catalog, please contact us or open a GitHub issue with the instrument and specific error indicated.

### Data over multiple nights

POTPyRI can stack data from multiple nights, although it will not use unique calibration frames for each stacked image.  To stack multiple nights, move all files to a single `raw` directory and run the pipeline.

## Documentation

- **Online:** Full API documentation (generated from docstrings) is at [https://ciera-transients.github.io/POTPyRI/](https://ciera-transients.github.io/POTPyRI/).
- **Local build:** Install with `pip install -e ".[docs]"` and run `cd docs && make html` (on Windows use `make.bat html`). Output is in `docs/build/html/`.
- **Deployment:** Docs are built and deployed to GitHub Pages on pushes to `main`; see `docs/DEPLOY.md` for setup details.

## Contributing and support

We welcome contributions and feedback. See [CONTRIBUTING.md](CONTRIBUTING.md) for:

- How to **report issues** or request features (GitHub Issues)
- How to **contribute** code or documentation (pull requests)
- Where to **seek support** (issues or developer contact)

Development often uses the `dev` branch; pull requests may target `main` or `dev` as noted in the repository. If you encounter an error in the pipeline, have a special data setup that does not run through the pipeline, or wish to add an instrument, please open an issue on [GitHub](https://github.com/CIERA-Transients/POTPyRI/issues) or contact the developers at `ckilpatrick@northwestern.edu`.

## License

This project is licensed under the GPL-3.0 license; see [LICENSE.txt](LICENSE.txt) for the full terms.
