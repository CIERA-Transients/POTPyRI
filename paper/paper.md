title: 'Pipeline for Optical/infrared Telescopes in Python for Reducing Images (POTPyRI)'

tags:
  - Python
  - astronomy
  - imaging
  - data reduction

authors:
  - name: K. Paterson^[corresponding author]
    orcid: 0000-0001-8340-3486
    affiliation: "1" #
  - name: C. D. Kilpatrick
    affiliation: 1
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

This pipeline was developed for the reduction and stacking of imaging data for a number of telescopes and instruments in the optical and near-infrared (NIR) bands. The purpose of the pipeline is to provide a semi-automated way to reduce imaging data, create a median stack with reliable astrometry and photometric calibration, and allow the user to perform manual aperture photometry in Python.

# Statement of need

The pipeline is written in Python (currently deployed and tested on Python 3.8) and uses packages from Astropy [@astropy:2013, @astropy:2018], ccdproc, astroquery [@Ginsburg2019] and photutils; with the additional use of SExtractor [@Bertin1996]. The code is available on GitHub at https://github.com/CIERA-Transients/POTPyRI, where instructions on the installation and detailed use can be found.

