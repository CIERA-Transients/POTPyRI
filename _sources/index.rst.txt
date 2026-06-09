*********************
POTPyRI Documentation
*********************

POTPyRI is a data reduction pipeline for imaging from large aperture telescopes.
It provides instrument-specific reduction workflows for bias, dark, and flat
calibration; image alignment and stacking; WCS solving with astrometry.net;
photometry with photutils; and flux calibration with astroquery.

The package includes nine instrument modules (BINOSPEC, DEIMOS, F2, FOURSTAR, GMOS, IMACS, LRIS, MMIRS, MOSFIRE), primitives for calibration, image processing, file sorting, WCS solving, photometry, and absolute photometry, plus CLI scripts for the main pipeline and archive data download. For installation, usage, and contributing, see the repository `README <https://github.com/CIERA-Transients/POTPyRI>`_ and `CONTRIBUTING.md <https://github.com/CIERA-Transients/POTPyRI/blob/main/CONTRIBUTING.md>`_.

Contents
========

.. toctree::
   :maxdepth: 2

   installation
   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
