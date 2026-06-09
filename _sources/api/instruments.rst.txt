.. _api-instruments:

Instruments
===========

Instrument implementations and factory for POTPyRI. Each instrument module
defines a subclass of :class:`~potpyri.instruments.instrument.Instrument` with
detector keywords, calibration behavior, and file-sorting rules.

Package overview
----------------
.. automodule:: potpyri.instruments
   :members:
   :imported-members:
   :no-undoc-members:

Instrument base class
---------------------
.. autosummary::
   :toctree: generated/
   :template: autosummary/module.rst

   potpyri.instruments.instrument

Instrument modules
------------------
.. autosummary::
   :toctree: generated/
   :template: autosummary/module.rst

   potpyri.instruments.BINOSPEC
   potpyri.instruments.DEIMOS
   potpyri.instruments.F2
   potpyri.instruments.FOURSTAR
   potpyri.instruments.GMOS
   potpyri.instruments.IMACS
   potpyri.instruments.LRIS
   potpyri.instruments.MMIRS
   potpyri.instruments.MOSFIRE
