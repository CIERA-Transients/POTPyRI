"""Photometry primitive, photloop worker, and photloop API."""
from __future__ import annotations

from potpyri.primitives.base_primitive import BasePrimitive

from .core import do_phot


def _photloop_worker(stack, phot_sn_min=3.0, phot_sn_max=40.0, fwhm_init=5.0, log=None):
    """Try PSF photometry with decreasing S/N threshold until do_phot succeeds."""
    signal_to_noise = phot_sn_max

    epsf = None
    fwhm = None
    while signal_to_noise > phot_sn_min:

        if log:
            log.info(f'Trying PSF generation with S/N={signal_to_noise}')
        star_param = {'snthresh_psf': signal_to_noise*2.0,
                      'fwhm_init': fwhm_init,
                      'snthresh_final': signal_to_noise}
        try:
            do_phot(stack, star_param=star_param)
        except Exception as e:
            if log:
                log.error(e)
            signal_to_noise = signal_to_noise / 2.0
            continue
        break


class PhotometryPrimitive(BasePrimitive):
    """PSF/aperture photometry on a stacked image (see :func:`photloop`)."""

    ARGUMENTS = {
        **BasePrimitive.ARGUMENTS,
        'phot_sn_min': 3.0,
        'phot_sn_max': 40.0,
        'fwhm_init': 5.0,
    }

    def _perform(self):
        _photloop_worker(
            self.input,
            phot_sn_min=self.phot_sn_min,
            phot_sn_max=self.phot_sn_max,
            fwhm_init=self.fwhm_init,
            log=self.log,
        )
        return {'output': None}


def photloop(stack, phot_sn_min=3.0, phot_sn_max=40.0, fwhm_init=5.0, log=None):
    """Try PSF photometry with decreasing S/N threshold until do_phot succeeds.

    Delegates to :class:`PhotometryPrimitive`.

    Parameters
    ----------
    stack : str
        Path to stack FITS (SCI, MASK, ERROR).
    phot_sn_min : float, optional
        Minimum S/N threshold to try. Default is 3.0.
    phot_sn_max : float, optional
        Initial S/N threshold. Default is 40.0.
    fwhm_init : float, optional
        Initial FWHM for star finding. Default is 5.0.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    None
        Stack FITS is updated in place when do_phot succeeds.
    """
    return PhotometryPrimitive(
        phot_sn_min=phot_sn_min,
        phot_sn_max=phot_sn_max,
        fwhm_init=fwhm_init,
    ).apply(input=stack, log=log)
