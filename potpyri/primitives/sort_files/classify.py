"""Header-based classification for raw FITS (bad, spec, flat, bias, dark, science)."""
from __future__ import annotations

import re

import numpy as np


def is_bad(hdr, tel):
    """Return True if the header matches bad_keywords/bad_values or has invalid binning.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.
    tel : Instrument
        Instrument instance (bad_keywords, bad_values, get_binning).

    Returns
    -------
    bool
        True if file should be excluded as bad.
    """
    keywords = tel.bad_keywords
    values = tel.bad_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    bad = np.any([bool(re.search(v, str(hdr[k]).lower()))
        for k,v in zip(keywords,values)])

    binn = str(tel.get_binning(hdr))
    if len(binn)>1:
        # Check if telescope is binned the same in all directions, we do not
        # want to reduce images with variable binning in different directions
        bad = not binn == len(binn) * binn[0]

    return(bad)

def is_spec(hdr, tel):
    """Return True if the header matches spectroscopic observation keywords.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.
    tel : Instrument
        Instrument instance (spec_keywords, spec_values).

    Returns
    -------
    bool
        True if file is spectroscopic.
    """
    keywords = tel.spec_keywords
    values = tel.spec_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    spec = np.all([bool(re.search(v, str(hdr[k]).lower()))
        for k,v in zip(keywords,values)])

    return(spec)

def is_flat(hdr, tel):
    """Return True if the header matches flat-field observation keywords.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.
    tel : Instrument
        Instrument instance (flat_keywords, flat_values).

    Returns
    -------
    bool
        True if file is a flat.
    """
    keywords = tel.flat_keywords
    values = tel.flat_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    flat = np.all([bool(re.search(v, str(hdr[k]).lower()))
        for k,v in zip(keywords,values)])

    return(flat)

def is_dark(hdr, tel):
    """Return True if the header matches dark observation keywords and valid binning.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.
    tel : Instrument
        Instrument instance (dark_keywords, dark_values, get_binning).

    Returns
    -------
    bool
        True if file is a dark.
    """
    keywords = tel.dark_keywords
    values = tel.dark_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    dark = np.all([bool(re.search(v, str(hdr[k]).lower()))
        for k,v in zip(keywords,values)])

    # Similar to bad, require that dark have equivalent binning in both dirs
    if dark:
        binn = str(tel.get_binning(hdr))
        if len(binn)>1:
            # Check if telescope is binned the same in all directions, we do not
            # want to reduce images with variable binning in different directions
            dark = binn == len(binn) * binn[0]

    return(dark)

def is_bias(hdr, tel):
    """Return True if the header matches bias observation keywords and valid binning.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.
    tel : Instrument
        Instrument instance (bias_keywords, bias_values, get_binning).

    Returns
    -------
    bool
        True if file is a bias.
    """
    keywords = tel.bias_keywords
    values = tel.bias_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    bias = np.all([bool(re.search(v, str(hdr[k]).lower()))
        for k,v in zip(keywords,values)])

    # Similar to bad, require that bias have equivalent binning in both dirs
    if bias:
        binn = str(tel.get_binning(hdr))
        if len(binn)>1:
            # Check if telescope is binned the same in all directions, we do not
            # want to reduce images with variable binning in different directions
            bias = binn == len(binn) * binn[0]

    return(bias)

def is_science(hdr, tel):
    """Return True if the header matches science observation keywords and min exptime.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.
    tel : Instrument
        Instrument instance (science_keywords, science_values, min_exptime, get_exptime).

    Returns
    -------
    bool
        True if file is science.
    """
    keywords = tel.science_keywords
    values = tel.science_values

    assert len(keywords)==len(values)
    if len(keywords)==0: return(False)

    science = np.all([bool(re.search(v, str(hdr[k]).lower()))
        for k,v in zip(keywords,values)])

    # Check minimum exposure time
    if tel.min_exptime:
        exptime = tel.get_exptime(hdr)
        if exptime < tel.min_exptime:
            return(False)

    return(science)
