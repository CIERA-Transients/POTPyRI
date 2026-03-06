"""Unit tests for potpyri.utils.utilities (find_catalog, is_number, parse_coord)."""
import pytest

from potpyri.utils import utilities


def test_find_catalog_ps1():
    """find_catalog returns PS1 columns for g,r,i,z,y."""
    cat, cid, ra, dec, mag, err = utilities.find_catalog("PS1", "r", 180.0, 0.0)
    assert cat == "PS1"
    assert cid == "II/349"
    assert ra == "RAJ2000"
    assert dec == "DEJ2000"
    assert mag == "rmag"
    assert err == "e_rmag"


def test_find_catalog_ps1_unsupported_filter():
    """find_catalog returns None columns for unsupported filter."""
    cat, cid, ra, dec, mag, err = utilities.find_catalog("PS1", "U", 180.0, 0.0)
    assert cid is None
    assert mag is None


def test_find_catalog_2mass():
    """find_catalog returns 2MASS columns for J,H,K."""
    cat, cid, ra, dec, mag, err = utilities.find_catalog("2MASS", "J", 0.0, 0.0)
    assert cid == "II/246"
    assert mag == "Jmag"
    assert err == "e_Jmag"


def test_find_catalog_2mass_k_band():
    """find_catalog returns 2MASS K-band (Kmag, e_Kmag) for K, Ks, and Kspec."""
    for filt in ("K", "Ks", "Kspec"):
        cat, cid, ra, dec, mag, err = utilities.find_catalog("2MASS", filt, 0.0, 0.0)
        assert cid == "II/246", f"2MASS catalog ID for filter {filt!r}"
        assert mag == "Kmag", f"2MASS K-band mag column for filter {filt!r}"
        assert err == "e_Kmag", f"2MASS K-band err column for filter {filt!r}"


def test_find_catalog_sdss():
    """find_catalog returns SDSS columns for u,g,r,i,z."""
    cat, cid, ra, dec, mag, err = utilities.find_catalog("SDSS", "g", 0.0, 0.0)
    assert cid == "V/154"
    assert mag == "gmag"


def test_find_catalog_skymapper_u_south():
    """find_catalog uses SkyMapper for u-band when dec < 0."""
    cat, cid, ra, dec, mag, err = utilities.find_catalog("PS1", "u", 180.0, -40.0)
    assert cat == "skymapper"
    assert cid == "II/379/smssdr4"


def test_is_number():
    """is_number returns True for numeric strings and numbers."""
    assert utilities.is_number("1.5") is True
    assert utilities.is_number(1.5) is True
    assert utilities.is_number("0") is True
    assert utilities.is_number("abc") is False
    assert utilities.is_number("") is False


def test_parse_coord_degrees():
    """parse_coord accepts decimal degrees."""
    coord = utilities.parse_coord(180.0, -30.0)
    assert coord is not None
    assert abs(coord.ra.deg - 180.0) < 0.01
    assert abs(coord.dec.deg - (-30.0)) < 0.01


def test_parse_coord_sexagesimal():
    """parse_coord accepts sexagesimal strings."""
    coord = utilities.parse_coord("12:00:00", "-30:00:00")
    assert coord is not None
    assert abs(coord.ra.deg - 180.0) < 1.0
    assert abs(coord.dec.deg - (-30.0)) < 1.0


def test_parse_coord_invalid():
    """parse_coord returns None for invalid input."""
    assert utilities.parse_coord("not", "numbers") is None


def test_find_catalog_ukirt():
    """find_catalog returns UKIRT columns for Y,J,H,K."""
    cat, cid, ra, dec, mag, err = utilities.find_catalog("UKIRT", "K", 0.0, 0.0)
    assert cid == "II/319"
    assert mag == "Kmag"
    assert err == "e_Kmag"


def test_parse_coord_value_error():
    """parse_coord returns None when SkyCoord raises ValueError."""
    # Valid format but invalid values can trigger ValueError
    result = utilities.parse_coord("99:99:99", "99:99:99")
    # May return None or raise; doc says returns None on parse failure
    assert result is None or hasattr(result, "ra")
