"""Unit tests for flux-calibration catalogs, DECaLS helpers, and point-source cuts."""
from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from potpyri.primitives import absphot
from potpyri.utils import catalogs


# ---------------------------------------------------------------------------
# find_catalog / aliases
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    'name,filt,expected_id',
    [
        ('PS1', 'r', 'II/349'),
        ('SDSS', 'g', 'V/154'),
        ('2MASS', 'J', 'II/246'),
        ('UKIRT', 'K', 'II/319'),
        ('SKYMAPPER', 'r', 'II/379/smssdr4'),
        ('DES', 'i', 'II/357'),
        ('DECALS', 'r', catalogs.DECALS_TRACTOR_TABLE),
        ('LEGACY', 'g', catalogs.DECALS_TRACTOR_TABLE),
        ('decals_dr10', 'z', catalogs.DECALS_TRACTOR_TABLE),
    ],
)
def test_find_catalog_flux_calibration_ids(name, filt, expected_id):
    cat, cid, ra, dec, mag, err = catalogs.find_catalog(name, filt, 180.0, 0.0)
    assert cid == expected_id
    assert ra and dec and mag and err
    assert cat == catalogs.normalize_flux_catalog_name(name)


def test_find_catalog_decals_unsupported_filter():
    cat, cid, ra, dec, mag, err = catalogs.find_catalog('DECALS', 'J', 0.0, 0.0)
    assert cid is None
    assert mag is None


def test_find_catalog_des_griz():
    for filt in ('g', 'r', 'i', 'z'):
        cat, cid, ra, dec, mag, err = catalogs.find_catalog('DES', filt, 10.0, -30.0)
        assert cat == 'DES'
        assert cid == 'II/357'
        assert mag == f'{filt}mag'


def test_find_catalog_skymapper_u_south_returns_canonical_name():
    cat, cid, ra, dec, mag, err = catalogs.find_catalog('PS1', 'u', 180.0, -40.0)
    assert cat == 'SKYMAPPER'
    assert cid == 'II/379/smssdr4'


def test_flux_calibration_catalogs_have_point_source_cuts():
    """Every flux ZP catalog ID from find_catalog must define a point-source cut."""
    for name in catalogs.FLUX_CALIBRATION_CATALOGS:
        # Pick a supported filter per catalog.
        filt = {
            'PS1': 'r', 'SDSS': 'r', '2MASS': 'J', 'UKIRT': 'J',
            'SKYMAPPER': 'r', 'DES': 'r', 'DECALS': 'r',
        }[name]
        cat, cid, *_ = catalogs.find_catalog(name, filt, 0.0, 0.0)
        assert cid is not None, name
        assert cid in catalogs.POINT_SOURCE_CUTS, (
            f'{name} catalog_id={cid!r} missing from POINT_SOURCE_CUTS'
        )


# ---------------------------------------------------------------------------
# Point-source cuts (synthetic tables; no network)
# ---------------------------------------------------------------------------

def test_apply_point_source_cut_ps1_psf_kron():
    tab = Table({
        'rmag': [18.0, 18.0, 19.0],
        'rKmag': [17.95, 17.5, 18.95],  # diffs: 0.05, 0.5, 0.05
        'e_rmag': [0.01, 0.01, 0.01],
    })
    out = catalogs.apply_point_source_cut(tab, 'II/349', 'r', mag_col='rmag')
    assert len(out) == 2
    np.testing.assert_allclose(out['rmag'], [18.0, 19.0])


def test_apply_point_source_cut_sdss_class_star():
    tab = Table({
        'class': [6, 3, 6, 6],
        'gmag': [18.0, 18.1, 19.0, 20.0],
    })
    out = catalogs.apply_point_source_cut(tab, 'V/154', 'g')
    assert len(out) == 3
    assert set(out['class']) == {6}


def test_apply_point_source_cut_2mass_xflg():
    tab = Table({'Xflg': [0, 1, 0], 'Jmag': [12.0, 12.1, 13.0]})
    out = catalogs.apply_point_source_cut(tab, 'II/246', 'J')
    assert len(out) == 2


def test_apply_point_source_cut_ukirt_merged_class():
    tab = Table({'mergedClass': [-1, 1, -1, 0], 'Kmag': [14.0, 14.1, 15.0, 16.0]})
    out = catalogs.apply_point_source_cut(tab, 'II/319', 'K')
    assert len(out) == 2
    assert np.all(out['mergedClass'] == -1)


def test_apply_point_source_cut_skymapper_classstar():
    tab = Table({'ClassStar': [0.95, 0.2, 0.99, 0.85], 'rPSF': [17.0, 17.1, 18.0, 19.0]})
    out = catalogs.apply_point_source_cut(tab, 'II/379/smssdr4', 'r')
    assert len(out) == 2
    assert np.all(out['ClassStar'] > 0.9)


def test_apply_point_source_cut_des_sg():
    tab = Table({
        'S/Gr': [0.8, 0.1, 0.6, 0.4],
        'rmag': [18.0, 18.1, 19.0, 20.0],
    })
    out = catalogs.apply_point_source_cut(tab, 'II/357', 'r')
    assert len(out) == 2
    assert np.all(out['S/Gr'] > 0.5)


def test_apply_point_source_cut_decals_type_psf():
    tab = Table({
        'type': ['PSF', 'REX', 'PSF', 'EXP'],
        'mag': [18.0, 18.1, 19.0, 20.0],
        'mag_err': [0.01, 0.01, 0.01, 0.01],
    })
    out = catalogs.apply_point_source_cut(tab, 'ls_dr10.tractor', 'r', mag_col='mag')
    assert len(out) == 2
    assert set(out['type']) == {'PSF'}


def test_apply_point_source_cut_missing_column_returns_unchanged():
    tab = Table({'rmag': [18.0, 19.0], 'e_rmag': [0.01, 0.01]})
    out = catalogs.apply_point_source_cut(tab, 'II/349', 'r', mag_col='rmag')
    assert len(out) == len(tab)


@pytest.mark.parametrize('catalog_id', list(catalogs.POINT_SOURCE_CUTS.keys()))
def test_point_source_cut_defined_for_every_rule(catalog_id):
    """Each POINT_SOURCE_CUTS entry has a runnable method and description."""
    rule = catalogs.POINT_SOURCE_CUTS[catalog_id]
    assert 'method' in rule
    assert 'description' in rule
    assert rule['method'] in {
        'psf_kron', 'class_equals', 'score_above', 'type_equals', 'des_sg',
    }


# ---------------------------------------------------------------------------
# DECaLS helpers + product schema consistency
# ---------------------------------------------------------------------------

def test_nanomaggy_to_ab_mag():
    flux = np.array([1.0, 10.0, 100.0])  # nanomaggy
    ivar = np.array([100.0, 100.0, 100.0])
    mag, magerr = catalogs.nanomaggy_to_ab_mag(flux, ivar)
    np.testing.assert_allclose(mag, 22.5 - 2.5 * np.log10(flux))
    assert np.all(np.isfinite(magerr))
    assert np.all(magerr > 0)


def test_query_decals_region_builds_consistent_product(monkeypatch):
    """Mock Data Lab TAP; ensure output matches absphot product columns."""

    class _FakeJob:
        def get_results(self):
            return Table({
                'ra': [150.0, 150.01, 150.02],
                'dec': [2.0, 2.01, 2.02],
                'type': ['PSF', 'REX', 'PSF'],
                'mag': [18.0, 18.5, 19.0],
                'flux_r': [10.0, 8.0, 5.0],
                'flux_ivar_r': [100.0, 100.0, 100.0],
            })

    class _FakeTap:
        def __init__(self, url):
            self.url = url

        def launch_job(self, query, maxrec=None):
            assert 'type = \'PSF\'' in query or "type = 'PSF'" in query
            assert 'BETWEEN' in query.upper()
            assert 'ls_dr10.tractor' in query or 'ls_dr9.tractor' in query
            return _FakeJob()

    monkeypatch.setattr(
        'astroquery.utils.tap.core.TapPlus', _FakeTap, raising=False)
    # Force import path used inside query_decals_region
    import astroquery.utils.tap.core as tap_core
    monkeypatch.setattr(tap_core, 'TapPlus', _FakeTap)

    center = SkyCoord(150.0, 2.0, unit='deg')
    out = catalogs.query_decals_region(center, 0.5 * u.deg, 'r', log=None)
    assert out is not None
    assert set(['ra', 'dec', 'mag', 'mag_err']).issubset(out.colnames)
    # SQL already requested type=PSF; fake job still has REX but product keeps rows
    # that pass flux cuts — type cut is also applied in get_catalog.
    assert len(out) >= 1
    assert np.all(np.isfinite(out['mag']))
    assert np.all(out['mag_err'] > 0)


def test_get_catalog_decals_product_matches_ps1_schema(monkeypatch):
    """absphot.get_catalog(DECALS) returns same core columns as a VizieR path."""

    def _fake_decals(center, width, filt, log=None, table=None):
        return Table({
            'ra': [10.0, 10.01],
            'dec': [-20.0, -20.01],
            'type': ['PSF', 'PSF'],
            'mag': [17.5, 18.2],
            'mag_err': [0.02, 0.03],
        })

    monkeypatch.setattr(catalogs, 'query_decals_region', _fake_decals)

    cal = absphot.absphot()
    coords = SkyCoord([10.0, 10.01], [-20.0, -20.01], unit='deg')
    cat, name, cid = cal.get_catalog(coords, 'DECALS', 'r', log=None)
    assert name == 'DECALS'
    assert cid == catalogs.DECALS_TRACTOR_TABLE
    assert list(cat.colnames) == ['ra', 'dec', 'mag', 'mag_err']
    assert len(cat) == 2


def test_get_catalog_ps1_applies_kron_cut(monkeypatch):
    """PS1 path still applies PSF–Kron cut and returns uniform product columns."""

    def _fake_vizier(center, width, catalog_id, columns, log=None):
        assert catalog_id == 'II/349'
        assert 'rKmag' in columns
        return Table({
            'RAJ2000': [180.0, 180.01, 180.02],
            'DEJ2000': [0.0, 0.01, 0.02],
            'rmag': [18.0, 18.0, 19.0],
            'e_rmag': [0.01, 0.01, 0.01],
            'rKmag': [17.95, 17.5, 18.95],
        })

    monkeypatch.setattr(catalogs, 'query_vizier_region', _fake_vizier)
    cal = absphot.absphot()
    coords = SkyCoord([180.0, 180.01], [0.0, 0.01], unit='deg')
    cat, name, cid = cal.get_catalog(coords, 'PS1', 'r', log=None)
    assert name == 'PS1'
    assert cid == 'II/349'
    assert {'ra', 'dec', 'mag', 'mag_err'}.issubset(set(cat.colnames))
    assert len(cat) == 2  # Kron cut removes the 0.5 mag outlier


def test_get_catalog_all_flux_catalogs_product_schema(monkeypatch):
    """Every flux-calibration catalog yields ra/dec/mag/mag_err via get_catalog."""

    def _fake_vizier(center, width, catalog_id, columns, log=None):
        # Build a minimal valid table for whichever columns were requested.
        data = {}
        for col in columns:
            if col.lower().startswith(('ra', 'de')):
                data[col] = [1.0, 1.01]
            elif 'flg' in col.lower() or col in ('Xflg', 'X'):
                data[col] = [0, 0]
            elif col in ('class', 'cl', 'Class'):
                data[col] = [6, 6]
            elif col in ('mergedClass', 'mergedclass'):
                data[col] = [-1, -1]
            elif 'ClassStar' in col or 'class_star' in col.lower():
                data[col] = [0.95, 0.96]
            elif col.startswith('S/G'):
                data[col] = [0.8, 0.9]
            elif 'Kmag' in col and col.startswith(('g', 'r', 'i', 'z', 'y')):
                # Kron mags slightly brighter/fainter within 0.1
                data[col] = [17.95, 18.95]
            elif col.startswith('e_') or col.endswith('err'):
                data[col] = [0.02, 0.03]
            else:
                data[col] = [18.0, 19.0]
        return Table(data)

    def _fake_decals(center, width, filt, log=None, table=None):
        return Table({
            'ra': [1.0, 1.01],
            'dec': [2.0, 2.01],
            'type': ['PSF', 'PSF'],
            'mag': [18.0, 19.0],
            'mag_err': [0.02, 0.03],
        })

    monkeypatch.setattr(catalogs, 'query_vizier_region', _fake_vizier)
    monkeypatch.setattr(catalogs, 'query_decals_region', _fake_decals)

    cal = absphot.absphot()
    coords = SkyCoord([1.0, 1.01], [2.0, 2.01], unit='deg')
    filt_by_cat = {
        'PS1': 'r', 'SDSS': 'r', '2MASS': 'J', 'UKIRT': 'J',
        'SKYMAPPER': 'r', 'DES': 'r', 'DECALS': 'r',
    }
    for name in catalogs.FLUX_CALIBRATION_CATALOGS:
        cat, cname, cid = cal.get_catalog(
            coords, name, filt_by_cat[name], log=None)
        assert cat is not None, name
        assert cname == name
        assert cid is not None
        assert {'ra', 'dec', 'mag', 'mag_err'}.issubset(set(cat.colnames)), name
        assert len(cat) >= 1, name
