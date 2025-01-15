import numpy as np
import pytest

import src.valladopy.astro.celestial.moon as moon

from ...conftest import DEFAULT_TOL, custom_isclose, custom_allclose


@pytest.fixture
def ttdb():
    # Julian centuries from J2000
    return -0.013641341546885694


def test_moon_latlon(ttdb):
    lon, lat = moon.get_moon_latlon(ttdb)
    assert np.isclose(np.degrees(lon), 138.13542521014608, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(lat), 358.86697486371054, rtol=DEFAULT_TOL)


def test_get_geodetic_dir_cosines(ttdb):
    lmn = moon.get_geodetic_dir_cosines(ttdb)
    lmn_exp = (-0.7445787078367155, 0.6200471059978293, 0.247273399661030)
    assert custom_allclose(lmn, lmn_exp, rtol=DEFAULT_TOL)


def test_position():
    # Vallado 2007, Example 5-3
    jd = 2449470.5
    rmoon, rtasc, decl = moon.position(jd)

    # Expected values
    rmoon_expected = [-134240.61111304158, -311571.55616382667, -126693.77078952147]
    assert np.allclose(rmoon, rmoon_expected, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(rtasc), -113.30879478775752, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(decl), -20.477717465404574, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "lon, moonrise_expected, moonset_expected, moonphaseang_expected",
    [
        # Example 5-4
        (0, 4.518799123473629, 18.51808524604512, 6.2177837312775415),
        # Non-zero longitude
        (np.radians(-74.3), 9.678940118281297, 23.593826737778475, 6.2615587982527865),
    ],
)
def test_rise_set(lon, moonrise_expected, moonset_expected, moonphaseang_expected):
    # Vallado 2007, Example 5-4
    jd = 2451046.5
    latgd = np.radians(40)
    moonrise, moonset, moonphaseang = moon.rise_set(jd, latgd, lon)

    # Expected values
    assert np.isclose(moonrise, moonrise_expected, rtol=DEFAULT_TOL)
    assert np.isclose(moonset, moonset_expected, rtol=DEFAULT_TOL)
    assert np.isclose(moonphaseang, moonphaseang_expected, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "elev, illum_expected",
    [
        (np.radians(100), 0.11663785663534358),
        (np.radians(15), 0.016132194140321483),
        (np.radians(3), 0.0024392692225373856),
        (np.radians(-1), 0),
    ],
)
def test_illumination(elev, illum_expected):
    f = np.radians(45)
    assert custom_isclose(moon.illumination(f, elev), illum_expected)
