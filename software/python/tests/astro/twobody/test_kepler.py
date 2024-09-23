import numpy as np
import pytest

import src.valladopy.astro.twobody.kepler as kepler
from ...conftest import custom_isclose


@pytest.mark.parametrize(
    "ecc, incl, expected_orbit_type",
    [
        (0.0, 0.0, kepler.OrbitType.CIR_EQUATORIAL),  # Circular equatorial
        (0.0, np.pi, kepler.OrbitType.CIR_EQUATORIAL),  # Circular equatorial
        (0.0, 0.1, kepler.OrbitType.CIR_INCLINED),  # Circular inclined
        (0.1, 0.0, kepler.OrbitType.EPH_EQUATORIAL),  # Elliptical equatorial
        (0.1, np.pi, kepler.OrbitType.EPH_EQUATORIAL),  # Elliptical equatorial
        (0.1, 0.1, kepler.OrbitType.EPH_INCLINED),  # Elliptical inclined
        (1.1, 0.0, kepler.OrbitType.EPH_EQUATORIAL),  # Hyperbolic equatorial
        (1.1, np.pi, kepler.OrbitType.EPH_EQUATORIAL),  # Hyperbolic equatorial
        (1.1, 0.1, kepler.OrbitType.EPH_INCLINED),  # Hyperbolic inclined
    ],
)
def test_determine_orbit_type(ecc, incl, expected_orbit_type):
    assert kepler.determine_orbit_type(ecc, incl) == expected_orbit_type


@pytest.mark.parametrize(
    "ecc, e0, expected_m, expected_nu",
    [
        # Test circular orbit (ecc = 0)
        (0, np.radians(42.42), np.radians(42.42), np.radians(42.42)),
        # Test elliptical orbit (0 < ecc < 1)
        (0.4, np.radians(334.566986), 6.011077700843401, -0.6639000611931205),
        # Test hyperbolic orbit (ecc > 1)
        (1.05, np.radians(123.45), 2.3123894903995055, -2.752304611108218),
        # Test parabolic orbit (ecc = 1)
        (1, np.radians(30), 0.5714479680061689, 0.9646958142020499),
    ],
)
def test_newtone(ecc, e0, expected_m, expected_nu):
    m, nu = kepler.newtone(ecc, e0)
    assert custom_isclose(m, expected_m)
    assert custom_isclose(nu, expected_nu)


def test_newtonnu():
    # TODO: test all cases
    ecc = 0.1
    nu = np.radians(45)
    e0, m = kepler.newtonnu(ecc, nu)
    assert abs(np.degrees(e0) - 41.078960346507934) < 1e-12
    assert abs(np.degrees(m)) - 37.31406335764441 < 1e-12


@pytest.mark.parametrize(
    "ecc, m, expected_e0, expected_nu",
    [
        # Test elliptical orbit (0 < ecc < 1)
        # Vallado Ex. 2-1
        (0.4, np.radians(235.4), 3.8486617450971616, -2.6674915592323702),
        # Test hyperbolic orbit (ecc > 1)
        (1.5, np.pi / 2, 1.449742150875263, 1.8916928725677562),
        # Test parabolic orbit (ecc = 1)
        (1.0, np.pi / 4, 0.6804016596037228, 1.1949025292471007),
        # Test circular orbit (ecc = 0)
        (0, np.pi / 3, np.pi / 3, np.pi / 3),  # e0 and nu should equal m
        # Test edge cases
        (0.9999, np.pi, np.pi, np.pi),  # Near-parabolic
        (0, 0, 0, 0),  # Circular and m = 0
    ],
)
def test_newtonm(ecc, m, expected_e0, expected_nu):
    e0, nu = kepler.newtonm(ecc, m)
    assert custom_isclose(e0, expected_e0, atol=1e-11)
    assert custom_isclose(nu, expected_nu)


@pytest.mark.parametrize(
    # TODO: find test case where np.abs(alpha) < SMALL
    "ro, vo, dtsec, r_expected, v_expected",
    [
        # alpha >= SMALL
        (
            [-6518.1083, -2403.8479, -22.1722],
            [2.604057, -7.105717, -0.263218],
            120,
            [-6150.8180772417, -3233.5647020192746, -53.47890911642744],
            [3.5085906394995146, -6.703227854017405, -0.2578170063354219],
        ),
        # alpha < SMALL
        (
            [15912.1061, 1352024060, 0],
            [-0.000783895533, 31889.3576, 0],
            100000,
            [15833.716546696181, 4540959819.99955, 0],
            [-0.000783895533049524, 31889.357599993506, 0],
        ),
        # alpha < SMALL, not converged
        (
            [15912.1061, 1352024060, 0],
            [-0.000783895533, 31889.3576, 0],
            5485549291.841488,
            [0, 0, 0],
            [0, 0, 0],
        ),
    ],
)
def test_kepler(ro, vo, dtsec, r_expected, v_expected):
    r, v = kepler.kepler(ro, vo, dtsec)
    assert np.allclose(r, r_expected, rtol=1e-12)
    assert np.allclose(v, v_expected, rtol=1e-12)
