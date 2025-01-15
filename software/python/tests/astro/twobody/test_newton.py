import numpy as np
import pytest

import src.valladopy.astro.twobody.newton as newton
from ...conftest import custom_isclose, DEFAULT_TOL


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
    m, nu = newton.newtone(ecc, e0)
    assert custom_isclose(m, expected_m)
    assert custom_isclose(nu, expected_nu)


def test_newtonnu():
    # TODO: test all cases
    ecc = 0.1
    nu = np.radians(45)
    e0, m = newton.newtonnu(ecc, nu)
    assert np.isclose(np.degrees(e0), 41.078960346507934, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(m), 37.31406335764441, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "ecc, m, expected_e0, expected_nu",
    [
        # Test elliptical orbit (0 < ecc < 1)
        # Vallado Ex. 2-1
        (0.4, np.radians(235.4), 3.8486617450971616, -2.6674915592323702),
        # Test hyperbolic orbit (ecc > 1)
        (1.5, np.pi / 2, 1.449742150875263, 1.8916928725677562),
        # Test parabolic orbit (ecc = 1)
        (1, np.pi / 4, 0.6804016596037228, 1.1949025292471007),
        # Test circular orbit (ecc = 0)
        (0, np.pi / 3, np.pi / 3, np.pi / 3),  # e0 and nu should equal m
        # Test edge cases
        (0.9999, np.pi, np.pi, np.pi),  # Near-parabolic
        (0, 0, 0, 0),  # Circular and m = 0
    ],
)
def test_newtonm(ecc, m, expected_e0, expected_nu):
    e0, nu = newton.newtonm(ecc, m)
    assert custom_isclose(e0, expected_e0, atol=1e-11)
    assert custom_isclose(nu, expected_nu)
