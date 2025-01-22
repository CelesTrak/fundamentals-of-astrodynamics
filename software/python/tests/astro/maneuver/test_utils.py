import numpy as np
import pytest

import src.valladopy.astro.maneuver.utils as utils

from ...conftest import DEFAULT_TOL


@pytest.fixture
def a():
    """Semi-major axis of the orbit in km"""
    return 7000


@pytest.fixture
def r():
    """Radius of the orbit in km"""
    return 6800


def test_time_of_flight(a):
    assert np.isclose(utils.time_of_flight(a), 2914.258319939692, rtol=DEFAULT_TOL)


def test_specific_mech_energy(a):
    sme = utils.specific_mech_energy(a)
    assert np.isclose(sme, -28.471460107142857, rtol=DEFAULT_TOL)


def test_velocity_mag(a, r):
    assert np.isclose(utils.velocity_mag(r, a), 7.7648247730584705, rtol=DEFAULT_TOL)


def test_angular_velocity(a):
    w = utils.angular_velocity(a)
    assert np.isclose(w, 0.0010780076124668337, rtol=DEFAULT_TOL)


def test_deltav():
    v1, v2 = 1.23, 4.56
    theta = np.radians(15)
    assert np.isclose(utils.deltav(v1, v2, theta), 3.386905734002608, rtol=DEFAULT_TOL)


def test_semimajor_axis(r):
    e = 0.1
    nu = np.radians(45)
    a = utils.semimajor_axis(r, e, nu)
    assert np.isclose(a, 7354.376374956417, rtol=DEFAULT_TOL)


def test_period(a):
    assert np.isclose(utils.period(a), 5828.516639879384, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "z, u_expected", [(0.1234, 0.26358469539982904), (0, 1e-6), (1, 0.999999)]
)
def test_lowuz(z, u_expected):
    u = utils.lowuz(z)
    assert np.isclose(u, u_expected, rtol=DEFAULT_TOL)
