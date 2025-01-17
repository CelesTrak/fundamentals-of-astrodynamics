import logging
import pytest

import numpy as np

import src.valladopy.astro.celestial.utils as utils
import src.valladopy.constants as const

from ...conftest import DEFAULT_TOL


@pytest.fixture
def t():
    # Julian centuries from J2000
    return -0.013641341546885694


@pytest.fixture
def coe_shadow():
    # Test case given in Vallado/Neta shadow paper
    # Neta, B., and Vallado, D. (1998) On Satellite Umbra/Penumbra Entry and Exit
    # Positions, Journal of the Astronautical Sciences, 46, No. 1, 91–104.
    e = 0.002
    a = (1.029 * const.RE) / (1 - e**2)
    raan, w = 0, 0
    i = np.radians(63.4)
    return [a, e, i, raan, w]


@pytest.mark.parametrize(
    "r2, earth_model, los, tmin",
    [
        (
            [0, 5740.323, 3189.068],
            utils.EarthModel.ELLIPSOIDAL,
            False,
            0.5082248650848982,
        ),
        (
            [0, 5740.323, 3189.068],
            utils.EarthModel.SPHERICAL,
            False,
            0.5082352992389487,
        ),
        (
            [122233179.72368076, -76150708.25425531, -33016373.913704105],
            utils.EarthModel.ELLIPSOIDAL,
            True,
            -2.334265707720442e-05,
        ),
        (
            [122233179.72368076, -76150708.25425531, -33016373.913704105],
            utils.EarthModel.SPHERICAL,
            True,
            -2.3290642130715177e-05,
        ),
    ],
)
def test_in_sight(r2, earth_model, los, tmin, caplog):
    # Vallado 2022, Example 5-6
    r1 = [0, -4464.696, -5102.509]

    # Call function with logging
    with caplog.at_level(logging.DEBUG):
        assert utils.in_sight(r1, r2, earth_model) == los
        assert f"Minimum parametric value (tmin): {tmin}" in caplog.messages[0]


def test_in_shadow_simple():
    # Test against values from Example 12.8 in Curtis
    r_sat = [2817.899, -14110.473, -7502.672]
    r_sun = [-11747041, 139486985, 60472278]
    assert utils.in_shadow_simple(r_sat, r_sun)


def test_in_shadow():
    r_eci = [-41260.1818237031, 8684.15782134066, 0]
    r_sun = [148470363.19330865, -9449738.11151353, -4096753.810182002]
    in_umbra, in_penumbra = utils.in_shadow(r_eci, r_sun)
    assert in_umbra
    assert in_penumbra


def test_cylindrical_shadow_roots(coe_shadow):
    # Test against paper values (use given temp params)
    a, e, *_ = coe_shadow
    beta_1 = 0.459588
    beta_2 = -0.6807135

    # Expected roots
    roots_expected = [
        0.9515384802192421,
        0.6383876664195322,
        -0.9573391650706946,
        -0.6284006781641529,
    ]

    # Call function
    roots = utils.cylindrical_shadow_roots(a, e, beta_1, beta_2)

    # Check results
    assert np.allclose(roots, roots_expected, rtol=DEFAULT_TOL)


def test_sun_ecliptic_parameters(t):
    mean_lon, mean_anomaly, ecliptic_lon = utils.sun_ecliptic_parameters(t)
    assert np.isclose(np.degrees(mean_lon), 149.36181814781156, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(mean_anomaly), 226.4523822485284, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(ecliptic_lon), 147.9940329397011, rtol=DEFAULT_TOL)


def test_obliquity_ecliptic(t):
    obliquity = utils.obliquity_ecliptic(t)
    assert np.isclose(np.degrees(obliquity), 23.439468394733744, rtol=DEFAULT_TOL)
