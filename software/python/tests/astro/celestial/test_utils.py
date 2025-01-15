import logging
import pytest

import numpy as np

import src.valladopy.astro.celestial.utils as utils

from ...conftest import DEFAULT_TOL


@pytest.fixture
def t():
    # Julian centuries from J2000
    return -0.013641341546885694


@pytest.mark.parametrize(
    "r2, earth_model, los, tmin",
    [
        (
            [0.0, 5740.323, 3189.068],
            utils.EarthModel.ELLIPSOIDAL,
            False,
            0.5082248650848982,
        ),
        (
            [0.0, 5740.323, 3189.068],
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


def test_sun_ecliptic_parameters(t):
    mean_lon, mean_anomaly, ecliptic_lon = utils.sun_ecliptic_parameters(t)
    assert np.isclose(np.degrees(mean_lon), 149.36181814781156, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(mean_anomaly), 226.4523822485284, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(ecliptic_lon), 147.9940329397011, rtol=DEFAULT_TOL)


def test_obliquity_ecliptic(t):
    obliquity = utils.obliquity_ecliptic(t)
    assert np.isclose(np.degrees(obliquity), 23.439468394733744, rtol=DEFAULT_TOL)
