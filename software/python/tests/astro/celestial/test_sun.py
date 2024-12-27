import logging
import numpy as np
import pytest

import src.valladopy.astro.celestial.sun as sun
from src.valladopy.astro.celestial.utils import EarthModel

from ...conftest import DEFAULT_TOL, custom_isclose


def test_position():
    # Vallado 2007, Example 5-1
    jd = 2453827.5
    rsun, rtasc, decl = sun.position(jd)

    # Expected values
    rsun_expected = [146185872.47656828, 28790431.768908378, 12481694.638482867]
    assert np.allclose(rsun, rsun_expected, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(rtasc), 11.14147300257444, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(decl), 4.788667312446681, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "event_type, lon, sunriset_expected, sunset_expected",
    [
        # Example 5-2
        (sun.SunEventType.SUNRISESET, 0, 5.972768319603184, 18.254933720588195),
        # Other options
        (sun.SunEventType.CIVIL_TWILIGHT, 0, 5.521700171944827, 18.706162691802458),
        (sun.SunEventType.NAUTICAL_TWILIGHT, 0, 4.992562953381093, 19.235701157332667),
        (
            sun.SunEventType.ASTRONOMICAL_TWILIGHT,
            0,
            4.452701865963909,
            19.776248578494837,
        ),
        # Non-zero longitude
        (
            sun.SunEventType.SUNRISESET,
            np.radians(-74.3),
            5.967174864254953,
            18.258450003347587,
        ),
    ],
)
def test_rise_set(event_type, lon, sunriset_expected, sunset_expected):
    # Vallado 2007, Example 5-2
    jd = 2450165.5
    latgd = np.radians(40)
    sunrise, sunset = sun.rise_set(jd, latgd, lon, event_type)

    # Expected values
    assert np.isclose(sunrise, sunriset_expected, rtol=DEFAULT_TOL)
    assert np.isclose(sunset, sunset_expected, rtol=DEFAULT_TOL)


def test_invalid_event_type():
    with pytest.raises(ValueError) as excinfo:
        sun.rise_set(2450165.5, np.radians(40), 0, "galactic")
    assert "Invalid event type" in str(excinfo.value)


@pytest.mark.parametrize(
    "earth_model, in_light, tmin",
    [
        (EarthModel.ELLIPSOIDAL, True, 1.0000233422875653),
        (EarthModel.SPHERICAL, True, 1.000023290273442),
    ],
)
def test_in_light(earth_model, in_light, tmin, caplog):
    # Vallado 2022, Example 5-6
    r = [0.0, -4464.696, -5102.509]
    jd = 2449763.5

    # Call function with logging
    with caplog.at_level(logging.DEBUG):
        assert sun.in_light(r, jd, earth_model) == in_light
        assert f"Minimum parametric value (tmin): {tmin}" in caplog.messages[0]


def test_illumination():
    # TODO: find examples to hit all sun elevation ranges
    jd = 2449763.5
    lat = np.radians(45)
    lon = np.radians(-75)
    assert custom_isclose(sun.illumination(jd, lat, lon), 0.0009451889253123211)
