import numpy as np
import pytest

import src.valladopy.astro.celestial.sun as sun

from ...conftest import DEFAULT_TOL


def test_position():
    # Vallado 2007, Example 5-1
    jd = 2453827.5
    rsun, rtasc, decl = sun.position(jd)

    # Expected values
    rsun_expected = [146186212.98684618, 28788976.311702874, 12481063.64508394]
    assert np.allclose(rsun, rsun_expected, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(rtasc), 11.140898551273013, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(decl), 4.788424663323541, rtol=DEFAULT_TOL)


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
def test_sunriset(event_type, lon, sunriset_expected, sunset_expected):
    # Vallado 2007, Example 5-2
    jd = 2450165.5
    latgd = np.radians(40)
    sunrise, sunset = sun.sunriset(jd, latgd, lon, event_type)

    # Expected values
    assert np.isclose(sunrise, sunriset_expected, rtol=DEFAULT_TOL)
    assert np.isclose(sunset, sunset_expected, rtol=DEFAULT_TOL)


def test_invalid_event_type():
    with pytest.raises(ValueError) as excinfo:
        sun.sunriset(2450165.5, np.radians(40), 0, "galactic")
    assert "Invalid event type" in str(excinfo.value)
