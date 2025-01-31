import logging
import numpy as np
import pytest

import src.valladopy.astro.celestial.sun as sun
from src.valladopy.astro.celestial.utils import EarthModel
import src.valladopy.constants as const

from ...conftest import DEFAULT_TOL, custom_isclose


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
    r = [0, -4464.696, -5102.509]
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


def test_in_shadow_simple():
    # Test against values from Example 12.8 in Curtis
    r_sat = [2817.899, -14110.473, -7502.672]
    r_sun = [-11747041, 139486985, 60472278]
    assert sun.in_shadow_simple(r_sat, r_sun)


def test_in_shadow():
    r_eci = [-41260.1818237031, 8684.15782134066, 0]
    r_sun = [148470363.19330865, -9449738.11151353, -4096753.810182002]
    in_umbra, in_penumbra = sun.in_shadow(r_eci, r_sun)
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
    roots = sun.cylindrical_shadow_roots(a, e, beta_1, beta_2)

    # Check results
    assert np.allclose(roots, roots_expected, rtol=DEFAULT_TOL)


def test_eclipse_entry_exit(coe_shadow):
    # From STK (Epoch = 18 Jan 2020 00:00:00.000 UTCG)
    r_sun = [67029379.328701, -120206236.987781, -52109013.709624]

    # Expected values in degrees
    # STK penumbra start and stop true anomalies: 50.061 and 196.692 deg
    # STK umbra start and stop true anomalies: 50.747 and 196.025 deg
    theta_en_exp = 50.043417929102404
    theta_ex_exp = -163.5273086027108  # 196.4726913972892 deg

    # Call function
    theta_en, theta_ex = sun.eclipse_entry_exit(r_sun, *coe_shadow)

    # Check results
    assert custom_isclose(float(np.degrees(theta_en)), theta_en_exp)
    assert custom_isclose(float(np.degrees(theta_ex)), theta_ex_exp)
