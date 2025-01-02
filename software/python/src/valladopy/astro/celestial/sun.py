# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import logging
from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const
from ..time.sidereal import lstime
from .utils import EarthModel, in_sight, sun_ecliptic_parameters, obliquity_ecliptic


# Set up logging
logger = logging.getLogger(__name__)


class SunEventType(Enum):
    SUNRISESET = "s"
    CIVIL_TWILIGHT = "c"
    NAUTICAL_TWILIGHT = "n"
    ASTRONOMICAL_TWILIGHT = "a"


def position(jd: float) -> Tuple[np.ndarray, float, float]:
    """Calculates the geocentric equatorial position vector of the Sun.

    This is the low precision formula and is valid for years from 1950 to 2050. The
    accuaracy of apparent coordinates is about 0.01 degrees.  notice many of the
    calculations are performed in degrees, and are not changed until later. This is due
    to the fact that the almanac uses degrees exclusively in their formulations.

    Sergey K (2022) has noted that improved results are found assuming the oputput is in
    a precessing frame (TEME) and converting to ICRF.

    References:
        Vallado: 2022, p. 285, Algorithm 29

    Args:
        jd (float): Julian date (days from 4713 BC)

    Returns:
        tuple: (rsun, rtasc, decl)
            rsun (np.ndarray): Inertial sun position vector in km
            rtasc (float): Right ascension of the sun in radians
            decl (float): Declination of the sun in radians
    """
    # Julian centuries from J2000
    tut1 = (jd - const.J2000) / const.CENT2DAY

    # Mean anomaly and ecliptic longitude of the sun in radians
    _, meananomaly, eclplong = sun_ecliptic_parameters(tut1)

    # Obliquity of the ecliptic in radians
    obliquity = obliquity_ecliptic(tut1)

    # Magnitude of the Sun vector in AU
    magr = (
        1.000140612
        - 0.016708617 * np.cos(meananomaly)
        - 0.000139589 * np.cos(2 * meananomaly)
    )

    # Sun position vector in geocentric equatorial coordinates
    rsun = np.array(
        [
            magr * np.cos(eclplong),
            magr * np.cos(obliquity) * np.sin(eclplong),
            magr * np.sin(obliquity) * np.sin(eclplong),
        ]
    )

    # Right ascension in radians
    rtasc = np.arctan(np.cos(obliquity) * np.tan(eclplong))

    # Ensure right ascension is in the same quadrant as ecliptic longitude
    if eclplong < 0:
        eclplong += const.TWOPI
    if abs(eclplong - rtasc) > np.pi / 2:
        rtasc += 0.5 * np.pi * round((eclplong - rtasc) / (0.5 * np.pi))

    # Declination (radians)
    decl = np.arcsin(np.sin(obliquity) * np.sin(eclplong))

    return rsun * const.AU2KM, rtasc, decl


def rise_set(
    jd: float,
    latgd: float,
    lon: float,
    event_type: SunEventType = SunEventType.SUNRISESET,
) -> Tuple[float, float]:
    """Finds the universal time for sunrise and sunset given the day and site location.

    References:
        Vallado: 2022, p. 289-290, Algorithm 30

    Args:
        jd (float): Julian date (days from 4713 BC)
        latgd (float): Geodetic latitude of the site in radians (-65 deg to 65 deg)
        lon (float): Longitude of the site in radians (-2pi to 2pi) (west is negative)
        event_type (SunEventType): Type of event to calculate
                                   (default is SunEventType.SUNRISESET)

    Returns:
        tuple: (sunrise, sunset)
            sunrise (float): Universal time of sunrise in hours
            sunset (float): Universal time of sunset in hours
    """
    # Normalize longitude to -π to π
    lon = (lon + np.pi) % const.TWOPI - np.pi

    # Select the sun angle based on the kind of event
    sunangle_map = {
        SunEventType.SUNRISESET: np.radians(90 + 50 / 60),
        SunEventType.CIVIL_TWILIGHT: np.radians(96),
        SunEventType.NAUTICAL_TWILIGHT: np.radians(102),
        SunEventType.ASTRONOMICAL_TWILIGHT: np.radians(108),
    }
    sunangle = sunangle_map.get(event_type, None)
    if sunangle is None:
        raise ValueError(f"Invalid event type: {event_type}")

    # Initialize results dictionary
    results = {"sunrise": np.nan, "sunset": np.nan}

    # Loop for sunrise and sunset
    initial_guess_times = {"sunrise": 6, "sunset": 18}
    for event, jd_offset in initial_guess_times.items():
        # Initialize Julian date for the day
        jdtemp = (
            jd
            + (np.degrees(-lon) / const.DEG2HR / const.DAY2HR)
            + jd_offset / const.DAY2HR
        )

        # Julian centuries from J2000.0
        tut1 = (jdtemp - const.J2000) / const.CENT2DAY

        # Ecliptic longitude of the Sun in radians
        *_, lonecliptic = sun_ecliptic_parameters(tut1)

        # Obliquity of the ecliptic in radians
        obliquity = obliquity_ecliptic(tut1)

        # Right ascension and declination in radians
        ra = np.arctan(np.cos(obliquity) * np.tan(lonecliptic))
        decl = np.arcsin(np.sin(obliquity) * np.sin(lonecliptic))

        # Local hour angle
        lha = (np.cos(sunangle) - np.sin(decl) * np.sin(latgd)) / (
            np.cos(decl) * np.cos(latgd)
        )
        if abs(lha) > 1:
            logger.error("Local hour angle out of range; sunrise/sunset not visible.")
            return results["sunrise"], results["sunset"]

        lha = np.arccos(lha)
        if event == "sunrise":
            lha = const.TWOPI - lha

        # GST and UT
        gst = (
            1.75336855923327
            + 628.331970688841 * tut1
            + 6.77071394490334e-06 * tut1**2
            - 4.50876723431868e-10 * tut1**3
        ) % const.TWOPI
        uttemp = (lha + ra - gst) % const.TWOPI
        uttemp = np.degrees(uttemp) / const.DEG2HR
        uttemp = uttemp % const.DAY2HR

        # Assign to sunrise or sunset
        results[event] = uttemp

    return results["sunrise"], results["sunset"]


def in_light(
    r: ArrayLike, jd: float, earth_model: EarthModel = EarthModel.ELLIPSOIDAL
) -> bool:
    """Determines if a spacecraft is in sunlight at a given time.

    References:
        Vallado: 2022, p. 312-315, Algorithm 35

    Args:
        r (array_like): Position vector of the spacecraft in km
        jd (float): Julian date (days from 4713 BC)
        earth_model (EarthModel, optional): Earth model to use (default is ELLIPSOIDAL)

    Returns:
        bool: True if the spacecraft is in sunlight, False otherwise
    """
    # Calculate the Sun's position vector
    rsun, *_ = position(jd)

    # Determine if the spacecraft is in sunlight
    return in_sight(rsun, r, earth_model)


def illumination(jd: float, lat: float, lon: float) -> float:
    """Calculates the illumination due to the sun at a given location and time.

    References:
        Vallado: 2022, p. 316-317, Eq. 5-10, Table 5-1

    Args:
        jd (float): Julian date (days from 4713 BC)
        lat (float): Latitude of the location in radians
        lon (float): Longitude of the location in radians (-2pi to 2pi) (West negative)

    Returns:
        float: Luminous emmittance, lux (lumen/m²)
    """
    # Sun right ascension and declination
    _, srtasc, sdecl = position(jd)

    # Local sidereal time
    lst, _ = lstime(lon, jd)

    # Local hour angle
    lha = lst - srtasc

    # Sun elevation
    sunel = np.arcsin(
        np.sin(sdecl) * np.sin(lat) + np.cos(sdecl) * np.cos(lat) * np.cos(lha)
    )

    # Convert sun elevation to degrees
    sunel_deg = np.degrees(sunel)

    # Compute illumination using ground illumination indices
    sunillum = 0
    if sunel_deg > -18.01:
        x = sunel_deg / 90

        # Determine coefficients based on sun elevation
        # See Table 5-1 in Vallado 2022 (p. 317)
        if sunel_deg >= 20:
            l0, l1, l2, l3 = 3.74, 3.97, -4.07, 1.47
        elif 5 <= sunel_deg < 20:
            l0, l1, l2, l3 = 3.05, 13.28, -45.98, 64.33
        elif -0.8 <= sunel_deg < 5:
            l0, l1, l2, l3 = 2.88, 22.26, -207.64, 1034.30
        elif -5 <= sunel_deg < -0.8:
            l0, l1, l2, l3 = 2.88, 21.81, -258.11, -858.36
        elif -12 <= sunel_deg < -5:
            l0, l1, l2, l3 = 2.70, 12.17, -431.69, -1899.83
        elif -18 <= sunel_deg < -12:
            l0, l1, l2, l3 = 13.84, 262.72, 1447.42, 2797.93
        else:
            l0, l1, l2, l3 = 0, 0, 0, 0

        # Compute illumination
        l1 = l0 + l1 * x + l2 * x**2 + l3 * x**3
        sunillum = 10**l1

        # Clamp sunillum to valid range
        if sunillum < 0 or sunillum >= 1e4:
            sunillum = 0

    return sunillum
