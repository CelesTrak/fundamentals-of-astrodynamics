# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import logging
from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const
from ...mathtime.vector import angle
from ..time.sidereal import lstime
from ..twobody.frame_conversions import perifocal_transform
from . import utils as utils


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
    _, meananomaly, eclplong = utils.sun_ecliptic_parameters(tut1)

    # Obliquity of the ecliptic in radians
    obliquity = utils.obliquity_ecliptic(tut1)

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
        *_, lonecliptic = utils.sun_ecliptic_parameters(tut1)

        # Obliquity of the ecliptic in radians
        obliquity = utils.obliquity_ecliptic(tut1)

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
    r: ArrayLike,
    jd: float,
    earth_model: utils.EarthModel = utils.EarthModel.ELLIPSOIDAL,
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
    return utils.in_sight(rsun, r, earth_model)


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


def in_shadow_simple(r_sat: ArrayLike, r_sun: ArrayLike) -> bool:
    """Check if the satellite is in Earth's shadow.

    References:
        Curtis, H.D.: Orbit Mechanics for Engineering Students, 2014, Algorithm 12.3

    Args:
        r_sat (array_like): Satellite position vector in km
        r_sun (array_like): Sun position vector in km

    Returns:
        bool: Whether satellite is in attracting body's shadow
    """
    # Calculate angles
    sun_sat_angle = angle(r_sun, r_sat)
    angle1 = np.arccos(const.RE / np.linalg.norm(r_sat))
    angle2 = np.arccos(const.RE / np.linalg.norm(r_sun))

    # Check line of sight (no LOS = eclipse)
    if (angle1 + angle2) <= sun_sat_angle:
        return True

    return False


def in_shadow(r_eci: ArrayLike, r_sun: ArrayLike):
    """Check if in Earth's shadow (umbra and penumbra).

    References:
        Vallado: 2022, p. 305-308, Algorithm 34

    Args:
        r_eci (array_like): ECI position vector in km or AU
        r_sun (array_like): Sun position vector in km

    Returns:
        dict: Dictionary containing the computed angles, horizon, vertical components,
              penumbra and umbra status, and distance parameters.
    """
    # Umbra/penumbra angles
    angumb = np.arctan((const.SUNRADIUS - const.RE) / const.AU2KM)
    angpen = np.arctan((const.SUNRADIUS + const.RE) / const.AU2KM)

    # Check if in umbra/penumbra
    in_umbra, in_penumbra = False, False

    if np.dot(r_eci, r_sun) < 0:
        # Get satellite's vertical and horizontal distances
        sun_sat_angle = angle(-np.array(r_sun), np.array(r_eci))
        sathoriz = np.linalg.norm(r_eci) * np.cos(sun_sat_angle)
        satvert = np.linalg.norm(r_eci) * np.sin(sun_sat_angle)

        # Calculte penumbra vertical distance
        x = const.RE / np.sin(angpen)
        penvert = np.tan(angpen) * (x + sathoriz)

        # Check if in penumbra
        if satvert <= penvert:
            in_penumbra = True
            y = const.RE / np.sin(angumb)

            # Calculate umbra vertical distance
            umbvert = np.tan(angumb) * (y - sathoriz)

            # Check if in umbra
            if satvert <= umbvert:
                in_umbra = True

    return in_umbra, in_penumbra


def cylindrical_shadow_roots(
    a: float, e: float, beta_1: float, beta_2: float
) -> np.ndarray:
    """Calculate roots of cylindrical shadow quartic equation.

    References:
        Vallado: 2022, p. 310, Equation 5-6

    Args:
        a (float): Semimajor axis in km
        e (float): Eccentricity
        beta_1 (float): First temporary parameter
        beta_2 (float): Second temporary parameter

    Returns:
        np.ndarray: Roots of cylindrical shadow model
    """
    alpha = const.RE / (a * (1 - e**2))

    # Shadow coefficients
    a0 = (
        alpha**4 * e**4
        - 2 * alpha**2 * (beta_2**2 - beta_1**2) * e**2
        + (beta_1**2 + beta_2**2) ** 2
    )
    a1 = 4 * alpha**4 * e**3 - 4 * alpha**2 * (beta_2**2 - beta_1**2) * e
    a2 = (
        6 * alpha**4 * e**2
        - 2 * alpha**2 * (beta_2**2 - beta_1**2)
        - 2 * alpha**2 * (1 - beta_2**2) * e**2
        + 2 * (beta_2**2 - beta_1**2) * (1 - beta_2**2)
        - 4 * beta_1**2 * beta_2**2
    )
    a3 = 4 * alpha**4 * e - 4 * alpha**2 * (1 - beta_2**2) * e
    a4 = alpha**4 - 2 * alpha**2 * (1 - beta_2**2) + (1 - beta_2**2) ** 2

    return np.roots([a0, a1, a2, a3, a4])


def eclipse_entry_exit(
    r_sun: ArrayLike,
    a: float,
    e: float,
    i: float,
    raan: float,
    w: float,
    adjust: bool = False,
    tol: float = 1e-6,
) -> Tuple[float, float]:
    """Calculate eclipse entry and exit angles using the cylindrical shadow model.

    References:
        Vallado: 2022, p. 309-311
        Falck, R., and Dankanich, J.: Optimization of Low-Thrust Spiral Trajectories
            by Collocation, AIAA 2012-4423

    Args:
        r_sun (array_like): Sun position vector in km
        a (float): Semimajor axis in km
        e (float): Eccentricity
        i (float): Inclination in radians
        raan (float): Right ascension of ascending node in radians
        w (float): Argument of periapsis in radians
        adjust (bool, optional): Whether to adjust angles so that entry > exit
                                 (default is False)
        tol (float, optional): Tolerance for root finding (default is 1e-6)

    Returns:
        tuple: (theta_en, theta_ex)
            theta_en (float): True anomaly of shadow entry in radians
            theta_ex (float): True anomaly of shadow exit in radians
    """
    # Initialize entry/exit angles
    theta_en, theta_ex = np.pi, -np.pi

    # PQW transformation vectors (just P & Q used here)
    p_, q_, _ = perifocal_transform(i, raan, w)

    # Semiparameter and alpha
    p = a * (1 - e**2)
    alpha = p**2 / const.RE**2

    # Temp parameters
    beta_1 = np.dot(r_sun, p_) / np.linalg.norm(r_sun)
    beta_2 = np.dot(r_sun, q_) / np.linalg.norm(r_sun)

    # Solve quartic equation roots
    roots = cylindrical_shadow_roots(a, e, beta_1, beta_2)

    # Check for complex roots and return zeros if encountered
    if any(isinstance(r, complex) for r in roots):
        logger.debug(
            f"Roots are complex for OE: \n"
            f"a = {a} km, e = {e}, i = {np.degrees(i)} deg, "
            f"raan = {np.degrees(raan)} deg, w = {np.degrees(w)} deg"
        )
        return theta_en, theta_ex

    # Check for validity of all roots to remove any false positives
    # There are 8 potential solutions from the roots (Falck and Dankanich)
    thetas_good = []
    for theta in np.arccos(roots):
        # Iterate over possible sine values to get other 4 roots
        for y in [np.sin(theta), -np.sin(theta)]:
            # Pass the roots through the shadow function
            # A valid root will return ~0 (there should be 4 of these total)
            sol = np.arctan2(y, np.cos(theta))
            s = (
                (1 + e * np.cos(sol)) ** 2
                + alpha * (beta_1 * np.cos(sol) + beta_2 * np.sin(sol)) ** 2
                - alpha
            )

            # Then check that the sun angular separation is >= 90 deg
            # i.e. the cosine of the separation angle is <= ~0
            check = beta_1 * np.cos(sol) + beta_2 * np.sin(sol)

            # Check both conditions
            if abs(s) < tol and check <= tol:
                thetas_good.append(sol)

    # Check if true anomaly is shadow entry or exit
    for theta in thetas_good:
        # Calculate shadow derivative
        ds_dtheta = 2 * (p**2 / const.RE**2) * (
            beta_1 * np.cos(theta) + beta_2 * np.sin(theta)
        ) * (-beta_1 * np.sin(theta) + beta_2 * np.cos(theta)) - 2 * e * np.sin(
            theta
        ) * (
            1 + e * np.cos(theta)
        )
        if ds_dtheta > 0:
            theta_en = theta
        elif ds_dtheta < 0:
            theta_ex = theta
        else:  # pragma: no cover
            logger.warning(
                f"Unable to determine shadow entry/exit for true anomaly of "
                f"{np.degrees(theta): .2f} deg"
            )

    # Adjust angles if exit angle > entry angle (for integration)
    if adjust and theta_ex > theta_en:
        theta_en += const.TWOPI

    return theta_en, theta_ex
