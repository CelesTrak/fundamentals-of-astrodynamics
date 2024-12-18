# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import logging

import numpy as np

from ... import constants as const
from ..time.sidereal import lstime
from .utils import sun_ecliptic_parameters


# Set up logging
logger = logging.getLogger(__name__)


def position(jd: float) -> tuple[np.ndarray, float, float]:
    """Calculates the geocentric equatorial position vector of the Moon.

    References:
        Vallado: 2007, p. 290, Algorithm 31

    Args:
        jd (float): Julian date (days from 4713 BC)

    Returns:
        tuple: (rmoon, rtasc, decl)
            rmoon (np.ndarray): Inertial moon position vector in km
            rtasc (float): Right ascension of the moon in radians
            decl (float): Declination of the moon in radians
    """
    # Julian centuries from J2000
    ttdb = (jd - const.J2000) / const.CENT2DAY

    # Ecliptic longitude (radians)
    eclplong = (
        np.radians(
            218.32
            + 481267.8813 * ttdb
            + 6.29 * np.sin(np.radians(134.9 + 477198.85 * ttdb))
            - 1.27 * np.sin(np.radians(259.2 - 413335.38 * ttdb))
            + 0.66 * np.sin(np.radians(235.7 + 890534.23 * ttdb))
            + 0.21 * np.sin(np.radians(269.9 + 954397.70 * ttdb))
            - 0.19 * np.sin(np.radians(357.5 + 35999.05 * ttdb))
            - 0.11 * np.sin(np.radians(186.6 + 966404.05 * ttdb))
        )
        % const.TWOPI
    )

    # Ecliptic latitude (radians)
    eclplat = (
        np.radians(
            5.13 * np.sin(np.radians(93.3 + 483202.03 * ttdb))
            + 0.28 * np.sin(np.radians(228.2 + 960400.87 * ttdb))
            - 0.28 * np.sin(np.radians(318.3 + 6003.18 * ttdb))
            - 0.17 * np.sin(np.radians(217.6 - 407332.20 * ttdb))
        )
        % const.TWOPI
    )

    # Horizontal parallax (radians)
    hzparal = (
        np.radians(
            0.9508
            + 0.0518 * np.cos(np.radians(134.9 + 477198.85 * ttdb))
            + 0.0095 * np.cos(np.radians(259.2 - 413335.38 * ttdb))
            + 0.0078 * np.cos(np.radians(235.7 + 890534.23 * ttdb))
            + 0.0028 * np.cos(np.radians(269.9 + 954397.70 * ttdb))
        )
        % const.TWOPI
    )

    # Obliquity of the ecliptic (radians)
    obliquity = np.radians(np.degrees(const.OBLIQUITYEARTH) - 0.0130042 * ttdb)

    # Intermediate values for debugging
    logger.debug(f"Ecliptic longitude: {np.degrees(eclplong): .2f} deg")
    logger.debug(f"Ecliptic latitude: {np.degrees(eclplat): .2f} deg")
    logger.debug(f"Horizontal parallax: {np.degrees(hzparal): .2f} deg")
    logger.debug(f"Obliquity: {np.degrees(obliquity): .2f} deg")

    # Geocentric direction cosines
    l = np.cos(eclplat) * np.cos(eclplong)  # noqa: E741
    m = np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong) - np.sin(
        obliquity
    ) * np.sin(eclplat)
    n = np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong) + np.cos(
        obliquity
    ) * np.sin(eclplat)

    # Moon's position vector
    magr = 1.0 / np.sin(hzparal)
    rmoon = np.array([magr * l, magr * m, magr * n])

    # Right ascension and declination
    rtasc = np.arctan2(m, l)
    decl = np.arcsin(n)

    return rmoon * const.RE, rtasc, decl


def moonriset(jd: float, latgd: float, lon: float, n_iters: int = 5, tol: float = 0.008) -> tuple[float, float, float]:
    """Finds the universal time for moonrise and moonset given the day and site
    location.

    References:
        Vallado: 2007, p. 292, Algorithm 32

    Args:
        jd (float): Julian date (days from 4713 BC)
        latgd (float): Geodetic latitude of the site in radians (-65 deg to 65 deg)
        lon (float): Longitude of the site in radians (-2pi to 2pi) (west is negative)
        n_iters (int, optional): Number of iterations to attempt to find the moon
                                 rise/set times (defaults to 5)
        tol (float, optional): Tolerance for the iteration (defaults to 0.008)

    Returns:
        tuple: (moonrise, moonset, moonphaseang)
            moonrise (float): Universal time of moonrise in hours
            moonset (float): Universal time of moonset in hours
            moonphaseang (float): Moon phase angle in degrees
    """
    # Normalize longitude to -π to π
    lon = (lon + np.pi) % const.TWOPI - np.pi

    # Initialize results and variables
    results = {"moonrise": np.nan, "moonset": np.nan}
    moongha, deltaut, ttdb, lha = 0.0, 0.0, 0.0, 0.0

    # Iteration parameters
    try1 = 1

    for event, jd_offset in [("moonrise", 6.0), ("moonset", 18.0)]:
        # Initial guess for UT
        sign = -1 if event == "moonrise" else 1
        uttemp = (jd_offset + sign * lon / const.DEG2HR) / const.DAY2HR
        # uttemp = (jd_offset - np.degrees(lon) / const.DEG2HR) / const.DAY2HR  # fix

        # Set if there's a problem
        if try1 == 2:
            uttemp = 0.5

        # Initialize variables for iteration
        i = 0
        tn = uttemp
        t = tn + 10.0
        jdtemp = jd + uttemp

        while abs(tn - t) >= tol and i <= n_iters:
            # Update the Julian date
            ttdb = (jdtemp - const.J2000) / const.CENT2DAY

            # Ecliptic longitude and latitude (radians)
            eclplong = (
                np.radians(
                    218.32
                    + 481267.8813 * ttdb
                    + 6.29 * np.sin(np.radians(134.9 + 477198.85 * ttdb))
                    - 1.27 * np.sin(np.radians(259.2 - 413335.38 * ttdb))
                    + 0.66 * np.sin(np.radians(235.7 + 890534.23 * ttdb))
                    + 0.21 * np.sin(np.radians(269.9 + 954397.70 * ttdb))
                    - 0.19 * np.sin(np.radians(357.5 + 35999.05 * ttdb))
                    - 0.11 * np.sin(np.radians(186.6 + 966404.05 * ttdb))
                )
                % const.TWOPI
            )
            eclplat = (
                np.radians(
                    5.13 * np.sin(np.radians(93.3 + 483202.03 * ttdb))
                    + 0.28 * np.sin(np.radians(228.2 + 960400.87 * ttdb))
                    - 0.28 * np.sin(np.radians(318.3 + 6003.18 * ttdb))
                    - 0.17 * np.sin(np.radians(217.6 - 407332.20 * ttdb))
                )
                % const.TWOPI
            )

            # Obliquity of the ecliptic (radians)
            obliquity = np.radians(np.degrees(const.OBLIQUITYEARTH) - 0.0130042 * ttdb)

            # Geocentric direction cosines
            l = np.cos(eclplat) * np.cos(eclplong)  # noqa: E741
            m = np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong) - np.sin(
                obliquity
            ) * np.sin(eclplat)
            n = np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong) + np.cos(
                obliquity
            ) * np.sin(eclplat)

            # Right ascension and declination
            rtasc = np.arctan2(m, l)
            decl = np.arcsin(n)

            # Correction for right ascension
            if abs(eclplong - rtasc) > np.pi * 0.5:
                rtasc += 0.5 * np.pi * round(0.5 + (eclplong - rtasc) / (0.5 * np.pi))

            # Local sidereal time
            lst, _ = lstime(lon, jdtemp)

            # Calculat hour angles
            moonghan = lst - lon - rtasc
            if i == 0:
                lha = moonghan + lon
                dgha = np.radians(347.81)
            else:
                dgha = (moonghan - moongha) / deltaut
            dgha = dgha + const.TWOPI / abs(deltaut) if dgha < 0 else dgha

            # Local hour angle at moonrise and moonset
            lhan = (0.00233 - np.sin(latgd) * np.sin(decl)) / (
                np.cos(latgd) * np.cos(decl)
            )
            lhan = np.clip(lhan, -1.0, 1.0)
            lhan = np.arccos(lhan)
            if event == "moonrise":
                lhan = const.TWOPI - lhan

            # Time adjustment
            if abs(dgha) > 1e-4:
                deltaut = (lhan - lha) / dgha
            else:
                deltaut = 1.0
                logger.warning("dgha is too small; setting deltaut to 1.0")

            # Adjust deltaut for convergence and handle wrap-around cases
            t = tn
            if abs(deltaut) > 0.5:
                if abs(dgha) > 1e-3:
                    deltaut += (
                        const.TWOPI / dgha if deltaut < 0 else -const.TWOPI / dgha
                    )
                    if abs(deltaut) > 0.51:
                        break

            # Update variables
            tn = uttemp + deltaut
            jdtemp = jdtemp - uttemp + tn
            moongha = moonghan
            i += 1

        # Convert UT to hours
        uttemp = tn * const.DAY2HR if i <= n_iters else 1e4
        uttemp = uttemp % const.DAY2HR if 0.0 <= uttemp < 1e4 else uttemp

        # Assign to results
        results[event] = uttemp

        # Update the iteration and check for solution
        try1 += 1
        if i > n_iters and try1 < 3:  # retry if the first attempt fails
            logger.debug(f"Retrying option {event} (attempt {try1})")
        else:
            if i > n_iters and try1 > 2:  # if all retries fail, set the error
                if event == "moonrise":
                    logger.error("Moonrise not found")
                    results[event] = np.inf
                elif event == "moonset":
                    logger.error("Moonset not found")
                    results[event] = np.inf
            try1 = 1

    # Mean longitude of the Moon
    meanlonmoon = (
        218.32
        + 481267.8813 * ttdb
        + 6.29 * np.sin(np.radians(134.9 + 477198.85 * ttdb))
        - 1.27 * np.sin(np.radians(259.2 - 413335.38 * ttdb))
        + 0.66 * np.sin(np.radians(235.7 + 890534.23 * ttdb))
        + 0.21 * np.sin(np.radians(269.9 + 954397.70 * ttdb))
        - 0.19 * np.sin(np.radians(357.5 + 35999.05 * ttdb))
        - 0.11 * np.sin(np.radians(186.6 + 966404.05 * ttdb))
    ) % np.degrees(const.TWOPI)

    # Mean longitude of the Sun in radians (TODO: put in utils)
    *_, loneclsun = sun_ecliptic_parameters(ttdb)

    # Moon phase angle
    moonphaseang = (meanlonmoon - np.degrees(loneclsun)) % np.degrees(const.TWOPI)

    return results["moonrise"], results["moonset"], moonphaseang


def illumination(f: float, elev: float) -> float:
    """Calculates the illumination due to the moon.

    References:
        Vallado: 2022, p. 316-317, Eq. 5-10, Table 5-1

    Args:
        f (float): Phase angle in radians
        elev (float): Moon elevation in radians

    Returns:
        float: Luminous emmittance, lux (lumen/m²)
    """
    # Convert inputs to degrees
    f = np.degrees(f)
    elev = np.degrees(elev)

    # Determine coefficients based on moon elevation
    # See Table 5-1 in Vallado 2022 (p. 317)
    if elev >= 20:
        l0, l1, l2, l3 = -1.95, 4.06, -4.24, 1.56
    elif 5.0 <= elev < 20.0:
        l0, l1, l2, l3 = -2.58, 12.58, -42.58, 59.06
    elif -0.8 < elev < 5.0:
        l0, l1, l2, l3 = -2.79, 24.27, -252.95, 1321.29
    else:
        l0, l1, l2, l3 = 0.0, 0.0, 0.0, 0.0
        f = 0.0

    # Compute l1 and l2 based on coefficients and phase angle
    x = elev / 90.0
    l1 = l0 + l1 * x + l2 * x**2 + l3 * x**3
    l2 = -0.00868 * f - 2.2e-9 * f**4

    # Compute moon illumination
    moonillum = 10.0 ** (l1 + l2)

    # Clamp moonillum to valid range
    if moonillum < 0 or moonillum >= 1:
        moonillum = 0.0

    return moonillum
