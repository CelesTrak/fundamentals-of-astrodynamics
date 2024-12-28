# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ...constants import RE, MU, ECCEARTHSQRD, SMALL, TWOPI, J2000_UTC


class OrbitType(Enum):
    CIR_EQUATORIAL = 1  # circular equatorial
    CIR_INCLINED = 2  # circular inclined
    EPH_EQUATORIAL = 3  # elliptical, parabolic, hyperbolic equatorial
    EPH_INCLINED = 4  # elliptical, parabolic, hyperbolic inclined


def determine_orbit_type(ecc: float, incl: float, tol: float = SMALL) -> OrbitType:
    """Determine the type of orbit based on eccentricity and inclination

    Args:
        ecc (float): The eccentricity of the orbit
        incl (float): The inclination of the orbit in radians
        tol (float, optional): Small value for tolerance

    Returns:
        OrbitType: The type of orbit categorized into one of the following:
                   - circular equatorial
                   - circular inclined
                   - elliptical, parabolic, hyperbolic equatorial
                   - elliptical, parabolic, hyperbolic inclined
    """
    if ecc < tol:
        if (incl < tol) or (abs(incl - np.pi) < tol):
            return OrbitType.CIR_EQUATORIAL
        else:
            return OrbitType.CIR_INCLINED
    elif (incl < tol) or (abs(incl - np.pi) < tol):
        return OrbitType.EPH_EQUATORIAL
    else:
        return OrbitType.EPH_INCLINED


def is_equatorial(inc: float) -> bool:
    """Equatorial check for inclinations.

    Args:
        inc (float): Inclination in radians

    Returns:
        (bool): True if the inclination is equatorial
    """
    return inc < SMALL or abs(inc - np.pi) < SMALL


def site(latgd: float, lon: float, alt: float) -> Tuple[np.ndarray, np.ndarray]:
    """Finds the position and velocity vectors for a site.

    The answer is returned in the geocentric equatorial (ECEF) coordinate system.
    Note that the velocity is zero because the coordinate system is fixed to the Earth.

    References:
        Vallado: 2001, p. 404-407, Algorithm 47

    Args:
        latgd (float): Geodetic latitude in radians
        lon (float): Longitude of the site in radians
        alt (float): Altitude in km

    Returns:
        tuple: (rsecef, vsecef)
            rsecef (np.ndarray): ECEF site position vector in km
            vsecef (np.ndarray): ECEF site velocity vector in km/s
    """
    # Compute site position vector
    sinlat = np.sin(latgd)
    cearth = RE / np.sqrt(1 - ECCEARTHSQRD * sinlat**2)
    rdel = (cearth + alt) * np.cos(latgd)
    rk = ((1 - ECCEARTHSQRD) * cearth + alt) * sinlat

    rsecef = np.array([rdel * np.cos(lon), rdel * np.sin(lon), rk])

    # Site velocity vector in ECEF frame is zero
    vsecef = np.zeros(3)

    return rsecef, vsecef


def findc2c3(znew: float) -> Tuple[float, float]:
    """Calculates the c2 and c3 functions for the universal variable z.

    References:
        Vallado: 2001, p. 70-71, Algorithm 1

    Args:
        znew (float): z variable in rad^2

    Returns:
        tuple: (c2new, c3new)
            c2new (float): c2 function value
            c3new (float): c3 function value
    """
    if znew > SMALL:
        sqrtz = np.sqrt(znew)
        c2new = (1 - np.cos(sqrtz)) / znew
        c3new = (sqrtz - np.sin(sqrtz)) / (sqrtz**3)
    elif znew < -SMALL:
        sqrtz = np.sqrt(-znew)
        c2new = (1 - np.cosh(sqrtz)) / znew
        c3new = (np.sinh(sqrtz) - sqrtz) / (sqrtz**3)
    else:
        c2new = 0.5
        c3new = 1 / 6

    return c2new, c3new


def lon2nu(jdut1: float, lon: float, incl: float, raan: float, argp: float) -> float:
    """Converts the longitude of the ascending node to the true anomaly.

    This function calculates the true anomaly (`nu`) of an object at a given Julian date
    using the Greenwich Mean Sidereal Time (GMST) and orbital elements.

    Args:
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lon (float): Longitude of the ascending node in radians
        incl (float): Orbital inclination in radians
        raan (float): Right ascension of the ascending node in radians
        argp (float): Argument of periapsis in radians

    Returns:
        float: True anomaly in radians (0 to 2pi)
    """
    # Calculate GMST
    ed = jdut1 - J2000_UTC
    gmst = 99.96779469 + 360.9856473662860 * ed + 0.29079e-12 * ed * ed  # deg
    gmst = np.remainder(np.radians(gmst), TWOPI)

    # Check quadrants
    if gmst < 0:
        gmst += TWOPI

    # Calculate lambdau
    lambdau = gmst + lon - raan

    # Ensure lambdau is within 0 to 2pi radians
    lambdau = np.mod(lambdau, TWOPI)

    # Calculate argument of latitude
    arglat = np.arctan(np.tan(lambdau) / np.cos(incl))

    # Adjust arglat for quadrants
    if (lambdau >= 0.5 * np.pi) and (lambdau < 1.5 * np.pi):
        arglat += np.pi

    return np.mod(arglat - argp, TWOPI)


def gc2gd(latgc: float) -> float:
    """Converts geocentric latitude to geodetic latitude for positions on the surface of
    the Earth.

    References:
        Vallado: 2001, p. 146, Eq. 3-11

    Args:
        latgc (float): Geocentric latitude in radians

    Returns:
        float: Geodetic latitude in radians (-pi/2 to pi/2)
    """
    return np.arctan(np.tan(latgc) / (1 - ECCEARTHSQRD))


def gd2gc(latgd: float) -> float:
    """Converts geodetic latitude to geocentric latitude for positions on the surface of
    the Earth.

    References:
        Vallado: 2001, p. 146, Eq. 3-11

    Args:
        latgd (float): Geodetic latitude in radians

    Returns:
        float: Geocentric latitude in radians (-pi/2 to pi/2)
    """
    return np.arctan((1 - ECCEARTHSQRD) * np.tan(latgd))


def checkhitearth(
    altpad: float, r1: ArrayLike, v1: ArrayLike, r2: ArrayLike, v2: ArrayLike, nrev: int
) -> Tuple[bool, str]:
    """Checks if the trajectory impacts Earth during the transfer.

    References:
        Vallado: 2013, p. 503, Algorithm 60

    Args:
        altpad (float): Altitude pad above the Earth's surface in km
        r1 (array_like): Initial position vector in km
        v1 (array_like): Initial velocity vector in km/s
        r2 (array_like): Final position vector in km
        v2 (array_like): Final velocity vector in km/s
        nrev (int): Number of revolutions (0, 1, 2, ...)

    Returns:
        tuple:
            bool: True if Earth is impacted (False otherwise)
            str: Explanation of the impact status
    """
    # Compute magnitudes of position vectors
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)

    # Define the padded radius (Earth's radius + altitude pad)
    rpad = RE + altpad

    # Check if the initial or final position vector is below the padded radius
    if magr1 < rpad or magr2 < rpad:
        return True, "Impact at initial/final radii"

    # Compute dot products of position and velocity vectors
    rdotv1, rdotv2 = np.dot(r1, v1), np.dot(r2, v2)

    # Solve for the reciprocal of the semi-major axis (1/a)
    ainv = 2 / magr1 - np.linalg.norm(v1) ** 2 / MU
    a = 1 / ainv

    # Find ecos(E)
    ecosea1, ecosea2 = 1 - magr1 * ainv, 1 - magr2 * ainv

    # Determine the radius of perigee for nrev > 0
    if nrev > 0:
        if a > 0:
            # Elliptical orbit
            esinea1 = rdotv1 / np.sqrt(MU * a)
            ecc = np.sqrt(ecosea1**2 + esinea1**2)
        else:
            # Hyperbolic orbit
            esinea1 = rdotv1 / np.sqrt(MU * abs(-a))
            ecc = np.sqrt(ecosea1**2 - esinea1**2)

        # Check if the radius of perigee is below the padded radius
        rp = a * (1 - ecc)
        if rp < rpad:
            return True, "Impact during nrev"

    # Check for special cases when nrev = 0
    else:
        if (
            (rdotv1 < 0 < rdotv2)
            or (rdotv1 > 0 < rdotv2 and ecosea1 < ecosea2)
            or (rdotv1 < 0 > rdotv2 and ecosea1 > ecosea2)
        ):

            # Check for parabolic impact
            if abs(ainv) <= SMALL:
                hbar = np.cross(r1, v1)
                magh = np.linalg.norm(hbar)
                rp = magh**2 * 0.5 / MU
                if rp < rpad:
                    return True, "Parabolic impact"

            else:
                esinea1 = rdotv1 / np.sqrt(MU * abs(a))
                if ainv > 0:
                    ecc = np.sqrt(ecosea1**2 + esinea1**2)
                else:
                    ecc = np.sqrt(ecosea1**2 - esinea1**2)

                # Check for elliptical impact
                rp = a * (1 - ecc)
                if ecc < 1 and rp < rpad:
                    return True, "Elliptical impact"

                # Check for hyperbolic impact
                elif rdotv1 < 0 < rdotv2 and rp < rpad:
                    return True, "Hyperbolic impact"

    return False, "No impact"


def findtof(ro: ArrayLike, r: ArrayLike, p: float) -> float:
    """Finds the time of flight for orbital transfer using p-iteration theory.

    References:
        Vallado: 2007, p. 134-135, Algorithm 11

    Args:
        ro (array_like): Interceptor position vector in km
        r (array_like): Target position vector in km
        p (float): Semiparameter in km

    Returns:
        float: Time of flight in seconds
    """
    # Magnitudes of position vectors
    magr, magro = np.linalg.norm(r), np.linalg.norm(ro)

    # Calculate cosine and sine of change in true anomaly
    cosdnu = np.dot(ro, r) / (magro * magr)
    rcrossr = np.cross(ro, r)
    sindnu = np.linalg.norm(rcrossr) / (magro * magr)

    # Intermediate calculations
    k = magro * magr * (1 - cosdnu)
    l_ = magro + magr
    m = magro * magr * (1 + cosdnu)
    a = (m * k * p) / ((2 * m - l_**2) * p**2 + 2 * k * l_ * p - k**2)

    # Compute f and g
    f = 1 - (magr / p) * (1 - cosdnu)
    g = magro * magr * sindnu / np.sqrt(MU * p)
    alpha = 1 / a

    # Find time of flight based on orbit type
    if alpha > SMALL:
        # Elliptical case
        dnu = np.arctan2(sindnu, cosdnu)
        fdot = (
            np.sqrt(MU / p)
            * np.tan(dnu * 0.5)
            * (((1 - cosdnu) / p) - (1 / magro) - (1 / magr))
        )
        cosdeltae = 1 - (magro / a) * (1 - f)
        sindeltae = (-magro * magr * fdot) / np.sqrt(MU * a)
        deltae = np.arctan2(sindeltae, cosdeltae)
        tof = g + np.sqrt(a**3 / MU) * (deltae - sindeltae)
    elif alpha < SMALL:
        # Hyperbolic case
        deltah = np.arccosh(1 - (magro / a) * (1 - f))
        tof = g + np.sqrt(-(a**3) / MU) * (np.sinh(deltah) - deltah)
    else:
        # Parabolic case
        c = np.sqrt(magr**2 + magro**2 - 2 * magr * magro * cosdnu)
        s = (magro + magr + c) * 0.5
        tof = (2 / 3) * np.sqrt((s**3) * 0.5 / MU) * (1 - ((s - c) / s) ** 1.5)

    return tof
