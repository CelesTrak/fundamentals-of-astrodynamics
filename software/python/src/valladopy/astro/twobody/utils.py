# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import numpy as np

from ...constants import RE, MU, ECCEARTHSQRD, SMALL, TWOPI


def is_equatorial(inc):
    """Equatorial check for inclinations.

    Args:
        inc (float): Inclination in radians

    Returns:
        (bool): True if the inclination is equatorial
    """
    return inc < SMALL or abs(inc - np.pi) < SMALL


def site(latgd, lon, alt):
    """Finds the position and velocity vectors for a site.

    The answer is returned in the geocentric equatorial (ECEF) coordinate
    system. Note that the velocity is zero because the coordinate system is
    fixed to the Earth.

    References:
        Vallado: 2001, p. 404-407, Algorithm 47

    Args:
        latgd (float): Geodetic latitude in radians
        lon (float): Longitude of the site in radians
        alt (float): Altitude in km

    Returns:
        rsecef (np.array): ECEF site position vector in km
        vsecef (np.array): ECEF site velocity vector in km/s
    """
    # Compute site position vector
    sinlat = np.sin(latgd)
    cearth = RE / np.sqrt(1.0 - ECCEARTHSQRD * sinlat**2)
    rdel = (cearth + alt) * np.cos(latgd)
    rk = ((1.0 - ECCEARTHSQRD) * cearth + alt) * sinlat

    rsecef = np.array([rdel * np.cos(lon), rdel * np.sin(lon), rk])

    # Site velocity vector in ECEF frame is zero
    vsecef = np.array([0.0, 0.0, 0.0])

    return rsecef, vsecef


def findc2c3(znew):
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
        c2new = (1.0 - np.cos(sqrtz)) / znew
        c3new = (sqrtz - np.sin(sqrtz)) / (sqrtz**3)
    elif znew < -SMALL:
        sqrtz = np.sqrt(-znew)
        c2new = (1.0 - np.cosh(sqrtz)) / znew
        c3new = (np.sinh(sqrtz) - sqrtz) / (sqrtz**3)
    else:
        c2new = 0.5
        c3new = 1.0 / 6.0

    return c2new, c3new


def lon2nu(jdut1, lon, incl, raan, argp):
    """Converts the longitude of the ascending node to the true anomaly.

    This function calculates the true anomaly (`nu`) of an object
    at a given Julian date (`jdut1`) using the Greenwich Mean Sidereal Time
    (GMST) and orbital elements.

    Args:
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lon (float): Longitude of the ascending node in radians
        incl (float): Orbital inclination in radians
        raan (float): Right ascension of the ascending node in radians
        argp (float): Argument of periapsis in radians

    Returns:
        float: True anomaly (`nu`) in radians (0 to 2pi)
    """
    # Calculate GMST
    ed = jdut1 - 2451544.5  # Elapsed days from 1 Jan 2000, 0 hrs
    gmst = 99.96779469 + 360.9856473662860 * ed + 0.29079e-12 * ed * ed  # deg
    gmst = np.remainder(np.radians(gmst), TWOPI)

    # Check quadrants
    if gmst < 0.0:
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


def gc2gd(latgc):
    """Converts geocentric latitude to geodetic latitude for positions on the
    surface of the Earth.

    References:
        Vallado: 2001, p. 146, Eq. 3-11

    Args:
        latgc (float): Geocentric latitude in radians

    Returns:
        (float): Geodetic latitude in radians (-pi/2 to pi/2)
    """
    return np.arctan(np.tan(latgc) / (1.0 - ECCEARTHSQRD))


def gd2gc(latgd):
    """Converts geodetic latitude to geocentric latitude for positions on the
    surface of the Earth.

    References:
        Vallado: 2001, p. 146, Eq. 3-11

    Args:
        latgd (float): Geodetic latitude in radians

    Returns:
        (float): Geocentric latitude in radians (-pi/2 to pi/2)
    """
    return np.arctan((1.0 - ECCEARTHSQRD) * np.tan(latgd))


def checkhitearth(altpad, r1, v1, r2, v2, nrev):
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
        tuple: (hitearth, hitearthstr)
            hitearth (bool): True if Earth is impacted (False otherwise)
            hitearthstr (str): Explanation of the impact status
    """
    # Initialize variables
    hitearth, hitearthstr = False, "No impact"

    # Compute magnitudes of position vectors
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)

    # Define the padded radius (Earth's radius + altitude pad)
    rpad = RE + altpad

    # Check if the initial or final position vector is below the padded radius
    if magr1 < rpad or magr2 < rpad:
        hitearth, hitearthstr = True, "Impact at initial/final radii"
    else:
        rdotv1 = np.dot(r1, v1)
        rdotv2 = np.dot(r2, v2)

        # Solve for the reciprocal of the semi-major axis (1/a)
        ainv = 2.0 / magr1 - np.linalg.norm(v1) ** 2 / MU

        # Find ecos(E)
        ecosea1 = 1.0 - magr1 * ainv
        ecosea2 = 1.0 - magr2 * ainv

        # Determine the radius of perigee for nrev > 0
        if nrev > 0:
            a = 1.0 / ainv
            if a > 0.0:
                # Elliptical orbit
                esinea1 = rdotv1 / np.sqrt(MU * a)
                ecc = np.sqrt(ecosea1**2 + esinea1**2)
            else:
                # Hyperbolic orbit
                esinea1 = rdotv1 / np.sqrt(MU * abs(-a))
                ecc = np.sqrt(ecosea1**2 - esinea1**2)

            # Check if the radius of perigee is below the padded radius
            rp = a * (1.0 - ecc)
            if rp < rpad:
                hitearth, hitearthstr = True, "Impact during nrev"

        # Check for special cases when nrev = 0
        else:
            if (
                (rdotv1 < 0.0 < rdotv2)
                or (rdotv1 > 0.0 < rdotv2 and ecosea1 < ecosea2)
                or (rdotv1 < 0.0 > rdotv2 and ecosea1 > ecosea2)
            ):

                # Check for parabolic impact
                if abs(ainv) <= SMALL:
                    hbar = np.cross(r1, v1)
                    magh = np.linalg.norm(hbar)
                    rp = magh**2 * 0.5 / MU
                    if rp < rpad:
                        hitearth, hitearthstr = True, "Parabolic impact"

                else:
                    a = 1.0 / ainv
                    esinea1 = rdotv1 / np.sqrt(MU * abs(a))
                    if ainv > 0.0:
                        ecc = np.sqrt(ecosea1**2 + esinea1**2)
                    else:
                        ecc = np.sqrt(ecosea1**2 - esinea1**2)

                    # Check for elliptical impact
                    if ecc < 1.0:
                        rp = a * (1.0 - ecc)
                        if rp < rpad:
                            hitearth, hitearthstr = True, "Elliptical impact"

                    # Check for hyperbolic impact
                    elif rdotv1 < 0.0 < rdotv2:
                        rp = a * (1.0 - ecc)
                        if rp < rpad:
                            hitearth, hitearthstr = True, "Hyperbolic impact"

    return hitearth, hitearthstr
