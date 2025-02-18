# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import numpy as np
from numpy.typing import ArrayLike
from typing import Tuple

from .utils import gibbs
from ... import constants as const


def calculate_time_intervals(
    jd1: float, jdf1: float, jd2: float, jdf2: float, jd3: float, jdf3: float
) -> Tuple[float, float, float]:
    """Calculate time intervals between three Julian dates.

    Args:
        jd1 (float): Julian date of first sighting (days from 4713 BC)
        jdf1 (float): Julian date fraction of first sighting (days from 4713 BC)
        jd2 (float): Julian date of second sighting (days from 4713 BC)
        jdf2 (float): Julian date fraction of second sighting (days from 4713 BC)
        jd3 (float): Julian date of third sighting (days from 4713 BC)
        jdf3 (float): Julian date fraction of third sighting (days from 4713 BC)

    Returns:
        tuple: (tau12, tau13, tau32)
            tau12 (float): Time interval between t1 and t2 in seconds
            tau13 (float): Time interval between t1 and t3 in seconds
            tau32 (float): Time interval between t3 and t2 in seconds
    """
    tau12 = ((jd1 - jd2) + (jdf1 - jdf2)) * const.DAY2SEC
    tau13 = ((jd1 - jd3) + (jdf1 - jdf3)) * const.DAY2SEC
    tau32 = ((jd3 - jd2) + (jdf3 - jdf2)) * const.DAY2SEC

    return tau12, tau13, tau32


def calculate_los_vectors(decl: list[float], rtasc: list[float]) -> list[np.ndarray]:
    """Calculate line-of-sight unit vectors from declination and right ascension.

    Args:
        decl (list[float]): List of declinations in radians
        rtasc (list[float]): List of right ascensions in radians

    Returns:
        list[np.ndarray]: List of line-of-sight unit vectors
    """
    return [
        np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
        for dec, ra in zip(decl, rtasc)
    ]


def _halley_iteration(
    poly: np.ndarray,
    bigr2c_default: float = 20000 / const.RE,
    tol: float = 8e-5,
    max_iter: int = 15,
):
    """Perform Halley's iteration to refine the value of bigr2."""
    # Select the max real root (use fallback of ~GPS altitude if necessary)
    roots = np.roots(poly)
    real_roots = roots[np.isreal(roots)].real
    bigr2c = max(real_roots) if real_roots.size > 0 else bigr2c_default

    # Halley iteration for refining bigr2
    bigr2, kk = 100, 0
    while abs(bigr2 - bigr2c) > tol and kk < max_iter:
        bigr2 = bigr2c
        deriv = bigr2**8 + poly[2] * bigr2**6 + poly[5] * bigr2**3 + poly[8]
        deriv1 = 8 * bigr2**7 + 6 * poly[2] * bigr2**5 + 3 * poly[5] * bigr2**2
        deriv2 = 56 * bigr2**6 + 30 * poly[2] * bigr2**4 + 6 * poly[5] * bigr2
        bigr2c = bigr2 - (2 * deriv * deriv1) / (2 * deriv1**2 - deriv * deriv2)
        kk += 1

    return bigr2c


def laplace(
    decl1: float,
    decl2: float,
    decl3: float,
    rtasc1: float,
    rtasc2: float,
    rtasc3: float,
    jd1: float,
    jdf1: float,
    jd2: float,
    jdf2: float,
    jd3: float,
    jdf3: float,
    diffsites: bool,
    rseci1: ArrayLike,
    rseci2: ArrayLike,
    rseci3: ArrayLike,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve orbit determination problem using three optical sightings and Laplace's
    method.

    References:
        Vallado: 2022, p. 441-445

    Args:
        decl1 (float): Declination of first sighting in radians
        decl2 (float): Declination of second sighting in radians
        decl3 (float): Declination of third sighting in radians
        rtasc1 (float): Right ascension of first sighting in radians
        rtasc2 (float): Right ascension of second sighting in radians
        rtasc3 (float): Right ascension of third sighting in radians
        jd1 (float): Julian date of first sighting (days from 4713 BC)
        jdf1 (float): Julian date fraction of first sighting (days from 4713 BC)
        jd2 (float): Julian date of second sighting (days from 4713 BC)
        jdf2 (float): Julian date fraction of second sighting (days from 4713 BC)
        jd3 (float): Julian date of third sighting (days from 4713 BC)
        jdf3 (float): Julian date fraction of third sighting (days from 4713 BC)
        diffsites (bool): True if different sites, False if same site
        rseci1 (array_like): ECI site position vector of first sighting in km
        rseci2 (array_like):  ECI site position vector of second sighting in km
        rseci3 (array_like):  ECI site position vector of third sighting in km

    Returns:
        tuple: (r2, v2)
            r2 (np.ndarray): ECI position vector in km
            v2 (np.ndarray): ECI velocity vector in km/s
    """
    # Convert Julian dates to time intervals
    tau12, tau13, tau32 = calculate_time_intervals(jd1, jdf1, jd2, jdf2, jd3, jdf3)

    # Line-of-sight unit vectors
    los1, los2, los3 = calculate_los_vectors(
        [decl1, decl2, decl3], [rtasc1, rtasc2, rtasc3]
    )

    # Check denominators for zero to avoid division errors
    if any(tau == 0 for tau in (tau12, tau13, tau32)):
        raise ValueError(
            "One or more time intervals (tau12, tau13, tau32) are zero, causing"
            "division by zero."
        )

    # Normalize to canonical units
    tau12c = tau12 / const.TUSEC
    tau13c = tau13 / const.TUSEC
    tau32c = tau32 / const.TUSEC
    rseci1c = np.array(rseci1) / const.RE
    rseci2c = np.array(rseci2) / const.RE
    rseci3c = np.array(rseci3) / const.RE

    # Lagrange coefficients
    s1 = -tau32c / (tau12c * tau13c)
    s2 = (tau12c + tau32c) / (tau12c * tau32c)
    s3 = -tau12c / (-tau13c * tau32c)
    s4 = 2 / (tau12c * tau13c)
    s5 = 2 / (tau12c * tau32c)
    s6 = 2 / (-tau13c * tau32c)

    # First and second derivatives of los vectors
    ldot = s1 * los1 + s2 * los2 + s3 * los3
    lddot = s4 * los1 + s5 * los2 + s6 * los3

    # Derivatives of site position vectors
    if not diffsites:
        # Same sites
        earth_rot_vec = [0, 0, const.EARTHROT * const.TUSEC]  # canonical Earth rotation
        rs2dot = np.cross(earth_rot_vec, rseci2c)
        rs2ddot = np.cross(earth_rot_vec, rs2dot)
    else:
        # Different sites
        rs2dot = s1 * rseci1c + s2 * rseci2c + s3 * rseci3c
        rs2ddot = s4 * rseci1c + s5 * rseci2c + s6 * rseci3c

    # Compute determinants
    dmat = np.column_stack((los2, ldot, lddot))
    dmat1 = np.column_stack((los2, ldot, rs2ddot))
    dmat2 = np.column_stack((los2, ldot, rseci2c))
    dmat3 = np.column_stack((los2, rs2ddot, lddot))
    dmat4 = np.column_stack((los2, rseci2c, lddot))

    d = 2 * np.linalg.det(dmat)
    d1c = np.linalg.det(dmat1)
    d2c = np.linalg.det(dmat2)
    d3c = np.linalg.det(dmat3)
    d4c = np.linalg.det(dmat4)

    # Check determinant value
    if abs(d) < const.SMALL:
        raise ValueError("Determinant is too small; system may be singular.")

    # Solve the 8th-order polynomial
    l2dotrs = np.dot(los2, rseci2c)
    poly = np.zeros(9)
    poly[0] = 1
    poly[2] = l2dotrs * 4 * d1c / d - 4 * d1c**2 / d**2 - np.linalg.norm(rseci2c) ** 2
    poly[5] = l2dotrs * 4 * d2c / d - 8 * d1c * d2c / d**2
    poly[8] = -4 * d2c**2 / d**2

    # Halley iteration for refining bigr2c
    bigr2c = _halley_iteration(poly)

    # Solve for rho and rho dot
    rho = -2 * d1c / d - 2 * d2c / (bigr2c**3 * d)
    rhodot = -d3c / d - d4c / (bigr2c**3 * d)

    # Position and velocity vectors at the middle
    r2c = rho * los2 + rseci2c
    v2c = rhodot * los2 + rho * ldot + rs2dot

    return r2c * const.RE, v2c * const.RE / const.TUSEC


def _update_fg_series(
    ll: int,
    r1: np.ndarray,
    r2: np.ndarray,
    r3: np.ndarray,
    v2: np.ndarray,
    tau12: float,
    tau32: float,
    bigr2: float,
    ll_iters: int = 8,
) -> Tuple[float, float, float, float]:
    """Helper function to compute or refine the f and g series.

    Args:
        ll (int): Current iteration count
        r1 (np.ndarray): ECI position vector at t1 in km
        r2 (np.ndarray): ECI position vector at t2 in km
        r3 (np.ndarray): ECI position vector at t3 in km
        v2 (np.ndarray): ECI velocity vector at t2 in km/s
        tau12 (float): Time interval between t1 and t2 in seconds
        tau32 (float): Time interval between t3 and t2 in seconds
        bigr2 (float): Magnitude of r2 in km
        ll_iters (int, optional): Maximum number of iterations for exact method
                                  (default 8)

    Returns:
        tuple: (f1, g1, f3, g3)
            f1 (float): Value of f1
            g1 (float): Value of g1
            f3 (float): Value of f3
            g3 (float): Value of g3
    """
    if ll <= ll_iters:
        # Approximate method for early iterations
        u = const.MU / np.linalg.norm(r2) ** 3
        rdot = np.dot(r2, v2) / np.linalg.norm(r2)
        udot = -3 * const.MU * rdot / np.linalg.norm(r2) ** 4

        f1 = 1 - 0.5 * u * tau12**2 - (1 / 6) * udot * tau12**3
        g1 = tau12 - (1 / 6) * u * tau12**3 - (1 / 12) * udot * tau12**4
        f3 = 1 - 0.5 * u * tau32**2 - (1 / 6) * udot * tau32**3
        g3 = tau32 - (1 / 6) * u * tau32**3 - (1 / 12) * udot * tau32**4
    else:
        # Exact method for later iterations
        theta12 = np.arccos(
            np.clip(np.dot(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2)), -1, 1)
        )
        theta23 = np.arccos(
            np.clip(np.dot(r2, r3) / (np.linalg.norm(r2) * np.linalg.norm(r3)), -1, 1)
        )

        f1 = 1 - np.linalg.norm(r1) * (1 - np.cos(theta12)) / bigr2
        g1 = np.linalg.norm(r1) * np.linalg.norm(r2) * np.sin(theta12) / np.sqrt(bigr2)
        f3 = 1 - np.linalg.norm(r3) * (1 - np.cos(theta23)) / bigr2
        g3 = np.linalg.norm(r3) * np.linalg.norm(r2) * np.sin(theta23) / np.sqrt(bigr2)

    return f1, g1, f3, g3


def gauss(
    decl1: float,
    decl2: float,
    decl3: float,
    rtasc1: float,
    rtasc2: float,
    rtasc3: float,
    jd1: float,
    jdf1: float,
    jd2: float,
    jdf2: float,
    jd3: float,
    jdf3: float,
    rseci1: ArrayLike,
    rseci2: ArrayLike,
    rseci3: ArrayLike,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve orbit determination problem using three optical sightings and the
    Gaussian method.

    References:
        Vallado: 2022, p. 448, Algorithm 52

    Args:
        decl1 (float): Declination of first sighting in radians
        decl2 (float): Declination of second sighting in radians
        decl3 (float): Declination of third sighting in radians
        rtasc1 (float): Right ascension of first sighting in radians
        rtasc2 (float): Right ascension of second sighting in radians
        rtasc3 (float): Right ascension of third sighting in radians
        jd1 (float): Julian date of first sighting (days from 4713 BC)
        jdf1 (float): Julian date fraction of first sighting (days from 4713 BC)
        jd2 (float): Julian date of second sighting (days from 4713 BC)
        jdf2 (float): Julian date fraction of second sighting (days from 4713 BC)
        jd3 (float): Julian date of third sighting (days from 4713 BC)
        jdf3 (float): Julian date fraction of third sighting (days from 4713 BC)
        rseci1 (array_like): ECI site position vector of first sighting in km
        rseci2 (array_like): ECI site position vector of second sighting in km
        rseci3 (array_like): ECI site position vector of third sighting in km

    Returns:
        tuple: (r2, v2)
            r2 (np.ndarray): ECI position vector in km
            v2 (np.ndarray): ECI velocity vector in km/s
    """
    # Time intervals in seconds
    tau12, tau13, tau32 = calculate_time_intervals(jd1, jdf1, jd2, jdf2, jd3, jdf3)

    # Line-of-sight unit vectors
    los1, los2, los3 = calculate_los_vectors(
        [decl1, decl2, decl3], [rtasc1, rtasc2, rtasc3]
    )

    # Normalize to canonical units
    tau12c = tau12 / const.TUSEC
    tau32c = tau32 / const.TUSEC
    rseci1c = np.array(rseci1) / const.RE
    rseci2c = np.array(rseci2) / const.RE
    rseci3c = np.array(rseci3) / const.RE

    # Construct l-matrix and determinant
    lmat = np.column_stack([los1, los2, los3])
    d = np.linalg.det(lmat)

    # Check determinant value
    if abs(d) < const.SMALL:
        raise ValueError("Determinant too small; unable to proceed with calculations.")

    # Inverse of l-matrix
    lmati = np.linalg.inv(lmat)

    # Range-site matrix
    rsmatc = np.column_stack([rseci1c, rseci2c, rseci3c])
    lir = lmati @ rsmatc

    # Calculate coefficients for polynomial
    a1 = tau32c / (tau32c - tau12c)
    a1u = (tau32c * ((tau32c - tau12c) ** 2 - tau32c**2)) / (6 * (tau32c - tau12c))
    a3 = -tau12c / (tau32c - tau12c)
    a3u = -(tau12c * ((tau32c - tau12c) ** 2 - tau12c**2)) / (6 * (tau32c - tau12c))

    # Polynomial terms
    d1c = lir[1, 0] * a1 - lir[1, 1] + lir[1, 2] * a3
    d2c = lir[1, 0] * a1u + lir[1, 2] * a3u
    magrs2 = np.linalg.norm(rseci2c)
    l2dotrs = np.dot(los2, rseci2c)

    # Solve the 8th-order polynomial
    poly = np.zeros(9)
    poly[0] = 1
    poly[2] = -(d1c**2 + 2 * d1c * l2dotrs + magrs2**2)
    poly[5] = -2 * (l2dotrs * d2c + d1c * d2c)
    poly[8] = -(d2c**2)

    # Halley iteration for refining bigr2c
    bigr2c = _halley_iteration(poly)

    # Ensure a valid radius solution
    if bigr2c < 0 or bigr2c * const.RE > 50000:
        bigr2c = 35000 / const.RE  # adjust to a GEO-like altitude

    # Convert some values back to non-canonical units
    bigr2 = bigr2c * const.RE
    a1u *= const.TUSEC**2
    a3u *= const.TUSEC**2

    # Solve for improved estimates of f and g series
    u = const.MU / bigr2**3
    c1 = a1 + a1u * u
    c2 = -1
    c3 = a3 + a3u * u

    # Compute range values
    rsmat = np.column_stack([rseci1, rseci2, rseci3])
    lir = lmati @ rsmat
    cmat = np.array([-c1, -c2, -c3]).reshape(-1, 1)
    rhomat = lir @ cmat

    # Form position vectors
    r1 = rhomat[0, 0] * los1 / c1 + rseci1
    r2 = rhomat[1, 0] * los2 / c2 + rseci2
    r3 = rhomat[2, 0] * los3 / c3 + rseci3

    # Find middle velocity vector with gibbs
    v2, *_ = gibbs(r1, r2, r3)

    return r2, v2


def doubler_iter(
    magr1in: float,
    magr2in: float,
    los1: ArrayLike,
    los2: ArrayLike,
    los3: ArrayLike,
    rsite1: ArrayLike,
    rsite2: ArrayLike,
    rsite3: ArrayLike,
    tau12: float,
    tau32: float,
    n12: int,
    n13: int,
    n23: int,
) -> Tuple[np.ndarray, np.ndarray, float, float, float, float, float, float, float]:
    """Perform the iterative work for the double-r angles-only routine.

    References:
        Vallado: 2022, p. 449-452

    Args:
        magr1in (float): Magnitude of the first sighting position vector
        magr2in (float): Magnitude of the second sighting position vector
        los1 (array_like): Line-of-sight unit vector for the first sighting
        los2 (array_like): Line-of-sight unit vector for the second sighting
        los3 (array_like): Line-of-sight unit vector for the third sighting
        rsite1 (array_like): ECI site position vector of the first sighting
        rsite2 (array_like): ECI site position vector of the second sighting
        rsite3 (array_like): ECI site position vector of the third sighting
        tau12 (float): Time interval between t1 and t2 in seconds
        tau32 (float): Time interval between t3 and t2 in seconds
        n12 (int): Number of days between the first and second sightings
        n13 (int): Number of days between the first and third sightings
        n23 (int): Number of days between the second and third sightings

        Returns:
            tuple: (r2, r3, f1, f2, q1, magr1, magr2, a, deltae32)
                r2 (np.ndarray): ECI position vector at t2 in km
                r3 (np.ndarray): ECI position vector at t3 in km
                f1 (float): Value of f1 coefficient
                f2 (float): Value of f2 coefficient
                q1 (float): Quality estimate of the solution
                magr1 (float): Magnitude of the first sighting position vector
                magr2 (float): Magnitude of the second sighting position vector
                a (float): Semi-major axis of the orbit in km
                deltae32 (float): Eccentric anomaly difference between obs 3 and 2 in
                                  radians

    """
    # Define default range value for when the square root is negative
    # Use this because hyperbolic likely at shorter times, lower alt
    default_range = 300

    # Range coefficients
    cc1 = 2 * np.dot(los1, rsite1)
    cc2 = 2 * np.dot(los2, rsite2)

    # Magnitude of the site position vectors
    magrsite1 = np.linalg.norm(rsite1)
    magrsite2 = np.linalg.norm(rsite2)

    # Compute rho1 and rho2
    tempsq1 = cc1**2 - 4 * (magrsite1**2 - magr1in**2)
    tempsq1 = default_range if tempsq1 < 0 else tempsq1
    rho1 = (-cc1 + np.sqrt(tempsq1)) * 0.5
    tempsq2 = cc2**2 - 4 * (magrsite2**2 - magr2in**2)
    tempsq2 = default_range if tempsq2 < 0 else tempsq2
    rho2 = (-cc2 + np.sqrt(tempsq2)) * 0.5

    # Compute r1 and r2
    r1 = rho1 * np.array(los1) + np.array(rsite1)
    r2 = rho2 * np.array(los2) + np.array(rsite2)
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)

    # Compute the cross product and determine rho3
    w = np.cross(r1, r2) / (magr1 * magr2)
    rho3 = -np.dot(rsite3, w) / np.dot(los3, w)
    r3 = rho3 * np.array(los3) + np.array(rsite3)
    magr3 = np.linalg.norm(r3)

    # Delta true anomaly between obs 2 and 1
    cosdv21 = np.dot(r2, r1) / (magr2 * magr1)
    sindv21 = np.linalg.norm(np.cross(r2, r1)) / (magr2 * magr1)
    dv21 = np.arctan2(sindv21, cosdv21) + const.TWOPI * n12

    # Delta true anomaly between obs 3 and 1
    cosdv31 = np.dot(r3, r1) / (magr3 * magr1)
    sindv31 = np.sqrt(1 - cosdv31**2)
    dv31 = np.arctan2(sindv31, cosdv31) + const.TWOPI * n13

    # Delta true anomaly between obs 3 and 2
    cosdv32 = np.dot(r3, r2) / (magr3 * magr2)
    sindv32 = np.linalg.norm(np.cross(r3, r2)) / (magr3 * magr2)

    # Compute the semi-parameter
    if dv31 > np.pi:
        p = (sindv21 - sindv31 + sindv32) / (
            -sindv31 / magr2 + sindv21 / magr3 + sindv32 / magr1
        )
    else:
        c1 = (magr1 * sindv31) / (magr2 * sindv32)
        c3 = (magr1 * sindv21) / (magr3 * sindv32)
        p = (c3 * magr3 - c1 * magr2 + magr1) / (-c1 + c3 + 1)

    # Compute the cosines of the eccentric anomalies
    ecosv1 = p / magr1 - 1
    ecosv2 = p / magr2 - 1
    ecosv3 = p / magr3 - 1

    # Compute the sines of the eccentric anomalies
    if dv21 != np.pi:
        esinv2 = (-cosdv21 * ecosv2 + ecosv1) / sindv21
    else:
        esinv2 = (cosdv32 * ecosv2 - ecosv3) / sindv32

    # Eccentricity and semi-major axis
    e = np.sqrt(ecosv2**2 + esinv2**2)
    a = p / (1 - e**2)

    # Compute the delta mean and eccentric anomalies
    deltam12 = 0
    if e**2 < 1:
        # Non-hyperbolic case
        n = np.sqrt(const.MU / a**3)
        s = magr2 / p * np.sqrt(1 - e**2) * esinv2
        c = magr2 / p * (e**2 + ecosv2)

        # Delta eccentric anomaly between obs 3 and 2
        sinde32 = magr3 / np.sqrt(a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s
        cosde32 = 1 - magr2 * magr3 / (a * p) * (1 - cosdv32)
        deltae32 = np.arctan2(sinde32, cosde32) + const.TWOPI * n23

        # Delta eccentric anomaly between obs 1 and 2
        sinde21 = magr1 / np.sqrt(a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s
        cosde21 = 1 - magr2 * magr1 / (a * p) * (1 - cosdv21)
        deltae21 = np.arctan2(sinde21, cosde21) + const.TWOPI * n12

        # Delta mean anomalies
        deltam32 = deltae32 + 2 * s * (np.sin(deltae32 / 2)) ** 2 - c * np.sin(deltae32)
        deltam12 = (
            -deltae21 + 2 * s * (np.sin(deltae21 / 2)) ** 2 + c * np.sin(deltae21)
        )
    else:
        # Hyperbolic case
        if a > 0:
            a = -a
            p = -p
        n = np.sqrt(const.MU / -(a**3))
        s = magr2 / p * np.sqrt(e**2 - 1) * esinv2
        c = magr2 / p * (e**2 + ecosv2)

        # Delta eccentric anomaly between obs 3 and 2
        sindh32 = magr3 / np.sqrt(-a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s
        deltah32 = np.log(sindh32 + np.sqrt(sindh32**2 + 1))
        deltam32 = (
            -deltah32 + 2 * s * (np.sinh(deltah32 / 2)) ** 2 + c * np.sinh(deltah32)
        )
        deltae32 = deltah32  # fix to match MATLAB

    # Calculate the f coefficients and quality estimate
    f1 = tau12 - deltam12 / n
    f2 = tau32 - deltam32 / n
    q1 = np.sqrt(f1**2 + f2**2)

    return r2, r3, f1, f2, q1, magr1, magr2, a, deltae32


def doubler(
    decl1: float,
    decl2: float,
    decl3: float,
    rtasc1: float,
    rtasc2: float,
    rtasc3: float,
    jd1: float,
    jdf1: float,
    jd2: float,
    jdf2: float,
    jd3: float,
    jdf3: float,
    rsite1: ArrayLike,
    rsite2: ArrayLike,
    rsite3: ArrayLike,
    rng1: float,
    rng2: float,
    pctchg: float = 0.005,
    tol_km: float = 0.1,
    max_iterations: int = 15,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve orbit determination problem using the double-r technique.

    References:
        Vallado: 2022, p. 449-452, Algorithm 53

    Args:
        decl1 (float): Declination of first sighting in radians
        decl2 (float): Declination of second sighting in radians
        decl3 (float): Declination of third sighting in radians
        rtasc1 (float): Right ascension of first sighting in radians
        rtasc2 (float): Right ascension of second sighting in radians
        rtasc3 (float): Right ascension of third sighting in radians
        jd1 (float): Julian date of first sighting (days from 4713 BC)
        jdf1 (float): Julian date fraction of first sighting (days from 4713 BC)
        jd2 (float): Julian date of second sighting (days from 4713 BC)
        jdf2 (float): Julian date fraction of second sighting (days from 4713 BC)
        jd3 (float): Julian date of third sighting (days from 4713 BC)
        jdf3 (float): Julian date fraction of third sighting (days from 4713 BC)
        rsite1 (array_like): ECI site position vector of first sighting in km
        rsite2 (array_like): ECI site position vector of second sighting in km
        rsite3 (array_like): ECI site position vector of third sighting in km
        rng1 (float): Range to first sighting in km
        rng2 (float): Range to second sighting in km
        pctchg (float, optional): Percentage change for iterative method (default 0.005)
        tol_km (float, optional): Position tolerance for convergence in km (default 0.1)
        max_iterations (int, optional): Maximum number of iterations (default 15)

    Returns:
        tuple: (r2, v2)
            r2 (np.ndarray): ECI position vector in km
            v2 (np.ndarray): ECI velocity vector in km/s
    """

    def calculate_fcoeffs(magr1in, magr2in, deltar1=None, deltar2=None):
        """Helper function to calculate partial derivatives and updated deltas."""
        if deltar1:
            magr1in += deltar1
        if deltar2:
            magr2in += deltar2

        # Call doubler to compute intermediate values
        _, _, f1, f2, q, *_ = doubler_iter(
            magr1in,
            magr2in,
            los1,
            los2,
            los3,
            rsite1,
            rsite2,
            rsite3,
            tau12,
            tau32,
            n12,
            n13,
            n23,
        )
        return f1, f2, q

    # Time differences in seconds
    tau12, tau13, tau32 = calculate_time_intervals(jd1, jdf1, jd2, jdf2, jd3, jdf3)

    # Period in seconds (assumed 1 day for Earth)
    n12 = np.floor(abs(tau12 / const.DAY2SEC))
    n13 = np.floor(abs(tau13 / const.DAY2SEC))
    n23 = np.floor(abs((tau12 + tau32) / const.DAY2SEC))

    # Line-of-sight unit vectors
    los1, los2, los3 = calculate_los_vectors(
        [decl1, decl2, decl3], [rtasc1, rtasc2, rtasc3]
    )

    # Iterative variables (make sure newqr is < oldqr to start)
    magr1in, magr2in = rng1, rng2
    magr1old, magr2old, oldqr, newqr = np.inf, np.inf, np.inf, 1e10

    # Main loop to get three values of the double-r for processing
    ktr = 0
    while (
        (abs(magr1in - magr1old) > tol_km or abs(magr2in - magr2old) > tol_km)
        and ktr < max_iterations
        and newqr < oldqr
    ):
        ktr += 1
        magr1o, magr2o, oldqr = magr1in, magr2in, newqr

        # Compute nominal values
        f1, f2, q1 = calculate_fcoeffs(magr1in, magr2in)

        # Compute perturbations
        # Re-calculate f1 and f2 with r1 = r1 + delta r1
        deltar1, deltar2 = pctchg * magr1o, pctchg * magr2o
        f1delr1, f2delr1, q2 = calculate_fcoeffs(magr1o, magr2o, deltar1=deltar1)
        pf1pr1, pf2pr1 = (f1delr1 - f1) / deltar1, (f2delr1 - f2) / deltar1

        # Re-calculate f1 and f2 with r2 = r2 + delta r2
        f1delr2, f2delr2, q3 = calculate_fcoeffs(magr1o, magr2o, deltar2=deltar2)
        pf1pr2, pf2pr2 = (f1delr2 - f1) / deltar2, (f2delr2 - f2) / deltar2

        # Compute delta updates
        magr1in, magr2in = magr1o, magr2o
        delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2
        delta1 = pf2pr2 * f1 - pf1pr2 * f2
        delta2 = pf1pr1 * f2 - pf2pr1 * f1
        deltar1 = -delta1 / delta if abs(delta) > const.SMALL else -delta1
        deltar2 = -delta2 / delta if abs(delta) > const.SMALL else -delta2

        # Limit corrections to avoid overshooting
        chkamt = 0.15
        deltar1 = np.clip(deltar1, -chkamt * magr1in, chkamt * magr1in)
        deltar2 = np.clip(deltar2, -chkamt * magr2in, chkamt * magr2in)

        # Update magnitudes and quality metric
        magr1old, magr2old = magr1in, magr2in
        magr1in += deltar1
        magr2in += deltar2
        newqr = np.sqrt(q1**2 + q2**2 + q3**2)

        # Reduce percentage change for next iteration
        pctchg *= 0.5

    # Final calculation for updated r2 and v2
    r2, r3, _, _, _, _, magr2, a, deltae32 = doubler_iter(
        magr1in,
        magr2in,
        los1,
        los2,
        los3,
        rsite1,
        rsite2,
        rsite3,
        tau12,
        tau32,
        n12,
        n13,
        n23,
    )

    f = 1 - a / magr2 * (1 - np.cos(deltae32))
    g = tau32 - np.sqrt(a**3 / const.MU) * (deltae32 - np.sin(deltae32))
    v2 = (r3 - f * r2) / g

    return r2, v2
