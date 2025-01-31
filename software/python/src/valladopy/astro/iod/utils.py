# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import warnings
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const
from ...mathtime.vector import unit


def gibbs(
    r1: ArrayLike,
    r2: ArrayLike,
    r3: ArrayLike,
    tol_angle: float = np.radians(1),
    orbit_tol: float = 1e-6,
) -> Tuple[np.ndarray, float, float, float]:
    """Determine the velocity at the middle point of the 3 given position vectors using
    the Gibbs method.

    References:
        Vallado: 2022, p. 460-467, Algorithm 54

    Args:
        r1 (array_like): ECI position vector #1 in km
        r2 (array_like): ECI position vector #2 in km
        r3 (array_like): ECI position vector #3 in km
        tol_angle (float, optional): Tolerance for angles (1 degree)
        orbit_tol (float, optional): Tolerance for orbit calculations (default 1e-6)

    Returns:
        tuple: (v2, theta12, theta23, copa)
            v2 (np.ndarray): ECI velocity vector at r2 in km/s
            theta12 (float): Angle between r1 and r2 in radians
            theta23 (float): Angle between r2 and r3 in radians
            copa (float): Co-planarity angle in radians
    """
    # Initialize variables
    theta12, theta23 = 0, 0

    # Magnitudes of position vectors
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)
    magr3 = np.linalg.norm(r3)

    # Initialize velocity vector
    v2 = np.zeros(3)

    # Cross products
    p = np.cross(r2, r3)
    q = np.cross(r3, r1)
    w = np.cross(r1, r2)

    # Co-planarity angle
    copa = np.arcsin(np.dot(unit(p), unit(r1)))

    # Check for coplanarity
    if abs(copa) > np.sin(tol_angle):
        warnings.warn(
            "Vectors are not coplanar - results might be inaccurate.", UserWarning
        )

    # Sum of cross products
    d = p + q + w
    magd = np.linalg.norm(d)

    # Weighted sum of position vectors
    n = magr1 * p + magr2 * q + magr3 * w
    magn = np.linalg.norm(n)

    # Check if orbit determination is possible
    # Both `d` and `n` must be in the same direction and non-zero
    if any(x < orbit_tol for x in [magd, magn, np.dot(unit(n), unit(d))]):
        warnings.warn("Orbit determination is not possible.", UserWarning)
        return v2, theta12, theta23, copa

    # Angles between position vectors
    theta12 = np.arccos(np.clip(np.dot(r1, r2) / (magr1 * magr2), -1, 1))
    theta23 = np.arccos(np.clip(np.dot(r2, r3) / (magr2 * magr3), -1, 1))

    # Differences in position vector magnitudes
    r1mr2 = magr1 - magr2
    r3mr1 = magr3 - magr1
    r2mr3 = magr2 - magr3

    # S vector
    s = r1mr2 * np.array(r3) + r3mr1 * np.array(r2) + r2mr3 * np.array(r1)

    # B vector
    b = np.cross(d, r2)

    # Scaling factor
    lg = np.sqrt(const.MU / (magd * magn))

    # Compute velocity at r2
    tover2 = lg / magr2
    v2 = tover2 * b + lg * s

    return v2, theta12, theta23, copa


def hgibbs(
    r1: ArrayLike,
    r2: ArrayLike,
    r3: ArrayLike,
    jd1: float,
    jd2: float,
    jd3: float,
    tol_angle: float = np.radians(1),
) -> Tuple[np.ndarray, float, float, float]:
    """Determines the velocity at the middle point of the 3 given position vectors
    using the Herrick-Gibbs method.

    References:
        Vallado: 2022, p. 467-472, Algorithm 55

    Args:
        r1 (array_like): ECI position vector #1 in km
        r2 (array_like): ECI position vector #2 in km
        r3 (array_likey): ECI position vector #3 in km
        jd1 (float): Julian date of the 1st sighting in days
        jd2 (float): Julian date of the 2nd sighting in days
        jd3 (float): Julian date of the 3rd sighting in days
        tol_angle (float, optional): Tolerance for angles in radians (default = 1 deg)

    Returns:
        tuple: (v2, theta12, theta23, copa)
            v2 (np.ndarray): Velocity vector at r2 in km/s
            theta12 (float): Angle between r1 and r2 in radians
            theta23 (float): Angle between r2 and r3 in radians
            copa (float): Co-planarity angle in radians
    """
    # Magnitudes of position vectors
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)
    magr3 = np.linalg.norm(r3)

    # Time differences in seconds
    dt21 = (jd2 - jd1) * const.DAY2SEC
    dt31 = (jd3 - jd1) * const.DAY2SEC
    dt32 = (jd3 - jd2) * const.DAY2SEC

    # Calculate coplanarity
    p = np.cross(r2, r3)
    copa = np.arcsin(np.dot(unit(p), unit(r1)))

    # Check for coplanarity
    if abs(copa) > np.sin(tol_angle):
        warnings.warn(
            "Vectors are not coplanar - results might be inaccurate.", UserWarning
        )

    # Check angle tolerance between position vectors
    theta12 = np.arccos(np.clip(np.dot(r1, r2) / (magr1 * magr2), -1, 1))
    theta23 = np.arccos(np.clip(np.dot(r2, r3) / (magr2 * magr3), -1, 1))

    # Warn if angles exceed tolerance
    if theta12 > tol_angle or theta23 > tol_angle:
        warnings.warn("Angles between vectors exceed tolerance.", UserWarning)

    # Herrick-Gibbs method
    term1 = -dt32 * (1 / (dt21 * dt31) + const.MU / (12 * magr1**3))
    term2 = (dt32 - dt21) * (1 / (dt21 * dt32) + const.MU / (12 * magr2**3))
    term3 = dt21 * (1 / (dt32 * dt31) + const.MU / (12 * magr3**3))

    # Calculate velocity at r2
    v2 = term1 * np.array(r1) + term2 * np.array(r2) + term3 * np.array(r3)

    return v2, theta12, theta23, copa
