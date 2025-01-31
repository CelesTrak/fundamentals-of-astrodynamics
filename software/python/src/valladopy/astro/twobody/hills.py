# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ...constants import RE, MU


def hillsr(
    r: ArrayLike, v: ArrayLike, alt: float, dts: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate position and velocity information for Hill's (Clohessy-Wiltshire)
    equations.

    References:
        Vallado: 2022, p. 397-401, Algorithm 48

    Args:
        r (array_like): Initial relative position of the interceptor in km
        v (array_like): Initial relative velocity of the interceptor in km/s
        alt (float): Altitude of the target satellite in km
        dts (float): Desired time in seconds

    Returns:
        tuple: (rint, vint)
            rint (np.ndarray): Final relative position of the interceptor in km
            vint (np.ndarray): Final relative velocity of the interceptor in km/s

    Notes:
        - Position and velocity vectors are in the RSW frame
        - Distance units for r and v are flexible, but must be consistent
    """
    # Calculate orbital parameters
    radius = RE + alt
    omega = np.sqrt(MU / (radius**3))
    nt = omega * dts
    cosnt = np.cos(nt)
    sinnt = np.sin(nt)

    # Determine new positions
    rint = np.zeros(3)
    rint[0] = (
        (v[0] / omega) * sinnt
        - ((2 * v[1] / omega) + 3 * r[0]) * cosnt
        + ((2 * v[1] / omega) + 4 * r[0])
    )
    rint[1] = (
        (2 * v[0] / omega) * cosnt
        + ((4 * v[1] / omega) + 6 * r[0]) * sinnt
        + (r[1] - (2 * v[0] / omega))
        - (3 * v[1] + 6 * omega * r[0]) * dts
    )
    rint[2] = r[2] * cosnt + (v[2] / omega) * sinnt

    # Determine new velocities
    vint = np.zeros(3)
    vint[0] = v[0] * cosnt + (2 * v[1] + 3 * omega * r[0]) * sinnt
    vint[1] = (
        -2 * v[0] * sinnt
        + (4 * v[1] + 6 * omega * r[0]) * cosnt
        - (3 * v[1] + 6 * omega * r[0])
    )
    vint[2] = -r[2] * omega * sinnt + v[2] * cosnt

    return rint, vint


def hillsv(r: ArrayLike, alt: float, dts: float, tol: float = 1e-6) -> np.ndarray:
    """Calculate the initial velocity for Hill's (Clohessy-Wiltshire) equations.

    References:
        Vallado: 2022, p. 410-414, Eq. 6-66

    Args:
        r (array_like): Initial position vector of the interceptor in km
        alt (float): Altitude of the target satellite in km
        dts (float): Desired time in seconds
        tol (float, optional): Tolerance for calculations (defaults to 1e-6)

    Returns:
        np.ndarray: Initial velocity vector of the interceptor in km/s

    Notes:
        - Position and velocity vectors are in the RSW frame
        - Distance units for r are flexible, and velocity units are consistent
    """
    # Calculate the orbital parameters
    radius = RE + alt
    omega = np.sqrt(MU / (radius**3))
    nt = omega * dts
    cosnt, sinnt = np.cos(nt), np.sin(nt)

    # Numerator and denominator for the initial velocity
    numkm = (6 * r[0] * (nt - sinnt) - r[1]) * omega * sinnt - 2 * omega * r[0] * (
        4 - 3 * cosnt
    ) * (1 - cosnt)
    denom = (4 * sinnt - 3 * nt) * sinnt + 4 * (1 - cosnt) ** 2

    # Determine initial velocity
    v = np.zeros(3)
    v[1] = numkm / denom if abs(denom) > tol else 0
    if abs(sinnt) > tol:
        v[0] = -(omega * r[0] * (4 - 3 * cosnt) + 2 * (1 - cosnt) * v[1]) / sinnt
    v[2] = -r[2] * omega / np.tan(nt)

    return v
