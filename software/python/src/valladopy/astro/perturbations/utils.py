# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import math
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const


def legpolyn(
    latgc: float, order: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Computes Legendre polynomials for the gravity field.

    References:
        Vallado: 2022, p. 600-601, Eq. 8-56

    Args:
        latgc: Geocentric latitude of the satellite in radians (-pi to pi)
        order: Size of the gravity field (1 to ~170)

    Returns:
        tuple: (legarr_mu, legarr_gu, legarr_mn, legarr_gn)
            legarr_mu (np.ndarray): Montenbruck approach Legendre polynomials
            legarr_gu (np.ndarray): GTDS approach Legendre polynomials
            legarr_mn (np.ndarray): Normalized Montenbruck polynomials
            legarr_gn (np.ndarray): Normalized GTDS polynomials

    Notes:
        - Some recursions at high degree tesseral terms experience error for resonant
          orbits - these are valid for normalized and unnormalized expressions, as long
          as the remaining equations are consistent.
        - For satellite operations, orders up to about 120 are valid.
    """
    legarr_mu = np.zeros((order + 1, order + 1))
    legarr_gu = np.zeros((order + 1, order + 1))
    legarr_mn = np.zeros((order + 1, order + 1))
    legarr_gn = np.zeros((order + 1, order + 1))

    # Perform recursions (Montenbruck approach)
    legarr_mu[:2, :2] = [[1, 0], [np.sin(latgc), np.cos(latgc)]]

    # Legendre functions, zonal
    for n in range(2, order + 1):
        legarr_mu[n, n] = (2 * n - 1) * legarr_mu[1, 1] * legarr_mu[n - 1, n - 1]

    # Associated Legendre functions
    for n in range(2, order + 1):
        for m in range(n):
            if n == m + 1:
                legarr_mu[n, m] = (2 * m + 1) * legarr_mu[1, 0] * legarr_mu[m, m]
            else:
                legarr_mu[n, m] = (1 / (n - m)) * (
                    (2 * n - 1) * legarr_mu[1, 0] * legarr_mu[n - 1, m]
                    - (n + m - 1) * legarr_mu[n - 2, m]
                )

    # Normalize the Legendre polynomials
    for n in range(order + 1):
        for m in range(n + 1):
            factor = 1 if m == 0 else 2
            conv = np.sqrt(
                (math.factorial(n - m) * factor * (2 * n + 1)) / math.factorial(n + m)
            )
            legarr_mn[n, m] = conv * legarr_mu[n, m]

    # Perform recursions (GTDS approach)
    legarr_gu[:2, :2] = [[1, 0], [np.sin(latgc), np.cos(latgc)]]

    for n in range(2, order + 1):
        for m in range(n + 1):
            legarr_gu[n, m] = 0

    for n in range(2, order + 1):
        for m in range(n + 1):
            # Legendre functions, zonal
            if m == 0:
                legarr_gu[n, m] = (
                    (2 * n - 1) * legarr_gu[1, 0] * legarr_gu[n - 1, m]
                    - (n - 1) * legarr_gu[n - 2, m]
                ) / n
            else:
                # Associated Legendre functions
                if m == n:
                    legarr_gu[n, m] = (
                        (2 * n - 1) * legarr_gu[1, 1] * legarr_gu[n - 1, m - 1]
                    )
                else:
                    legarr_gu[n, m] = (
                        legarr_gu[n - 2, m]
                        + (2 * n - 1) * legarr_gu[1, 1] * legarr_gu[n - 1, m - 1]
                    )

    # Normalize the Legendre polynomials
    for n in range(order + 1):
        for m in range(n + 1):
            factor = 1 if m == 0 else 2
            conv1 = np.sqrt(
                (math.factorial(n - m) * factor * (2 * n + 1)) / math.factorial(n + m)
            )
            legarr_gn[n, m] = conv1 * legarr_gu[n, m]

    return legarr_mu, legarr_gu, legarr_mn, legarr_gn


def trigpoly(
    recef: ArrayLike, latgc: float, lon: float, order: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes accumulated Legendre polynomials and trigonometric terms.

    References:
        Vallado: 2022, p. 600-602

    Args:
        recef (array_like): ECEF satellite position vector in km
        latgc (float): Geocentric latitude of the satellite in radians
        lon (float): Longitude of the satellite in radians
        order (int): Size of the gravity field (1.. 2160..)

    Returns:
        tuple: (trig_arr, v_arr, w_arr)
            trig_arr (np.ndarray): Array of trigonometric terms
            v_arr (np.ndarray): V array of trigonometric terms
            w_arr (np.ndarray): W array of trigonometric terms
    """
    magr = np.linalg.norm(recef)
    l_ = 0

    # Initialize arrays
    trig_arr = np.zeros((order + 1, 3))
    v_arr = np.zeros((order + 2, order + 2))
    w_arr = np.zeros((order + 2, order + 2))

    # Trigonometric terms (GTDS approach)
    trig_arr[0, 0] = 0  # sin terms
    trig_arr[0, 1] = 1  # cos terms
    tlon = np.tan(latgc)
    trig_arr[1, 0] = np.sin(lon)
    trig_arr[1, 1] = np.cos(lon)
    clon = np.cos(lon)

    for m in range(2, order + 1):
        # Sine terms
        trig_arr[m, 0] = 2 * clon * trig_arr[m - 1, 0] - trig_arr[m - 2, 0]
        # Cosine terms
        trig_arr[m, 1] = 2 * clon * trig_arr[m - 1, 1] - trig_arr[m - 2, 1]
        # Tangent terms
        trig_arr[m, 2] = (m - 1) * tlon + tlon

    # Montenbruck approach for V and W arrays
    temp = const.RE / (magr * magr)
    v_arr[0, 0] = const.RE / magr
    v_arr[1, 0] = v_arr[0, 0] ** 2 * np.sin(latgc)

    for l_ in range(2, order + 2):
        x1 = ((2 * l_ - 1) / l_) * recef[1] * temp
        x2 = ((l_ - 1) / l_) * temp * const.RE
        v_arr[l_, 0] = x1 * v_arr[l_ - 1, 0] - x2 * v_arr[l_ - 2, 0]

    # Tesseral and sectoral values for L = m
    for l_ in range(1, order + 2):
        m = l_
        x1 = (2 * m - 1) * recef[0] * temp
        x2 = recef[1] * temp
        v_arr[l_, m] = x1 * v_arr[l_ - 1, m - 1] - x2 * w_arr[l_ - 1, m - 1]
        w_arr[l_, m] = x1 * w_arr[l_ - 1, m - 1] - x2 * v_arr[l_ - 1, m - 1]

    for m in range(l_ + 1, order + 1):
        if m <= order:
            x = (2 * l_ - 1) / (l_ - m) * recef[1] * temp
            v_arr[l_ + 1, m] = x * v_arr[l_, m]
            w_arr[l_ + 1, m] = x * w_arr[l_, m]

        for l2 in range(m + 2, order + 2):
            x1 = ((2 * l2 - 1) / (l2 - m)) * recef[1] * temp
            x2 = ((l2 + m - 1) / (l2 - m)) * temp * const.RE
            v_arr[l2, m] = x1 * v_arr[l2 - 1, m] - x2 * v_arr[l2 - 2, m]
            w_arr[l2, m] = x1 * w_arr[l2 - 1, m] - x2 * w_arr[l2 - 2, m]

    return trig_arr, v_arr, w_arr
