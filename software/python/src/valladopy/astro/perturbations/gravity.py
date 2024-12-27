# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const


def trigpoly(
    recef: ArrayLike, latgc: float, lon: float, order: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes accumulated Legendre polynomials and trigonometric terms.

    References:
        Vallado: 2013, p. 597, Eq. 8-57

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
