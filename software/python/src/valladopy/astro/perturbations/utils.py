# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass
class GravityFieldData:
    c: np.ndarray = None
    s: np.ndarray = None
    normalized: bool = False


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


def read_gravity_field(filename: str, normalized: bool) -> GravityFieldData:
    """Reads and stores gravity field coefficients.

    Args:
        filename (str): The filename of the gravity field data
        normalized (bool): True if the gravity field data is normalized

    Returns:
        GravityFieldData: A dataclass containing gravity field data:
            - c (np.ndarray): Cosine coefficients
            - s (np.ndarray): Sine coefficients
            - normalized (bool): True if the gravity field data is normalized
    """
    # Load gravity field data
    file_data = np.loadtxt(filename)

    # Get the maximum degree of the gravity field
    max_degree = int(np.max(file_data[:, 0]))

    # Initialize gravity field data
    gravarr = GravityFieldData()
    gravarr.c = np.zeros((max_degree + 1, max_degree + 1))
    gravarr.s = np.zeros((max_degree + 1, max_degree + 1))
    gravarr.normalized = normalized

    # Store gravity field coefficients
    for row in file_data:
        n, m = int(row[0]), int(row[1])
        c_value, s_value = row[2], row[3]
        gravarr.c[n, m] = c_value
        gravarr.s[n, m] = s_value

    return gravarr
