# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass
class GravityFieldData:
    c: np.ndarray = None
    s: np.ndarray = None
    normalized: bool = False


def leg_poly(x: float | ArrayLike, order: int) -> np.ndarray:
    """Calculates the Legendre polynomial of the specified order for the input x.
    x can be a scalar, a vector, or a matrix.

    Args:
        x: Input values. Can be a scalar, vector, or matrix.
        order: The order of the Legendre polynomial.

    Returns:
        np.ndarray: The computed Legendre polynomial values of the specified order.
    """
    # Initialize base cases for Legendre polynomials
    x = np.array(x, dtype=float)
    pn_p1 = 1
    pn = np.ones_like(x)

    if order > 0:
        pn_1 = x * pn  # first-order polynomial
        if order == 1:
            pn = pn_1
        else:
            for n in range(2, order + 1):
                pn_p1 = (x * (2 * n - 1) * pn_1 - (n - 1) * pn) / n
                pn = pn_1
                pn_1 = pn_p1
            pn = pn_p1

    return pn


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
