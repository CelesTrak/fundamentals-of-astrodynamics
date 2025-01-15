# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import numpy as np
from numpy.typing import ArrayLike


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
