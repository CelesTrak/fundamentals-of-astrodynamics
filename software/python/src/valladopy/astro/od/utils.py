# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 15 Jan 2008
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike


def finite_diff(
    pertelem: int,
    percentchg: float,
    deltaamtchg: float,
    xnom: ArrayLike,
    growth_factor: float = 1.4,
) -> Tuple[float, np.ndarray]:
    """Perturbs the components of the state vector for finite differencing.

    References:
        Vallado: 2022, p. 779-780

    Args:
        pertelem (int): Index of the element to perturb (0-based)
        percentchg (float): Amount to modify the vector by in finite differencing
        deltaamtchg (float): Tolerance for small value in finite differencing
        xnom (np.ndarray): State vector to perturb
        growth_factor (float, optional): Factor to increase perturbation by
                                         (defaults to 1.4)

    Returns:
        tuple (deltaamt, xnomp):
            deltaamt (float): Amount of perturbation
            xnomp (np.ndarray): Perturbed state vector
    """
    deltaamt = 0
    xnomp = xnom.copy()

    for _ in range(len(xnom)):
        deltaamt = xnom[pertelem] * percentchg
        xnomp[pertelem] += deltaamt

        if abs(deltaamt) >= deltaamtchg:
            break  # perturbation is sufficient
        else:
            percentchg *= growth_factor  # increase perturbation by growth factor

    return deltaamt, xnomp
