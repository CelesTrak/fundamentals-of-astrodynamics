# --------------------------------------------------------------------------------------
# Authors: Sal Alfano, David Vallado
# Date: 31 March 2011
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import numpy as np
from numpy.typing import ArrayLike


def posvelcov2pts(reci: ArrayLike, veci: ArrayLike, cov: ArrayLike) -> np.ndarray:
    """Generates 12 sigma points from position, velocity, and covariance using the
    Cholesky method.

    Args:
        reci (array_like): 3x1 ECI position vector
        veci (array_like): 3x1 ECI velocity vector
        cov (array_like): 6x6 ECI covariance matrix

    Returns:
        np.ndarray: 6x12 matrix of sigma points (position and velocity)

    Notes:
        - Units can be in km or m, but must be consistent
    """
    # Initialize the sigma points matrix
    sigmapts = np.zeros((6, 12))

    # Compute matrix square root using Cholesky decomposition
    s = np.sqrt(6) * np.linalg.cholesky(cov)

    # Generate sigma points
    for i in range(6):
        offset = s[:, i]
        jj = i * 2  # index for positive/negative perturbations

        # Positive perturbation
        sigmapts[:3, jj] = np.array(reci) + offset[:3]
        sigmapts[3:, jj] = np.array(veci) + offset[3:]

        # Negative perturbation
        sigmapts[:3, jj + 1] = np.array(reci) - offset[:3]
        sigmapts[3:, jj + 1] = np.array(veci) - offset[3:]

    return sigmapts
