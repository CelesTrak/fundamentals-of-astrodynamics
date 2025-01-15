# --------------------------------------------------------------------------------------
# Authors: Sal Alfano, David Vallado
# Date: 31 March 2011
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike


def poscov2pts(reci: ArrayLike, cov: ArrayLike) -> np.ndarray:
    """Generates 6 sigma points from position and covariance using the Cholesky method.

    Args:
        reci (array_like): 3x1 ECI position vector
        cov (array_like): 3x3 ECI covariance matrix

    Returns:
        np.ndarray: 3x6 matrix of sigma points (position only)

    Notes:
        - Units can be in km or m, but must be consistent
    """
    # Initialize the sigma points matrix
    sigmapts = np.zeros((3, 6))

    # Compute matrix square root using Cholesky decomposition
    s = np.sqrt(3) * np.linalg.cholesky(cov)

    # Generate sigma points
    for i in range(3):
        offset = s[:, i]
        jj = i * 2  # index for positive/negative perturbations

        # Positive perturbation
        sigmapts[:, jj] = np.array(reci) + offset[:3]

        # Negative perturbation
        sigmapts[:, jj + 1] = np.array(reci) - offset[:3]

    return sigmapts


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


def remakecov(sigmapts: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
    """Calculates the mean vector and covariance matrix from propagated sigma points.

    Args:
        sigmapts (array_like): n_dim x n_pts matrix of propagated points from the
        square root algorithm

    Returns:
        tuple: (yu, cov)
            yu (np.ndarray): n_dim x 1 mean vector
            cov (np.ndarray): n_dim x n_dim covariance matrix
    """
    sigmapts = np.array(sigmapts)
    n_dim, n_pts = sigmapts.shape

    # Compute the mean vector
    yu = np.mean(sigmapts, axis=1, keepdims=True)

    # Compute the deviation matrix
    y = sigmapts - yu

    # Compute the covariance matrix
    cov = (y @ y.T) / n_pts

    return yu, cov
