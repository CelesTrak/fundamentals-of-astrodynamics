# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import numpy as np
from numpy.typing import ArrayLike

from ..constants import SMALL


###############################################################################
# Axes Rotations
###############################################################################


def rot1(vec, xval):
    """Rotation about the 1st axis (x-axis)

    Args:
        vec (array_like): Input vector
        xval (float): Angle of rotation in radians

    Returns:
        array_like: Rotated vector
    """
    c, s = np.cos(xval), np.sin(xval)
    outvec = np.array([vec[0], c * vec[1] + s * vec[2], c * vec[2] - s * vec[1]])
    return outvec


def rot2(vec, xval):
    """Rotation about the 2nd axis (y-axis)

    Args:
        vec (array_like): Input vector
        xval (float): Angle of rotation in radians

    Returns:
        array_like: Rotated vector
    """
    c, s = np.cos(xval), np.sin(xval)
    outvec = np.array([c * vec[0] - s * vec[2], vec[1], c * vec[2] + s * vec[0]])
    return outvec


def rot3(vec, xval):
    """Rotation about the 3rd axis (z-axis)

    Args:
        vec (array_like): Input vector
        xval (float): Angle of rotation in radians

    Returns:
        array_like: Rotated vector
    """
    c, s = np.cos(xval), np.sin(xval)
    outvec = np.array([c * vec[0] + s * vec[1], c * vec[1] - s * vec[0], vec[2]])
    return outvec


###############################################################################
# Rotation Matrices
###############################################################################


def rot1mat(xval: float) -> ArrayLike:
    """Rotation matrix for an input angle about the first axis.
    Assume: Use of "column" vectors

    Args:
        xval (float): Angle of rotation in radians

    Returns:
        outmat (array_like): Rotation matrix
    """
    c, s = np.cos(xval), np.sin(xval)

    outmat = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

    return outmat


def rot2mat(xval: float) -> ArrayLike:
    """Rotation matrix for an input angle about the second axis.
    Assume: Use of "column" vectors

    Args:
        xval (float): Angle of rotation in radians

    Returns:
        outmat (array_like): Rotation matrix
    """
    c, s = np.cos(xval), np.sin(xval)

    outmat = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])

    return outmat


def rot3mat(xval: float) -> ArrayLike:
    """Rotation matrix for an input angle about the third axis.
    Assume: Use of "column" vectors

    Args:
        xval (float): Angle of rotation in radians

    Returns:
        outmat (array_like): Rotation matrix
    """
    c, s = np.cos(xval), np.sin(xval)

    outmat = np.array(
        [
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    return outmat


###############################################################################
# Vector Math
###############################################################################


def angle(v1, v2):
    """Calculates the angle between two vectors in radians

    This function computes the angle between two vectors using the dot product
    and the magnitudes of the vectors. The function handles cases where the dot
    product might slightly exceed the interval [-1, 1] due to numerical
    precision issues by clipping the value. If either vector has zero
    magnitude, the function returns `np.nan`to indicate that the angle is not
    computable.

    Args:
        v1 (numpy.ndarray): The first vector.
        v2 (numpy.ndarray): The second vector.

    Returns:
        float: The angle between the two vectors in radians. Returns `np.nan`
               if either vector is zero.
    """
    mag_v1 = np.linalg.norm(v1)
    mag_v2 = np.linalg.norm(v2)

    if mag_v1 * mag_v2 > SMALL**2:
        cos_angle = np.dot(v1, v2) / (mag_v1 * mag_v2)
        cos_angle = np.clip(cos_angle, -1, 1)  # Keep cosine within domain
        return np.arccos(cos_angle)
    else:
        return np.nan


def unit(v):
    """Returns the unit vector of a given vector

    Args:
        v (array-like): The input vector

    Returns:
        numpy.ndarray: The unit vector of the input vector
                       (v / ||v|| if ||v|| > 0, 0 otherwise)
    """
    mag = np.linalg.norm(v)
    return v / mag if mag > SMALL else np.zeros_like(v)
