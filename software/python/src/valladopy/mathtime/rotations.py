# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 20 Jan 2025
#
# Copyright (c) 2025
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import numpy as np
from numpy.typing import ArrayLike

from .vector import unit


def quat_multiply(
    qa: ArrayLike, qb: ArrayLike, dir_a: int = 1, dir_b: int = 1
) -> np.ndarray:
    """Multiplies two quaternions with optional direction signs.

    Args:
        qa (array_like): First quaternion as a 4-element array [x, y, z, w]
        qb (array_like): Second quaternion as a 4-element array [x, y, z, w]
        dir_a (int, optional): Direction of first quaternion (1=direct, -1=inverse)
        dir_b (int, optional): Direction of second quaternion (1=direct, -1=inverse)

    Returns:
        np.ndarray: Resulting quaternion as a 4-element array [x, y, z, w]

    Raises:
        ValueError: If direction signs are not -1 or 1.
    """
    # Validate input direction signs
    if dir_a not in (-1, 1) or dir_b not in (-1, 1):
        raise ValueError("Direction signs must be -1 or 1.")

    # Copy input quaternions and apply direction signs
    qa = np.asarray(qa, dtype=float).copy()
    qb = np.asarray(qb, dtype=float).copy()
    qa[3] *= dir_a
    qb[3] *= dir_b

    # Multiply the quaternions
    q = np.zeros(4)
    q[0] = qb[3] * qa[0] + qb[2] * qa[1] - qb[1] * qa[2] + qb[0] * qa[3]
    q[1] = qb[3] * qa[1] - qb[2] * qa[0] + qb[1] * qa[3] + qb[0] * qa[2]
    q[2] = qb[3] * qa[2] + qb[2] * qa[3] + qb[1] * qa[0] - qb[0] * qa[1]
    q[3] = qb[3] * qa[3] - qb[2] * qa[2] - qb[1] * qa[1] - qb[0] * qa[0]

    if q[3] < 0:
        q = -q

    return q


def quat_transform(qi: ArrayLike, qf: ArrayLike) -> np.ndarray:
    """Computes the transformation quaternion qt such that qf = qi * qt.

    Args:
        qi (array_like): Initial quaternion as a 4-element array [x, y, z, w]
        qf (array_like): Final quaternion as a 4-element array [x, y, z, w]

    Returns:
        np.ndarray: Transformation quaternion qt as a 4-element array [x, y, z, w]
    """
    dq = np.zeros(4)
    dq[0] = qi[3] * qf[0] - qi[0] * qf[3] - qi[1] * qf[2] + qi[2] * qf[1]
    dq[1] = qi[3] * qf[1] - qi[1] * qf[3] - qi[2] * qf[0] + qi[0] * qf[2]
    dq[2] = qi[3] * qf[2] - qi[2] * qf[3] - qi[0] * qf[1] + qi[1] * qf[0]
    dq[3] = qi[0] * qf[0] + qi[1] * qf[1] + qi[2] * qf[2] + qi[3] * qf[3]

    return unit(dq)


def quat2body(q: ArrayLike) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes body-frame unit vectors from a quaternion.

    Args:
        q (array_like): Quaternion as a 4-element array [x, y, z, w]

    Returns:
        tuple: (x_axis, y_axis, z_axis)
            x_axis (np.ndarray): Unit vector in the x direction
            y_axis (np.ndarray): Unit vector in the y direction
            z_axis (np.ndarray): Unit vector in the z direction
    """
    x, y, z, w = q
    x_axis = np.array([1 - 2 * (y**2 + z**2), 2 * (x * y + z * w), 2 * (x * z - y * w)])
    y_axis = np.array([2 * (x * y - z * w), 1 - 2 * (x**2 + z**2), 2 * (y * z + x * w)])
    z_axis = np.array([2 * (x * z + y * w), 2 * (y * z - x * w), 1 - 2 * (x**2 + y**2)])

    return x_axis, y_axis, z_axis
