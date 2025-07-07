# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 20 Jan 2025
#
# Copyright (c) 2025
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import numpy as np
from numpy.typing import ArrayLike


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
