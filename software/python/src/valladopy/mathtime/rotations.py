# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 20 Jan 2025
#
# Copyright (c) 2025
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import numpy as np
from numpy.typing import ArrayLike
from scipy.spatial.transform import Rotation as Rot

from .vector import unit


def quat_multiply(
    qa: ArrayLike, qb: ArrayLike, dir_a: int = 1, dir_b: int = 1
) -> np.ndarray:
    """Multiplies two quaternions with optional direction signs.

    Args:
        qa (array_like): First quaternion as a 4-element array [x, y, z, w]
        qb (array_like): Second quaternion as a 4-element array [x, y, z, w]
        dir_a (int, optional): Direction of first quaternion (+1 or -1)
        dir_b (int, optional): Direction of second quaternion (+1 or -1)

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


def vec_by_quat(q: ArrayLike, vec: ArrayLike, direction: int = 1) -> np.ndarray:
    """Rotates a 3D vector by a quaternion.

    Args:
        q (array_like): Quaternion as a 4-element array [x, y, z, w]
        vec (array_like): 3D vector to rotate as a 3-element array
        direction (int, optional): Direction of rotation (+1 or -1)

    Returns:
        np.ndarray: Rotated vector as a 3-element array
    """
    # Validate input quaternion and vector
    if direction not in (-1, 1):
        raise ValueError("Direction must be +1 or -1.")

    x, y, z, w = q
    w *= direction

    t1 = z * vec[1] - y * vec[2]
    t2 = x * vec[2] - z * vec[0]
    t3 = y * vec[0] - x * vec[1]

    vec_out = np.zeros(3)
    vec_out[0] = vec[0] + 2 * (t1 * w + t2 * z - t3 * y)
    vec_out[1] = vec[1] + 2 * (t2 * w + t3 * x - t1 * z)
    vec_out[2] = vec[2] + 2 * (t3 * w + t1 * y - t2 * x)

    return vec_out


def quat2rv(q: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Converts a 7-element orbit state quaternion to position and velocity vectors.

    Args:
        q (array_like): Orbit state quaternion as a 7-element array
                        [x, y, z, w, rmag, dot, omega], where:
                        - x, y, z, w (floats): Quaternion components
                        - rmag (float): Magnitude of the position vector
                        - dot (float): Radial velocity component
                        - omega (float): Angular velocity component

    Returns:
        tuple: (r, v)
            r (np.ndarray): Position vector in the inertial frame as a 3-element array
            v (np.ndarray): Velocity vector in the inertial frame as a 3-element array
    """
    x, y, z, w, rmag, dot, omega = q

    # LVLH-frame vectors
    vel_local = np.array([rmag * omega, 0, -dot])
    r_local = np.array([0, 0, -rmag])  # +Z nadir (toward Earth)

    # Inverse transform from LVLH to inertial
    rot = Rot.from_quat([x, y, z, w])
    r = rot.apply(r_local)
    v = rot.apply(vel_local)

    return r, v


def rv2quat(r: ArrayLike, v: ArrayLike) -> np.ndarray:
    """Converts a position and velocity vector to a 7-element orbit quaternion state.

    Args:
        r (array_like): Position vector in the inertial frame as a 3-element array
        v (array_like): Velocity vector in the inertial frame as a 3-element array

    Returns:
        np.ndarray: 7-element array representing the orbit quaternion state
                    [x, y, z, w, rmag, dot, omega], where:
                    - x, y, z, w (floats): Quaternion components
                    - rmag (float): Magnitude of the position vector
                    - dot (float): Radial velocity component
                    - omega (float): Angular velocity component
    """
    # Calculate orbit quantities
    r2 = np.dot(r, r)
    magr = np.sqrt(r2)
    dot = np.dot(r, v) / magr
    omega_vec = np.cross(r, v) / r2
    omega = np.linalg.norm(omega_vec)

    # Construct LVLH basis vectors
    z_hat = -np.array(r) / magr
    y_hat = -omega_vec / omega
    x_hat = np.cross(y_hat, z_hat)

    # Rotation matrix from LVLH to inertial
    dcm = np.column_stack([x_hat, y_hat, z_hat])

    # Convert rotation matrix to quaternion
    q = Rot.from_matrix(dcm).as_quat()  # [x, y, z, w]

    return np.array([*q, magr, dot, omega])
