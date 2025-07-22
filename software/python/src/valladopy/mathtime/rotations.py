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
from .. import constants as const


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
    q = (Rot.from_quat(qa) * Rot.from_quat(qb)).as_quat()

    return -q if q[3] < 0 else q


def quat_transform(qi: ArrayLike, qf: ArrayLike) -> np.ndarray:
    """Computes the transformation quaternion qt such that qf = qi * qt.

    Args:
        qi (array_like): Initial quaternion as a 4-element array [x, y, z, w]
        qf (array_like): Final quaternion as a 4-element array [x, y, z, w]

    Returns:
        np.ndarray: Transformation quaternion qt as a 4-element array [x, y, z, w]
    """
    return (Rot.from_quat(qi).inv() * Rot.from_quat(qf)).as_quat()


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
    dcm = quat2dcm(q)
    return dcm[0, :], dcm[1, :], dcm[2, :]


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


def vec_by_dcm(dcm: ArrayLike, vec: ArrayLike) -> np.ndarray:
    """Applies a DCM rotation to a vector.

    Args:
        dcm (array_like): 3×3 direction cosine matrix
        vec (array_like): 3-element vector

    Returns:
        np.ndarray: Rotated 3-element vector
    """
    return np.asarray(dcm, dtype=float) @ np.asarray(vec, dtype=float)


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


def quat2dcm(q: ArrayLike) -> np.ndarray:
    """Converts a quaternion to a direction cosine matrix (DCM).

    Args:
        q (array_like): Quaternion as a 4-element array [x, y, z, w]

    Returns:
        np.ndarray: 3×3 direction cosine matrix
    """
    return Rot.from_quat(q).as_matrix().T


def dcm2quat(dcm: ArrayLike) -> np.ndarray:
    """Converts a direction cosine matrix (DCM) to a quaternion.

    Args:
        dcm (array_like): 3×3 direction cosine matrix

    Returns:
        np.ndarray: Quaternion as a 4-element array [x, y, z, w]
    """
    dcm = np.asarray(dcm, dtype=float)
    return Rot.from_matrix(dcm.T).as_quat()


def quat2euler(q: ArrayLike) -> tuple[float, float, float]:
    """Converts a quaternion to Euler angles (y–x–z sequence).

    Args:
        q (array_like): Quaternion as a 4-element array [x, y, z, w]

    Returns:
        tuple: (theta, phi, psi)
            theta (float): Angle around the y-axis in radians (0 to 2π)
            phi (float): Angle around the x-axis in radians (-π/2 to π/2)
            psi (float): Angle around the z-axis in radians (0 to 2π)
    """
    psi, phi, theta = Rot.from_quat(q).as_euler("zxy")  # intrinsic rotations
    return theta % const.TWOPI, phi, psi % const.TWOPI


def euler2quat(theta: float, phi: float, psi: float) -> np.ndarray:
    """Converts Euler angles (y–x–z rotation order) to a quaternion.

    Args:
        theta (float): Rotation about y in radians
        phi (float): Rotation about x in radians
        psi (float): Rotation about z in radians

    Returns:
        np.ndarray: Quaternion as a 4-element array [x, y, z, w]
    """
    q = Rot.from_euler("zxy", [psi, phi, theta]).as_quat()  # intrinsic rotations
    return -q if q[3] < 0 else q


def quat2eigen(q: ArrayLike) -> tuple[np.ndarray, float]:
    """Converts a quaternion to a rotation angle and eigen-axis.

    Args:
        q (array_like): Quaternion as a 4-element array [x, y, z, w]

    Returns:
        tuple: (axis, angle)
            axis (np.ndarray): 3-element unit vector
            angle (float): Rotation angle in radians

    Raises:
        ValueError: If the rotation angle is zero (undefined axis)
    """
    rot = Rot.from_quat(q)
    angle = rot.magnitude()
    axis = rot.as_rotvec() / angle if angle != 0 else np.zeros(3)

    if np.isclose(angle, 0):
        raise ValueError(
            "Eigen-axis is not well defined because the rotation angle is zero."
        )

    return axis, angle


def eigen2quat(axis: ArrayLike, angle: float) -> np.ndarray:
    """Converts a rotation axis and angle to a quaternion.

    Args:
        axis (array_like): Eigen-axis as a 3-element vector
        angle (float): Rotation angle in radians

    Returns:
        np.ndarray: Quaternion as a 4-element array [x, y, z, w]
    """
    rotvec = unit(axis) * angle
    return Rot.from_rotvec(rotvec).as_quat()


def dcm2euler(dcm: ArrayLike) -> tuple[float, float, float]:
    """Converts a direct cosine matrix (DCM) to Euler angles in y–x–z rotation order.

    Args:
        dcm (array_like): 3×3 direction cosine matrix

    Returns:
        tuple: (theta, phi, psi)
            theta (float): Angle around the y-axis in radians (0 to 2π)
            phi (float): Angle around the x-axis in radians (-π/2 to π/2)
            psi (float): Angle around the z-axis in radians (0 to 2π)
    """
    dcm = np.asarray(dcm, dtype=float)
    psi, phi, theta = Rot.from_matrix(dcm.T).as_euler("zxy", degrees=False)

    return theta % const.TWOPI, phi, psi % const.TWOPI


def euler2dcm(theta: float, phi: float, psi: float) -> np.ndarray:
    """Converts Euler angles to a direction cosine matrix (DCM).

    Args:
        theta (float): Rotation about y in radians
        phi (float): Rotation about x in radians
        psi (float): Rotation about z in radians

    Returns:
        np.ndarray: 3×3 direction cosine matrix
    """
    return Rot.from_euler("zxy", [psi, phi, theta]).as_matrix().T
