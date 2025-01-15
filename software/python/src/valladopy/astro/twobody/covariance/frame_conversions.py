# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 21 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ...time.data import IAU80Array
from .. import frame_conversions as fc
from ..newton import newtonnu
from ....constants import KM2M, MUM, SMALL, TWOPI
from ....mathtime.vector import unit


class AnomalyType(Enum):
    MEAN_A = "meana"
    TRUE_A = "truea"
    MEAN_N = "meann"
    TRUE_N = "truen"


def is_mean_anomaly(anom_type: AnomalyType) -> bool:
    """Checks if the anomaly type is a mean anomaly.

    Args:
        anom_type (AnomalyType): Anomaly type

    Returns:
        bool: True if the anomaly type is a mean anomaly
    """
    return anom_type in [AnomalyType.MEAN_A, AnomalyType.MEAN_N]


########################################################################################
# Classical <-> Cartesian Elements
# TODO: This section could use some cleanup
########################################################################################


def _compute_partials_a(a, n, reci_m, veci_m):
    """Compute partial derivatives of a w.r.t. (rx, ry, rz, vx, vy, vz)."""
    tm_a = np.zeros(6)
    p0 = 2.0 * a**2 / np.linalg.norm(reci_m) ** 3
    p1 = 2.0 / (n**2 * a)
    tm_a[:3] = p0 * reci_m
    tm_a[3:] = p1 * veci_m
    return tm_a


def _compute_partials_ecc(reci_m, veci_m, ecc_vec, ecc):
    """Compute partial derivatives of ecc w.r.t. (rx, ry, rz, vx, vy, vz)."""
    rx, ry, rz = reci_m
    vx, vy, vz = veci_m
    tm_ecc = np.zeros(6)
    p0 = 1.0 / (MUM * ecc)

    magr = np.linalg.norm(reci_m)
    magr3 = magr**3

    tm_ecc[0] = -p0 * (
        ((vx * vy - MUM * rx * ry / magr3) * ecc_vec[1])
        + ((vx * vz - MUM * rx * rz / magr3) * ecc_vec[2])
        - (vy**2 + vz**2 - MUM / magr + MUM * rx**2 / magr3) * ecc_vec[0]
    )
    tm_ecc[1] = -p0 * (
        ((vx * vy - MUM * rx * ry / magr3) * ecc_vec[0])  # TODO: check
        + ((vy * vz - MUM * ry * rz / magr3) * ecc_vec[2])
        - (vx**2 + vz**2 - MUM / magr + MUM * ry**2 / magr3) * ecc_vec[1]
    )
    tm_ecc[2] = -p0 * (
        ((vx * vz - MUM * rx * rz / magr3) * ecc_vec[0])
        + ((vy * vz - MUM * ry * rz / magr3) * ecc_vec[1])
        - (vx**2 + vy**2 - MUM / magr + MUM * rz**2 / magr3) * ecc_vec[2]
    )
    tm_ecc[3] = -p0 * (
        ((rx * vy - 2 * ry * vx) * ecc_vec[1])
        + ((ry * vy + rz * vz) * ecc_vec[0])
        + ((rx * vz - 2 * rz * vx) * ecc_vec[2])
    )
    tm_ecc[4] = -p0 * (
        ((ry * vx - 2 * rx * vy) * ecc_vec[0])
        + ((rx * vx + rz * vz) * ecc_vec[1])
        + ((ry * vz - 2 * rz * vy) * ecc_vec[2])
    )
    tm_ecc[5] = -p0 * (
        ((rx * vx + ry * vy) * ecc_vec[2])
        + ((rz * vx - 2 * rx * vz) * ecc_vec[0])
        + ((rz * vy - 2 * ry * vz) * ecc_vec[1])
    )

    return tm_ecc


def _compute_partials_incl(reci_m, veci_m, h_vec, node):
    """Compute partial derivatives of inclination w.r.t. (rx, ry, rz, vx, vy, vz)."""
    rx, ry, rz = reci_m
    vx, vy, vz = veci_m
    h = np.linalg.norm(h_vec)
    tm_incl = np.zeros(6)
    p3 = 1.0 / node
    tm_incl[0] = -p3 * (vy - h_vec[2] * (vy * h_vec[2] - vz * h_vec[1]) / h**2)
    tm_incl[1] = p3 * (vx - h_vec[2] * (vx * h_vec[2] - vz * h_vec[0]) / h**2)
    tm_incl[2] = -p3 * (h_vec[2] * (vy * h_vec[0] - vx * h_vec[1]) / h**2)
    tm_incl[3] = p3 * (ry - h_vec[2] * (ry * h_vec[2] - rz * h_vec[1]) / h**2)
    tm_incl[4] = -p3 * (rx - h_vec[2] * (rx * h_vec[2] - rz * h_vec[0]) / h**2)
    tm_incl[5] = p3 * (h_vec[2] * (ry * h_vec[0] - rx * h_vec[1]) / h**2)

    return tm_incl


def _compute_partials_node(reci_m, veci_m, node_vec, node):
    """Compute partial derivatives of node w.r.t. (rx, ry, rz, vx, vy, vz)."""
    rx, ry, rz = reci_m
    vx, vy, vz = veci_m
    tm_node = np.zeros(6)
    p4 = 1.0 / (node**2)
    tm_node[0] = -p4 * vz * node_vec[1]
    tm_node[1] = p4 * vz * node_vec[0]
    tm_node[2] = p4 * (vx * node_vec[1] - vy * node_vec[0])
    tm_node[3] = p4 * rz * node_vec[1]
    tm_node[4] = -p4 * rz * node_vec[0]
    tm_node[5] = p4 * (ry * node_vec[0] - rx * node_vec[1])

    return tm_node


def _compute_partials_argp(
    tm_ecc, reci_m, veci_m, ecc_vec, ecc, h_vec, node, n_dot_e, w_scale
):
    """Compute partial derivatives of argp w.r.t. (rx, ry, rz, vx, vy, vz)."""
    rx, ry, rz = reci_m
    vx, vy, vz = veci_m
    magr = np.linalg.norm(reci_m)
    tm_argp = np.zeros(6)
    tm_argp[0] = (
        -h_vec[1] * (vy**2 + vz**2 - MUM / magr + MUM * rx**2 / magr**3)
        - h_vec[0] * (vx * vy - MUM * rx * ry / magr**3)
        + vz * MUM * ecc_vec[0]
    )
    tm_argp[0] = (
        tm_argp[0] / (MUM * node * ecc)
        + vz * h_vec[1] * n_dot_e / (node**3 * ecc)
        - tm_ecc[0] * n_dot_e / (node * ecc**2)
    ) * w_scale

    tm_argp[1] = (
        h_vec[0] * (vx**2 + vz**2 - MUM / magr + MUM * ry**2 / magr**3)
        + h_vec[1] * (vx * vy - MUM * rx * ry / magr**3)
        + vz * MUM * ecc_vec[1]
    )
    tm_argp[1] = (
        tm_argp[1] / (MUM * node * ecc)
        - vz * h_vec[0] * n_dot_e / (node**3 * ecc)
        - tm_ecc[1] * n_dot_e / (node * ecc**2)
    ) * w_scale

    tm_argp[2] = (
        -h_vec[1] * (vx * vz - MUM * rx * rz / magr**3)
        + h_vec[0] * (vy * vz - MUM * ry * rz / magr**3)
        + vx * MUM * ecc_vec[0]
        + vy * MUM * ecc_vec[1]
    )
    tm_argp[2] = (
        -tm_argp[2] / (MUM * node * ecc)
        + (vy * h_vec[0] - vx * h_vec[1]) * n_dot_e / (node**3 * ecc)
        - tm_ecc[2] * n_dot_e / (node * ecc**2)
    ) * w_scale

    tm_argp[3] = (
        (rx * vy - 2 * ry * vx) * h_vec[0]
        - h_vec[1] * (ry * vy + rz * vz)
        + rz * MUM * ecc_vec[0]
    )
    tm_argp[3] = (
        -tm_argp[3] / (MUM * node * ecc)
        - rz * h_vec[1] * n_dot_e / (node**3 * ecc)
        - tm_ecc[3] * n_dot_e / (node * ecc**2)
    ) * w_scale

    tm_argp[4] = (
        -(ry * vx - 2 * rx * vy) * h_vec[1]
        + h_vec[0] * (rx * vx + rz * vz)
        + rz * MUM * ecc_vec[1]
    )
    tm_argp[4] = (
        -tm_argp[4] / (MUM * node * ecc)
        + rz * h_vec[0] * n_dot_e / (node**3 * ecc)
        - tm_ecc[4] * n_dot_e / (node * ecc**2)
    ) * w_scale

    tm_argp[5] = (
        -(rz * vx - 2 * rx * vz) * h_vec[1]
        + h_vec[0] * (rz * vy - 2 * ry * vz)
        - rx * MUM * ecc_vec[0]
        - ry * MUM * ecc_vec[1]
    )
    tm_argp[5] = (
        -tm_argp[5] / (MUM * node * ecc)
        + (rx * h_vec[1] - ry * h_vec[0]) * n_dot_e / (node**3 * ecc)
        - tm_ecc[5] * n_dot_e / (node * ecc**2)
    ) * w_scale

    return tm_argp


def _compute_partials_nu(
    tm_ecc, reci_m, veci_m, ecc, ecc_term, r_dot_v, r_dot_e, nu_scale
):
    """Compute partial derivatives of nu/M w.r.t. (rx, ry, rz, vx, vy, vz)."""
    rx, ry, rz = reci_m
    vx, vy, vz = veci_m
    magr = np.linalg.norm(reci_m)
    tm_nu = np.zeros(6)
    tm_nu[0] = (
        ry * (vx * vy - MUM * rx * ry / magr**3)
        - rx * ecc_term
        + rz * (vx * vz - MUM * rx * rz / magr**3)
        - rx * (vy**2 + vz**2 - MUM / magr + MUM * rx**2 / magr**3)
        + vx * r_dot_v
    )
    tm_nu[0] = (
        -tm_nu[0] / (MUM * magr * ecc)
        - rx * r_dot_e / (magr**3 * ecc)
        - tm_ecc[0] * r_dot_e / (magr * ecc**2)
    ) * nu_scale

    tm_nu[1] = (
        rx * (vx * vy - MUM * rx * ry / magr**3)
        - ry * ecc_term
        + rz * (vy * vz - MUM * ry * rz / magr**3)
        - ry * (vx**2 + vz**2 - MUM / magr + MUM * ry**2 / magr**3)
        + vy * r_dot_v
    )
    tm_nu[1] = (
        -tm_nu[1] / (MUM * magr * ecc)
        - ry * r_dot_e / (magr**3 * ecc)
        - tm_ecc[1] * r_dot_e / (magr * ecc**2)
    ) * nu_scale

    tm_nu[2] = (
        rx * (vx * vz - MUM * rx * rz / magr**3)
        - rz * ecc_term
        + ry * (vy * vz - MUM * ry * rz / magr**3)
        - rz * (vx**2 + vy**2 - MUM / magr + MUM * rz**2 / magr**3)
        + vz * r_dot_v
    )
    tm_nu[2] = (
        -tm_nu[2] / (MUM * magr * ecc)
        - rz * r_dot_e / (magr**3 * ecc)
        - tm_ecc[2] * r_dot_e / (magr * ecc**2)
    ) * nu_scale

    tm_nu[3] = (
        ry * (rx * vy - 2 * ry * vx)
        + rx * (ry * vy + rz * vz)
        + rz * (rx * vz - 2 * rz * vx)
    )
    tm_nu[3] = (
        -tm_nu[3] / (MUM * magr * ecc) - tm_ecc[3] * r_dot_e / (magr * ecc**2)
    ) * nu_scale

    tm_nu[4] = (
        rx * (ry * vx - 2 * rx * vy)
        + ry * (rx * vx + rz * vz)
        + rz * (ry * vz - 2 * rz * vy)
    )
    tm_nu[4] = (
        -tm_nu[4] / (MUM * magr * ecc) - tm_ecc[4] * r_dot_e / (magr * ecc**2)
    ) * nu_scale

    tm_nu[5] = (
        rz * (rx * vx + ry * vy)
        + rx * (rz * vx - 2 * rx * vz)
        + ry * (rz * vy - 2 * ry * vz)
    )
    tm_nu[5] = (
        -tm_nu[5] / (MUM * magr * ecc) - tm_ecc[5] * r_dot_e / (magr * ecc**2)
    ) * nu_scale

    return tm_nu


def covct2cl(
    cartcov: ArrayLike, cartstate: ArrayLike, use_mean_anom: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from Cartesian elements to classical orbital
    elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        cartcov (array_like): 6x6 Cartesian covariance matrix in m and m/s
        cartstate (array_like): 6x1 Cartesian orbit state in km and km/s
                                (rx, ry, rz, vx, vy, vz)
        use_mean_anom (bool): Flag to use mean anomaly instead of true anomaly
                              (defaults to False)

    Returns:
        tuple: (classcov, tm)
            classcov (np.ndarray): 6x6 Classical orbital elements covariance matrix
                                   in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Parse the input state vector
    reci_m = np.array(cartstate[:3]) * KM2M
    veci_m = np.array(cartstate[3:]) * KM2M

    # Convert to classical orbital elements
    _, a, ecc, incl, omega, argp, nu, *_ = fc.rv2coe(reci_m / KM2M, veci_m / KM2M)
    a *= KM2M
    n = np.sqrt(MUM / a**3)

    # Common quantities
    sqrt1me2 = np.sqrt(1 - ecc**2)
    magr = np.linalg.norm(reci_m)
    magv = np.linalg.norm(veci_m)

    # Eccentricity vector
    r_dot_v = np.dot(reci_m, veci_m)
    ecc_term = magv**2 - MUM / magr
    ecc_vec = (ecc_term * reci_m - r_dot_v * veci_m) / MUM

    # Node vector
    h_vec = np.cross(reci_m, veci_m)
    node_vec = np.cross([0, 0, 1], h_vec)
    node = np.linalg.norm(node_vec)

    # Additional terms for argument of periapsis and true anomaly
    n_dot_e = np.dot(node_vec, ecc_vec)
    sign_w = np.sign((magv**2 - MUM / magr) * reci_m[2] - r_dot_v * veci_m[2])
    cos_w = n_dot_e / (ecc * node)
    w_scale = -sign_w / np.sqrt(1 - cos_w**2)

    # Additional terms for true anomaly
    r_dot_e = np.dot(reci_m, ecc_vec)
    cos_nu = r_dot_e / (magr * ecc)
    sign_nu = np.sign(r_dot_v)
    nu_scale = -sign_nu / np.sqrt(1 - cos_nu**2)

    # Compute partial derivatives
    tm = np.zeros((6, 6))

    # Compute partials
    tm[0, :] = _compute_partials_a(a, n, reci_m, veci_m)
    tm[1, :] = _compute_partials_ecc(reci_m, veci_m, ecc_vec, ecc)
    tm[2, :] = _compute_partials_incl(reci_m, veci_m, h_vec, node)
    tm[3, :] = _compute_partials_node(reci_m, veci_m, node_vec, node)
    tm[4, :] = _compute_partials_argp(
        tm[1, :], reci_m, veci_m, ecc_vec, ecc, h_vec, node, n_dot_e, w_scale
    )
    tm[5, :] = _compute_partials_nu(
        tm[1, :], reci_m, veci_m, ecc, ecc_term, r_dot_v, r_dot_e, nu_scale
    )

    # Update partials for mean anomaly, if specified
    if use_mean_anom:
        dmdnu = (sqrt1me2**2) ** 1.5 / (1 + ecc * np.cos(nu)) ** 2
        dmde = (
            -np.sin(nu)
            * (
                (ecc * np.cos(nu) + 1)
                * (ecc + np.cos(nu))
                / np.sqrt((ecc + np.cos(nu)) ** 2)
                + 1
                - 2 * ecc**2
                - ecc**3 * np.cos(nu)
            )
            / ((ecc * np.cos(nu) + 1) ** 2 * sqrt1me2)
        )
        tm[5, :] = tm[5, :] * dmdnu + tm[1, :] * dmde

    # Calculate the classical covariance matrix
    classcov = tm @ cartcov @ tm.T

    return classcov, tm


def covcl2ct(
    classcov: ArrayLike, classstate: ArrayLike, use_mean_anom: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from classical elements to Cartesian elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        classcov (array_like): 6x6 Classical orbital elements covariance matrix
                               in m and m/s
        classstate (array_like): 6x1 Classical orbital elements in km and radians
                                (a, ecc, incl, node, argp, nu or m)
        use_mean_anom (bool): Flag to use mean anomaly instead of true anomaly
                              (defaults to False)

    Returns:
        tuple: (cartcov, tm)
            cartcov (np.ndarray): 6x6 Cartesian covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """

    def update_row_with_mean_anomaly(mat, row):
        if use_mean_anom:
            mat[row, 5] /= dmdnu
            mat[row, 1] -= mat[row, 5] * dmde
        return mat

    # Parse the classical elements
    a, ecc, incl, raan, argp, anom = classstate
    a *= KM2M

    # Convert anomaly as needed
    if use_mean_anom:
        e, nu = fc.newtonm(ecc, anom)
    else:
        nu = anom
        e, _ = fc.newtonnu(ecc, nu)

    # Compute trigonometric values
    sin_inc, cos_inc = np.sin(incl), np.cos(incl)
    sin_raan, cos_raan = np.sin(raan), np.cos(raan)
    sin_w, cos_w = np.sin(argp), np.cos(argp)
    sin_nu, cos_nu = np.sin(nu), np.cos(nu)

    # Define PQW to ECI transformation elements (p. 168)
    p11 = cos_raan * cos_w - sin_raan * sin_w * cos_inc
    p12 = -cos_raan * sin_w - sin_raan * cos_w * cos_inc
    p13 = sin_raan * sin_inc
    p21 = sin_raan * cos_w + cos_raan * sin_w * cos_inc
    p22 = -sin_raan * sin_w + cos_raan * cos_w * cos_inc
    p23 = -cos_raan * sin_inc
    p31 = sin_w * sin_inc
    p32 = cos_w * sin_inc

    # Define constants for efficiency
    p0 = np.sqrt(MUM / (a * (1 - ecc**2)))
    p1 = (1 - ecc**2) / (1 + ecc * cos_nu)
    p2 = 1 / (2 * a) * p0
    p3 = (2 * a * ecc + a * cos_nu + a * cos_nu * ecc**2) / ((1 + ecc * cos_nu) ** 2)
    p4 = ecc * MUM / (a * (1 - ecc**2) ** 2 * p0)
    p5 = a * p1
    p6 = a * (1 - ecc**2) / ((1 + ecc * cos_nu) ** 2)

    dmdnu = (1 - ecc**2) ** 1.5 / (1 + ecc * cos_nu) ** 2
    dmde = (
        -sin_nu
        * (
            (ecc * cos_nu + 1) * (ecc + cos_nu) / np.sqrt((ecc + cos_nu) ** 2)
            + 1
            - 2 * ecc**2
            - ecc**3 * cos_nu
        )
        / ((ecc * cos_nu + 1) ** 2 * np.sqrt(1 - ecc**2))
    )

    # Compute partial derivatives
    tm = np.zeros((6, 6))

    # Partials of (a, ecc, incl, node, argp, nu) w.r.t. rx
    tm[0, 0] = p1 * (p11 * cos_nu + p12 * sin_nu)
    tm[0, 1] = -p3 * (p11 * cos_nu + p12 * sin_nu)
    tm[0, 2] = p5 * p13 * (sin_w * cos_nu + cos_w * sin_nu)
    tm[0, 3] = -p5 * (p21 * cos_nu + p22 * sin_nu)
    tm[0, 4] = p5 * (p12 * cos_nu - p11 * sin_nu)
    tm[0, 5] = p6 * (-p11 * sin_nu + p12 * (ecc + cos_nu))
    tm = update_row_with_mean_anomaly(tm, 0)

    # Partials of (a, ecc, incl, node, argp, nu) w.r.t. ry
    tm[1, 0] = p1 * (p21 * cos_nu + p22 * sin_nu)
    tm[1, 1] = -p3 * (p21 * cos_nu + p22 * sin_nu)
    tm[1, 2] = p5 * p23 * (sin_w * cos_nu + cos_w * sin_nu)
    tm[1, 3] = p5 * (p11 * cos_nu + p12 * sin_nu)
    tm[1, 4] = p5 * (p22 * cos_nu - p21 * sin_nu)
    tm[1, 5] = p6 * (-p21 * sin_nu + p22 * (ecc + cos_nu))
    tm = update_row_with_mean_anomaly(tm, 1)

    # Partials of (a, ecc, incl, node, argp, nu) w.r.t. rz
    tm[2, 0] = p1 * (p31 * cos_nu + p32 * sin_nu)
    tm[2, 1] = -p3 * sin_inc * (cos_w * sin_nu + sin_w * cos_nu)
    tm[2, 2] = p5 * cos_inc * (cos_w * sin_nu + sin_w * cos_nu)
    tm[2, 3] = 0.0
    tm[2, 4] = p5 * sin_inc * (cos_w * cos_nu - sin_w * sin_nu)
    tm[2, 5] = p6 * (-p31 * sin_nu + p32 * (ecc + cos_nu))
    tm = update_row_with_mean_anomaly(tm, 2)

    # Partials of (a, ecc, incl, node, argp, nu) w.r.t. vx
    tm[3, 0] = p2 * (p11 * sin_nu - p12 * (ecc + cos_nu))
    tm[3, 1] = -p4 * (p11 * sin_nu - p12 * (ecc + cos_nu)) + p12 * p0
    tm[3, 2] = -p0 * sin_raan * (p31 * sin_nu - p32 * (ecc + cos_nu))
    tm[3, 3] = p0 * (p21 * sin_nu - p22 * (ecc + cos_nu))
    tm[3, 4] = -p0 * (p12 * sin_nu + p11 * (ecc + cos_nu))
    tm[3, 5] = -p0 * (p11 * cos_nu + p12 * sin_nu)
    tm = update_row_with_mean_anomaly(tm, 3)

    # Partials of (a, ecc, incl, node, argp, nu) w.r.t. vy
    tm[4, 0] = p2 * (p21 * sin_nu - p22 * (ecc + cos_nu))
    tm[4, 1] = -p4 * (p21 * sin_nu - p22 * (ecc + cos_nu)) + p22 * p0
    tm[4, 2] = p0 * cos_raan * (p31 * sin_nu - p32 * (ecc + cos_nu))
    tm[4, 3] = p0 * (-p11 * sin_nu + p12 * (ecc + cos_nu))
    tm[4, 4] = -p0 * (p22 * sin_nu + p21 * (ecc + cos_nu))
    tm[4, 5] = -p0 * (p21 * cos_nu + p22 * sin_nu)
    tm = update_row_with_mean_anomaly(tm, 4)

    # Partials of (a, ecc, incl, node, argp, nu) w.r.t. vz
    tm[5, 0] = p2 * (p31 * sin_nu - p32 * (ecc + cos_nu))
    tm[5, 1] = -p4 * (p31 * sin_nu - p32 * (ecc + cos_nu)) + p32 * p0
    tm[5, 2] = p0 * cos_inc * (cos_w * cos_nu - sin_w * sin_nu + ecc * cos_w)
    tm[5, 3] = 0.0
    tm[5, 4] = -p0 * (p32 * sin_nu + p31 * (ecc + cos_nu))
    tm[5, 5] = -p0 * (p31 * cos_nu + p32 * sin_nu)
    tm = update_row_with_mean_anomaly(tm, 5)

    # Calculate the Cartesian covariance matrix
    cartcov = tm @ classcov @ tm.T

    return cartcov, tm


########################################################################################
# Equinoctial <-> Cartesian Elements
########################################################################################


def newton_mean_anomaly(
    meanlon_m: float, af: float, ag: float, tol: float = 1e-8, max_iter: int = 25
) -> float:
    """Solves for the equinoctial anomaly (F) using Newton's method.

    Args:
        meanlon_m (float): Mean longitude (mean anomaly + mean argument of latitude)
        af (float): Component of the eccentricity vector in equinoctial elements
        ag (float): Component of the eccentricity vector in equinoctial elements
        tol (float, optional): Convergence tolerance (defaults to 1e-8)
        max_iter (int, optional): Maximum number of iterations (default is 25)

    Returns:
        float: The solved equinoctial anomaly (F)
    """
    # Initial guess
    f0 = meanlon_m
    for ktr in range(max_iter):
        # Compute the function value and its derivative
        f1 = f0 - (f0 + ag * np.cos(f0) - af * np.sin(f0) - meanlon_m) / (
            1 - ag * np.sin(f0) - af * np.cos(f0)
        )

        # Check for convergence
        if abs(f1 - f0) < tol:
            return f1

        # Update the guess
        f0 = f1

    # If the loop completes without convergence, raise an error
    raise RuntimeError(f"Newton method failed to converge after {max_iter} iterations.")


def covct2eq(
    cartcov: ArrayLike, cartstate: ArrayLike, fr: int, use_mean_anom: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from Cartesian to equinoctial elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        cartcov (array_like): 6x6 Cartesian covariance matrix in m and m/s
        cartstate (array_like): 6x1 Cartesian state vector in km and km/s
                                (rx, ry, rz, vx, vy, vz)
        fr (int): Retrograde factor (+1 or -1)
        use_mean_anom (bool): Flag to use mean anomaly instead of true anomaly
                              (defaults to False)

    Returns:
        tuple: (eqcov, tm)
            eqcov (np.ndarray): 6x6 Equinoctial covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Parse input vectors
    reci_m = np.array(cartstate[:3]) * KM2M
    veci_m = np.array(cartstate[3:]) * KM2M
    magr = np.linalg.norm(reci_m)
    magv = np.linalg.norm(veci_m)

    # Classical orbital elements
    a = 1 / (2 / magr - magv**2 / MUM)
    n = np.sqrt(MUM / a**3)
    h_vec = np.cross(reci_m, veci_m)
    w_vec = h_vec / np.linalg.norm(h_vec)
    chi = w_vec[0] / (1 + fr * w_vec[2])
    psi = -w_vec[1] / (1 + fr * w_vec[2])

    # Equinoctial components
    p0 = 1 / (1 + chi**2 + psi**2)
    f_vec = p0 * np.array([1 - chi**2 + psi**2, 2 * chi * psi, -2 * fr * chi])
    g_vec = p0 * np.array([2 * fr * chi * psi, fr * (1 + chi**2 - psi**2), 2 * psi])

    # Compute eccentricity vector
    r_dot_v = np.dot(reci_m, veci_m)
    p1 = magv**2 - MUM / magr
    ecc_vec = (p1 * reci_m - r_dot_v * veci_m) / MUM

    # Get eccentricity components
    af = np.dot(ecc_vec, f_vec)
    ag = np.dot(ecc_vec, g_vec)

    # Intermediate terms
    x = np.dot(reci_m, f_vec)
    y = np.dot(reci_m, g_vec)

    b = 1 / (1 + np.sqrt(1 - af**2 - ag**2))
    p0 = 1 / (a * np.sqrt(1 - af**2 - ag**2))
    sinf = ag + p0 * ((1 - ag**2 * b) * y - ag * af * b * x)
    cosf = af + p0 * ((1 - af**2 * b) * x - ag * af * b * y)
    f = np.arctan2(sinf, cosf) % TWOPI

    xd = n * a**2 / magr * (af * ag * b * np.cos(f) - (1 - ag**2 * b) * np.sin(f))
    yd = n * a**2 / magr * ((1 - af**2 * b) * np.cos(f) - af * ag * b * np.sin(f))

    # Compute partial derivatives
    a_ = np.sqrt(MUM * a)
    b_ = np.sqrt(1 - ag**2 - af**2)
    c_ = 1 + chi**2 + psi**2

    a_term = a_ / magr**3
    term1_b = a / (1 + b_)
    term2_b = 1 / b_

    partxdaf = a * xd * yd / (a_ * b_) - a_term * (term1_b * ag * x + term2_b * x * y)
    partydaf = -a * xd**2 / (a_ * b_) - a_term * (term1_b * ag * y - term2_b * x**2)
    partxdag = a * yd**2 / (a_ * b_) + a_term * (term1_b * af * x - term2_b * y**2)
    partydag = -a * xd * yd / (a_ * b_) + a_term * (term1_b * af * y + term2_b * x * y)

    # Initialize transformation matrix
    tm = np.zeros((6, 6))

    # Partials of sma w.r.t. (rx ry rz vx vy vz)
    if use_mean_anom:
        p0 = 2 * a**2 / magr**3
        p1 = 2 / (n**2 * a)
    else:
        p0 = -3 * n * a / magr**3
        p1 = -3 / (n * a**2)

    tm[0, :3] = p0 * reci_m
    tm[0, 3:] = p1 * veci_m

    # Partials of v w.r.t. ag
    tm_v_ag_parts = partxdag * f_vec + partydag * g_vec

    # Partials of af w.r.t. (rx ry rz vx vy vz)
    p0 = 1 / MUM
    af_term = -a * b * af * b_ / (magr**3)
    chi_psi_term = ag * (chi * xd - psi * fr * yd) / (a_ * b_)
    tm[1, :3] = af_term * reci_m - chi_psi_term * w_vec + (b_ / a_) * tm_v_ag_parts

    xy_term = p0 * ((2 * x * yd - xd * y) * g_vec - y * yd * f_vec)
    psi_chi_term = ag * (psi * fr * y - chi * x) / (a_ * b_)
    tm[1, 3:] = xy_term - psi_chi_term * w_vec

    # Partials of v w.r.t. af
    tm_v_af_parts = partxdaf * f_vec + partydaf * g_vec

    # Partials of af w.r.t. (rx ry rz vx vy vz)
    ag_term = -a * b * ag * b_ / (magr**3)
    chi_psi_term_af = af * (chi * xd - psi * fr * yd) / (a_ * b_)
    tm[2, :3] = ag_term * reci_m + chi_psi_term_af * w_vec - (b_ / a_) * tm_v_af_parts

    xy_term = p0 * ((2 * xd * y - x * yd) * f_vec - x * xd * g_vec)
    psi_chi_term_af = af * (psi * fr * y - chi * x) / (a_ * b_)
    tm[2, 3:] = xy_term + psi_chi_term_af * w_vec

    # Partials of chi w.r.t. (rx ry rz vx vy vz)
    den = 2 * a_ * b_
    tm[3, :3] = -c_ * yd * w_vec / den
    tm[3, 3:] = c_ * y * w_vec / den

    # Partials of psi w.r.t. (rx ry rz vx vy vz)
    tm[4, :3] = -fr * c_ * xd * w_vec / den
    tm[4, 3:] = fr * c_ * x * w_vec / den

    # Partials of true/mean anomaly w.r.t. (rx ry rz vx vy vz)
    tm[5, :] = 0
    if use_mean_anom:
        tm[5, :3] = (
            -veci_m / a_
            + (chi * xd - psi * fr * yd) * w_vec / (a_ * b_)
            - (b * b_ / a_) * (ag * tm_v_ag_parts + af * tm_v_af_parts)
        )
        tm[5, 3:] = (
            -2.0 * reci_m / a_
            + (af * tm[2, 3:6] - ag * tm[1, 3:6]) / (1 + b_)
            + (fr * psi * y - chi * x) * w_vec / a_
        )

    # Compute the equinoctial covariance matrix
    eqcov = tm @ cartcov @ tm.T

    return eqcov, tm


def coveq2ct(
    eqcov: ArrayLike, eqstate: ArrayLike, fr: int, anom_type: AnomalyType
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from equinoctial to cartesian elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        eqcov (array_like): 6x6 equinoctial covariance matrix in m and m/s
        eqstate (array_like): 6x1 equinoctial orbit state in km and km/s
                              (a/n, af, ag, chi, psi, lm/ln)
        fr (int): Retrograde factor (+1 or -1)
        anom_type (AnomalyType): Anomaly type (MEAN_A, TRUE_A, MEAN_N, TRUE_N)

    Returns:
        tuple: (cartcov, tm)
            cartcov (np.ndarray): 6x6 Cartesian covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix

    TODO:
        - Return cartesian covariance matrix not as expected?
    """
    # Parse eqstate and anomaly-dependent calculations
    use_anom_a = anom_type in {AnomalyType.TRUE_A, AnomalyType.MEAN_A}
    if use_anom_a:
        a = eqstate[0] * KM2M  # in meters
        n = np.sqrt(MUM / a**3)
    else:
        n = eqstate[0]  # rad/s
        a = (MUM / n**2) ** (1 / 3)  # in meters

    # Get orbital elements
    af, ag, chi, psi = eqstate[1:5]
    use_mean_anom = anom_type in {AnomalyType.MEAN_A, AnomalyType.MEAN_N}
    meanlon_m = eqstate[5] if use_mean_anom else None
    meanlon_nu = eqstate[5] if not use_mean_anom else None
    omega = np.arctan2(chi, psi)
    argp = np.arctan2(ag, af) - fr * omega
    ecc = np.sqrt(af**2 + ag**2)

    # Update for true anomaly
    if not use_mean_anom:
        nu = np.mod(meanlon_nu - fr * omega - argp, TWOPI)
        _, m = newtonnu(ecc, nu)
        meanlon_m = np.mod(fr * omega + argp + m, TWOPI)

    # Convert equinoctial to cartesian state
    reci, veci = fc.eq2rv(a / KM2M, af, ag, chi, psi, meanlon_m, fr)
    reci_m = reci * KM2M
    veci_m = veci * KM2M
    magr = np.linalg.norm(reci_m)

    # Constants for transformation
    a_ = n * a**2
    b_ = np.sqrt(1 - ag**2 - af**2)
    c_ = 1 + chi**2 + psi**2
    b = 1 / (1 + b_)
    g_ = a_ * b_

    # Compute partial derivatives
    f_ = newton_mean_anomaly(meanlon_m, af, ag)
    x = a * ((1 - ag**2 * b) * np.cos(f_) + af * ag * b * np.sin(f_) - af)
    y = a * ((1 - af**2 * b) * np.sin(f_) + af * ag * b * np.cos(f_) - ag)
    xd = n * a**2 / magr * (af * ag * b * np.cos(f_) - (1 - ag**2 * b) * np.sin(f_))
    yd = n * a**2 / magr * ((1 - af**2 * b) * np.cos(f_) - af * ag * b * np.sin(f_))

    # Equinoctial system components
    p0 = 1 / (1 + chi**2 + psi**2)
    f_vec = p0 * np.array([1 - chi**2 + psi**2, 2 * chi * psi, -2 * fr * chi])
    g_vec = p0 * np.array([2 * fr * chi * psi, fr * (1 + chi**2 - psi**2), 2 * psi])
    w_vec = p0 * np.array([2 * chi, -2 * psi, fr * (1 - chi**2 - psi**2)])

    # Partial derivatives wrt af and ag
    partxaf = ag * b * xd / n + a * y * xd / g_ - a
    partyaf = ag * b * yd / n - a * x * xd / g_
    partxag = -af * b * xd / n + a * y * yd / g_
    partyag = -af * b * yd / n - a * x * yd / g_ - a

    partxdaf = a * xd * yd / g_ - a_ / magr**3 * (a * ag * x / (1 + b_) + x * y / b_)
    partydaf = -a * xd**2 / g_ - a_ / magr**3 * (a * ag * y / (1 + b_) - x**2 / b_)
    partxdag = a * yd**2 / g_ + a_ / magr**3 * (a * af * x / (1 + b_) - y**2 / b_)
    partydag = -a * xd * yd / g_ + a_ / magr**3 * (a * af * y / (1 + b_) + x * y / b_)

    # Initialize transformation matrix
    tm = np.zeros((6, 6))

    # Partials of (rx ry rz vx vy vz) w.r.t. sma
    if use_anom_a:
        p0 = 1 / a
        p1 = -1 / (2 * a)
    else:
        p0 = -2 / (3 * n)
        p1 = 1 / (3 * n)

    tm[:3, 0] = p0 * reci_m
    tm[3:, 0] = p1 * veci_m

    # Partials of (rx ry rz vx vy vz) w.r.t. af and ag
    tm[:3, 1] = partxaf * f_vec + partyaf * g_vec
    tm[3:, 1] = partxdaf * f_vec + partydaf * g_vec
    tm[:3, 2] = partxag * f_vec + partyag * g_vec
    tm[3:, 2] = partxdag * f_vec + partydag * g_vec

    # Partials of (rx ry rz vx vy vz) w.r.t. chi and psi
    p0 = 2 * fr / c_
    tm[:3, 3] = p0 * (psi * (y * f_vec - x * g_vec) - x * w_vec)
    tm[3:, 3] = p0 * (psi * (yd * f_vec - xd * g_vec) - xd * w_vec)
    tm[:3, 4] = p0 * (chi * (x * g_vec - y * f_vec) + y * w_vec)
    tm[3:, 4] = p0 * (chi * (xd * g_vec - yd * f_vec) + yd * w_vec)

    # Partials of (rx ry rz vx vy vz) w.r.t. mean longitude
    p0 = 1 / n
    p1 = -n * a**3 / magr**3
    tm[:3, 5] = p0 * veci_m
    tm[3:, 5] = p1 * reci_m

    # Final Cartesian covariance
    cartcov = tm @ eqcov @ tm.T

    return cartcov, tm


########################################################################################
# Equinoctial <-> Classical Elements
########################################################################################


def covcl2eq(
    classcov: ArrayLike, classstate: ArrayLike, fr: int, anom: AnomalyType
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from classical to equinoctial elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        classcov (array_like): 6x6 Classical orbital elements covariance matrix
                               in m and m/s
        classstate (array_like): 6x1 Classical orbital elements in km and radians
                                (a, ecc, incl, node, argp, nu or m)
        fr (int): Retrograde factor (+1 or -1)
        anom (AnomalyType): Anomaly type (MEAN_A, TRUE_A, MEAN_N, TRUE_N)

    Returns:
        tuple: (eqcov, tm)
            eqcov (np.ndarray): 6x6 Equinoctial covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Parse the orbit state
    a, ecc, incl, omega, argp = classstate[:5]
    a *= KM2M

    # Initialize transformation matrix
    tm = np.zeros((6, 6))

    # Partials of a/n wrt (a, ecc, incl, node, argp, nu/M)
    if anom in {AnomalyType.TRUE_A, AnomalyType.MEAN_A}:
        tm[0, 0] = 1
    else:
        tm[0, 0] = -(3 * np.sqrt(MUM / a**3)) / (2 * a)

    # Partials of af wrt (a, ecc, incl, node, argp, nu/M)
    tm[1, 1] = np.cos(fr * omega + argp)
    tm[1, 3] = -ecc * fr * np.sin(fr * omega + argp)
    tm[1, 4] = -ecc * np.sin(fr * omega + argp)

    # Partials of ag wrt (a, ecc, incl, node, argp, nu/M)
    tm[2, 1] = np.sin(fr * omega + argp)
    tm[2, 3] = ecc * fr * np.cos(fr * omega + argp)
    tm[2, 4] = ecc * np.cos(fr * omega + argp)

    # Partials of chi wrt (a, ecc, incl, node, argp, nu/M)
    mult = (0.5 * np.tan(incl * 0.5) ** 2 + 0.5) * fr * np.tan(incl * 0.5) ** (fr - 1)
    tm[3, 2] = np.sin(omega) * mult
    tm[3, 3] = np.tan(incl * 0.5) ** fr * np.cos(omega)

    # Partials of psi wrt (a, ecc, incl, node, argp, nu/M)
    tm[4, 2] = np.cos(omega) * mult
    tm[4, 3] = -np.tan(incl * 0.5) ** fr * np.sin(omega)

    # Partials of meanlonM/meanlonNu wrt (a, ecc, incl, node, argp, nu/M)
    tm[5, 3] = fr
    tm[5, 4] = tm[5, 5] = 1

    # Calculate the output covariance matrix
    eqcov = tm @ classcov @ tm.T

    return eqcov, tm


def coveq2cl(
    eqcov: ArrayLike, eqstate: ArrayLike, fr: int, anom: AnomalyType
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from equinoctial to classical elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        eqcov (array_like): 6x6 Equinoctial covariance matrix in m and m/s
        eqstate (array_like): 6x1 Equinoctial orbital state in km and radians
                              (a/n, af, ag, chi, psi, lm/ln)
        fr (int): Retrograde factor (+1 or -1)
        anom (AnomalyType): Anomaly type (MEAN_A, TRUE_A, MEAN_N, TRUE_N)

    Returns:
        tuple: (classcov, tm)
            classcov (np.ndarray): 6x6 Classical covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Parse the equinoctial state
    if anom in {AnomalyType.TRUE_A, AnomalyType.MEAN_A}:
        a = eqstate[0]
        n = np.sqrt(MUM / a**3)
    else:
        n = eqstate[0]

    # Get the equinoctial elements
    af, ag, chi, psi = eqstate[1:5]

    # Initialize transformation matrix
    tm = np.zeros((6, 6))

    # Partials for semi-major axis or mean motion
    if anom in {AnomalyType.TRUE_A, AnomalyType.MEAN_A}:
        tm[0, 0] = 1
    else:
        tm[0, 0] = -2 / (3 * n) * (MUM / n**2) ** (1 / 3)

    # Partials for eccentricity
    p0 = 1 / np.sqrt(af**2 + ag**2)
    tm[1, 1] = p0 * af
    tm[1, 2] = p0 * ag

    # Partials for inclination
    p1 = 2 * fr / ((1 + chi**2 + psi**2) * np.sqrt(chi**2 + psi**2))
    tm[2, 3] = p1 * chi
    tm[2, 4] = p1 * psi

    # Partials for RAAN (node)
    p2 = 1 / (chi**2 + psi**2)
    tm[3, 3] = p2 * psi
    tm[3, 4] = -p2 * chi

    # Partials for argument of perigee
    p3 = 1 / (af**2 + ag**2)
    tm[4, 1] = -p3 * ag
    tm[4, 2] = p3 * af
    tm[4, 3] = -fr * p2 * psi
    tm[4, 4] = fr * p2 * chi

    # Partials for anomaly
    if anom in {AnomalyType.TRUE_A, AnomalyType.TRUE_N}:
        p4 = 1 / (af**2 + ag**2)
        tm[5, 1] = p4 * ag
        tm[5, 2] = -p4 * af
    else:
        tm[5, 1] = p3 * ag
        tm[5, 2] = -p3 * af

    tm[5, 5] = 1

    # Compute the transformed covariance matrix
    classcov = tm @ eqcov @ tm.T

    return classcov, tm


########################################################################################
# Cartesian <-> Flight Parameters
########################################################################################


def _compute_partials_az(reci_m, veci_m):
    """Compute partial derivatives of azimuth w.r.t. (rx, ry, rz, vx, vy, vz)."""
    rx, ry, rz = reci_m
    vx, vy, vz = veci_m
    magr = np.linalg.norm(reci_m)
    magv = np.linalg.norm(veci_m)
    rdotv = np.dot(reci_m, veci_m)

    # Sal from mathcad methoc
    p2 = 1 / ((magv**2 - (rdotv / magr) ** 2) * (rx**2 + ry**2))
    k1 = np.linalg.norm(reci_m) * np.cross(reci_m[:2], veci_m[:2])
    k2 = ry * (ry * vz - rz * vy) + rx * (rx * vz - rz * vx)
    k12_sq = k1**2 + k2**2

    # Construct the transformation matrix
    tm_az = np.zeros(6)
    tm_az[0] = p2 * (
        vy * (magr * vz - rz * rdotv / magr)
        - (rx * vy - ry * vx) / magr * (rx * vz - rz * vx + rx * rdotv / magr)
    )
    p2 = 1 / (magr * k12_sq)
    tm_az[0] = p2 * (
        k1 * magr * (rz * vx - 2 * rx * vz)
        + k2 * (-ry * vx * rx + vy * rx**2 + vy * magr**2)
    )
    tm_az[1] = p2 * (
        k1 * magr * (rz * vy - 2 * ry * vz)
        + k2 * (rx * vy * ry - vx * ry**2 - vx * magr**2)
    )
    p2 = k1 / (magr**2 * k12_sq)
    tm_az[2] = p2 * (k2 * rz + (rx * vx + ry * vy) * magr**2)
    p2 = 1 / k12_sq
    tm_az[3] = p2 * (k1 * rx * rz - k2 * ry * magr)
    tm_az[4] = p2 * (k1 * ry * rz + k2 * rx * magr)
    tm_az[5] = -p2 * (k1 * (rx**2 + ry**2))

    return tm_az


def covct2fl(
    cartcov: ArrayLike,
    cartstate: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
    use_latlon: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from Cartesian to flight parameters.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        cartcov (array_like): 6x6 Cartesian covariance matrix in m and m/s
        cartstate (array_like): 6x1 Cartesian state vector in km and km/s
                                (rx, ry, rz, vx, vy, vz)
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 nutation parameters
        eqeterms (bool, optional): Add terms for ast calculation (default True)
        use_latlon (bool, optional): Flag to use lat/lon instead of ra/dec
                                     (default True)

    Returns:
        tuple: (flcov, tm)
            flcov (np.ndarray): 6x6 Flight parameters covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Parse input state into cartesian components (in meters and m/s)
    reci_m = np.array(cartstate[:3]) * KM2M
    veci_m = np.array(cartstate[3:]) * KM2M

    # Convert to ECEF coordinates if using lat/lon
    r = reci_m
    if use_latlon:
        aeci = np.array([0, 0, 0])
        recef, *_ = fc.eci2ecef(
            reci_m / KM2M,
            veci_m / KM2M,
            aeci,
            ttt,
            jdut1,
            lod,
            xp,
            yp,
            ddpsi,
            ddeps,
            iau80arr,
            eqeterms,
        )
        r = recef * KM2M

    # Calculate common quantities
    magr = np.linalg.norm(reci_m)
    magv = np.linalg.norm(veci_m)
    h = np.linalg.norm(np.cross(reci_m, veci_m))

    # Initialize transformation matrix
    tm = np.zeros((6, 6))

    # Transformation matrix components for latlon or radec
    r_xy2 = np.sum(r[:2] ** 2)
    r_xy = np.sqrt(r_xy2)

    p0 = 1 / r_xy2
    tm[0, 0:2] = [-p0 * r[1], p0 * r[0]]

    p0 = 1 / (magr**2 * r_xy)
    tm[1, 0:2] = -p0 * r[:2] * r[2]
    tm[1, 2] = r_xy / magr**2

    # Partial of flight path angle (fpa) wrt (x, y, z, vx, vy, vz)
    p0 = 1 / (magr**2 * h)
    p1 = 1 / (magv**2 * h)
    tm[2, :3] = p0 * (veci_m * np.sum(reci_m**2) - reci_m * np.sum(reci_m * veci_m))
    tm[2, 3:] = p1 * (reci_m * np.sum(veci_m**2) - veci_m * np.sum(reci_m * veci_m))

    # Partial of azimuth (az) wrt (x, y, z, vx, vy, vz)
    tm[3, :] = _compute_partials_az(reci_m, veci_m)

    # Partial of r and v wrt (x, y, z, vx, vy, vz)
    tm[4, :3] = reci_m / magr
    tm[5, 3:] = veci_m / magv

    # Calculate the output covariance matrix
    flcov = tm @ cartcov @ tm.T

    return flcov, tm


def covfl2ct(
    flcov: ArrayLike,
    flstate: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
    use_latlon: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 covariance matrix from flight parameters to Cartesian elements.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        flcov (array_like): 6x6 Flight parameters covariance matrix in m and m/s
        flstate (array_like): 6x1 Flight parameters state vector in km, km/s and radians
                              (lon or ra, latgc or dec, fpa, az, r, v)
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 nutation parameters
        eqeterms (bool, optional): Add terms for ast calculation (default True)
        use_latlon (bool, optional): Flag to use lat/lon instead of ra/dec
                                     (default True)

    Returns:
        tuple: (cartcov, tm)
            cartcov (np.ndarray): 6x6 Cartesian covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix

    TODO:
        - The returned covariance for lat/lon seems incorrect - investigate
    """
    # Parse input flight state
    lon, latgc, fpa, az, magr, magv = flstate
    magr *= KM2M  # convert to meters
    magv *= KM2M

    # fpa/az trigonometric terms
    cfpa, sfpa = np.cos(fpa), np.sin(fpa)
    caz, saz = np.cos(az), np.sin(az)

    # Initialize transformation matrix
    tm = np.zeros((6, 6))

    # Decide which coordinates to use
    craf, sraf, cdf, sdf = 0, 0, 0, 0
    if use_latlon:
        # Compute ECEF coordinates
        craf, sraf = np.cos(lon), np.sin(lon)
        cdf, sdf = np.cos(latgc), np.sin(latgc)

        # Get ECEF vectors in km and km/s
        recef = np.array([magr * cdf * craf, magr * cdf * sraf, magr * sdf]) / KM2M
        vecef = (
            magv
            / KM2M
            * np.array(
                [
                    -craf * sdf * caz * cfpa - sraf * saz * cfpa + craf * cdf * sfpa,
                    -sraf * sdf * caz * cfpa + craf * saz * cfpa + sraf * cdf * sfpa,
                    sdf * sfpa + cdf * caz * cfpa,
                ]
            )
        )

        # Convert to ECI
        aecef = np.zeros(3)
        reci, veci, _ = fc.ecef2eci(
            recef,
            vecef,
            aecef,
            ttt,
            jdut1,
            lod,
            xp,
            yp,
            ddpsi,
            ddeps,
            iau80arr,
            eqeterms,
        )
        reci *= KM2M  # now in meters
        veci *= KM2M

        # Compute RA/dec from ECI vectors
        temp = np.sqrt(reci[0] ** 2 + reci[1] ** 2)
        rtasc = (
            np.arctan2(reci[1], reci[0])
            if temp >= SMALL
            else np.arctan2(veci[1], veci[0])
        )
        decl = np.arcsin(reci[2] / magr)
    else:
        # Use RA/dec directly
        rtasc, decl = lon, latgc

    # ra/dec trigonometric terms
    cra, sra = np.cos(rtasc), np.sin(rtasc)
    cd, sd = np.cos(decl), np.sin(decl)

    # Reusable variables
    if use_latlon:
        s1, c1 = sraf, craf
        s2, c2 = sdf, cdf
    else:
        s1, c1 = sra, cra
        s2, c2 = sd, cd

    # Position partials w.r.t. (lon latgc fpa az r v)
    tm[0, :2] = [-magr * c2 * s1, -magr * s2 * c1]  # rx
    tm[0, 4] = cd * cra
    tm[1, :2] = [magr * c2 * c1, -magr * s2 * s1]  # ry
    tm[1, 4] = cd * sra
    tm[2, 1] = magr * c2  # rz
    tm[2, 4] = sd

    # Partial of vx wrt (lon, latgc, fpa, az, r, v)
    tm[3, 0] = -magv * (-s1 * caz * s2 * cfpa + c1 * saz * cfpa + c2 * s1 * sfpa)
    tm[3, 1] = -c1 * magv * (s2 * sfpa + c2 * caz * cfpa)
    tm[3, 2] = magv * (cra * caz * sd * sfpa + sra * saz * sfpa + cd * cra * cfpa)
    tm[3, 3] = magv * (cra * saz * sd * cfpa - sra * caz * cfpa)
    tm[3, 5] = -cra * caz * sd * cfpa - sra * saz * cfpa + cd * cra * sfpa

    # Partial of vy wrt (lon, latgc, fpa, az, r, v)
    tm[4, 0] = magv * (-c1 * caz * s2 * cfpa - s1 * saz * cfpa + c2 * c1 * sfpa)
    tm[4, 1] = -s1 * magv * (s2 * sfpa + c2 * caz * cfpa)
    tm[4, 2] = magv * (sra * caz * sd * sfpa - cra * saz * sfpa + cd * sra * cfpa)
    tm[4, 3] = magv * (sra * saz * sd * cfpa + cra * caz * cfpa)
    tm[4, 5] = -sra * caz * sd * cfpa + cra * saz * cfpa + cd * sra * sfpa

    # Partial of vz wrt (lon, latgc, fpa, az, r, v)
    tm[5, 1] = magv * (c2 * sfpa - s2 * caz * cfpa)
    tm[5, 2] = magv * (sd * cfpa - cd * caz * sfpa)
    tm[5, 3] = -magv * cd * saz * cfpa
    tm[5, 5] = sd * sfpa + cd * caz * cfpa

    # Compute Cartesian covariance
    cartcov = tm @ flcov @ tm.T

    return cartcov, tm


########################################################################################
# Satellite Coordinate Systems
########################################################################################


def covct2rsw(
    cartcov: ArrayLike, cartstate: ArrayLike
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 Cartesian covariance matrix to an orbit plane RSW frame.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        cartcov (array_like): 6x6 Cartesian covariance matrix in (m, m/s) or (km, km/s)
        cartstate (array_like): 6x1 Cartesian state vector in km and km/s
                                (rx, ry, rz, vx, vy, vz)

    Returns:
        tuple: (covrsw, tm)
            covrsw (np.ndarray): 6x6 RSW covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Extract position and velocity vectors
    r = np.array(cartstate[:3])
    v = np.array(cartstate[3:])

    # Define RSW unit vectors
    rv = unit(r)  # along the position vector (radial direction)
    temv = np.cross(r, v)
    wv = unit(temv)  # along the angular momentum vector (out of plane)
    sv = np.cross(wv, rv)  # along the direction perpendicular to rv and wv (in-plane)

    # Initialize the transformation matrix
    tm = np.zeros((6, 6))

    # Populate the position transformation submatrix
    tm[0, 0:3] = rv
    tm[1, 0:3] = sv
    tm[2, 0:3] = wv

    # Populate the velocity transformation submatrix
    tm[3, 3:6] = rv
    tm[4, 3:6] = sv
    tm[5, 3:6] = wv

    # Transform the covariance matrix
    covrsw = tm @ cartcov @ tm.T

    return covrsw, tm


def covct2ntw(
    cartcov: ArrayLike, cartstate: ArrayLike
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms a 6x6 Cartesian covariance matrix to an orbit plane NTW frame.

    References:
        Vallado and Alfano 2015, AAS 15-537

    Args:
        cartcov (array_like): 6x6 Cartesian covariance matrix in (m, m/s) or (km, km/s)
        cartstate (array_like): 6x1 Cartesian state vector in km and km/s
                                (rx, ry, rz, vx, vy, vz)

    Returns:
        tuple: (covntw, tm)
            covntw (np.ndarray): 6x6 NTW covariance matrix in m and m/s
            tm (np.ndarray): 6x6 Transformation matrix
    """
    # Extract position and velocity vectors
    r_eci_m = np.array(cartstate[:3])
    v_eci_m = np.array(cartstate[3:])

    # Define NTW unit vectors
    tv = unit(v_eci_m)  # along the velocity vector
    temv = np.cross(r_eci_m, v_eci_m)
    wv = unit(temv)  # along the angular momentum vector
    nv = np.cross(tv, wv)  # normal to both tv and wv

    # Initialize the transformation matrix
    tm = np.zeros((6, 6))

    # Position transformation submatrix
    tm[0, 0:3] = nv
    tm[1, 0:3] = tv
    tm[2, 0:3] = wv

    # Velocity transformation submatrix
    tm[3, 3:6] = nv
    tm[4, 3:6] = tv
    tm[5, 3:6] = wv

    # Transform the covariance matrix
    covntw = tm @ cartcov @ tm.T

    return covntw, tm
