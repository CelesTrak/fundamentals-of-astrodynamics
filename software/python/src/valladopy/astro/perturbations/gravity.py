# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from .utils import legpolyn, trigpoly
from ..twobody.frame_conversions import ecef2ll
from ... import constants as const


@dataclass
class GravityFieldData:
    c: np.ndarray = None
    s: np.ndarray = None
    c_unc: np.ndarray = None
    s_unc: np.ndarray = None
    normalized: bool = False


def read_gravity_field(filename: str, normalized: bool) -> GravityFieldData:
    """Reads and stores gravity field coefficients.

    References:
        Vallado: 2022, p. 550-551

    Args:
        filename (str): The filename of the gravity field data
        normalized (bool): True if the gravity field data is normalized

    Returns:
        GravityFieldData: A dataclass containing gravity field data:
            - c (np.ndarray): Cosine coefficients
            - s (np.ndarray): Sine coefficients
            - normalized (bool): True if the gravity field data is normalized
    """
    # Load gravity field data
    file_data = np.loadtxt(filename)

    # Get the maximum degree of the gravity field
    max_degree = int(np.max(file_data[:, 0]))
    size = max_degree + 1

    # Initialize gravity field data
    gravarr = GravityFieldData(
        c=np.zeros((size, size)), s=np.zeros((size, size)), normalized=normalized
    )

    # Check if uncertainties are included in the data (columns 5 and 6)
    has_uncertainty = file_data.shape[1] >= 6
    if has_uncertainty:
        gravarr.c_unc = np.zeros((size, size))
        gravarr.s_unc = np.zeros((size, size))

    # Store gravity field coefficients
    for row in file_data:
        n, m = int(row[0]), int(row[1])
        c_value, s_value = row[2], row[3]
        gravarr.c[n, m] = c_value
        gravarr.s[n, m] = s_value

        if has_uncertainty:
            gravarr.c_unc[n, m] = row[4]
            gravarr.s_unc[n, m] = row[5]

    return gravarr


def get_norm(degree: int) -> np.ndarray:
    """Computes normalization constants for the gravity field.

    This normalization is useful for GTDS and Montenbruck-based gravity models.

    References:
        Vallado: 2022, p. 550

    Args:
        degree (int): Maximum degree of the gravity field (2 to 120)

    Returns:
        norm_arr (np.ndarray): Normalization array of shape (degree + 1, degree + 1)

    Notes:
        - Above degree 170, the factorial will return 0, thus affecting the results.
    """
    size = degree + 1
    norm_arr = np.zeros((size, size))

    for n in range(degree + 1):
        for m in range(n + 1):
            if m == 0:
                norm_arr[n, m] = np.sqrt(
                    (math.factorial(n) * (2 * n + 1)) / math.factorial(n)
                )
            else:
                norm_arr[n, m] = np.sqrt(
                    (math.factorial(n - m) * 2 * (2 * n + 1)) / math.factorial(n + m)
                )

    return norm_arr


def get_norm_gott(
    degree: int,
) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    """Get normalization arrays for the Gottlieb approach.

    References:
        Vallado: 2022, p. 600, Eq. 8-56
        Eckman, Brown, Adamo 2016 NASA report

    Args:
        degree (int): Maximum degree of the gravity field (zonals)

    Returns:
        tuple: (norm1, norm2, norm11, normn10, norm1m, norm2m, normn1)
            norm1 (np.ndarray): Normalization Legendre polynomial (1 x degree)
            norm2 (np.ndarray): Normalization Legendre polynomial (1 x degree)
            norm11 (np.ndarray): Normalization Legendre polynomial (1 x degree)
            normn10 (np.ndarray): Normalization Legendre polynomial (1 x degree)
            norm1m (np.ndarray): Normalization Legendre polynomial (degree x degree)
            norm2m (np.ndarray): Normalization Legendre polynomial (degree x degree)
            normn1 (np.ndarray): Normalization Legendre polynomial (degree x degree)
    """
    # Normalization arrays
    size = degree + 1
    norm1, norm2, norm11, normn10 = (np.zeros(size) for _ in range(4))
    norm1m, norm2m, normn1 = (np.zeros((size, size)) for _ in range(3))

    for n in range(2, size + 1):
        norm1[n - 1] = np.sqrt((2 * n + 1) / (2 * n - 1))
        norm2[n - 1] = np.sqrt((2 * n + 1) / (2 * n - 3))
        norm11[n - 1] = np.sqrt((2 * n + 1) / (2 * n)) / (2 * n - 1)
        normn10[n - 1] = np.sqrt((n + 1) * n * 0.5)
        for m in range(1, n + 1):
            norm1m[n - 1, m - 1] = np.sqrt(
                (n - m) * (2 * n + 1) / ((n + m) * (2 * n - 1))
            )
            norm2m[n - 1, m - 1] = np.sqrt(
                (n - m)
                * (n - m - 1)
                * (2 * n + 1)
                / ((n + m) * (n + m - 1) * (2 * n - 3))
            )
            normn1[n - 1, m - 1] = np.sqrt((n + m + 1) * (n - m))

    return norm1, norm2, norm11, normn10, norm1m, norm2m, normn1


def accel_gott(
    recef: ArrayLike, gravarr: GravityFieldData, degree: int, order: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute gravity acceleration using the normalized Gottlieb approach.

    References:
        Eckman, Brown, Adamo 2016 NASA report

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field
        order (int): Maximum order of the gravity field

    Returns:
        tuple: (leg_gott_n, accel)
            leg_gott_n (np.ndarray): Legendre terms ([degree + 1 x degree + 1] array)
            accel (np.ndarray): ECEF acceleration vector in km/s^2 (1 x 3 array)

    Notes:
        - This function is able to handle degree and order terms larger than 170 due to
          the formulation.
        - Includes two-body contribution

    TODO:
        - Add support for partials?
    """
    # Check to make sure gravity field data is normalized
    if not gravarr.normalized:
        raise ValueError("Gravity field data must be normalized")

    # Definitions
    ri = 1 / np.linalg.norm(recef)
    xor, yor, zor = recef * ri
    sinlat = zor
    reor = const.RE * ri
    reorn = reor
    muor2 = const.MU * ri**2

    # Normalization arrays
    norm1, norm2, norm11, normn10, norm1m, norm2m, normn1 = get_norm_gott(degree)

    # Legendre terms initialization
    size = degree + 1
    leg_gott_n = np.zeros((size, size))
    leg_gott_n[0, 0] = 1
    leg_gott_n[1, 1] = np.sqrt(3)
    leg_gott_n[1, 0] = np.sqrt(3) * sinlat

    for n in range(2, size):
        leg_gott_n[n, n] = norm11[n - 1] * leg_gott_n[n - 1, n - 1] * (2 * n - 1)

    ctil, stil = np.zeros(size), np.zeros(size)
    ctil[0], ctil[1] = 1, xor
    stil[1] = yor

    sumh, sumgm, sumj, sumk = 0, 1, 0, 0
    for n in range(2, size):
        reorn *= reor
        n2m1 = 2 * n - 1
        nm1 = n - 1
        np1 = n + 1

        # Tesseral (n, m=ni-1) initial value
        leg_gott_n[n, nm1] = normn1[nm1, nm1 - 1] * sinlat * leg_gott_n[n, n]

        # Zonal (n, m=0)
        leg_gott_n[n, 0] = (
            n2m1 * sinlat * norm1[nm1] * leg_gott_n[nm1, 0]
            - nm1 * norm2[nm1] * leg_gott_n[nm1 - 1, 0]
        ) / n

        # Tesseral (n, m=1) initial value
        leg_gott_n[n, 1] = (
            n2m1 * sinlat * norm1m[nm1, 0] * leg_gott_n[nm1, 1]
            - n * norm2m[nm1, 0] * leg_gott_n[nm1 - 1, 1]
        ) / nm1

        sumhn = normn10[nm1] * leg_gott_n[n, 1] * gravarr.c[n, 0]
        sumgmn = leg_gott_n[n, 0] * gravarr.c[n, 0] * np1

        if order > 0:
            for m in range(2, nm1):
                leg_gott_n[n, m] = (
                    n2m1 * sinlat * norm1m[nm1, m - 1] * leg_gott_n[nm1, m]
                    - (nm1 + m) * norm2m[nm1, m - 1] * leg_gott_n[nm1 - 1, m]
                ) / (n - m)

            sumjn = sumkn = 0
            ctil[n] = ctil[1] * ctil[nm1] - stil[1] * stil[nm1]
            stil[n] = stil[1] * ctil[nm1] + ctil[1] * stil[nm1]

            lim = min(n, order)
            for m in range(1, lim + 1):
                mxpnm = m * leg_gott_n[n, m]
                bnmtil = gravarr.c[n, m] * ctil[m] + gravarr.s[n, m] * stil[m]

                if m + 1 < leg_gott_n.shape[1]:
                    sumhn += normn1[nm1, m - 1] * leg_gott_n[n, m + 1] * bnmtil
                sumgmn += (n + m + 1) * leg_gott_n[n, m] * bnmtil

                bnmtm1 = gravarr.c[n, m] * ctil[m - 1] + gravarr.s[n, m] * stil[m - 1]
                anmtm1 = gravarr.c[n, m] * stil[m - 1] - gravarr.s[n, m] * ctil[m - 1]
                sumjn += mxpnm * bnmtm1
                sumkn -= mxpnm * anmtm1

            sumj += reorn * sumjn
            sumk += reorn * sumkn

        sumh += reorn * sumhn
        sumgm += reorn * sumgmn

    lambda_val = sumgm + sinlat * sumh
    accel = -muor2 * np.array(
        [lambda_val * xor - sumj, lambda_val * yor - sumk, lambda_val * zor - sumh]
    )

    return leg_gott_n, accel


def accel_gtds(recef: ArrayLike, gravarr, degree: int) -> np.ndarray:
    """Compute gravity acceleration using the GTDS approach.

    References:
        Vallado: 2022, p. 600-602

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field (1 to ~85)

    Returns:
        np.ndarray: ECEF acceleration vector in km/s² (1 x 3 array)
    """
    # Check to make sure gravity field data is normalized
    if not gravarr.normalized:
        raise ValueError("Gravity field data must be normalized")

    # Get normalization coefficients
    norm_arr = get_norm(degree)

    # Find latitude and longitude
    latgc, _, lon, _ = ecef2ll(recef)

    # Find Legendre and trigonometric polynomials
    _, legarr_gu, *_ = legpolyn(latgc, degree + 2)
    trig_arr, *_ = trigpoly(recef, latgc, lon, degree + 2)

    # Intermediate variables
    r_mag = np.linalg.norm(recef)
    oor = 1 / r_mag
    reor = const.RE * oor
    d_r_dr = d_r_dlat = d_r_dlon = 0

    for n in range(2, degree + 1):
        temp = reor**2
        sum1 = sum2 = sum3 = 0

        for m in range(0, n + 1):
            cnm = gravarr.c[n, m]
            snm = gravarr.s[n, m]
            norm = np.array(norm_arr)[n, m]

            # Take normalized coefficients and revert to unnormalized
            temparg = norm * (cnm * trig_arr[m, 1] + snm * trig_arr[m, 0])
            sum1 += legarr_gu[n, m] * temparg
            sum2 += (legarr_gu[n, m + 1] - trig_arr[m, 2] * legarr_gu[n, m]) * temparg
            sum3 += (
                m
                * legarr_gu[n, m]
                * norm
                * (snm * trig_arr[m, 1] - cnm * trig_arr[m, 0])
            )

        d_r_dr += temp * (n + 1) * sum1
        d_r_dlat += temp * sum2
        d_r_dlon += temp * sum3

    muor = const.MU * oor
    d_r_dr *= -muor * oor
    d_r_dlat *= muor
    d_r_dlon *= muor

    # Non-spherical perturbative acceleration
    x, y, z = recef
    r_delta = np.sqrt(x**2 + y**2)
    oor_delta = 1 / r_delta
    temp1 = oor * d_r_dr - z * oor**2 * oor_delta * d_r_dlat

    ax = temp1 * x - oor_delta**2 * d_r_dlon * y
    ay = temp1 * y + oor_delta**2 * d_r_dlon * x
    az = oor * d_r_dr * z + oor**2 * r_delta * d_r_dlat

    return np.array([ax, ay, az])
