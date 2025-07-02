# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import math
from dataclasses import dataclass
from enum import Enum
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


class GravityAccelMethod(Enum):
    GOTT = "gott"
    LEAR = "lear"
    GTDS = "gtds"
    MONTENBRUCK = "mont"
    PINES = "pines"


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
    """Get normalization arrays for the Gottlieb (and Lear) approach.

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
) -> np.ndarray:
    """Compute gravity acceleration using the normalized Gottlieb approach.

    This returns the full acceleration that contains the two-body contribution.

    References:
        Eckman, Brown, Adamo 2016 NASA report

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field
        order (int): Maximum order of the gravity field

    Returns:
        np.ndarray: ECEF acceleration vector in km/s² (1 x 3 array)

    Raises:
        ValueError: If the gravity field data is not normalized

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

    return accel


def accel_lear(
    recef: ArrayLike, gravarr: GravityFieldData, degree: int, order: int
) -> np.ndarray:
    """Compute gravity acceleration using the normalized Lear approach.

    This returns the full acceleration that contains the two-body contribution.

    References:
        Eckman, Brown, Adamo 2016 NASA report

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field
        order (int): Maximum order of the gravity field

    Returns:
        np.ndarray: ECEF acceleration vector in km/s² (1 x 3 array)

    Raises:
        ValueError: If the gravity field data is not normalized
    """
    # Check to make sure gravity field data is normalized
    if not gravarr.normalized:
        raise ValueError("Gravity field data must be normalized")

    # Normalization arrays
    norm1, norm2, norm11, _, norm1m, norm2m, _ = get_norm_gott(degree)

    # Definitions
    size = degree + 1
    pnm = np.zeros((size, size))
    ppnm = np.zeros((size, size))

    e1 = recef[0] ** 2 + recef[1] ** 2
    magr2 = e1 + recef[2] ** 2
    magr = np.sqrt(magr2)
    r1 = np.sqrt(e1)
    sphi = recef[2] / magr
    cphi = r1 / magr

    sm = np.zeros(size + 1)
    cm = np.zeros(size + 1)
    sm[0] = recef[1] / r1 if r1 != 0 else 0
    cm[0] = recef[0] / r1 if r1 != 0 else 1
    sm[1] = 2 * cm[0] * sm[0]
    cm[1] = 2 * cm[0] ** 2 - 1

    reor = np.zeros(size + 1)
    reor[0] = const.RE / magr
    reor[1] = reor[0] ** 2

    root3, root5 = np.sqrt(3), np.sqrt(5)
    pn = np.zeros(size + 1)
    ppn = np.zeros(size + 1)
    pn[0] = root3 * sphi
    pn[1] = root5 * (3 * sphi**2 - 1) * 0.5
    ppn[0] = root3
    ppn[1] = root5 * 3 * sphi

    pnm[0, 0] = root3
    pnm[1, 1] = root5 * root3 * cphi * 0.5
    pnm[1, 0] = root5 * root3 * sphi
    ppnm[0, 0] = -root3 * sphi
    ppnm[1, 1] = -root3 * root5 * sphi * cphi
    ppnm[1, 0] = root5 * root3 * (1 - 2 * sphi**2)

    if degree >= 3:
        for n in range(3, size):
            nm1, nm2 = n - 1, n - 2
            reor[nm1] = reor[nm2] * reor[0]
            sm[n - 1] = 2 * cm[0] * sm[nm1 - 1] - sm[nm2 - 1]
            cm[n - 1] = 2 * cm[0] * cm[nm1 - 1] - cm[nm2 - 1]
            e1 = 2 * n - 1
            pn[n - 1] = (
                e1 * sphi * norm1[n - 1] * pn[nm1 - 1]
                - nm1 * norm2[n - 1] * pn[nm2 - 1]
            ) / n
            ppn[n - 1] = norm1[n - 1] * (sphi * ppn[nm1 - 1] + n * pn[nm1 - 1])
            pnm[n - 1, n - 1] = e1 * cphi * norm11[n - 1] * pnm[nm1 - 1, nm1 - 1]
            ppnm[n - 1, n - 1] = -n * sphi * pnm[n - 1, n - 1]

        for n in range(3, size):
            nm1 = n - 1
            e1 = (2 * n - 1) * sphi
            e2 = -n * sphi
            for m in range(1, nm1 + 1):
                e3 = norm1m[n - 1, m - 1] * pnm[nm1 - 1, m - 1]
                e4 = n + m
                e5 = (e1 * e3 - (e4 - 1) * norm2m[n - 1, m - 1] * pnm[n - 3, m - 1]) / (
                    n - m
                )
                pnm[n - 1, m - 1] = e5
                ppnm[n - 1, m - 1] = e2 * e5 + e4 * e3

    asph = np.zeros(3)
    asph[0] = -1
    for n in range(2, size):
        e1 = gravarr.c[n, 0] * reor[n - 1]
        asph[0] -= (n + 1) * e1 * pn[n - 1]
        asph[2] += e1 * ppn[n - 1]
    asph[2] *= cphi

    t1 = t3 = 0
    for n in range(2, size):
        e1 = e2 = e3 = 0
        nmodel = min(n, order)
        for m in range(1, nmodel + 1):
            tsnm = gravarr.s[n, m]
            tcnm = gravarr.c[n, m]
            tsm = sm[m - 1]
            tcm = cm[m - 1]
            tpnm = pnm[n - 1, m - 1]
            e4 = tsnm * tsm + tcnm * tcm
            e1 += e4 * tpnm
            e2 += m * (tsnm * tcm - tcnm * tsm) * tpnm
            e3 += e4 * ppnm[n - 1, m - 1]

        t1 += (n + 1) * reor[n - 1] * e1
        asph[1] += reor[n - 1] * e2
        t3 += reor[n - 1] * e3

    e4 = const.MU / magr2
    asph[0] = e4 * (asph[0] - cphi * t1)
    asph[1] = e4 * asph[1]
    asph[2] = e4 * (asph[2] + t3)
    e5 = asph[0] * cphi - asph[2] * sphi

    # Compute acceleration vector
    accel = np.array(
        [
            e5 * cm[0] - asph[1] * sm[0],
            e5 * sm[0] + asph[1] * cm[0],
            asph[0] * sphi + asph[2] * cphi,
        ]
    )

    return accel


def accel_gtds(recef: ArrayLike, gravarr: GravityFieldData, degree: int) -> np.ndarray:
    """Compute gravity acceleration perturbation using the GTDS approach.

    This returns the acceleration perturbation only (no two-body contribution).

    References:
        Vallado: 2022, p. 600-602

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field (1 to ~85)

    Returns:
        np.ndarray: ECEF acceleration perturbation vector in km/s² (1 x 3 array)

    Raises:
        ValueError: If the gravity field data is not normalized
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


def accel_mont(
    recef: ArrayLike, gravarr: GravityFieldData, degree: int, order: int
) -> np.ndarray:
    """Compute gravity acceleration perturbation using the Montenbruck approach.

    This returns the acceleration perturbation only (no two-body contribution).

    References:
        Vallado: 2022, p. 600-602

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field
        order (int): Maximum order of the gravity field

    Returns:
        np.ndarray: ECEF acceleration vector in km/s² (1 x 3 array)

    Raises:
        ValueError: If the gravity field data is not normalized
    """
    # Check to make sure gravity field data is normalized
    if not gravarr.normalized:
        raise ValueError("Gravity field data must be normalized")

    # Initialize acceleration vector
    apert = np.zeros(3)

    # Get normalization coefficients
    norm_arr = get_norm(degree + 1)

    # Body-fixed position and auxiliary quantities
    r2 = np.dot(recef, recef)
    rho = (const.RE**2) / r2
    r0 = const.RE * np.array(recef) / r2
    x0, y0, z0 = r0

    # Zonal terms
    v = np.zeros((degree + 3, order + 3))
    w = np.zeros_like(v)
    v[0, 0] = const.RE / np.sqrt(r2)
    w[0, 0] = 0
    v[1, 0] = z0 * v[0, 0]
    w[1, 0] = 0

    for n in range(2, degree + 2):
        v[n, 0] = (
            (2 * (n - 1) + 1) * z0 * v[n - 1, 0] - (n - 1) * rho * v[n - 2, 0]
        ) / n
        w[n, 0] = 0

    # Tesseral and sectoral terms
    for m in range(1, order + 2):
        v[m, m] = (2 * m - 1) * (x0 * v[m - 1, m - 1] - y0 * w[m - 1, m - 1])
        w[m, m] = (2 * m - 1) * (x0 * w[m - 1, m - 1] + y0 * v[m - 1, m - 1])
        if m <= degree:
            v[m + 1, m] = (2 * m + 1) * z0 * v[m, m]
            w[m + 1, m] = (2 * m + 1) * z0 * w[m, m]
        for n in range(m + 2, degree + 2):
            v[n, m] = (
                (2 * (n - 1) + 1) * z0 * v[n - 1, m] - (n + m - 1) * rho * v[n - 2, m]
            ) / (n - m)
            w[n, m] = (
                (2 * (n - 1) + 1) * z0 * w[n - 1, m] - (n + m - 1) * rho * w[n - 2, m]
            ) / (n - m)

    # Calculate acceleration contributions
    for m in range(order + 1):
        for n in range(m, degree + 1):
            c = gravarr.c[n, m] * norm_arr[n, m]
            s = gravarr.s[n, m] * norm_arr[n, m]
            if m == 0:
                apert[0] -= c * v[n + 1, 1]
                apert[1] -= c * w[n + 1, 1]
                apert[2] -= (n + 1) * c * v[n + 1, 0]
            else:
                fac = 0.5 * (n - m + 1) * (n - m + 2)
                apert[0] += 0.5 * (-c * v[n + 1, m + 1] - s * w[n + 1, m + 1])
                apert[0] += fac * (c * v[n + 1, m - 1] + s * w[n + 1, m - 1])
                apert[1] += 0.5 * (-c * w[n + 1, m + 1] + s * v[n + 1, m + 1])
                apert[1] += fac * (-c * w[n + 1, m - 1] + s * v[n + 1, m - 1])
                apert[2] += (n - m + 1) * (-c * v[n + 1, m] - s * w[n + 1, m])

    return const.MU / (const.RE**2) * apert


def accel_pines(
    recef: ArrayLike, gravarr: GravityFieldData, degree: int, order: int
) -> np.ndarray:
    """Compute gravity acceleration perturbation using the normalized Pines approach.

    This returns the acceleration perturbation only (no two-body contribution).

    References:
        Eckman, Brown, Adamo 2016 NASA report

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Normalized gravity field data
        degree (int): Maximum degree of the gravity field
        order (int): Maximum order of the gravity field

    Returns:
        np.ndarray: ECEF acceleration vector in km/s² (1 x 3 array)

    Raises:
        ValueError: If the gravity field data is not normalized
    """
    # Check to make sure gravity field data is normalized
    if not gravarr.normalized:
        raise ValueError("Gravity field data must be normalized")

    # Definitions
    size = degree + 3
    magr = np.linalg.norm(recef)
    s, t, u = recef / magr
    leg_pines_n = np.zeros((size, size))

    leg_pines_n[0, 0] = np.sqrt(2)
    for m in range(size):
        if m != 0:  # diagonal recursion
            leg_pines_n[m, m] = np.sqrt(1 + 1 / (2 * m)) * leg_pines_n[m - 1, m - 1]
        if m != degree + 2:  # first off-diagonal recursion
            leg_pines_n[m + 1, m] = np.sqrt(2 * m + 3) * u * leg_pines_n[m, m]
        if m < degree + 1:
            for n in range(m + 2, size):
                alpha = np.sqrt((2 * n + 1) * (2 * n - 1) / ((n - m) * (n + m)))
                beta = np.sqrt(
                    (2 * n + 1)
                    * (n - m - 1)
                    * (n + m - 1)
                    / ((2 * n - 3) * (n + m) * (n - m))
                )
                leg_pines_n[n, m] = (
                    alpha * u * leg_pines_n[n - 1, m] - beta * leg_pines_n[n - 2, m]
                )

    leg_pines_n[:size, 0] *= np.sqrt(0.5)
    rm = np.zeros(order + 2)
    im = np.zeros(order + 2)
    rm[1] = 1
    for m in range(1, order + 1):
        rm[m + 1] = s * rm[m] - t * im[m]
        im[m + 1] = s * im[m] + t * rm[m]

    rho = const.MU / (const.RE * magr)
    reor = const.RE / magr
    g1 = g2 = g3 = g4 = 0

    for n in range(degree + 1):
        g1temp = g2temp = g3temp = g4temp = 0
        sm = 0.5
        nmodel = min(order, n)

        for m in range(nmodel + 1):
            dnm = gravarr.c[n, m] * rm[m + 1] + gravarr.s[n, m] * im[m + 1]
            enm = gravarr.c[n, m] * rm[m] + gravarr.s[n, m] * im[m]
            fnm = gravarr.s[n, m] * rm[m] - gravarr.c[n, m] * im[m]
            alpha = np.sqrt(sm * (n - m) * (n + m + 1))

            g1temp += leg_pines_n[n, m] * m * enm
            g2temp += leg_pines_n[n, m] * m * fnm
            g3temp += alpha * leg_pines_n[n, m + 1] * dnm
            g4temp += (
                (n + m + 1) * leg_pines_n[n, m] + alpha * u * leg_pines_n[n, m + 1]
            ) * dnm

            if m == 0:
                sm = 1

        rho *= reor
        g1 += rho * g1temp
        g2 += rho * g2temp
        g3 += rho * g3temp
        g4 += rho * g4temp

    return np.array([g1 - g4 * s, g2 - g4 * t, g3 - g4 * u])


def _parse_method(method: str | GravityAccelMethod) -> GravityAccelMethod:
    if isinstance(method, GravityAccelMethod):
        return method
    try:
        return GravityAccelMethod(method.lower())
    except ValueError as e:
        raise ValueError(f"Unknown gravity acceleration method: {method}") from e


def accel(
    recef: ArrayLike,
    gravarr: GravityFieldData,
    degree: int,
    order: int,
    method: GravityAccelMethod = GravityAccelMethod.LEAR,
) -> np.ndarray:
    """Compute gravity acceleration using the specified method.

    Args:
        recef (array_like): ECEF position vector in km
        gravarr (GravityFieldData): Gravity field data
        degree (int): Maximum degree of the gravity field
        order (int): Maximum order of the gravity field
        method (GravityAccelMethod | str, optional): Method to use for acceleration
                                                     calculation (default is LEAR)

    Returns:
        np.ndarray: ECEF acceleration vector in km/s² (1 x 3 array)

    Raises:
        ValueError: If an unknown method is specified.
    """
    # Parse the method to handle strings
    method = _parse_method(method)

    # Full-body acceleration methods (includes two-body contribution)
    if method == GravityAccelMethod.GOTT:
        return accel_gott(recef, gravarr, degree, order)
    elif method == GravityAccelMethod.LEAR:
        return accel_lear(recef, gravarr, degree, order)

    # Perturbation-only acceleration methods (no two-body contribution)
    elif method == GravityAccelMethod.GTDS:
        accel = accel_gtds(recef, gravarr, degree)
    elif method == GravityAccelMethod.MONTENBRUCK:
        accel = accel_mont(recef, gravarr, degree, order)
    elif method == GravityAccelMethod.PINES:
        accel = accel_pines(recef, gravarr, degree, order)

    # Un-recognized method
    else:
        raise ValueError(f"Unknown gravity acceleration method: {method}")

    # Add two-body contribution to perturbation acceleration
    accel -= const.MU / np.linalg.norm(recef) ** 3 * np.array(recef)

    return accel
