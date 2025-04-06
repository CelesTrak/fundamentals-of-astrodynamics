# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import logging
import math
import numpy as np
from enum import Enum
from numpy.typing import ArrayLike
from pydantic import BaseModel, ConfigDict
from typing import Optional, Tuple

from ..twobody.utils import findc2c3
from ...constants import MU, SMALL, TWOPI
from ...mathtime.utils import safe_sqrt
from ...mathtime.vector import unit


logger = logging.getLogger(__name__)


# Constants
OOMU = 1 / np.sqrt(MU)


class DirectionOfMotion(Enum):
    LONG = "L"  # Long way
    SHORT = "S"  # Short way


class DirectionOfEnergy(Enum):
    LOW = "L"  # Low
    HIGH = "H"  # High


class DirectionOfFlight(Enum):
    DIRECT = "D"  # Direct motion
    RETROGRADE = "R"  # Retrograde motion


class LambertParams(BaseModel):
    # TODO: use this model directly?
    r1: ArrayLike
    r2: ArrayLike
    v1: Optional[ArrayLike] = None
    dm: Optional[DirectionOfMotion] = None
    de: Optional[DirectionOfEnergy] = None
    df: Optional[DirectionOfFlight] = None
    nrev: Optional[int] = None
    dtsec: Optional[float] = None

    model_config = ConfigDict(arbitrary_types_allowed=True, use_enum_values=True)


########################################################################################
# Supporting Functions
########################################################################################


def calculate_mag_and_angle(r1: ArrayLike, r2: ArrayLike) -> Tuple[float, float, float]:
    """Calculate the magnitudes of two position vectors and the cosine of the
    angle between them.

    Args:
        r1 (array_like): Initial position vector
        r2 (array_like): Final position vector

    Returns:
        tuple: (magr1, magr2, cosdeltanu)
            magr1 (float): Magnitude of the initial position vector
            magr2 (float): Magnitude of the final position vector
            cosdeltanu (float): Cosine of the angle between the two position
    """
    magr1 = float(np.linalg.norm(r1))
    magr2 = float(np.linalg.norm(r2))
    cosdeltanu = np.dot(r1, r2) / (magr1 * magr2)
    cosdeltanu = np.clip(cosdeltanu, -1, 1)  # ensure within bounds [-1, 1]
    return magr1, magr2, float(cosdeltanu)


########################################################################################
# Lambert's Problem
########################################################################################


def min_energy(
    r1: ArrayLike, r2: ArrayLike, dm: DirectionOfMotion, nrev: int
) -> Tuple[np.ndarray, float, float, float]:
    """Solves the Lambert minimum energy problem.

    References:
        Vallado: 2022, p. 479-481, Algorithm 56

    Args:
        r1 (array_like): Initial ECI position vector in km
        r2 (array_like): Final ECI position vector in km
        dm (DirectionOfMotion): Direction of motion (LONG or SHORT)
        nrev (int): Number of revolutions (0, 1, 2, ...)

    Returns:
        tuple: (v, aminenergy, tminenergy, tminabs)
            v (np.ndarray): Minimum energy velocity vector in km/s
            aminenergy (float): Minimum energy semi-major axis in km
            tminenergy (float): Minimum energy time of flight in seconds
            tminabs (float): Minimum time of flight (parabolic) in seconds
    """
    # Validate the Pydantic model
    _ = LambertParams(r1=r1, r2=r2, dm=dm, nrev=nrev)

    # Calculate r1 and r2 mag and the cosine of the angle between them
    magr1, magr2, cosdeltanu = calculate_mag_and_angle(r1, r2)

    # Compute the minimum energy semi-major axis
    c = np.sqrt(magr1**2 + magr2**2 - 2 * magr1 * magr2 * cosdeltanu)
    s = 0.5 * (magr1 + magr2 + c)
    aminenergy = 0.5 * s

    # Define alphae and betae
    alphae = np.pi
    betae = 2 * np.arcsin(np.sqrt((s - c) / s))

    # Compute the minimum energy time of flight
    # Use multiplier based on direction of motion
    sign = 1 if dm == DirectionOfMotion.SHORT else -1
    tminenergy = np.sqrt(aminenergy**3 / MU) * (
        2 * nrev * np.pi + alphae + sign * (betae - np.sin(betae))
    )

    # Calculate the parabolic time of flight, which serves as the minimum limit
    tminabs = (1 / 3) * np.sqrt(2 / MU) * (s**1.5 - (s - c) ** 1.5)

    # Compute intermediate values
    rcrossr = np.cross(r1, r2)
    magrcrossr = np.linalg.norm(rcrossr)
    pmin = magr1 * magr2 / c * (1 - cosdeltanu)
    sindeltanu = magrcrossr / (magr1 * magr2) * sign

    # Compute the minimum energy velocity vector
    v = (np.sqrt(MU * pmin) / (magr1 * magr2 * sindeltanu)) * (
        r2 - (1 - magr2 / pmin * (1 - cosdeltanu)) * r1
    )

    return v, aminenergy, tminenergy, tminabs


def min_time(
    r1: ArrayLike,
    r2: ArrayLike,
    dm: DirectionOfMotion,
    nrev: int,
    fa_tol: float = 1e-5,
    fa_iter: int = 20,
) -> Tuple[float, float, float]:
    """Solves Lambert's problem to find the minimum time of flight for the
    multi-revolution cases.

    References:
        Vallado: 2022, p. 481-482, Algorithm 57
        Prussing: JAS 2000

    Args:
        r1 (array_like): Initial ECI position vector in km
        r2 (array_like): Final ECI position vector in km
        dm (DirectionOfMotion): Direction of motion (LONG or SHORT)
        nrev (int): Number of revolutions (0, 1, 2, ...)
        fa_tol (float, optional): Tolerance for the Prussing method min TOF
                                  (defaults to 1e-5)
        fa_iter (int, optional): Maximum number of iterations for the Prussing method
                                 min TOF (defaults to 20)

    Returns:
        tuple: (tmin, tminp, tminenergy)
            tmin (float): Minimum time of flight in seconds
            tminp (float): Minimum time of flight (parabolic) in seconds
            tminenergy (float): Minimum energy time of flight in seconds
    """
    # Validate the Pydantic model
    _ = LambertParams(r1=r1, r2=r2, dm=dm, nrev=nrev)

    # Create numpy arrays and compute magnitudes of r1 and r2
    magr1, magr2, cosdeltanu = calculate_mag_and_angle(r1, r2)

    # Calculate chord and semiperimeter
    chord = np.sqrt(magr1**2 + magr2**2 - 2 * magr1 * magr2 * cosdeltanu)
    s = (magr1 + magr2 + chord) * 0.5

    # Multipliers based on direction of motion and energy
    sign_dm = 1 if dm == DirectionOfMotion.SHORT else -1

    # Calculate minimum parabolic time of flight to see if orbit is possible
    sindeltanu = sign_dm * np.linalg.norm(np.cross(r1, r2)) / (magr1 * magr2)
    sign = -1 if sindeltanu > 0 else 1
    tminp = (1 / 3) * np.sqrt(2 / MU) * ((s**1.5) + sign * (s - chord) ** 1.5)

    # Calculate minimum energy ellipse time of flight
    amin = 0.5 * s
    beta = 2 * np.arcsin(np.sqrt((s - chord) / s))
    tminenergy = (
        (amin**1.5)
        * ((2 * nrev + 1) * np.pi + sign_dm * (beta - np.sin(beta)))
        / np.sqrt(MU)
    )

    # Iteratively calculate the minimum time of flight (ellipse)
    # using Prussing method (Prussing 1992 AAS, 2000 JAS, Stern 1964 p. 230)
    an = 1.001 * amin
    i = 1
    fa, xi, eta = 10, 0, 0
    while abs(fa) > fa_tol and i <= fa_iter:
        a = an
        alp = 1 / a
        alpha = 2 * np.arcsin(np.sqrt(0.5 * s * alp))
        beta = sign_dm * 2 * np.arcsin(np.sqrt(0.5 * (s - chord) * alp))
        xi = alpha - beta
        eta = np.sin(alpha) - np.sin(beta)
        fa = (6 * nrev * np.pi + 3 * xi - eta) * (np.sin(xi) + eta) - 8 * (
            1 - np.cos(xi)
        )
        fadot = (
            (6 * nrev * np.pi + 3 * xi - eta) * (np.cos(xi) + np.cos(alpha))
            + (3 - np.cos(alpha)) * (np.sin(xi) + eta)
            - 8 * np.sin(xi)
        ) * (-alp * np.tan(0.5 * alpha)) + (
            (6 * nrev * np.pi + 3 * xi - eta) * (-np.cos(xi) - np.cos(beta))
            + (-3 - np.cos(beta)) * (np.sin(xi) + eta)
            + 8 * np.sin(xi)
        ) * (
            -alp * np.tan(0.5 * beta)
        )
        an = a - fa / fadot
        i += 1

    # Calculate the minimum time of flight
    tmin = (an**1.5) * (TWOPI * nrev + xi - eta) / np.sqrt(MU)

    return tmin, tminp, tminenergy


def tmax_rp(
    r1: ArrayLike, r2: ArrayLike, dm: DirectionOfMotion, nrev: int, tol: float = 1e-3
) -> Tuple[float, np.ndarray]:
    """Solves Lambert's problem and finds the TOF for maximum perigee radius.

    References:
        Thompson: 2019

    Args:
        r1 (array_like): Initial ECI position vector in km
        r2 (array_like): Final ECI position vector in km
        dm (DirectionOfMotion): Direction of motion (LONG or SHORT)
        nrev (int): Number of revolutions (0, 1, 2, ...)
        tol (float, optional): Distance tolerance in km (defaults to 1e-3, or 1m)

    Returns:
        tuple: (tmaxrp, v1t)
            tmaxrp (float): Maximum perigee time of flight in seconds
            v1t (np.ndarray): Initial velocity vector at r1 in km/s
    """
    v1t = np.zeros(3)
    magr1, magr2 = np.linalg.norm(r1), np.linalg.norm(r2)

    # Calculate the cosine of the angle between the two position vectors
    cos_deltanu = np.dot(r1, r2) / (magr1 * magr2)
    cos_deltanu = np.clip(cos_deltanu, -1, 1)  # ensure within [-1, 1]

    # Calculate the sine of the angle
    rcrossr = np.cross(r1, r2)
    sindeltanu = np.linalg.norm(rcrossr) / (magr1 * magr2)
    if dm == DirectionOfMotion.LONG:
        sindeltanu *= -1  # flip sign for LONG way

    # Calculate the time of flight for maximum perigee radius
    y1, y2 = MU / magr1, MU / magr2
    tempvec = r1 + r2
    if np.linalg.norm(tempvec) < tol:  # nearly circular endpoints
        c = np.sqrt(y1)
        tmaxrp = (TWOPI * nrev + np.arctan2(sindeltanu, cos_deltanu)) * (MU / c**3)
    else:
        if magr1 < magr2:
            c = np.sqrt((y2 - y1 * cos_deltanu) / (1 - cos_deltanu))
            x1, x2 = 0, (y1 - c**2) * sindeltanu
        else:
            c = np.sqrt((y1 - y2 * cos_deltanu) / (1 - cos_deltanu))
            x1, x2 = (-y2 + c**2) * sindeltanu, 0

        r = np.sqrt(x1**2 + (y1 - c**2) ** 2) / c

        # Check if acos is larger than 1
        temp = c * (r**2 + y1 - c**2) / (r * y1)
        e1 = np.arccos(np.clip(temp, -1, 1))
        if x1 < 0:
            e1 = TWOPI - e1

        temp = c * (r**2 + y2 - c**2) / (r * y2)
        e2 = np.arccos(np.clip(temp, -1, 1))
        if x2 < 0:
            e2 = TWOPI - e2
        if e2 < e1:
            e2 += TWOPI

        k = (e2 - e1) - np.sin(e2 - e1)

        tmaxrp = MU * (
            (TWOPI * nrev + k) / np.abs(c**2 - r**2) ** 1.5
            + (c * sindeltanu) / (y1 * y2)
        )

        # Close to 180 deg transfer case
        if magr2 * sindeltanu > tol:
            nunit = unit(rcrossr)
            if sindeltanu < 0:
                nunit *= -1

            nr1 = np.cross(nunit, r1)
            v1t = (x1 / c) * r1 / magr1 + (y1 / c) * (nr1 / magr1)

    return tmaxrp, v1t


########################################################################################
# Battin's Method
########################################################################################


def seebatt(v: float) -> float:
    """Recursively calculates a value used in the Lambert Battin problem using
    pre-defined coefficients.

    Args:
        v (float): Input value for the recursive calculations (v > -1)

    Returns:
        float: The computed value

    Raises:
        ValueError: If `v` is less than -1
    """
    # Check that v is greater than -1
    if v <= -1:
        raise ValueError("Input value v must be greater than -1.")

    # Coefficients derived from Battin's recursive series
    c = [
        9 / 7,
        16 / 63,
        25 / 99,
        36 / 143,
        49 / 195,
        64 / 255,
        81 / 323,
        100 / 399,
        121 / 483,
        144 / 575,
        169 / 675,
        196 / 783,
        225 / 899,
        256 / 1023,
        289 / 1155,
        324 / 1295,
        361 / 1443,
        400 / 1599,
        441 / 1763,
        484 / 1935,
    ]

    # Recursive formulaiton for the Lambert problem
    sqrtopv = np.sqrt(1 + v)
    eta = v / (1 + sqrtopv) ** 2
    ktr = 20
    term2 = 1 + c[ktr - 1] * eta
    for j in range(ktr - 2, -1, -1):
        term2 = 1 + (c[j] * eta / term2)

    return 8 * (1 + sqrtopv) / (3 + (1 / (5 + eta + ((9 / 7) * eta / term2))))


def kbatt(v: float) -> float:
    """Recursively calculates a value used in the Lambert Battin problem using
    pre-defined coefficients.

    Args:
        v (float): Input value for the recursive calculations

    Returns:
        float: The computed value
    """
    # Coefficients derived from Battin's recursive series
    d = [
        1 / 3,
        4 / 27,
        8 / 27,
        2 / 9,
        22 / 81,
        208 / 891,
        340 / 1287,
        418 / 1755,
        598 / 2295,
        700 / 2907,
        928 / 3591,
        1054 / 4347,
        1330 / 5175,
        1480 / 6075,
        1804 / 7047,
        1978 / 8091,
        2350 / 9207,
        2548 / 10395,
        2968 / 11655,
        3190 / 12987,
        3658 / 14391,
    ]

    # Initial values
    sum1, delold, termold = d[0], 1, d[0]

    # Process forwards
    i, ktr = 1, 21
    while i < ktr and abs(termold) > 1e-8:
        del_ = 1 / (1 + d[i] * v * delold)
        term = termold * (del_ - 1)
        sum1 += term
        i += 1
        delold = del_
        termold = term

    # Process backwards
    term2 = 1 + d[-1] * v
    for i in range(ktr - 2):
        sum2 = d[ktr - i - 2] * v / term2
        term2 = 1 + sum2

    return d[0] / term2


def hodograph(
    r1: ArrayLike,
    r2: ArrayLike,
    v1: ArrayLike,
    p: float,
    ecc: float,
    dnu: float,
    dtsec: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Accomplishes a 180-degree (and 360-degree) transfer for the Lambert problem.

    References:
        Thompson JGCD 2013 v34 n6 1925
        Thompson AAS GNC 2018

    Args:
        r1 (array_like): Initial position vector in km
        r2 (array_like): Final position vector in km
        v1 (array_like): Initial velocity vector in km/s
        p (float): Semi-parameter of transfer orbit in km
        ecc (float): Eccentricity of the transfer orbit
        dnu (float): Change in true anomaly in radians
        dtsec (float): Time between r1 and r2 in seconds

    Returns:
        tuple: (v1t, v2t)
            v1t (np.ndarray): Transfer velocity vector at r1 in km/s
            v2t (np.ndarray): Transfer velocity vector at r2 in km/s
    """
    # Create numpy arrays and compute magnitudes of r1 and r2
    r1, r2 = np.array(r1), np.array(r2)
    magr1, magr2 = np.linalg.norm(r1), np.linalg.norm(r2)

    # Compute parameters a and b
    a = MU * (1 / magr1 - 1 / p)
    b = (MU * ecc / p) ** 2 - a**2

    # Calculate x1 based on b
    x1 = 0 if b <= 0 else -np.sqrt(b)

    # 180-degree or multiple 180-degree transfers
    if abs(np.sin(dnu)) < SMALL:
        # Check that the cross product norm is not zero
        cross_product_r1_v1 = np.cross(r1, v1)
        norm_cross_r1_v1 = np.linalg.norm(cross_product_r1_v1)
        if norm_cross_r1_v1 < SMALL:
            raise ValueError(
                "Vectors r1 and v1 are parallel or nearly parallel;"
                " the vector normal is undefined."
            )

        # Normal vector
        nvec = cross_product_r1_v1 / norm_cross_r1_v1

        # Adjust the direction of x1 based on the time of flight
        if ecc < 1:
            ptx = TWOPI * np.sqrt(p**3 / (MU * (1 - ecc**2) ** 3))
            if dtsec % ptx > ptx * 0.5:
                x1 = -x1
    else:
        # Common path
        y2a = MU / p - x1 * np.sin(dnu) + a * np.cos(dnu)
        y2b = MU / p + x1 * np.sin(dnu) + a * np.cos(dnu)
        if abs(MU / magr2 - y2b) < abs(MU / magr2 - y2a):
            x1 = -x1

        # Check that the cross product norm is not zero
        cross_product_r1_r2 = np.cross(r1, r2)
        norm_cross_r1_r2 = np.linalg.norm(cross_product_r1_r2)
        if norm_cross_r1_r2 < SMALL:
            raise ValueError(
                "Vectors r1 and r2 are parallel or nearly parallel;"
                " the vector normal is undefined."
            )

        # Normal vector
        # Depending on the cross product, this will be normal, in plane, or even a fan
        nvec = cross_product_r1_r2 / norm_cross_r1_r2
        if dnu % TWOPI > np.pi:
            nvec = -nvec

    # Compute transfer velocity vectors
    v1t = (np.sqrt(MU * p) / magr1) * ((x1 / MU) * r1 + np.cross(nvec, r1) / magr1)
    x2 = x1 * np.cos(dnu) + a * np.sin(dnu)
    v2t = (np.sqrt(MU * p) / magr2) * ((x2 / MU) * r2 + np.cross(nvec, r2) / magr2)

    return v1t, v2t


def battin(
    r1: ArrayLike,
    r2: ArrayLike,
    v1: ArrayLike,
    dm: DirectionOfMotion,
    de: DirectionOfEnergy,
    nrev: int,
    dtsec: float,
    n_loops_he: int = 20,
    n_loops_le: int = 30,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solves Lambert's problem using Battin's method.

    This method is developed in battin (1987) and explained by Thompson 2018. It uses
    continued fractions to speed the solution and has several parameters that are
    defined differently than the traditional Gaussian technique.

    References:
        Vallado: 2022, p. 505-510, Algorithm 61
        Battin: 1987, p. 325-342
        Thompson: AAS GNC 2018

    Args:
        r1 (array_like): Initial ECI position vector in km
        r2 (array_like): Final ECI position vector in km
        v1 (array_like): Initial ECI velocity vector in km/s
                         (needed for 180-degree transfer)
        dm (DirectionOfMotion): Direction of motion (LONG or SHORT)
        de (DirectionOfEnergy): Direction of energy (LOW or HIGH)
        nrev (int): Number of revolutions (0, 1, 2, ...)
        dtsec (float): Time between r1 and r2 in seconds
        n_loops_he (int, optional): Number of loops for high energy case
                                    (defaults to 20)
        n_loops_le (int, optional): Number of loops for low energy case
                                    (defaults to 30)

    Returns:
        tuple: (v1t, v2t)
            v1t (np.ndarray): Transfer velocity vector at r1 in km/s
            v2t (np.ndarray): Transfer velocity vector at r2 in km/s

    Notes:
        - The performance over time varies depending on the input TOF (dtsec)
          for the transfer, and Some values of TOF simply have no solutions.
          See plot of time vs. psi (Figure 7-16) in Vallado for more details.
    """
    # Validate the Pydantic model
    _ = LambertParams(r1=r1, v1=v1, r2=r2, dm=dm, de=de, nrev=nrev, dtsec=dtsec)

    # Initialize values
    v1t, v2t = np.array([np.NAN] * 3), np.array([np.NAN] * 3)
    y = 0

    # Create numpy arrays and compute magnitudes of r1 and r2
    magr1, magr2, cosdeltanu = calculate_mag_and_angle(r1, r2)

    # Determine direction of flight
    magrcrossr = np.linalg.norm(np.cross(r1, r2))
    sign = 1 if dm == DirectionOfMotion.SHORT else -1
    sindeltanu = sign * magrcrossr / (magr1 * magr2)

    # Compute delta nu
    dnu = np.arctan2(sindeltanu, cosdeltanu)
    dnu = dnu if dnu >= 0 else TWOPI + dnu  # Ensure positive angle

    # Calculate chord and semiperimeter
    chord = np.sqrt(magr1**2 + magr2**2 - 2 * magr1 * magr2 * cosdeltanu)
    s = (magr1 + magr2 + chord) * 0.5
    ror = magr2 / magr1
    eps = ror - 1

    # Calculate lambda, L, and m
    lam = 1 / s * np.sqrt(magr1 * magr2) * np.cos(dnu * 0.5)
    l_ = ((1 - lam) / (1 + lam)) ** 2
    m = 8 * MU * dtsec**2 / (s**3 * (1 + lam) ** 6)

    # Initial guess for x
    xn = 1 + 4 * l_ if nrev > 0 else l_

    # Context for the safe square root function errors
    con = "Battin's method intermediate calculations: please adjust `dtsec`"

    # Default values for y and xn
    y_default, xn_default = 75, 1
    warn_msg = f"Failed to calculate y, setting to {y_default} and xn to {xn_default}"

    # High energy case adjustments for long way, retrograde multi-rev
    if de == DirectionOfEnergy.HIGH and nrev > 0:
        xn, x = 1e-20, 10
        loops = 1
        while abs(xn - x) >= SMALL and loops <= n_loops_he:
            # Calculate h1 and h2
            x = xn
            temp = 1 / (2 * (l_ - x**2))
            temp1 = safe_sqrt(x, con)
            temp2 = (nrev * np.pi * 0.5 + np.arctan(temp1)) / temp1
            h1 = temp * (l_ + x) * (1 + 2 * x + l_)
            h2 = temp * m * temp1 * ((l_ - x**2) * temp2 - (l_ + x))

            # Calculate b and f
            b = 0.25 * 27 * h2 / ((temp1 * (1 + h1)) ** 3)
            if b < 0:
                f = 2 * np.cos(1 / 3 * np.arccos(safe_sqrt(b + 1), con))
            else:
                a_ = (safe_sqrt(b, con) + safe_sqrt(b + 1, con)) ** (1 / 3)
                f = a_ + 1 / a_

            # Calculate y and xn
            y = 2 / 3 * temp1 * (1 + h1) * (safe_sqrt(b + 1, con) / f + 1)
            xn = 0.5 * (
                (m / (y**2) - (1 + l_))
                - safe_sqrt((m / (y**2) - (1 + l_)) ** 2 - 4 * l_, con)
            )

            # Check for NaN values
            if np.isnan(y):
                y, xn = y_default, xn_default
                logger.warning(warn_msg)

            loops += 1

        # Determine transfer velocity vectors for high energy case
        x = xn
        a = s * (1 + lam) ** 2 * (1 + x) * (l_ + x) / (8 * x)
        p = (2 * magr1 * magr2 * (1 + x) * np.sin(dnu * 0.5) ** 2) / (
            s * (1 + lam) ** 2 * (l_ + x)
        )
        ecc = safe_sqrt(1 - p / a, con)
        v1t, v2t = hodograph(r1, r2, v1, p, ecc, dnu, dtsec)
    else:
        # Standard processing, low energy case
        loops, x = 1, 10
        while abs(xn - x) >= SMALL and loops <= n_loops_le:
            # Calculate h1 and h2
            x = xn
            if nrev > 0:
                temp = 1 / ((1 + 2 * x + l_) * (4 * x**2))
                temp1 = (nrev * np.pi * 0.5 + np.arctan(safe_sqrt(x, con))) / safe_sqrt(
                    x, con
                )
                h1 = temp * (l_ + x) ** 2 * (3 * (1 + x) ** 2 * temp1 - (3 + 5 * x))
                h2 = temp * m * ((x**2 - x * (1 + l_) - 3 * l_) * temp1 + (3 * l_ + x))
            else:
                tempx = seebatt(x)
                denom = 1 / ((1 + 2 * x + l_) * (4 * x + tempx * (3 + x)))
                h1 = (l_ + x) ** 2 * (1 + 3 * x + tempx) * denom
                h2 = m * (x - l_ + tempx) * denom

            # Calculate y and xn
            b = 0.25 * 27 * h2 / ((1 + h1) ** 3)
            u = 0.5 * b / (1 + safe_sqrt(1 + b, con))
            k2 = kbatt(u)
            y = ((1 + h1) / 3) * (2 + safe_sqrt(1 + b, con) / (1 + 2 * u * k2 * k2))
            xn = safe_sqrt(((1 - l_) * 0.5) ** 2 + m / (y**2), con) - (1 + l_) * 0.5

            # Check for NaN values
            if np.isnan(y):
                y, xn = y_default, xn_default
                logger.warning(warn_msg)

            loops += 1

        # Determine transfer velocity vectors for standard case
        if loops < n_loops_le:
            p = (2 * magr1 * magr2 * y**2 * (1 + x) ** 2 * np.sin(dnu * 0.5) ** 2) / (
                m * s * (1 + lam) ** 2
            )
            ecc = safe_sqrt(
                (
                    eps**2
                    + 4
                    * magr2
                    / magr1
                    * np.sin(dnu * 0.5) ** 2
                    * ((l_ - x) / (l_ + x)) ** 2
                )
                / (eps**2 + 4 * magr2 / magr1 * np.sin(dnu * 0.5) ** 2),
                con,
            )
            v1t, v2t = hodograph(r1, r2, v1, p, ecc, dnu, dtsec)

    return v1t, v2t


########################################################################################
# Universal Variable Lambert Problem
########################################################################################


def _calculate_c2dot_c3dot(
    psi: float, c2: float, c3: float, tol: float = 1e-5
) -> Tuple[float, float, float, float]:
    """Calculate the derivatives of c2 and c3 with respect to psi.

    Args:
        psi (float): Psi value
        c2 (float): Coefficient c2
        c3 (float): Coefficient c3
        tol (float, optional): Tolerance for small psi (defaults to 1e-5)

    Returns:
        tuple: (c2dot, c3dot, c2ddot, c3ddot)
            c2dot (float): Derivative of c2 with respect to psi
            c3dot (float): Derivative of c3 with respect to psi
            c2ddot (float): Second derivative of c2 with respect to psi
            c3ddot (float): Second derivative of c3 with respect to psi
    """
    if abs(psi) > tol:
        c2dot = 0.5 / psi * (1 - psi * c3 - 2 * c2)
        c3dot = 0.5 / psi * (c2 - 3 * c3)
        c2ddot = 1 / (4 * psi**2) * ((8 - psi) * c2 + 5 * psi * c3 - 4)
        c3ddot = 1 / (4 * psi**2) * ((15 - psi) * c3 - 7 * c2 + 1)
    else:
        # Taylor series expansion for small psi
        c2dot = sum(
            (-1) ** (i + 1) * (i + 1) * psi**i / math.factorial(2 * i + 4)
            for i in range(5)
        )
        c3dot = sum(
            (-1) ** (i + 1) * (i + 1) * psi**i / math.factorial(2 * i + 5)
            for i in range(5)
        )
        c2ddot, c3ddot = 0, 0

    return c2dot, c3dot, c2ddot, c3ddot


def _calculate_dtdpsi(
    x: float, c2: float, c3: float, c2dot: float, c3dot: float, vara: float, y: float
) -> float:
    """Calculate the derivative of time of flight with respect to psi.

    Args:
        x (float): Value of x
        c2 (float): Coefficient c2
        c3 (float): Coefficient c3
        c2dot (float): Derivative of c2 with respect to psi
        c3dot (float): Derivative of c3 with respect to psi
        vara (float): Value of vara
        y (float): Value of y

    Returns:
        float: Derivative of time of flight with respect to psi
    """
    return (
        x**3 * (c3dot - 3 * c3 * c2dot / (2 * c2))
        + 0.125 * vara * (3 * c3 * np.sqrt(y) / c2 + vara / x)
    ) * OOMU


def get_kbiu(r1: ArrayLike, r2: ArrayLike, order: int) -> Tuple[np.ndarray, np.ndarray]:
    """Form the minimum time and universal variable matrix for multi-rev cases.

    Args:
        r1 (array_like): Initial ECI position vector in km
        r2 (array_like): Final ECI position vector in km
        order (int): The number of revolutions to consider

    Returns:
        tuple: (kbi_arr, kbil_arr)
            kbi_arr (np.ndarray): Matrix of kbi and time of flight for SHORT direction
            kbil_arr (np.ndarray): Matrix of kbi and time of flight for LONG direction
    """
    kbi_arr = np.zeros((order, 2))
    for i in range(order):
        kbi, tof = universal_min(r1, r2, DirectionOfMotion.SHORT, i + 1)
        kbi_arr[i, 0] = kbi
        kbi_arr[i, 1] = tof

    kbil_arr = np.zeros((order, 2))
    for i in range(order):
        kbil, tofl = universal_min(r1, r2, DirectionOfMotion.LONG, i + 1)
        kbil_arr[i, 0] = kbil
        kbil_arr[i, 1] = tofl

    return kbi_arr, kbil_arr


def universal_min(
    r1: ArrayLike, r2: ArrayLike, dm: DirectionOfMotion, nrev: int, n_iter: int = 20
) -> Tuple[float, float]:
    """Find the minimum kbi value for the universal variable Lambert problem for the
    multi-rev cases.

    References:
        Arora and Russell: AAS 10-198

    Args:
        r1 (array_like): Initial ECI position vector in km
        r2 (array_like): Final ECI position vector in km
        dm (DirectionOfMotion): Direction of motion (LONG or SHORT)
        nrev (int): Number of revolutions (0, 1, 2, ...)
        n_iter (int): Number of iterations to perform (defaults to 20)

    Returns:
        tuple: (kbi, tof)
            kbi (float): K value for min TOF
            tof (float): Time of flight in seconds

    TODO:
        - Identify and capture any exceptions that may occur due to bad inputs
    """
    # Validate the Pydantic model
    _ = LambertParams(r1=r1, r2=r2, dm=dm, nrev=nrev)

    # Create numpy arrays and compute magnitudes of r1 and r2
    magr1, magr2, cosdeltanu = calculate_mag_and_angle(r1, r2)

    # Determine vara based on direction of motion
    sign_dm = -1 if dm == DirectionOfMotion.LONG else 1
    vara = sign_dm * np.sqrt(magr1 * magr2 * (1.0 + cosdeltanu))

    # Outer bounds for the nrev case
    lower = 4 * nrev**2 * np.pi**2
    upper = 4 * (nrev + 1) ** 2 * np.pi**2

    # Streamline by narrowing down the bounds (since it's near the center)
    upper = lower + (upper - lower) * 0.6
    lower = lower + (upper - lower) * 0.3

    # Initial psi guess (put in center of bounds)
    psiold = (upper + lower) * 0.5

    # Get initial values of c2 and c3
    c2, c3 = findc2c3(psiold)

    # Iterative loop to find minimum psi
    x, sqrty, loops = 0, 0, 0
    dtdpsi, psinew = 200, psiold
    while abs(dtdpsi) >= 0.1 and loops < n_iter:
        # Calculate y and x
        if abs(c2) > SMALL:
            y = magr1 + magr2 - (vara * (1 - psiold * c3) / np.sqrt(c2))
        else:
            y = magr1 + magr2
        x = np.sqrt(y / c2) if abs(c2) > SMALL else 0

        # Calculate derivatives of c2 and c3
        c2dot, c3dot, c2ddot, c3ddot = _calculate_c2dot_c3dot(psiold, c2, c3)

        # Solve for dt = 0.0
        dtdpsi = _calculate_dtdpsi(x, c2, c3, c2dot, c3dot, vara, y)

        # Solve for second derivative of dt with respect to psi
        sqrty = np.sqrt(y)
        q = 0.25 * vara * np.sqrt(c2) - x**2 * c2dot
        s1 = -24 * q * x**3 * c2 * sqrty * c3dot
        s2 = 36 * q * x**3 * sqrty * c3 * c2dot - 16 * x**5 * sqrty * c3ddot * c2**2
        s3 = (
            24 * x**5 * sqrty * (c3dot * c2dot * c2 + c3 * c2ddot * c2 - c3 * c2dot**2)
            - 6 * vara * c3dot * y * c2 * x**2
        )
        s4 = (
            -0.75 * vara**2 * c3 * c2**1.5 * x**2
            + 6 * vara * c3 * y * c2dot * x**2
            + (vara**2 * c2 * (0.25 * vara * np.sqrt(c2) - x**2 * c2)) * sqrty / x
        )
        dtdpsi2 = -(s1 + s2 + s3 + s4) / (16 * np.sqrt(MU) * (c2**2 * sqrty * x**2))

        # Newton-Raphson update for psi
        psinew = psiold - dtdpsi / dtdpsi2

        # Update psi and c2, c3 for the next iteration
        psiold = psinew
        c2, c3 = findc2c3(psiold)
        loops += 1

    # Calculate time of flight and final psi
    dtnew = (x**3 * c3 + vara * sqrty) * OOMU
    tof = dtnew
    kbi = psinew

    return kbi, tof


def universal(
    r1: ArrayLike,
    v1: ArrayLike,
    r2: ArrayLike,
    dtsec: float,
    dm: DirectionOfMotion,
    de: DirectionOfEnergy,
    nrev: int,
    kbi: float,
    tol: float = 1e-05,
    n_iter: int = 20,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solves Lambert's problem for orbit determination, returns the velocity vectors at
    each of two given position vectors.

    The solution uses universal variables for calculation and a bissection technique
    updating psi.

    References:
        Vallado: 2022, p. 499-505, Algorithm 60

    Args:
        r1 (array_like): Initial position vector in km
        v1 (array_like): Initial velocity vector in km/s
        r2 (array_like): Final position vector in km
        dtsec (float): Time of flight in seconds
        dm (DirectionOfMotion): Direction of motion (LONG or SHORT period)
        de (DirectionOfEnergy): Direction of energy (HIGH or LOW)
        nrev (int): Number of revolutions (0, 1, 2, ...)
        kbi (float): K value for the minimum time of flight
        tol (float): Tolerance for the Lambert problem (defaults to 1e-05)
                     (can affect cases where znew is multiples of 2pi^2)
        n_iter (int): Maximum number of iterations to perform (defaults to 20)

    Returns:
        tuple: (v1dv, v2dv)
            v1dv (np.ndarray): Transfer ECI velocity vector at r1 in km/s
            v2dv (np.ndarray): Transfer ECI velocity vector at r2 in km/s

    Notes:
        - If the orbit is not possible, the method will log an error and return the
          Battin method or a Hohmann transfer
        - This method is sensitive to inputs, specifically `dtsec` and `psi` values in
          `psi_vec` (for multi-rev cases) - a bad combination can lead to no solutions
          (see Vallado 2022, Figure 7-16)
    """
    # Validate the Pydantic model
    _ = LambertParams(r1=r1, v1=v1, r2=r2, dtsec=dtsec, dm=dm, de=de, nrev=nrev)

    # Definitions and initialization
    max_ynegktr_iters = 10  # maximum number of iterations for y < 0
    v1dv, v2dv = np.zeros(3), np.zeros(3)

    # Compute magnitudes of r1 and r2
    magr1, magr2, cosdeltanu = calculate_mag_and_angle(r1, r2)

    # Determine vara based on direction of motion
    sign = -1 if dm == DirectionOfMotion.LONG else 1
    vara = sign * np.sqrt(magr1 * magr2 * (1 + cosdeltanu))

    # Set up initial bounds for bisection
    if nrev == 0:
        # Hyperbolic and parabolic solutions
        lower = -16 * np.pi**2
        upper = 4 * np.pi**2
    else:
        lower = 4 * nrev**2 * np.pi**2
        upper = 4 * (nrev + 1) ** 2 * np.pi**2
        if de == DirectionOfEnergy.HIGH:
            upper = kbi
        else:
            lower = kbi

    # Form initial guess for psi
    if nrev == 0:
        # From empirical data (Arora and Russell, 2010)
        psiold = (np.log(dtsec) - 9.61202327) / 0.10918231
        psiold = min(psiold, upper - np.pi)
    else:
        psiold = lower + (upper - lower) * 0.5

    # Initial values for c2 and c3
    c2new, c3new = findc2c3(psiold)

    # Compute initial y and dtold
    if abs(c2new) > tol:
        y = magr1 + magr2 - (vara * (1 - psiold * c3new) / np.sqrt(c2new))
    else:
        y = magr1 + magr2

    xold = np.sqrt(y / c2new) if abs(c2new) > tol else 0
    dtold = (xold**3 * c3new + vara * np.sqrt(y)) * OOMU

    # Check if orbit is not possible
    if abs(vara) < 0.2:  # not exactly zero
        # Can't do bissection because w series is not accurate
        logger.error("Orbit is not possible")

        # Call Battin method or use a Hohmann transfer in 3D
        atx = (MU * (dtsec / np.pi) ** 2) ** (1 / 3)  # half period
        v1tmag = np.sqrt(2 * MU / magr1 - MU / atx)
        v2tmag = np.sqrt(2 * MU / magr2 - MU / atx)

        # Compute the direction of the velocity vectors
        wxu = unit(np.cross(r1, v1))
        v1diru = unit(np.cross(r1, wxu))
        v2diru = unit(np.cross(r2, wxu))

        # Compute the velocities using the direction and magnitude
        v1dv = -v1tmag * v1diru
        v2dv = -v2tmag * v2diru

        return v1dv, v2dv

    # Loop for iteration
    loops, ynegktr, dtnew = 0, 1, -10

    # Loop to find psi
    while abs(dtnew - dtsec) >= tol and loops < n_iter and ynegktr <= max_ynegktr_iters:
        if abs(c2new) > tol:
            y = magr1 + magr2 - (vara * (1 - psiold * c3new) / np.sqrt(c2new))
        else:
            y = magr1 + magr2

        # Check for negative y values
        if y < 0 < vara:
            while y < 0 and ynegktr < max_ynegktr_iters:
                psinew = (
                    0.8 * (1 / c3new) * (1 - (magr1 + magr2) * np.sqrt(c2new) / vara)
                )
                c2new, c3new = findc2c3(psinew)
                psiold = psinew
                lower = psiold
                y = magr1 + magr2 - (vara * (1 - psiold * c3new) / np.sqrt(c2new))
                ynegktr += 1

        loops += 1

        # Check for convergence
        if ynegktr < max_ynegktr_iters:
            xold = np.sqrt(y / c2new) if abs(c2new) > tol else 0
            dtnew = (xold**3 * c3new + vara * np.sqrt(y)) * OOMU

            # Newton-Raphson iteration for psi update
            c2dot, c3dot, _, _ = _calculate_c2dot_c3dot(psiold, c2new, c3new)

            # Calculate new psi
            dtdpsi = _calculate_dtdpsi(xold, c2new, c3new, c2dot, c3dot, vara, y)
            psinew = psiold - (dtnew - dtsec) / dtdpsi

            # Newton iteration test to see if it keeps within the bounds
            if psinew > upper or psinew < lower:
                # psi is outside bounds (too steep a slope)
                # Re-adjust upper and lower bounds
                if de == DirectionOfEnergy.LOW or nrev == 0:
                    if dtold < dtsec:
                        lower = psiold
                    else:
                        upper = psiold
                else:
                    if dtold < dtsec:
                        upper = psiold
                    else:
                        lower = psiold

                psinew = (upper + lower) * 0.5

            # Update c2 and c3
            c2new, c3new = findc2c3(psinew)
            psiold = psinew
            dtold = dtnew

            # Make sure the first guess isn't too close
            if abs(dtnew - dtsec) < tol and loops == 1:
                dtnew = dtsec - 1

    # Check for non-convergence
    if loops >= n_iter or ynegktr >= max_ynegktr_iters:
        logger.error("Lambert did not converge! Try different values of dtsec or psi")
        if ynegktr >= max_ynegktr_iters:
            logger.error("y is negative")
    else:
        # Compute velocity vectors using f and g series
        f = 1 - y / magr1
        gdot = 1 - y / magr2
        g = 1 / (vara * np.sqrt(y / MU))

        for i in range(3):
            v1dv[i] = (r2[i] - f * r1[i]) * g
            v2dv[i] = (gdot * r2[i] - r1[i]) * g

    return v1dv, v2dv
