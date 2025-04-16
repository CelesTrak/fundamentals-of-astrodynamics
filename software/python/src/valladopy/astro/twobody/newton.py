# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import numpy as np
from typing import Tuple

from ...constants import SMALL, TWOPI


def newtone(ecc: float, e0: float) -> Tuple[float, float]:
    """Solves for the mean anomaly and true anomaly given the eccentric, parabolic, or
    hyperbolic anomalies.

    References:
        Vallado: 2022, p. 78, Algorithm 6

    Args:
        ecc (float): Eccentricity
        e0 (float): Eccentric anomaly in radians (-2pi to 2pi)

    Returns:
        tuple: (m, nu)
            m (float): Mean anomaly in radians (0 to 2pi)
            nu (float): True anomaly in radians (0 to 2pi)
    """
    # Circular orbit case - values are same as eccentric anomaly
    if abs(ecc) < SMALL:
        return e0, e0

    # Non-circular cases
    if ecc < 0.999:
        # Elliptical orbit
        m = e0 - ecc * np.sin(e0)
        sinv = (np.sqrt(1 - ecc**2) * np.sin(e0)) / (1 - ecc * np.cos(e0))
        cosv = (np.cos(e0) - ecc) / (1 - ecc * np.cos(e0))
        nu = np.arctan2(sinv, cosv)
    elif ecc > 1.0001:
        # Hyperbolic orbit
        m = ecc * np.sinh(e0) - e0
        sinv = (np.sqrt(ecc**2 - 1) * np.sinh(e0)) / (1 - ecc * np.cosh(e0))
        cosv = (np.cosh(e0) - ecc) / (1 - ecc * np.cosh(e0))
        nu = np.arctan2(sinv, cosv)
    else:
        # Parabolic orbit
        m = e0 + (1 / 3) * e0**3
        nu = 2 * np.arctan(e0)

    return m, nu


def newtonnu(
    ecc: float, nu: float, parabolic_lim_deg: float = 168
) -> Tuple[float, float]:
    """Solves for the eccentric anomaly and mean anomaly given the true anomaly.

    This function solves Kepler's equation when the true anomaly is known. The mean and
    eccentric, parabolic, or hyperbolic anomaly is also found. The default parabolic
    limit at 168 deg is arbitrary. The hyperbolic anomaly is also limited. The
    hyperbolic sine is used because it's not double-valued.

    References:
        Vallado: 2022, p. 78, Algorithm 5

    Args:
        ecc (float): Eccentricity of the orbit
        nu (float): True anomaly in radians
        parabolic_lim_deg (float, optional): The paraboloic limit in degrees
                                             (default is 168)

    Returns:
        tuple: (e0, m)
            e0 (float): Eccentric anomaly in radians
            m (float): Mean anomaly in radians
    """
    e0, m = np.inf, np.inf

    # Circular case
    if abs(ecc) < SMALL:
        e0, m = nu, nu
    # Elliptical case
    elif ecc < 1 - SMALL:
        sine = (np.sqrt(1 - ecc**2) * np.sin(nu)) / (1 + ecc * np.cos(nu))
        cose = (ecc + np.cos(nu)) / (1 + ecc * np.cos(nu))
        e0 = np.arctan2(sine, cose)
        m = e0 - ecc * np.sin(e0)
    # Hyperbolic case
    elif ecc > 1 + SMALL:
        if ecc > 1 and abs(nu) < np.pi - np.arccos(1 / ecc):
            sine = (np.sqrt(ecc**2 - 1) * np.sin(nu)) / (1 + ecc * np.cos(nu))
            e0 = np.arcsinh(sine)
            m = ecc * np.sinh(e0) - e0
    # Parabolic case
    else:
        if abs(nu) < np.radians(parabolic_lim_deg):
            e0 = np.tan(nu / 2)
            m = e0 + (e0**3) / 3

    # Update eccentric and mean anomaly to be within (0, 2pi) range
    if ecc < 1:
        m = np.fmod(m, TWOPI)
        if m < 0:
            m += TWOPI
        e0 = np.fmod(e0, TWOPI)

    return e0, m


def newtonm(ecc: float, m: float, n_iter: int = 50) -> Tuple[float, float]:
    """Solves for the eccentric anomaly and true anomaly given the mean anomaly using
    Newton-Raphson iteration.

    References:
        Vallado: 2022, p. 65, Algorithm 2
        Oltrogge: JAS 2015

    Args:
        ecc (float): Eccentricity of the orbit
        m (float): Mean anomaly in radians
        n_iter (int, optional): Number of iterations for eccentric anomaly solving
                                (default is 50)

    Returns:
        tuple: (e0, nu)
            e0 (float): Eccentric anomaly in radians
            nu (float): True anomaly in radians
    """
    # Define eccentricity thresholds
    # TODO: better definition/notes
    ecc_thresh_mid, ecc_thresh_high = 1.6, 3.6

    # Hyperbolic orbit
    if (ecc - 1) > SMALL:
        if ecc < ecc_thresh_mid:
            if (0 > m > -np.pi) or (m > np.pi):
                e0 = m - ecc
            else:
                e0 = m + ecc
        else:
            if ecc < ecc_thresh_high and abs(m) > np.pi:
                e0 = m - np.sign(m) * ecc
            else:
                e0 = m / (ecc - 1)

        e1 = e0 + ((m - ecc * np.sinh(e0) + e0) / (ecc * np.cosh(e0) - 1))
        ktr = 1
        while abs(e1 - e0) > SMALL and ktr <= n_iter:
            e0 = e1
            e1 = e0 + (m - ecc * np.sinh(e0) + e0) / (ecc * np.cosh(e0) - 1)
            ktr += 1

        sinv = -(np.sqrt(ecc**2 - 1) * np.sinh(e1)) / (1 - ecc * np.cosh(e1))
        cosv = (np.cosh(e1) - ecc) / (1 - ecc * np.cosh(e1))
        nu = np.arctan2(sinv, cosv)

    # Parabolic orbit
    elif abs(ecc - 1) < SMALL:
        s = 0.5 * (np.pi * 0.5 - np.arctan(1.5 * m))
        w = np.arctan(np.tan(s) ** (1 / 3))
        e0 = 2 / np.tan(2 * w)
        nu = 2 * np.arctan(e0)

    # Elliptical orbit
    elif ecc > SMALL:
        if (0 > m > -np.pi) or (m > np.pi):
            e0 = m - ecc
        else:
            e0 = m + ecc
        e1 = e0 + (m - e0 + ecc * np.sin(e0)) / (1 - ecc * np.cos(e0))
        ktr = 1
        while abs(e1 - e0) > SMALL and ktr <= n_iter:
            e0 = e1
            e1 = e0 + (m - e0 + ecc * np.sin(e0)) / (1 - ecc * np.cos(e0))
            ktr += 1

        sinv = (np.sqrt(1 - ecc**2) * np.sin(e1)) / (1 - ecc * np.cos(e1))
        cosv = (np.cos(e1) - ecc) / (1 - ecc * np.cos(e1))
        nu = np.arctan2(sinv, cosv)

    # Circular orbit
    else:
        nu, e0 = m, m

    return e0, nu
