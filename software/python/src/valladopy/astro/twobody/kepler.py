# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import logging
from enum import Enum

import numpy as np

from ...constants import SMALL, MU, TWOPI
from .utils import findc2c3


logger = logging.getLogger(__name__)


class OrbitType(Enum):
    CIR_EQUATORIAL = 1  # circular equatorial
    CIR_INCLINED = 2  # circular inclined
    EPH_EQUATORIAL = 3  # elliptical, parabolic, hyperbolic equatorial
    EPH_INCLINED = 4  # elliptical, parabolic, hyperbolic inclined


def determine_orbit_type(ecc, incl, tol=SMALL):
    """Determine the type of orbit based on eccentricity and inclination

    Args:
        ecc (float): The eccentricity of the orbit
        incl (float): The inclination of the orbit in radians
        tol (float, optional): Small value for tolerance

    Returns:
        OrbitType: The type of orbit categorized into one of the following:
                   - circular equatorial
                   - circular inclined
                   - elliptical, parabolic, hyperbolic equatorial
                   - elliptical, parabolic, hyperbolic inclined
    """
    if ecc < tol:
        if (incl < tol) or (abs(incl - np.pi) < tol):
            return OrbitType.CIR_EQUATORIAL
        else:
            return OrbitType.CIR_INCLINED
    elif (incl < tol) or (abs(incl - np.pi) < tol):
        return OrbitType.EPH_EQUATORIAL
    else:
        return OrbitType.EPH_INCLINED


def newtone(ecc, e0):
    """Solves for the mean anomaly and true anomaly given the eccentric,
    parabolic, or hyperbolic anomalies.

    References:
        vallado: 2001, p. 85, Algorithm 6

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
        sinv = (np.sqrt(1.0 - ecc * ecc) * np.sin(e0)) / (1.0 - ecc * np.cos(e0))
        cosv = (np.cos(e0) - ecc) / (1.0 - ecc * np.cos(e0))
        nu = np.arctan2(sinv, cosv)
    elif ecc > 1.0001:
        # Hyperbolic orbit
        m = ecc * np.sinh(e0) - e0
        sinv = (np.sqrt(ecc * ecc - 1.0) * np.sinh(e0)) / (1.0 - ecc * np.cosh(e0))
        cosv = (np.cosh(e0) - ecc) / (1.0 - ecc * np.cosh(e0))
        nu = np.arctan2(sinv, cosv)
    else:
        # Parabolic orbit
        m = e0 + (1.0 / 3.0) * e0 * e0 * e0
        nu = 2.0 * np.arctan(e0)

    return m, nu


def newtonnu(ecc, nu, parabolic_lim_deg=168):
    """Solves for the eccentric anomaly and mean anomaly given the true
    anomaly.

    This function solves Kepler's equation when the true anomaly is known.
    The mean and eccentric, parabolic, or hyperbolic anomaly is also found.
    The default parabolic limit at 168 deg is arbitrary. The hyperbolic anomaly
    is also limited. The hyperbolic sine is used because it's not double
    valued.

    References:
        vallado: 2007, p. 85, Algorithm 5

    Args:
        ecc (float): Eccentricity of the orbit
        nu (float): True anomaly in radians
        parabolic_lim_deg (float, optional): The paraboloic limit in degrees

    Returns:
        e0 (float): Eccentric anomaly in radians
        m (float): Mean anomaly in radians
    """
    e0, m = np.inf, np.inf

    # Circular case
    if abs(ecc) < SMALL:
        e0 = nu
        m = nu
    # Elliptical case
    elif ecc < 1.0 - SMALL:
        sine = (np.sqrt(1.0 - ecc**2) * np.sin(nu)) / (1.0 + ecc * np.cos(nu))
        cose = (ecc + np.cos(nu)) / (1.0 + ecc * np.cos(nu))
        e0 = np.arctan2(sine, cose)
        m = e0 - ecc * np.sin(e0)
    # Hyperbolic case
    elif ecc > 1.0 + SMALL:
        if ecc > 1.0 and abs(nu) < np.pi - np.arccos(1.0 / ecc):
            sine = (np.sqrt(ecc**2 - 1.0) * np.sin(nu)) / (1.0 + ecc * np.cos(nu))
            e0 = np.arcsinh(sine)
            m = ecc * np.sinh(e0) - e0
    # Parabolic case
    else:
        if abs(nu) < np.radians(parabolic_lim_deg):
            e0 = np.tan(nu / 2.0)
            m = e0 + (e0**3) / 3.0

    # Update eccentric and mean anomaly to be within (0, 2pi) range
    if ecc < 1.0:
        m = np.fmod(m, 2.0 * np.pi)
        if m < 0.0:
            m += 2.0 * np.pi
        e0 = np.fmod(e0, 2.0 * np.pi)

    return e0, m


def newtonm(ecc, m, n_iter=50):
    """Solves for the eccentric anomaly and true anomaly given the mean anomaly
    using Newton-Raphson iteration.

    References:
        vallado: 2001, p. 72-75, Algorithm 2, Ex. 2-1

    Args:
        ecc (float): Eccentricity of the orbit
        m (float): Mean anomaly in radians
        n_iter (int): Number of iterations for eccentric anomaly solving

    Returns:
        e0 (float): Eccentric anomaly in radians
        nu (float): True anomaly in radians
    """
    # Define eccentricity thresholds
    # TODO: better definition/notes
    ecc_thresh_mid = 1.6
    ecc_thresh_high = 3.6

    # Hyperbolic orbit
    if (ecc - 1.0) > SMALL:
        if ecc < ecc_thresh_mid:
            if (0.0 > m > -np.pi) or (m > np.pi):
                e0 = m - ecc
            else:
                e0 = m + ecc
        else:
            if ecc < ecc_thresh_high and abs(m) > np.pi:
                e0 = m - np.sign(m) * ecc
            else:
                e0 = m / (ecc - 1.0)

        e1 = e0 + ((m - ecc * np.sinh(e0) + e0) / (ecc * np.cosh(e0) - 1.0))
        ktr = 1
        while abs(e1 - e0) > SMALL and ktr <= n_iter:
            e0 = e1
            e1 = e0 + (m - ecc * np.sinh(e0) + e0) / (ecc * np.cosh(e0) - 1.0)
            ktr += 1

        sinv = -(np.sqrt(ecc**2 - 1.0) * np.sinh(e1)) / (1.0 - ecc * np.cosh(e1))
        cosv = (np.cosh(e1) - ecc) / (1.0 - ecc * np.cosh(e1))
        nu = np.arctan2(sinv, cosv)

    # Parabolic orbit
    elif abs(ecc - 1.0) < SMALL:
        s = 0.5 * (np.pi * 0.5 - np.arctan(1.5 * m))
        w = np.arctan(np.tan(s) ** (1.0 / 3.0))
        e0 = 2.0 / np.tan(2.0 * w)
        nu = 2.0 * np.arctan(e0)

    # Elliptical orbit
    elif ecc > SMALL:
        if (0.0 > m > -np.pi) or (m > np.pi):
            e0 = m - ecc
        else:
            e0 = m + ecc
        e1 = e0 + (m - e0 + ecc * np.sin(e0)) / (1.0 - ecc * np.cos(e0))
        ktr = 1
        while abs(e1 - e0) > SMALL and ktr <= n_iter:
            e0 = e1
            e1 = e0 + (m - e0 + ecc * np.sin(e0)) / (1.0 - ecc * np.cos(e0))
            ktr += 1

        sinv = (np.sqrt(1.0 - ecc**2) * np.sin(e1)) / (1.0 - ecc * np.cos(e1))
        cosv = (np.cos(e1) - ecc) / (1.0 - ecc * np.cos(e1))
        nu = np.arctan2(sinv, cosv)

    # Circular orbit
    else:
        nu, e0 = m, m

    return e0, nu


def kepler(ro, vo, dtsec, n_iters=50):
    """Solves Kepler's problem for orbit determination and returns a future
    geocentric equatorial (ECI) position and velocity vector using universal
    variables.

    Args:
        ro (array_like): Initial ECI position vector in km
        vo (array_like): Initial ECI velocity vector in km/s
        dtsec (float): Time interval to propagate in seconds
        n_iters (int, optional): Number of iterations for Newton-Raphson
                                 method

    Returns:
        tuple: (r, v)
            r (np.array): Propagated ECI position vector in km
            v (np.array): Propagated ECI velocity vector in km/s
    """
    # Convert to numpy arrays
    ro, vo = np.array(ro), np.array(vo)

    # Initialize values
    ktr = 0
    xnew, znew, c2new, c3new = 0.0, 0.0, 0.0, 0.0
    dtnew = -10.0
    smu = np.sqrt(MU)
    magro = np.linalg.norm(ro)
    magvo = np.linalg.norm(vo)
    rdotv = np.dot(ro, vo)

    # Find specific mechanical energy, alpha, and semi-major axis
    sme = (magvo**2 / 2) - (MU / magro)
    alpha = -2 * sme / MU
    a = -MU / (2 * sme) if np.abs(sme) > SMALL else np.inf
    alpha = 0.0 if np.abs(alpha) < SMALL else alpha

    # Setup initial guess for x
    if alpha >= SMALL:
        # Circular and elliptical orbits
        period = TWOPI * np.sqrt(np.abs(a) ** 3 / MU)
        if np.abs(dtsec) > np.abs(period):
            dtsec = dtsec % period
        xold = smu * dtsec * alpha
    elif np.abs(alpha) < SMALL:
        # Parabolic orbit
        h = np.cross(ro, vo)
        magh = np.linalg.norm(h)
        p = magh**2 / MU
        s = 0.5 * (np.pi / 2 - np.arctan(3.0 * np.sqrt(MU / (p**3)) * dtsec))
        w = np.arctan(np.tan(s) ** (1 / 3))
        xold = np.sqrt(p) * (2.0 / np.tan(2.0 * w))
        alpha = 0.0
    else:
        # Hyperbolic orbit
        temp = (
            -2.0
            * MU
            * dtsec
            / (a * (rdotv + np.sign(dtsec) * np.sqrt(-MU * a) * (1.0 - magro * alpha)))
        )
        xold = np.sign(dtsec) * np.sqrt(-a) * np.log(temp)

    # Newton-Raphson iteration to find x
    tmp = 1.0 / smu
    while (np.abs(dtnew * tmp - dtsec) >= SMALL) and (ktr < n_iters):
        xoldsqrd = xold * xold
        znew = xoldsqrd * alpha

        # Find c2 and c3 functions
        c2new, c3new = findc2c3(znew)

        # Use a newton iteration for new values
        rval = (
            xoldsqrd * c2new
            + rdotv * tmp * xold * (1.0 - znew * c3new)
            + magro * (1.0 - znew * c2new)
        )
        dtnew = (
            xoldsqrd * xold * c3new
            + rdotv * tmp * xoldsqrd * c2new
            + magro * xold * (1.0 - znew * c3new)
        )

        # Calculate new value for x
        temp1 = (dtsec * smu - dtnew) / rval
        xnew = xold + temp1

        # Check if the univ param goes negative; if so, use bissection
        if (xnew < 0.0) and (dtsec > 0.0):
            xnew = xold * 0.5

        ktr += 1
        xold = xnew

    # Check for convergence
    if ktr >= n_iters:
        logger.error(
            f"Kepler not converged in {n_iters} iterations for dtsec = {dtsec}"
        )
        return np.zeros(3), np.zeros(3)

    # Find f and g values
    xnewsqrd = xnew * xnew
    f = 1.0 - (xnewsqrd * c2new / magro)
    g = dtsec - xnewsqrd * xnew * c3new / smu

    # Find position and velocity vectors at new time
    r = f * ro + g * vo
    magr = np.linalg.norm(r)
    gdot = 1.0 - (xnewsqrd * c2new / magr)
    fdot = (smu * xnew / (magro * magr)) * (znew * c3new - 1.0)
    v = fdot * ro + gdot * vo

    # Check if f and g values are consistent
    temp = f * gdot - fdot * g
    if np.abs(temp - 1.0) > 0.00001:
        logger.warning("f and g values are inconsistent")

    return r, v
