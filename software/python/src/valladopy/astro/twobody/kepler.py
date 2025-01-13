# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import logging
from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ...constants import SMALL, RE, J2, MU, TWOPI
from .frame_conversions import rv2coe, coe2rv, rv2pqw
from .newton import newtonm
from .utils import is_equatorial, findc2c3


logger = logging.getLogger(__name__)


class FGMethod(Enum):
    PQW = "pqw"
    SERIES = "series"
    C2C3 = "c2c3"


def kepler(
    ro: ArrayLike, vo: ArrayLike, dtsec: float, n_iters: int = 50
) -> Tuple[np.ndarray, np.ndarray]:
    """Solves Kepler's problem for orbit determination and returns a future geocentric
    equatorial (ECI) position and velocity vector using universal variables.

    References:
        Vallado: 2022, p. 94-96, Algorithm 8

    Args:
        ro (array_like): Initial ECI position vector in km
        vo (array_like): Initial ECI velocity vector in km/s
        dtsec (float): Time interval to propagate in seconds
        n_iters (int, optional): Number of iterations for Newton-Raphson
                                 method

    Returns:
        tuple: (r, v)
            r (np.ndarray): Propagated ECI position vector in km
            v (np.ndarray): Propagated ECI velocity vector in km/s
    """
    # Convert to numpy arrays
    ro, vo = np.array(ro), np.array(vo)

    # Initialize values
    ktr = 0
    xnew, znew, c2new, c3new = 0, 0, 0, 0
    dtnew = -10
    smu = np.sqrt(MU)
    magro = np.linalg.norm(ro)
    magvo = np.linalg.norm(vo)
    rdotv = np.dot(ro, vo)

    # Find specific mechanical energy, alpha, and semi-major axis
    sme = (magvo**2 / 2) - (MU / magro)
    alpha = -2 * sme / MU
    a = -MU / (2 * sme) if np.abs(sme) > SMALL else np.inf
    alpha = 0 if np.abs(alpha) < SMALL else alpha

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
        s = 0.5 * (np.pi / 2 - np.arctan(3 * np.sqrt(MU / (p**3)) * dtsec))
        w = np.arctan(np.tan(s) ** (1 / 3))
        xold = np.sqrt(p) * (2 / np.tan(2 * w))
        alpha = 0
    else:
        # Hyperbolic orbit
        temp = (
            -2
            * MU
            * dtsec
            / (a * (rdotv + np.sign(dtsec) * np.sqrt(-MU * a) * (1 - magro * alpha)))
        )
        xold = np.sign(dtsec) * np.sqrt(-a) * np.log(temp)

    # Newton-Raphson iteration to find x
    tmp = 1 / smu
    while (np.abs(dtnew * tmp - dtsec) >= SMALL) and (ktr < n_iters):
        xoldsqrd = xold * xold
        znew = xoldsqrd * alpha

        # Find c2 and c3 functions
        c2new, c3new = findc2c3(znew)

        # Use a newton iteration for new values
        rval = (
            xoldsqrd * c2new
            + rdotv * tmp * xold * (1 - znew * c3new)
            + magro * (1 - znew * c2new)
        )
        dtnew = (
            xoldsqrd * xold * c3new
            + rdotv * tmp * xoldsqrd * c2new
            + magro * xold * (1 - znew * c3new)
        )

        # Calculate new value for x
        temp1 = (dtsec * smu - dtnew) / rval
        xnew = xold + temp1

        # Check if the univ param goes negative; if so, use bissection
        if (xnew < 0) and (dtsec > 0):
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
    f = 1 - (xnewsqrd * c2new / magro)
    g = dtsec - xnewsqrd * xnew * c3new / smu

    # Find position and velocity vectors at new time
    r = f * ro + g * vo
    magr = np.linalg.norm(r)
    gdot = 1 - (xnewsqrd * c2new / magr)
    fdot = (smu * xnew / (magro * magr)) * (znew * c3new - 1)
    v = fdot * ro + gdot * vo

    # Check if f and g values are consistent
    temp = f * gdot - fdot * g
    if np.abs(temp - 1) > 1e-5:
        logger.warning("f and g values are inconsistent")

    return r, v


def pkepler(
    ro: ArrayLike, vo: ArrayLike, dtsec: float, ndot: float, nddot: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Propagates a satellite's position and velocity vectors over a given time period
    accounting for J2 perturbations.

    References:
        Vallado: 2022, p. 715-717, Algorithm 66

    Args:
        ro (array_like): Initial ECI position vector in km
        vo (array_like): Initial ECI velocity vector in km/s
        dtsec (float): Time interval to propagate in seconds
        ndot (float): First time derivative of mean motion in rad/s²
        nddot (float): Second time derivative of mean motion in rad/s³

    Returns:
        tuple: (r, v)
            r (np.ndarray): Propagated ECI position vector in km
            v (np.ndarray): Propagated ECI velocity vector in km/s

    TODO:
        - Move to perturbations?
    """
    # Convert position and velocity to orbital elements
    output = rv2coe(ro, vo)
    processed_output = tuple(
        0 if isinstance(x, float) and np.isnan(x) else x for x in output
    )
    p, a, ecc, incl, raan, _, nu, m, arglat, truelon, lonper, _ = processed_output

    # Check for negative semi-major axis
    if a < 0:
        logger.error("Negative semi-major axis encountered")
        return np.zeros(3), np.zeros(3)

    # Mean motion
    n = np.sqrt(MU / (a**3))

    # J2 perturbation effects
    j2op2 = (n * 1.5 * RE**2 * J2) / (p**2)
    raandot = -j2op2 * np.cos(incl)
    argpdot = j2op2 * (2 - 2.5 * np.sin(incl) ** 2)
    mdot = n

    # Update semi-major axis and eccentricity
    a -= 2 * ndot * dtsec * a / (3 * n)
    ecc -= 2 * (1 - ecc) * ndot * dtsec / (3 * n)
    p = a * (1 - ecc**2)

    # Update orbital elements
    if ecc < SMALL:
        # Circular orbit
        if is_equatorial(incl):
            # Circular equatorial
            truelondot = raandot + argpdot + mdot
            truelon = truelon + truelondot * dtsec
            truelon = np.mod(truelon, TWOPI)
        else:
            # Circular inclined
            raan = raan + raandot * dtsec
            raan = np.mod(raan, TWOPI)
            arglatdot = argpdot + mdot
            arglat = arglat + arglatdot * dtsec
            arglat = np.mod(arglat, TWOPI)
    else:
        # Elliptical orbit
        if is_equatorial(incl):
            # Elliptical equatorial
            lonperdot = raandot + argpdot
            lonper = lonper + lonperdot * dtsec
            lonper = np.mod(lonper, TWOPI)
        else:
            # Elliptical inclined
            raan = raan + raandot * dtsec
            raan = np.mod(raan, TWOPI)

        m = m + mdot * dtsec + ndot * dtsec**2 + nddot * dtsec**3
        m = np.mod(m, TWOPI)
        e0, nu = newtonm(ecc, m)

    # Convert updated orbital elements back to position and velocity
    r, v = coe2rv(p, ecc, incl, raan, nu, arglat, truelon, lonper)

    return r, v


def findfandg(
    r1: ArrayLike,
    v1: ArrayLike,
    r2: ArrayLike,
    v2: ArrayLike,
    dtsec: float,
    x: float,
    z: float,
    c2: float,
    c3: float,
    method: FGMethod,
) -> Tuple[float, float, float, float]:
    """Calculates the f and g functions for use in various applications.

    References:
        Vallado: 2022, p. 83-88

    Args:
        r1 (array_like): First position vector in km
        v1 (array_like): First velocity vector in km/s
        r2 (array_like): Second position vector in km
        v2 (array_like): Second velocity vector in km/s
        dtsec (float): Step size in seconds
        x (float): Universal variable x
        z (float): Universal variable z
        c2 (float): c2 function value
        c3 (float): c3 function value
        method (FGMethod): Method to use for calculating f and g values

    Returns:
        tuple: (f, g, fdot, gdot)
            f (float): f function value
            g (float): g function value
            fdot (float): fdot function value
            gdot (float): gdot function value

    Notes:
        - The step size `dtsec` should be small (on the order of 60-120 seconds)
    """
    f = g = fdot = gdot = 0
    magr1, magv1 = np.linalg.norm(r1), np.linalg.norm(v1)

    if method == FGMethod.PQW:
        hbar = np.cross(r1, v1)
        h = np.linalg.norm(hbar)
        rpqw1, vpqw1 = rv2pqw(r1, v1)
        rpqw2, vpqw2 = rv2pqw(r2, v2)

        f = (rpqw2[0] * vpqw1[1] - vpqw2[0] * rpqw1[1]) / h
        g = (rpqw1[0] * rpqw2[1] - rpqw2[0] * rpqw1[1]) / h
        gdot = (rpqw1[0] * vpqw2[1] - vpqw2[0] * rpqw1[1]) / h
        fdot = (vpqw2[0] * vpqw1[1] - vpqw2[1] * vpqw1[0]) / h

    elif method == FGMethod.SERIES:
        u = MU / (magr1**3)
        p = np.dot(r1, v1) / (magr1**2)
        q = (magv1**2 - u * magr1**2) / (magr1**2)

        p2, p4, p6 = p**2, p**4, p**6
        u2, u3 = u**2, u**3
        q2, q3 = q**2, q**3
        dt2, dt3, dt4, dt5, dt6, dt7, dt8 = [dtsec**i for i in range(2, 9)]

        # fmt: off
        f = (
            1 - 0.5 * u * dt2 + 0.5 * u * p * dt3
            + u / 24 * (-15 * p2 + 3 * q + u) * dt4
            + p * u / 8 * (7 * p2 - 3 * q - u) * dt5
            + u / 720
            * (-945 * p4 + 630 * p2 * q + 210 * u * p2 - 45 * q2 - 24 * u * q - u2)
            * dt6
            + p * u / 80
            * (165 * p4 - 150 * p2 * q - 50 * u * p2 + 25 * q2 + 14 * u * q + u2) * dt7
            + u / 40320
            * (
                -135135 * p6
                + 155925 * p4 * q
                + 51975 * u * p4
                - 42525 * p2 * q2
                - 24570 * u * p2 * q
                - 2205 * u2 * p2
                + 1575 * q3
                + 1107 * u * q2
                + 117 * u2 * q
                + u3
            )
            * dt8
        )

        g = (
            dtsec - 1 / 6 * u * dt3 + 0.25 * u * p * dt4
            + u / 120 * (-45 * p2 + 9 * q + u) * dt5
            + p * u / 24 * (14 * p2 - 6 * q - u) * dt6
            + u / 5040
            * (-4725 * p4 + 3150 * p2 * q + 630 * u * p2 - 225 * q2 - 54 * u * q - u2)
            * dt7
            + p * u / 320
            * (495 * p4 - 450 * p2 * q - 100 * u * p2 + 75 * q2 + 24.0 * u * q + u2)
            * dt8
        )

        fdot = (
            -u * dtsec + 1.5 * u * p * dt2
            + u / 6 * (-15 * p2 + 3 * q + u) * dt3
            + 5 * p * u / 8 * (7 * p2 - 3 * q - u) * dt4
            + u / 120
            * (-945 * p4 + 630 * p2 * q + 210 * u * p2 - 45 * q2 - 24 * u * q - u2)
            * dt5
            + 7 * p * u / 80
            * (165 * p4 - 150 * p2 * q - 50 * u * p2 + 25 * q2 + 14 * u * q + u2) * dt6
            + u / 5040
            * (-135135 * p6 + 155925 * p4 * q + 51975 * u * p4 - 42525 * p2 * q2
               - 24570 * u * p2 * q - 2205 * u2 * p2 + 1575 * q3 + 1107 * u * q2
               + 117 * u2 * q + u3) * dt7
        )

        gdot = (
            1 - 0.5 * u * dt2 + u * p * dt3
            + u / 24 * (-45 * p2 + 9 * q + u) * dt4
            + p * u / 4 * (14 * p2 - 6 * q - u) * dt5
            + u / 720
            * (-4725 * p4 + 3150 * p2 * q + 630 * u * p2 - 225 * q2 - 54 * u * q - u2)
            * dt6
            + p * u / 40
            * (495 * p4 - 450 * p2 * q - 100 * u * p2 + 75 * q2 + 24 * u * q + u2) * dt7
        )
        # fmt: on

    elif method == FGMethod.C2C3:
        xsqrd = x**2
        magr2 = np.linalg.norm(r2)
        f = 1 - (xsqrd * c2 / magr1)
        g = dtsec - xsqrd * x * c3 / np.sqrt(MU)
        gdot = 1 - (xsqrd * c2 / magr2)
        fdot = (np.sqrt(MU) * x / (magr1 * magr2)) * (z * c3 - 1)

    return f, g, fdot, gdot
