# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 1 March 2001
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np

from ... import constants as const
from ...mathtime.vector import unit
from ..twobody.frame_conversions import rv2ntw
from ..twobody.kepler import kepler
from . import utils as utils


########################################################################################
# Coplanar Transfers
########################################################################################


def hohmann(
    rinit: float,
    rfinal: float,
    einit: float,
    efinal: float,
    nuinit: float,
    nufinal: float,
) -> Tuple[float, float, float]:
    """Calculates the delta-v values for a Hohmann transfer, either circle-to-circle
    or ellipse-to-ellipse.

    References:
        Vallado 2022, p. 326-329, Algorithm 36

    Args:
        rinit (float): Initial position magnitude in km
        rfinal (float): Final position magnitude in km
        einit (float): Eccentricity of the initial orbit
        efinal (float): Eccentricity of the final orbit
        nuinit (float): True anomaly of the initial orbit in radians (0 or pi)
        nufinal (float): True anomaly of the final orbit in radians (0 or pi)

    Returns:
        tuple: (deltava, deltavb, dtsec)
            deltava (float): Change in velocity at point A in km/s
            deltavb (float): Change in velocity at point B in km/s
            dtsec (float): Time of flight for the transfer in seconds
    """
    # Semi-major axes of initial, transfer, and final orbits
    ainit = utils.semimajor_axis(rinit, einit, nuinit)
    atran = (rinit + rfinal) / 2
    afinal = utils.semimajor_axis(rfinal, efinal, nufinal)

    # Initialize outputs
    deltava, deltavb, dtsec = 0, 0, 0

    if einit < 1 or efinal < 1:
        # Calculate delta-v at point A
        vinit = utils.velocity_mag(rinit, ainit)
        vtrana = utils.velocity_mag(rinit, atran)
        deltava = np.abs(vtrana - vinit)

        # Calculate delta-v at point B
        vfinal = utils.velocity_mag(rfinal, afinal)
        vtranb = utils.velocity_mag(rfinal, atran)
        deltavb = np.abs(vfinal - vtranb)

        # Calculate transfer time of flight
        dtsec = utils.time_of_flight(atran)

    return deltava, deltavb, dtsec


def bielliptic(
    rinit: float,
    rb: float,
    rfinal: float,
    einit: float,
    efinal: float,
    nuinit: float,
    nufinal: float,
) -> Tuple[float, float, float, float]:
    """Calculates the delta-v values for a bi-elliptic transfer, either circle-to-circle
    or ellipse-to-ellipse.

    References:
        Vallado 2022, pp. 326-330, Algorithm 37

    Args:
        rinit (float): Initial position magnitude in km
        rb (float): Intermediate position magnitude in km
        rfinal (float): Final position magnitude in km
        einit (float): Eccentricity of the initial orbit
        efinal (float): Eccentricity of the final orbit
        nuinit (float): True anomaly of the initial orbit in radians (0 or pi)
        nufinal (float): True anomaly of the final orbit in radians (0 or pi)

    Returns:
        tuple: (deltava, deltavb, deltavc, dtsec)
            deltava (float): Change in velocity at point A in km/s
            deltavb (float): Change in velocity at point B in km/s
            deltavc (float): Change in velocity at point C in km/s
            dtsec (float): Time of flight for the transfer in seconds
    """
    # Semi-major axes of initial, transfer, and final orbits
    ainit = utils.semimajor_axis(rinit, einit, nuinit)
    atran1 = (rinit + rb) * 0.5
    atran2 = (rb + rfinal) * 0.5
    afinal = utils.semimajor_axis(rfinal, efinal, nufinal)

    # Initialize outputs
    deltava, deltavb, deltavc, dtsec = 0, 0, 0, 0

    # Check if inputs represent elliptical orbits
    if einit < 1 and efinal < 1:
        # Calculate delta-v at point A
        vinit = utils.velocity_mag(rinit, ainit)
        vtran1a = utils.velocity_mag(rinit, atran1)
        deltava = np.abs(vtran1a - vinit)

        # Calculate delta-v at point B
        vtran1b = utils.velocity_mag(rb, atran1)
        vtran2b = utils.velocity_mag(rb, atran2)
        deltavb = np.abs(vtran1b - vtran2b)

        # Calculate delta-v at point C
        vtran2c = utils.velocity_mag(rfinal, atran2)
        vfinal = utils.velocity_mag(rfinal, afinal)
        deltavc = np.abs(vfinal - vtran2c)

        # Calculate total time of flight
        dtsec = utils.time_of_flight(atran1) + utils.time_of_flight(atran2)

    return deltava, deltavb, deltavc, dtsec


def one_tangent(
    rinit: float,
    rfinal: float,
    efinal: float,
    nuinit: float,
    nutran: float,
    tol: float = 1e-6,
) -> Tuple[float, float, float, float, float, float, float]:
    """Calculates the delta-v values for a one tangent transfer, either circle-to-circle
    or ellipse-to-ellipse.

    References:
        Vallado 2022, p. 335-338, Algorithm 38

    Args:
        rinit (float): Initial position magnitude in km
        rfinal (float): Final position magnitude in km
        efinal (float): Eccentricity of the final orbit
        nuinit (float): True anomaly of the initial orbit in radians (0 or pi)
        nutran (float): True anomaly of the transfer orbit in radians
                        (same quadrant as `nuinit`)
        tol (float): Tolerance for checking if transfer orbit is valid (default 1e-6)

    Returns:
        tuple: (deltava, deltavb, dtsec, etran, atran, vtrana, vtranb)
            deltava (float): Change in velocity at point A in km/s
            deltavb (float): Change in velocity at point B in km/s
            dtsec (float): Time of flight for the transfer in seconds
            etran (float): Eccentricity of the transfer orbit
            atran (float): Semi-major axis of the transfer orbit in km
            vtrana (float): Velocity of the transfer orbit at point A in km/s
            vtranb (float): Velocity of the transfer orbit at point B in km/s

    Raises:
        ValueError: If the one-tangent burn is not possible for the given inputs
    """
    # Initialize transfer time
    dtsec = 0

    # Ratio of initial to final orbit radii
    ratio = rinit / rfinal

    # Determine eccentricity of transfer orbit
    if abs(nuinit) < 0.01:  # near 0 or 180 degrees
        etran = (ratio - 1) / (np.cos(nutran) - ratio)
    else:
        etran = (ratio - 1) / (np.cos(nutran) + ratio)

    # Check if transfer orbit is valid
    if etran >= 0:
        # Semi-major axes of initial, final, and transfer orbits
        afinal = utils.semimajor_axis(rfinal, efinal, nutran)

        if abs(etran - 1) > tol:
            atran = utils.semimajor_axis(rinit, etran, nuinit)
        else:
            atran = np.inf  # parabolic orbit (infinite semi-major axis)

        # Calculate delta-V at point A
        vinit = np.sqrt(const.MU / rinit)
        vtrana = utils.velocity_mag(rinit, atran)
        deltava = np.abs(vtrana - vinit)

        # Calculate delta-V at point B
        vfinal = utils.velocity_mag(rfinal, afinal)
        vtranb = utils.velocity_mag(rfinal, atran)
        fpatranb = np.arctan2(etran * np.sin(nutran), 1 + etran * np.cos(nutran))
        fpafinal = np.arctan2(efinal * np.sin(nutran), 1 + efinal * np.cos(nutran))
        deltavb = utils.deltav(vfinal, vtranb, fpatranb - fpafinal)

        # Calculate time of flight
        if etran < 1:
            sinv = (np.sqrt(1 - etran**2) * np.sin(nutran)) / (
                1 + etran * np.cos(nutran)
            )
            cosv = (etran + np.cos(nutran)) / (1 + etran * np.cos(nutran))
            e = np.arctan2(sinv, cosv)
            eainit = 0 if abs(nuinit) < 0.01 else np.pi
            dtsec = np.sqrt(atran**3 / const.MU) * (
                e - etran * np.sin(e) - (eainit - etran * np.sin(eainit))
            )
    else:
        raise ValueError("The one-tangent burn is not possible for this case.")

    return deltava, deltavb, dtsec, etran, atran, vtrana, vtranb


########################################################################################
# Non-Coplanar Transfers
########################################################################################


def incl_only(deltai: float, vinit: float, fpa: float) -> float:
    """Calculates the delta-v for a change in inclination only.

    References:
        Vallado 2022, p. 346-348, Algorithm 39

    Args:
        deltai (float): Change in inclination in radians
        vinit (float): Initial velocity in km/s
        fpa (float): Flight path angle in radians

    Returns:
        float: Delta-v required for inclination change in km/s

    Notes:
        - Units are flexible for `vinit` and the output will match its units
    """
    return 2 * vinit * np.cos(fpa) * np.sin(0.5 * deltai)


def node_only(
    iinit: float,
    ecc: float,
    deltaraan: float,
    vinit: float,
    fpa: float,
    incl: float,
    tol: float = 1e-7,
) -> Tuple[float, float, float, float]:
    """Calculates the delta-v for a change in longitude of the ascending node.

    References:
        Vallado 2022, pp. 349-351, Algorithm 40

    Args:
        iinit (float): Initial inclination in radians
        ecc (float): Eccentricity of the initial orbit
        deltaraan (float): Change in right ascension of the ascending node in radians
        vinit (float): Initial velocity in km/s
        fpa (float): Flight path angle in radians
        incl (float): Inclination in radians
        tol (float): Tolerance for checking if transfer is elliptical (default 1e-7)

    Returns:
        tuple: (ifinal, deltav, arglat_init, arglat_final)
            ifinal (float): Final inclination in radians
            deltav (float): Change in velocity in km/s
            arglat_init (float): Initial argument of latitude in radians
            arglat_final (float): Final argument of latitude in radians

    Notes:
        - Units are flexible for `vinit` and `deltav` will match its units
    """
    if ecc > tol:
        # Elliptical orbit
        theta = np.arctan(np.sin(iinit) * np.tan(deltaraan))
        ifinal = np.arcsin(np.sin(theta) / np.sin(deltaraan))
        deltav = 2 * vinit * np.cos(fpa) * np.sin(0.5 * theta)

        # Initial argument of latitude
        arglat_init = np.pi / 2  # set at 90 degrees

    else:
        # Circular orbit
        ifinal = incl
        theta = np.arccos(np.cos(iinit) ** 2 + np.sin(iinit) ** 2 * np.cos(deltaraan))
        deltav = 2 * vinit * np.sin(0.5 * theta)

        # Initial argument of latitude
        arglat_init = np.arccos(
            (np.tan(iinit) * (np.cos(deltaraan) - np.cos(theta))) / np.sin(theta)
        )

    # Final argument of latitude
    arglat_final = np.arccos(
        (np.cos(incl) * np.sin(incl) * (1 - np.cos(deltaraan))) / np.sin(theta)
    )

    return ifinal, deltav, arglat_init, arglat_final


def incl_and_node(
    iinit: float, ifinal: float, deltaraan: float, vinit: float, fpa: float
) -> Tuple[float, float, float]:
    """Calculates the delta-v for a change in inclination and right ascension of the
    ascending node.

    References:
        Vallado 2022, p. 352, Algorithm 41

    Args:
        iinit (float): Initial inclination in radians
        ifinal (float): Final inclination in radians
        deltaraan (float): Change in right ascension of the ascending node in radians
        vinit (float): Initial velocity in km/s
        fpa (float): Flight path angle in radians

    Returns:
        tuple: (deltav, arglat_init, arglat_final)
            deltav (float): Change in velocity in km/s
            arglat_init (float): Initial argument of latitude in radians
            arglat_final (float): Final argument of latitude in radians

    Notes:
        - Units are flexible for `vinit` and `deltav` will match its units
    """
    # Pre-compute trigonometric values for efficiency
    cosdraan = np.cos(deltaraan)
    sinii, cosii = np.sin(iinit), np.cos(iinit)
    sinif, cosif = np.sin(ifinal), np.cos(ifinal)

    # Calculate theta
    cost = cosii * cosif + sinii * sinif * cosdraan
    theta = np.arccos(cost)
    sint = np.sin(theta)

    # Calculate delta-v
    deltav = incl_only(theta, vinit, fpa)

    # Calculate argument of latitude changes
    arglat_init = np.arccos((sinif * cosdraan - cost * sinii) / (sint * cosii))
    arglat_final = np.arccos((cosii * sinif - sinii * cosif * cosdraan) / sint)

    return deltav, arglat_init, arglat_final


########################################################################################
# Combined Transfers
########################################################################################


def min_combined(
    rinit: float,
    rfinal: float,
    einit: float,
    efinal: float,
    nuinit: float,
    nufinal: float,
    iinit: float,
    ifinal: float,
    use_optimal: bool = True,
    tol: float = 1e-6,
) -> Tuple[float, float, float, float, float]:
    """Calculates the delta-v and inclination change for the minimum velocity change
    between two non-coplanar orbits.

    References:
        Vallado 2022, p. 355-357, Algorithm 42

    Args:
        rinit (float): Initial position magnitude in km
        rfinal (float): Final position magnitude in km
        einit (float): Eccentricity of the initial orbit
        efinal (float): Eccentricity of the final orbit
        nuinit (float): True anomaly of the initial orbit in radians (0 or pi)
        nufinal (float): True anomaly of the final orbit in radians (0 or pi)
        iinit (float): Initial inclination in radians
        ifinal (float): Final inclination in radians
        use_optimal (bool): Use iterative optimization for inclination change
                            (default True)
        tol (float): Tolerance for inclination iteration (default 1e-6)

    Returns:
        tuple: (deltai_init, deltai_final, deltava, deltavb, dtsec)
            deltai_init (float): Inclination change at point A in radians
            deltai_final (float): Inclination change at point B in radians
            deltava (float): Delta-v at point A in km/s
            deltavb (float): Delta-v at point B in km/s
            dtsec (float): Time of flight for the transfer in seconds
    """
    # Compute semi-major axes
    a1 = (rinit * (1 + einit * np.cos(nuinit))) / (1 - einit**2)
    a2 = 0.5 * (rinit + rfinal)
    a3 = (rfinal * (1 + efinal * np.cos(nufinal))) / (1 - efinal**2)

    # Compute velocities
    vinit = utils.velocity_mag(rinit, a1)
    v1t = utils.velocity_mag(rinit, a2)
    vfinal = utils.velocity_mag(rfinal, a3)
    v3t = utils.velocity_mag(rfinal, a2)

    # Delta inclination
    tdi = ifinal - iinit

    # Delta-Vs
    temp = (1 / tdi) * np.arctan(np.sin(tdi) / ((rfinal / rinit) ** 1.5 + np.cos(tdi)))
    deltava = utils.deltav(v1t, vinit, temp * tdi)
    deltavb = utils.deltav(v3t, vfinal, tdi * (1 - temp))

    # Inclination change
    deltai_init = temp * tdi
    deltai_final = tdi * (1 - temp)

    # Compute transfer time of flight
    dtsec = utils.time_of_flight(a2)

    if not use_optimal:
        return deltai_init, deltai_final, deltava, deltavb, dtsec

    # Iterative optimization
    deltai_final_iter, n_iter = 100, 0
    while abs(deltai_init - deltai_final_iter) > tol:
        deltai_final_iter = deltai_init
        deltava = utils.deltav(v1t, vinit, deltai_final_iter)
        deltavb = utils.deltav(v3t, vfinal, tdi - deltai_final_iter)
        deltai_init = np.arcsin(
            (deltava * vfinal * v3t * np.sin(tdi - deltai_final_iter))
            / (vinit * v1t * deltavb)
        )
        n_iter += 1

    return deltai_init, tdi - deltai_init, deltava, deltavb, dtsec


def combined(
    rinit: float, rfinal: float, einit: float, nuinit: float, deltai: float
) -> Tuple[float, float, float, float, float, float, float]:
    """Calculates the delta-v for a Hohmann transfer between two orbits, considering
    inclination changes. The final orbit is assumed to be circular.

    References:
        Vallado 2007, p. 360-361, Example 6-7

    Args:
        rinit (float): Initial position magnitude in km
        rfinal (float): Final position magnitude in km
        einit (float): Eccentricity of the initial orbit
        nuinit (float): True anomaly of the initial orbit in radians (0 or pi)
        deltai (float): Inclination change in radians (final - initial)

    Returns:
        tuple: (deltai1, deltai2, deltava, deltavb, dtsec, gama, gamb)
            deltai1 (float): Inclination change at point A in radians
            deltai2 (float): Inclination change at point B in radians
            deltava (float): Delta-v at point A in km/s
            deltavb (float): Delta-v at point B in km/s
            dtsec (float): Time of flight for the transfer in seconds
            gama (float): Firing angle at point A in radians
            gamb (float): Firing angle at point B in radians

    TODO:
        - Support non-circular final orbits?
    """
    # Semi-major axes
    ainit = utils.semimajor_axis(rinit, einit, nuinit)
    atran = (rinit + rfinal) * 0.5

    # Velocities
    vinit = utils.velocity_mag(rinit, ainit)
    vtransa = utils.velocity_mag(rinit, atran)
    vfinal = np.sqrt(const.MU / rfinal)  # assumes circular orbit
    vtransb = utils.velocity_mag(rfinal, atran)

    # Proportions of inclination change
    ratio = rfinal / rinit
    s = (1 / deltai) * np.arctan(np.sin(deltai) / (ratio**1.5 + np.cos(deltai)))
    deltai1 = s * deltai
    deltai2 = (1 - s) * deltai

    # Delta-v calculations
    deltava = utils.deltav(vinit, vtransa, deltai1)
    deltavb = utils.deltav(vfinal, vtransb, deltai2)

    # Time of flight for the transfer
    dtsec = utils.time_of_flight(atran)

    # Firing angles
    gama = np.arccos(-(vinit**2 + deltava**2 - vtransa**2) / (2 * vinit * deltava))
    gamb = np.arccos(-(vtransb**2 + deltavb**2 - vfinal**2) / (2 * vtransb * deltavb))

    return deltai1, deltai2, deltava, deltavb, dtsec, gama, gamb


########################################################################################
# Rendezvous
########################################################################################


def rendezvous_coplanar(
    rcsint: float,
    rcstgt: float,
    phasei: float,
    einit: float,
    efinal: float,
    nuinit: float,
    nufinal: float,
    kint: int,
    ktgt: int,
    tol: float = 1e-6,
) -> Tuple[float, float, float]:
    """Calculates parameters for a Hohmann transfer rendezvous.

    References:
        Vallado 2022, pp. 361-367, Algorithms 44 and 45

    Args:
        rcsint (float): Radius of the interceptor's circular orbit in km
        rcstgt (float): Radius of the target's circular orbit in km
        phasei (float): Initial phase angle (target - interceptor) in radians
        einit (float): Eccentricity of the initial orbit
        efinal (float): Eccentricity of the final orbit
        nuinit (float): True anomaly of the initial orbit in radians (0 or pi)
        nufinal (float): True anomaly of the final orbit in radians (0 or pi)
        kint (int): Number of interceptor orbits
        ktgt (int): Number of target orbits to wait
        tol (float): Tolerance for checking if orbit is the same (default 1e-6 rad/s)

    Returns:
        tuple: (phasef, waittime, deltav)
            phasef (float): Final phase angle in radians
            waittime (float): Wait time before next intercept opportunity in seconds
            deltav (float): Total change in velocity in km/s
    """
    # Angular velocities
    angvelint = utils.angular_velocity(rcsint)
    angveltgt = utils.angular_velocity(rcstgt)
    vint = np.sqrt(const.MU / rcsint)

    # Same orbit case
    if abs(angvelint - angveltgt) < tol:
        periodtrans = (ktgt * const.TWOPI + phasei) / angveltgt
        atrans = (const.MU * (periodtrans / (const.TWOPI * kint)) ** 2) ** (1 / 3)

        # Check for intersection with Earth
        rp = 2 * atrans - rcsint
        if rp < 1:
            raise ValueError("Error: the transfer orbit intersects the Earth.")

        # Calculate delta-V
        vtrans = np.sqrt((2 * const.MU / rcsint) - (const.MU / atrans))
        deltav = 2 * (vtrans - vint)

        # Calculate final phase angle and wait time
        phasef = phasei
        waittime = periodtrans

    else:
        # Different orbits
        atrans = (rcsint + rcstgt) / 2
        dttutrans = np.pi * np.sqrt(atrans**3 / const.MU)

        # Calculate final phase angle
        leadang = angveltgt * dttutrans
        phasef = np.pi - leadang
        if phasef < 0:
            phasef += np.pi

        # Calculate wait time
        waittime = (phasef - phasei + const.TWOPI * ktgt) / (angvelint - angveltgt)

        # Semi-major axes
        a1 = utils.semimajor_axis(rcsint, einit, nuinit)
        a2 = (rcsint + rcstgt) / 2
        a3 = utils.semimajor_axis(rcstgt, efinal, nufinal)

        # Delta-V at point A
        vinit = utils.velocity_mag(rcsint, a1)
        vtransa = utils.velocity_mag(rcsint, a2)
        deltava = abs(vtransa - vinit)

        # Delta-V at point B
        vfinal = utils.velocity_mag(rcstgt, a3)
        vtransb = utils.velocity_mag(rcstgt, a2)
        deltavb = abs(vfinal - vtransb)

        # Total delta-V
        deltav = deltava + deltavb

    return phasef, waittime, deltav


def rendezvous_noncoplanar(
    phasei: float,
    aint: float,
    atgt: float,
    kint: int,
    ktgt: int,
    nodeint: float,
    truelon: float,
    deltai: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculates parameters for a non-coplanar Hohmann transfer maneuver.

    References:
        Vallado 2022, pp. 369-373, Algorithm 46

    Args:
        phasei (float): Initial phase angle in radians
        aint (float): Semi-major axis of interceptor orbit in km
        atgt (float): Semi-major axis of target orbit in km
        kint (int): Number of interceptor orbits
        ktgt (int): Number of target orbits
        nodeint (float): Longitude of the ascending node of the interceptor in radians
        truelon (float): True longitude of the target in radians
        deltai (float): Change in inclination in radians (final - initial)

    Returns:
        tuple: (ttrans, tphase, dvphase, dvtrans1, dvtrans2, aphase)
            ttrans (float): Transfer time in seconds
            tphase (float): Phase time in seconds
            dvphase (float): Delta-v for phasing in km/s
            dvtrans1 (float): Delta-v for first transfer in km/s
            dvtrans2 (float): Delta-v for second transfer in km/s
            aphase (float): Semi-major axis for phasing orbit in km/s
    """
    # Angular velocities
    angvelint = utils.angular_velocity(aint)
    angveltgt = utils.angular_velocity(atgt)

    # Calculate transfer time
    atrans = (aint + atgt) / 2
    ttrans = utils.time_of_flight(atrans)

    # Calculate phase time
    deltatnode = phasei / angvelint
    lead = angveltgt * ttrans
    omeganode = angveltgt * deltatnode
    phasenew = nodeint + np.pi - (truelon + omeganode)
    leadnew = np.pi + phasenew
    tphase = (leadnew - lead + const.TWOPI * ktgt) / angveltgt

    # Semi-major axis of phasing orbit
    aphase = (const.MU * (tphase / (const.TWOPI * kint)) ** 2) ** (1 / 3)

    # Calculate phasing delta-V
    vphase = utils.velocity_mag(aint, aphase)
    dvphase = vphase - np.sqrt(const.MU / aint)

    # Calculate delta-V for first transfer
    vtrans1 = utils.velocity_mag(aint, atrans)
    dvtrans1 = vtrans1 - vphase

    # Calculate delta-V for second transfer
    vtrans2 = utils.velocity_mag(atgt, atrans)
    vtgt = np.sqrt(const.MU / atgt)
    dvtrans2 = utils.deltav(vtgt, vtrans2, deltai)

    return ttrans, tphase, dvphase, dvtrans1, dvtrans2, aphase


########################################################################################
# Low Thrust
########################################################################################


def low_thrust(
    ainit: float,
    afinal: float,
    iinit: float,
    mdot: float,
    accelinit: float,
    lambda_: float,
    steps: int = 360,
) -> Tuple[float, float]:
    """Calculate the delta-V and time of flight for a low-thrust, circular transfer.

    References:
        Vallado 2022, p. 382-392, Algorithm 47

    Args:
        ainit (float): Initial semi-major axis in km
        afinal (float): Final semi-major axis in km
        iinit (float): Initial inclination in radians
        mdot (float): Mass flow rate in kg/s
        accelinit (float): Initial acceleration in km/s^2
        lambda_ (float): Control parameter for optimization
        steps (int, optional): Number of steps per orbit (defaults to 360)

    Returns:
        tuple: (deltav, tof)
            deltav (float): Total change in velocity in km/s
            tof (float): Time of flight in seconds
    """
    # Calculate ratio and transition SMA
    ratio = afinal / ainit
    atrans = (afinal + ainit) * 0.5  # km

    # Set up canonical units
    du = ainit
    duptu = np.sqrt(const.MU / ainit)
    atrans /= du  # DU

    # Initial state vector (assuming circular orbit)
    r1 = np.array([ainit, 0, 0])
    magr = np.linalg.norm(r1)
    v = np.sqrt(const.MU / magr)
    v1 = np.array([0, v * np.cos(iinit), v * np.sin(iinit)])

    # Get step size
    dtsec = utils.period(magr) / steps  # steps per orbit

    # Initialize variables
    # cosf is the angle from the node (perigee) - start at 0
    a1 = 1  # normalized to initial radius DU
    tof = rev = cosfold = 0
    nodevec = unit(r1)

    while a1 < (afinal / ainit):
        # Calculate steering angle and acceleration
        cv = 1 / (4 * lambda_**2 * a1**2 + 1)
        cosf = utils.lowuz(cv)
        steering = np.arctan(cosf / np.sqrt(1 / (1 / cv - 1)))
        accel = accelinit / (1 + mdot * tof)  # km/s^2

        # Calculate acceleration vector
        avec = accel * np.array([0, np.cos(steering), np.sin(steering)])

        # Convert to NTW frame
        rntw1, vntw1, tmntw = rv2ntw(r1, v1)
        acc1 = np.dot(tmntw.T, avec)  # NTW to ECI
        vtot = v1 + acc1 * dtsec

        # Propagate orbit to next step
        r1n, v1n = kepler(r1, vtot, dtsec)

        # Calculate the new node angle
        magr, magv = np.linalg.norm(r1n), np.linalg.norm(v1n)
        a1 = utils.specific_mech_energy((magv**2 * 0.5) - (const.MU / magr))
        cosf = np.dot(r1n, nodevec) / magr
        if cosf * cosfold < 0:
            rev += 0.5

        # Update variables
        tof += dtsec
        r1, v1 = r1n, v1n
        dtsec = utils.period(a1) / steps
        cosfold = cosf
        a1 /= du

    # Calculate delta-V
    deltav = (1 - np.sqrt(1 / ratio)) * duptu  # km/s

    return deltav, tof
