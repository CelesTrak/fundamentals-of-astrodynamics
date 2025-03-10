# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 6 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import ellipeinc

from ... import constants as const
from ...mathtime.vector import rot1, rot2, rot3, angle, unit
from ..time.data import IAU80Array
from ..time.frame_conversions import ecef2eci, eci2ecef
from .newton import newtonnu, newtonm, newtone
from . import utils


########################################################################################
# Spherical Elements
########################################################################################


def adbar2rv(
    rmag: float, vmag: float, rtasc: float, decl: float, fpav: float, az: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Conversion from spherical elements to ECI position & velocity vectors.

    References:
        Vallado: 2001, XX
        Chobotov: 2006, p. 70

    Args:
        rmag (float): ECI position vector magnitude in km
        vmag: (float): ECI velocity vector magnitude in km/sec
        rtasc (float): Right ascension of satellite in radians
        decl (float): Declination of satellite in radians
        fpav (float): Satellite flight path angle from vertical in radians
        az (float): Satellite flight path azimuth in radians

    Returns:
        tuple: (r, v)
            r (np.ndarray): ECI position vector
            v (np.ndarray): ECI velocity vector
    """
    # Form position vector
    r = np.array(
        [
            rmag * np.cos(decl) * np.cos(rtasc),
            rmag * np.cos(decl) * np.sin(rtasc),
            rmag * np.sin(decl),
        ]
    )

    # Form velocity vector
    v = np.array(
        [
            vmag
            * (
                np.cos(rtasc)
                * (
                    -np.cos(az) * np.sin(fpav) * np.sin(decl)
                    + np.cos(fpav) * np.cos(decl)
                )
                - np.sin(az) * np.sin(fpav) * np.sin(rtasc)
            ),
            vmag
            * (
                np.sin(rtasc)
                * (
                    -np.cos(az) * np.sin(fpav) * np.sin(decl)
                    + np.cos(fpav) * np.cos(decl)
                )
                + np.sin(az) * np.sin(fpav) * np.cos(rtasc)
            ),
            vmag
            * (np.cos(az) * np.cos(decl) * np.sin(fpav) + np.cos(fpav) * np.sin(decl)),
        ]
    )

    return r, v


def rv2adbar(
    r: ArrayLike, v: ArrayLike
) -> Tuple[float, float, float, float, float, float]:
    """Conversion from position & velocity vectors to spherical elements.

    References:
        Vallado: 2001, xx
        Chobotov: 2006, p. 70

    Args:
        r (array_like): ECI position vector
        v (array_like): ECI velocity vector

    Returns:
        tuple: (rmag, vmag, rtasc, decl, fpav, az)
            rmag (float): ECI position vector magnitude in km
            vmag: (float): ECI velocity vector magnitude in km/sec
            rtasc (float): Right ascension of satellite in radians
            decl (float): Declination of satellite in radians
            fpav (float): Satellite flight path angle from vertical in radians
            az (float): Satellite flight path azimuth in radians
    """
    rmag, vmag = np.linalg.norm(r), np.linalg.norm(v)
    rtemp = np.sqrt(r[0] ** 2 + r[1] ** 2)
    vtemp = np.sqrt(v[0] ** 2 + v[1] ** 2)

    # Right ascension of sateillite
    if rtemp < const.SMALL:
        rtasc = np.arctan2(v[1], v[0]) if vtemp > const.SMALL else 0
    else:
        rtasc = np.arctan2(r[1], r[0])

    # Declination of satellite
    decl = np.arcsin(r[2] / rmag)

    # Flight path angle from vertical
    h = np.cross(r, v)
    fpav = np.arctan2(np.linalg.norm(h), np.dot(r, v))

    # Flight path azimuth
    hcrossr = np.cross(h, r)
    az = np.arctan2(r[0] * hcrossr[1] - r[1] * hcrossr[0], hcrossr[2] * rmag)

    return rmag, vmag, rtasc, decl, fpav, az


########################################################################################
# Classical Elements
########################################################################################


def coe2rv(
    p: float,
    ecc: float,
    incl: float,
    raan: float,
    argp: float,
    nu: float = 0.0,
    arglat: float = 0.0,
    truelon: float = 0.0,
    lonper: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert from classical elements to position & velocity vectors.

    References:
        Vallado: 2022, p. 120-121, Algorithm 10
        Chobotov: 2006, p. 70

    Args:
        p (float): Semi-latus rectum of the orbit in km
        ecc (float): Eccentricity of the orbit
        incl (float): Inclination of the orbit in radians
        raan (float): Right ascension of the ascending node (RAAN) in radians
        argp (float): Argument of perigee in radians
        nu (float, optional): True anomaly in radians
        arglat (float, optional): Argument of latitude in radians
        truelon (float, optional): True longitude in radians
        lonper (float, optional): Longitude of periapsis in radians

    Returns:
        tuple: (r, v)
            r (np.ndarray): ECI position vector in km
            v (np.ndarray): ECI velocity vector in km/s
    """
    # Handle special cases for orbit type
    if ecc < const.SMALL:
        if utils.is_equatorial(incl):
            # Circular equatorial
            argp, raan = 0, 0
            nu = truelon
        else:
            # Circular inclined
            argp = 0
            nu = arglat
    else:
        # Elliptical equatorial
        if utils.is_equatorial(incl):
            argp = lonper
            raan = 0

    # Compute position and velocity in the perifocal coordinate system
    cosnu, sinnu = np.cos(nu), np.sin(nu)
    r_pqw = np.array([cosnu, sinnu, 0]) * (p / (1 + ecc * cosnu))
    v_pqw = np.array([-sinnu, ecc + cosnu, 0]) * (np.sqrt(const.MU / p))

    # Transform from PQW to IJK (GEC)
    r = rot3(rot1(rot3(r_pqw, -argp), -incl), -raan)
    v = rot3(rot1(rot3(v_pqw, -argp), -incl), -raan)

    return r, v


def rv2coe(
    r: ArrayLike, v: ArrayLike
) -> Tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    utils.OrbitType | None,
]:
    """Converts position and velocity vectors into classical orbital elements.

    Args:
        r (array_like): Position vector in km
        v (array_like): Velocity vector in km/s

    References:
        Vallado: 2022, p. 115-116, Algorithm 9

    Returns:
        tuple: (p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper,
                orbit_type)
            p (float): Semilatus rectum in km
            a (float): Semimajor axis in km
            ecc (float): Eccentricity
            incl (float): Inclination in radians (0 to 2pi)
            raan (float): Right ascension of the ascending node in radians
                          (0 to 2pi)
            argp (float): Argument of perigee in radians (0 to 2pi)
            nu (float): True anomaly in radians (0 to 2pi)
            m (float): Mean anomaly in radians (0 to 2pi)
            arglat (float): Argument of latitude in radians (0 to 2pi)
            truelon (float): True longitude in radians (0 to 2pi)
            lonper (float): Longitude of periapsis in radians (0 to 2pi)
            orbit_type (enum): Type of orbit as defined in the OrbitType enum
    """

    def adjust_angle(ang):
        """Adjust angle by subtracting it from 2pi"""
        return const.TWOPI - ang

    # Initialize variables
    p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper = (np.nan,) * 11
    orbit_type = None

    # Make sure position and velocity vectors are numpy arrays
    r, v = np.array(r), np.array(v)

    # Get magnitude of position and velocity vectors
    r_mag, v_mag = np.linalg.norm(r), np.linalg.norm(v)

    # Get angular momentum
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    # Elements are undefined for negative angular momentum
    if h_mag < 0:
        return p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper, orbit_type

    # Define line of nodes vector
    n_vec = np.array([-h[1], h[0], 0])
    n_mag = np.linalg.norm(n_vec)

    # Get eccentricity vector
    e_vec = ((v_mag**2 - const.MU / r_mag) * r - np.dot(r, v) * v) / const.MU
    ecc = np.linalg.norm(e_vec)

    # find a, e, and p (semi-latus rectum)
    sme = (v_mag**2 / 2) - (const.MU / r_mag)
    a = -const.MU / (2 * sme) if abs(sme) > const.SMALL else np.inf

    # Semi-latus rectum
    p = h_mag**2 / const.MU

    # Find inclination
    incl = np.arccos(h[2] / h_mag)

    # Determine orbit type
    orbit_type = utils.determine_orbit_type(ecc, incl, tol=const.SMALL)

    # Find right ascension of ascending node
    if n_mag > const.SMALL:
        raan = np.arccos(np.clip(n_vec[0] / n_mag, -1, 1))
        if n_vec[1] < 0:
            raan = adjust_angle(raan)

    # Find argument of periapsis
    if orbit_type is utils.OrbitType.EPH_INCLINED:
        argp = angle(n_vec, e_vec)
        if e_vec[2] < 0:
            argp = adjust_angle(argp)

    # Find true anomaly at epoch
    if orbit_type in [utils.OrbitType.EPH_INCLINED, utils.OrbitType.EPH_EQUATORIAL]:
        nu = angle(e_vec, r)
        if np.dot(r, v) < 0:
            nu = adjust_angle(nu)

    # Find argument of latitude (inclined cases)
    if orbit_type in [utils.OrbitType.CIR_INCLINED, utils.OrbitType.EPH_INCLINED]:
        arglat = angle(n_vec, r)
        if r[2] < 0:
            arglat = adjust_angle(arglat)

    # Find longitude of periapsis
    if ecc > const.SMALL and orbit_type is utils.OrbitType.EPH_EQUATORIAL:
        lonper = np.arccos(np.clip(e_vec[0] / ecc, -1, 1))
        if e_vec[1] < 0:
            lonper = adjust_angle(lonper)
        if incl > np.pi / 2:
            lonper = adjust_angle(lonper)

    # Find true longitude
    if r_mag > const.SMALL and orbit_type is utils.OrbitType.CIR_EQUATORIAL:
        truelon = np.arccos(np.clip(r[0] / r_mag, -1, 1))
        if r[1] < 0:
            truelon = adjust_angle(truelon)
        if incl > np.pi / 2:
            truelon = adjust_angle(truelon)

    # Find mean anomaly for eccentric orbits
    if np.isreal(ecc):
        _, m = newtonnu(ecc, nu)

    return p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper, orbit_type


########################################################################################
# Equinoctial Elements
########################################################################################


def eq2rv(
    a: float, af: float, ag: float, chi: float, psi: float, meanlon: float, fr: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert from equinoctial elements to position & velocity vectors.

    This function finds the position and velocity vectors in geocentric equatorial (ijk)
    system given the equinoctial orbit elements.

    References:
        Vallado: 2022, p. 110-111

    Args:
        a (float): Semi-major axis in km
        af (float): Component of eccentricity vector (also called k)
        ag (float): Component of eccentricity vector (also called h)
        chi (float): Component of node vector in eqw (also called p)
        psi (float): Component of node vector in eqw (also called q)
        meanlon (float): Mean longitude in radians
        fr (int): Retrograde factor (+1 for prograde, -1 for retrograde)

    Returns:
        tuple: (r, v)
            r (np.ndarray): ECI position vector in km
            v (np.ndarray): ECI velocity vector in km/s

    TODO:
        - Add vector option for conversion
    """
    # Initialize variables
    arglat, truelon, lonper = (0,) * 3

    # Compute eccentricity
    ecc = np.sqrt(af**2 + ag**2)
    p = a * (1 - ecc**2)
    incl = np.pi * ((1 - fr) * 0.5) + 2 * fr * np.arctan(np.sqrt(chi**2 + psi**2))
    raan = np.arctan2(chi, psi)
    argp = np.arctan2(ag, af) - fr * raan

    if ecc < const.SMALL:
        # Circular orbits
        if utils.is_equatorial(incl):
            # Circular equatorial
            truelon = raan
            raan = 0
        else:
            # Circular inclined
            arglat = argp
    else:
        # Elliptical equatorial
        if utils.is_equatorial(incl):
            lonper = argp
            raan = 0

    # Mean anomaly
    m = meanlon - fr * raan - argp
    m = np.mod(m, const.TWOPI)

    # Solve for eccentric anomaly and true anomaly
    e0, nu = newtonm(ecc, m)

    if ecc < const.SMALL:
        # Circular orbits
        if utils.is_equatorial(incl):
            # Circular equatorial
            truelon = nu
        else:
            # Circular inclined
            arglat = nu - fr * raan

    # Convert back to position and velocity vectors
    return coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper)


def rv2eq(
    r: ArrayLike, v: ArrayLike
) -> Tuple[float, float, float, float, float, float, float, float, int]:
    """Convert from position & velocity vectors to equinoctial elements.

    References:
        Vallado: 2022, p. 110-111
        Chobotov: 2006, p. 30

    Args:
        r (array_like): ECI position vector in km
        v (array_like): ECI velocity vector in km/s

    Returns:
        tuple: (a, n, af, ag, chi, psi, meanlon, truelon, fr)
            a (float): Semi-major axis in km
            n (float): Mean motion in rad/s
            af (float): Component of eccentricity vector
            ag (float): Component of eccentricity vector
            chi (float): Component of node vector in eqw
            psi (float): Component of node vector in eqw
            meanlon (float): Mean longitude in radians
            meanlonNu (float): True longitude in radians
            fr (int): Retrograde factor (+1 for prograde, -1 for retrograde)
    """
    # Convert to classical orbital elements
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper, _ = rv2coe(r, v)

    # Determine retrograde factor
    fr = -1 if abs(incl - np.pi) < const.SMALL else 1

    if ecc < const.SMALL:
        # Circular orbits
        if utils.is_equatorial(incl):
            # Circular Equatorial
            argp, omega = 0, 0
            nu, m = truelon, truelon
        else:
            # Circular inclined
            argp = 0
            nu, m = arglat
    else:
        # Elliptical equatorial
        if utils.is_equatorial(incl):
            argp = lonper
            omega = 0

    # Calculate mean motion
    # TODO: put in separate utility function
    n = np.sqrt(const.MU / (a**3))

    # Get eccentricity vector components
    af = ecc * np.cos(fr * omega + argp)
    ag = ecc * np.sin(fr * omega + argp)

    # Get EQW node vector components
    if fr > 0:
        chi = np.tan(incl * 0.5) * np.sin(omega)
        psi = np.tan(incl * 0.5) * np.cos(omega)
    else:
        chi = 1 / np.tan(incl * 0.5) * np.sin(omega)
        psi = 1 / np.tan(incl * 0.5) * np.cos(omega)

    # Determine mean longitude
    meanlon = fr * omega + argp + m
    meanlon = np.mod(meanlon + const.TWOPI, const.TWOPI)

    # Determine true longitude
    truelon = fr * omega + argp + nu
    truelon = np.mod(truelon + const.TWOPI, const.TWOPI)

    return a, n, af, ag, chi, psi, meanlon, truelon, fr


########################################################################################
# Topocentric Elements
########################################################################################


def tradec2rv(
    trr: float,
    trtasc: float,
    tdecl: float,
    dtrr: float,
    tdrtasc: float,
    tddecl: float,
    rseci: ArrayLike,
    vseci: ArrayLike,
) -> Tuple[np.ndarray, np.ndarray]:
    """Converts topocentric coordinates (range, right ascension, declination, and their
    rates) into geocentric equatorial (ECI) position and velocity vectors.

    References:
        Vallado: 2022, p. 257, Algorithm 26

    Args:
        trr (float): Satellite range from site in km
        trtasc (float): Topocentric right ascension in radians
        tdecl (float): Topocentric declination in radians
        dtrr (float): Range rate in km/s
        tdrtasc (float): Topocentric right ascension rate in rad/s
        tddecl (float): Topocentric declination rate in rad/s
        rseci (array_like): ECI site position vector in km
        vseci (array_like): ECI site velocity vector in km/s

    Returns:
        tuple: (reci, veci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
    """
    # Calculate topocentric slant range vectors
    rhov = np.array(
        [
            trr * np.cos(tdecl) * np.cos(trtasc),
            trr * np.cos(tdecl) * np.sin(trtasc),
            trr * np.sin(tdecl),
        ]
    )

    # Slant range rate vectors
    drhov = np.array(
        [
            dtrr * np.cos(tdecl) * np.cos(trtasc)
            - trr * np.sin(tdecl) * np.cos(trtasc) * tddecl
            - trr * np.cos(tdecl) * np.sin(trtasc) * tdrtasc,
            dtrr * np.cos(tdecl) * np.sin(trtasc)
            - trr * np.sin(tdecl) * np.sin(trtasc) * tddecl
            + trr * np.cos(tdecl) * np.cos(trtasc) * tdrtasc,
            dtrr * np.sin(tdecl) + trr * np.cos(tdecl) * tddecl,
        ]
    )

    # ECI position and velocity vectors
    reci = rhov + rseci
    veci = drhov + vseci

    return reci, veci


def rv2tradec(
    reci: ArrayLike, veci: ArrayLike, rseci: ArrayLike, vseci: ArrayLike
) -> Tuple[float, float, float, float, float, float]:
    """Converts geocentric equatorial (ECI) position and velocity vectors into range,
    topocentric right ascension, declination, and rates.

    References:
        Vallado: 2022, p. 257, Algorithm 26

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like)): ECI velocity vector in km/s
        rseci (array_like)): ECI site position vector in km
        vseci (array_like)): ECI site velocity vector in km/s

    Returns:
        tuple: (rho, trtasc, tdecl, drho, dtrtasc, dtdecl)
            rho (float): Satellite range from site in km
            trtasc (float): Topocentric right ascension in radians (0 to 2pi)
            tdecl (float): Topocentric declination in radians (-pi/2 to pi/2)
            drho (float): Range rate in km/s
            dtrtasc (float): Topocentric right ascension rate in rad/s
            dtdecl (float): Topocentric declination rate in rad/s
    """
    # Find ECI slant range vector from site to satellite
    rhoveci = np.array(reci) - np.array(rseci)
    drhoveci = np.array(veci) - np.array(vseci)
    rho = np.linalg.norm(rhoveci)

    # Calculate topocentric right ascension and declination
    temp = np.sqrt(rhoveci[0] ** 2 + rhoveci[1] ** 2)
    if temp < const.SMALL:
        trtasc = np.arctan2(drhoveci[1], drhoveci[0])
    else:
        trtasc = np.arctan2(rhoveci[1], rhoveci[0])

    # Directly over the North Pole
    if temp < const.SMALL:
        tdecl = np.sign(rhoveci[2]) * np.pi / 2  # +- 90 deg
    else:
        tdecl = np.arcsin(rhoveci[2] / rho)

    if trtasc < 0:
        trtasc += const.TWOPI

    # Calculate topocentric right ascension and declination rates
    temp1 = -rhoveci[1] ** 2 - rhoveci[0] ** 2
    drho = np.dot(rhoveci, drhoveci) / rho
    dtrtasc, dtdecl = 0, 0
    if abs(temp1) > const.SMALL:
        dtrtasc = (drhoveci[0] * rhoveci[1] - drhoveci[1] * rhoveci[0]) / temp1
    if abs(temp) > const.SMALL:
        dtdecl = (drhoveci[2] - drho * np.sin(tdecl)) / temp

    return rho, trtasc, tdecl, drho, dtrtasc, dtdecl


########################################################################################
# Flight Elements
########################################################################################


def flt2rv(
    rmag: float,
    vmag: float,
    latgc: float,
    lon: float,
    fpa: float,
    az: float,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Converts flight elements into ECI position and velocity vectors.

    References:
        Vallado: 2022, p. 111-112
        Escobal: 1985, p. 397
        Chobotov: 2006, p. 67

    Args:
        rmag (float): Position vector magnitude in km
        vmag (float): Velocity vector magnitude in km/s
        latgc (float): Geocentric latitude in radians
        lon (float): Longitude in radians
        fpa (float): Flight path angle in radians
        az (float): Flight path azimuth in radians
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (reci, veci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
    """
    # Form position vector
    recef = np.array(
        [
            rmag * np.cos(latgc) * np.cos(lon),
            rmag * np.cos(latgc) * np.sin(lon),
            rmag * np.sin(latgc),
        ]
    )

    # Convert r to ECI
    vecef = np.zeros(3)  # this is a dummy for now
    aecef = np.zeros(3)
    reci, veci, _ = ecef2eci(
        recef, vecef, aecef, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms
    )

    # Calculate right ascension and declination
    if np.sqrt(reci[0] ** 2 + reci[1] ** 2) < const.SMALL:
        rtasc = np.arctan2(veci[1], veci[0])
    else:
        rtasc = np.arctan2(reci[1], reci[0])
    decl = np.arcsin(reci[2] / rmag)

    # Form velocity vector
    fpav = const.HALFPI - fpa
    veci = vmag * np.array(
        [
            # First element
            (
                -np.cos(rtasc)
                * np.sin(decl)
                * (
                    np.cos(az) * np.cos(fpav)
                    - np.sin(rtasc) * np.sin(az) * np.cos(fpav)
                )
                + np.cos(rtasc) * np.sin(decl) * np.sin(fpav)
            ),
            # Second element
            (
                -np.sin(rtasc)
                * np.sin(decl)
                * (
                    np.cos(az) * np.cos(fpav)
                    + np.cos(rtasc) * np.sin(az) * np.cos(fpav)
                )
                + np.sin(rtasc) * np.cos(decl) * np.sin(fpav)
            ),
            # Third element
            (np.sin(decl) * np.sin(fpav) + np.cos(decl) * np.cos(az) * np.cos(fpav)),
        ]
    )

    return reci, veci


def rv2flt(
    reci: ArrayLike,
    veci: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """Transforms a position and velocity vector to flight elements.

    References:
        Vallado: 2022, p. 111-112

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (lon, latgc, rtasc, decl, fpa, az, rmag, vmag)
            lon (float): Longitude in radians
            latgc (float): Geocentric latitude in radians
            rtasc (float): Right ascension angle in radians
            decl (float): Declination angle in radians
            fpa (float): Flight path angle in radians
            az (float): Flight path azimuth in radians
            rmag (float): Position vector magnitude in km
            vmag (float): Velocity vector magnitude in km/s
    """
    # Get magnitude of position and velocity vectors
    rmag, vmag = np.linalg.norm(reci), np.linalg.norm(veci)

    # Convert r to ECEF for lat/lon calculations
    aecef = np.zeros(3)
    recef, vecef, aecef = eci2ecef(
        reci, veci, aecef, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms
    )

    # Calculate longitude
    if np.sqrt(recef[0] ** 2 + recef[1] ** 2) < const.SMALL:
        lon = np.arctan2(vecef[1], vecef[0])
    else:
        lon = np.arctan2(recef[1], recef[0])

    latgc = np.arcsin(recef[2] / rmag)

    # Calculate right ascension and declination
    if np.sqrt(reci[0] ** 2 + reci[1] ** 2) < const.SMALL:
        rtasc = np.arctan2(veci[1], veci[0])
    else:
        rtasc = np.arctan2(reci[1], reci[0])

    decl = np.arcsin(reci[2] / rmag)

    # Calculate flight path angle
    h = np.cross(reci, veci)
    hmag = np.linalg.norm(h)
    rdotv = np.dot(reci, veci)
    fpav = np.arctan2(hmag, rdotv)
    fpa = const.HALFPI - fpav

    # Calculte azimuth
    hcrossr = np.cross(h, reci)
    az = np.arctan2(reci[0] * hcrossr[1] - reci[1] * hcrossr[0], hcrossr[2] * rmag)

    return lon, latgc, rtasc, decl, fpa, az, rmag, vmag


########################################################################################
# Ecliptic Elements
########################################################################################


def ell2rv(
    rr: float, ecllon: float, ecllat: float, drr: float, decllon: float, decllat: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms ecliptic latitude and longitude to position and velocity vectors.

    References:
        Vallado: 2022, p. 265-267

    Args:
        rr (float): Radius of the satellite in km
        ecllon (float): Ecliptic longitude in radians
        ecllat (float): Ecliptic latitude in radians
        drr (float): Radius of the satellite rate in km/s
        decllon (float): Ecliptic longitude rate of change in rad/s
        decllat (float): Ecliptic latitude rate of change in rad/s

    Returns:
        tuple: (reci, veci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
    """
    # Calculate position vector in ecliptic coordinates
    r = np.array(
        [
            rr * np.cos(ecllat) * np.cos(ecllon),
            rr * np.cos(ecllat) * np.sin(ecllon),
            rr * np.sin(ecllat),
        ]
    )

    # Calculate velocity vector in ecliptic coordinates
    v = np.array(
        [
            # X component
            drr * np.cos(ecllat) * np.cos(ecllon)
            - rr * np.sin(ecllat) * np.cos(ecllon) * decllat
            - rr * np.cos(ecllat) * np.sin(ecllon) * decllon,
            # Y component
            drr * np.cos(ecllat) * np.sin(ecllon)
            - rr * np.sin(ecllat) * np.sin(ecllon) * decllat
            + rr * np.cos(ecllat) * np.cos(ecllon) * decllon,
            # Z component
            drr * np.sin(ecllat) + rr * np.cos(ecllat) * decllat,
        ]
    )

    # Rotate position and velocity vectors to the ECI frame
    reci = rot1(r, -const.OBLIQUITYEARTH)
    veci = rot1(v, -const.OBLIQUITYEARTH)

    return reci, veci


def rv2ell(
    reci: ArrayLike, veci: ArrayLike
) -> Tuple[float, float, float, float, float, float]:
    """Transforms position and velocity vectors to ecliptic latitude and longitude.

    References:
        Vallado: 2022, p. 265-267

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s

    Returns:
        tuple: (rr, ecllon, ecllat, drr, decllon, decllat)
            rr (float): Radius of the satellite in km
            ecllon (float): Ecliptic longitude in radians
            ecllat (float): Ecliptic latitude in radians
            drr (float): Radius of the satellite rate in km/s
            decllon (float): Ecliptic longitude rate of change in rad/s
            decllat (float): Ecliptic latitude rate of change in rad/s
    """
    # Perform rotation about the x-axis by the obliquity
    r = rot1(reci, const.OBLIQUITYEARTH)
    v = rot1(veci, const.OBLIQUITYEARTH)

    # Calculate magnitudes
    rr = np.linalg.norm(r)
    temp = np.sqrt(r[0] ** 2 + r[1] ** 2)

    # Calculate ecliptic longitude
    if temp < const.SMALL:
        temp1 = np.sqrt(v[0] ** 2 + v[1] ** 2)
        ecllon = np.arctan2(v[1], v[0]) if abs(temp1) > const.SMALL else 0
    else:
        ecllon = np.arctan2(r[1], r[0])

    # Calculate ecliptic latitude
    ecllat = np.arcsin(r[2] / rr)

    # Calculate rates
    temp1 = -r[1] ** 2 - r[0] ** 2  # Different now
    drr = np.dot(r, v) / rr
    decllon = (v[0] * r[1] - v[1] * r[0]) / temp1 if abs(temp1) > const.SMALL else 0
    decllat = (v[2] - drr * np.sin(ecllat)) / temp if abs(temp) > const.SMALL else 0

    return rr, ecllon, ecllat, drr, decllon, decllat


########################################################################################
# Celestial Elements
########################################################################################


def radec2rv(
    rr: float, rtasc: float, decl: float, drr: float, drtasc: float, ddecl: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms celestial (right ascension and declination) elements to position and
    velocity vectors.

    References:
        Vallado: 2022, p. 254-256, Algorithm 25

    Args:
        rr (float): Radius of the satellite in km
        rtasc (float): Right ascension in radians
        decl (float): Declination in radians
        drr (float): Radius of the satellite rate in km/s
        drtasc (float): Right ascension rate in rad/s
        ddecl (float): Declination rate in rad/s

    Returns:
        tuple: (r, v)
            r (np.ndarray): ECI position vector in km
            v (np.ndarray): ECI velocity vector in km/s
    """
    # Position vector
    r = np.array(
        [
            rr * np.cos(decl) * np.cos(rtasc),
            rr * np.cos(decl) * np.sin(rtasc),
            rr * np.sin(decl),
        ]
    )

    # Velocity vector
    v = np.array(
        [
            # X component
            drr * np.cos(decl) * np.cos(rtasc)
            - rr * np.sin(decl) * np.cos(rtasc) * ddecl
            - rr * np.cos(decl) * np.sin(rtasc) * drtasc,
            # Y component
            drr * np.cos(decl) * np.sin(rtasc)
            - rr * np.sin(decl) * np.sin(rtasc) * ddecl
            + rr * np.cos(decl) * np.cos(rtasc) * drtasc,
            # Z component
            drr * np.sin(decl) + rr * np.cos(decl) * ddecl,
        ]
    )

    return r, v


def rv2radec(
    r: ArrayLike, v: ArrayLike
) -> Tuple[float, float, float, float, float, float]:
    """Transforms position and velocity vectors to celestial (right ascension and
    declination) elements.

    References:
        Vallado: 2022, p. 254-256, Algorithm 25

    Args:
        r (array_like): ECI position vector in km
        v (array_like): ECI velocity vector in km/s

    Returns:
        tuple: (rr, rtasc, decl, drr, drtasc, ddecl)
            rr (float): Radius of the satellite in km
            rtasc (float): Right ascension in radians
            decl (float): Declination in radians
            drr (float): Radius of the satellite rate in km/s
            drtasc (float): Right ascension rate in rad/s
            ddecl (float): Declination rate in rad/s
    """
    # Calculate the magnitude of the position vector
    rr = np.linalg.norm(r)
    temp = np.sqrt(r[0] ** 2 + r[1] ** 2)

    # Calculate right ascension
    rtasc = np.arctan2(v[1], v[0]) if temp < const.SMALL else np.arctan2(r[1], r[0])
    rtasc += const.TWOPI if rtasc < 0 else rtasc

    # Calculate declination
    decl = np.arcsin(r[2] / rr)

    # Calculate radius rate
    drr = np.dot(r, v) / rr
    temp1 = -r[1] * r[1] - r[0] * r[0]

    # Calculate right ascension rate
    drtasc = (v[0] * r[1] - v[1] * r[0]) / temp1 if abs(temp1) > const.SMALL else 0

    # Calculate declination rate
    ddecl = (v[2] - drr * np.sin(decl)) / temp if abs(temp) > const.SMALL else 0

    return rr, rtasc, decl, drr, drtasc, ddecl


########################################################################################
# Azimuth-Elevation Elements
########################################################################################


def razel2rvsez(
    rho: float, az: float, el: float, drho: float, daz: float, del_el: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Converts range, azimuth, and elevation values with slant range and velocity
    vectors for a satellite from a radar site in the topocentric horizon (SEZ) system.

    References:
        Vallado: 2022, p. 258-259, Eqs. 4-4 and 4-5

    Args:
        rho (float): Satellite range from site in km
        az (float): Azimuth in radians (0 to 2pi)
        el (float): Elevation in radians (-pi/2 to pi/2)
        drho (float): Range rate in km/s
        daz (float): Azimuth rate in rad/s
        del_el (float): Elevation rate in rad/s

    Returns:
        tuple: (rhosez, drhosez)
            rhosez (np.ndarray): SEZ range vector in km
            drhosez (np.ndarray): SEZ velocity vector in km/s
    """
    # Initialize values
    sinel, cosel = np.sin(el), np.cos(el)
    sinaz, cosaz = np.sin(az), np.cos(az)

    # Form SEZ range vector
    rhosez = np.array([-rho * cosel * cosaz, rho * cosel * sinaz, rho * sinel])

    # Form SEZ velocity vector
    drhosez = np.array(
        [
            -drho * cosel * cosaz + rhosez[2] * del_el * cosaz + rhosez[1] * daz,
            drho * cosel * sinaz - rhosez[2] * del_el * sinaz - rhosez[0] * daz,
            drho * sinel + rho * del_el * cosel,
        ]
    )

    return rhosez, drhosez


def rvsez2razel(
    rhosez: ArrayLike, drhosez: ArrayLike
) -> Tuple[float, float, float, float, float, float]:
    """Transforms SEZ range and velocity vectors to range, azimuth, and elevation values
    and their rates.

    References:
        Vallado: 2022, p. 259-263, Algorithm 27

    Args:
        rhosez (array_like): SEZ range vector in km
        drhosez (array_like): SEZ velocity vector in km/s

    Returns:
        tuple: (rho, az, el, drho, daz, del_el)
            rho (float): Satellite range from site in km
            az (float): Azimuth in radians (0 to 2pi)
            el (float): Elevation in radians (-pi/2 to pi/2)
            drho (float): Range rate in km/s
            daz (float): Azimuth rate in rad/s
            del_el (float): Elevation rate in rad/s
    """
    # Range magnitude
    rho = np.linalg.norm(rhosez)

    # Calculate azimuth
    temp = np.sqrt(rhosez[0] ** 2 + rhosez[1] ** 2)
    if abs(rhosez[1]) < const.SMALL:
        if temp < const.SMALL:
            az = np.arctan2(drhosez[1], -drhosez[0])
        else:
            az = np.pi if rhosez[0] > 0 else 0
    else:
        az = np.arctan2(rhosez[1], -rhosez[0])

    # Calculate elevation
    el = (
        np.sign(rhosez[2]) * const.HALFPI
        if temp < const.SMALL
        else np.arcsin(rhosez[2] / rho)
    )

    # Range rate
    drho = np.dot(rhosez, drhosez) / rho

    # Azimuth rate
    daz = (
        (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) / (temp**2)
        if abs(temp**2) > const.SMALL
        else 0
    )

    # Elevation rate
    del_el = (drhosez[2] - drho * np.sin(el)) / temp if abs(temp) > const.SMALL else 0

    return rho, az, el, drho, daz, del_el


def razel2rv(
    rho: float,
    az: float,
    el: float,
    drho: float,
    daz: float,
    del_el: float,
    latgd: float,
    lon: float,
    alt: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms range, azimuth, elevation, and their rates to the geocentric
    equatorial (ECI) position and velocity vectors.

    References:
        Vallado: 2022, p. 259-263, Algorithm 27

    Args:
        rho (float): Satellite range from site in km
        az (float): Azimuth in radians
        el (float): Elevation in radians
        drho (float): Range rate in km/s
        daz (float): Azimuth rate in rad/s
        del_el (float): Elevation rate in rad/s
        latgd (float): Geodetic latitude of site in radians
        lon (float): Longitude of site in radians
        alt (float): Altitude of site in km

    Returns:
        tuple: (recef, vecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
    """
    # Find SEZ range and velocity vectors
    rhosez, drhosez = razel2rvsez(rho, az, el, drho, daz, del_el)

    # Perform SEZ to ECEF transformation
    rhoecef = rot3(rot2(rhosez, latgd - const.HALFPI), -lon).T
    drhoecef = rot3(rot2(drhosez, latgd - const.HALFPI), -lon).T

    # Find ECEF range and velocity vectors
    rs, vs = utils.site(latgd, lon, alt)
    recef = rhoecef + rs
    vecef = drhoecef

    return recef, vecef


def rv2razel(
    recef: ArrayLike, vecef: ArrayLike, latgd: float, lon: float, alt: float
) -> Tuple[float, float, float, float, float, float]:
    """Transforms ECEF position and velocity vectors to range, azimuth, elevation, and
    their rates.

    The value of `SMALL` can affect the rate term calculations. the solution uses the
    velocity vector to find the singular cases. also, the elevation and azimuth rate
    terms are not observable unless the acceleration vector is available.

    References:
        Vallado: 2022, p. 259-263, Algorithm 27

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        latgd (float): Geodetic latitude of site in radians
        lon (float): Longitude of site in radians
        alt (float): Altitude of site in km

    Returns:
        tuple: (rho, az, el, drho, daz, del_el)
            rho (float): Satellite range from site in km
            az (float): Azimuth in radians (0 to 2pi)
            el (float): Elevation in radians (-pi/2 to pi/2)
            drho (float): Range rate in km/s
            daz (float): Azimuth rate in rad/s
            del_el (float): Elevation rate in rad/s
    """
    # Get site vector in ECEF
    rsecef, _ = utils.site(latgd, lon, alt)

    # Find ECEF range vector from site to satellite
    rhoecef = recef - rsecef
    drhoecef = vecef
    rho = np.linalg.norm(rhoecef)

    # Convert to SEZ for calculations
    tempvec = rot3(rhoecef, lon)
    rhosez = rot2(tempvec, const.HALFPI - latgd)

    tempvec = rot3(drhoecef, lon)
    drhosez = rot2(tempvec, const.HALFPI - latgd)

    # Calculate azimuth and elevation
    temp = np.sqrt(rhosez[0] ** 2 + rhosez[1] ** 2)
    if temp < const.SMALL:
        el = np.sign(rhosez[2]) * const.HALFPI
        az = np.arctan2(drhosez[1], -drhosez[0])
    else:
        magrhosez = np.linalg.norm(rhosez)
        el = np.arcsin(rhosez[2] / magrhosez)
        az = np.arctan2(rhosez[1] / temp, -rhosez[0] / temp)

    # Calculate range, azimuth, and elevation rates
    drho = np.dot(rhosez, drhosez) / rho
    daz = (
        (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) / (temp * temp)
        if temp > const.SMALL
        else 0
    )
    del_el = (drhosez[2] - drho * np.sin(el)) / temp if abs(temp) > const.SMALL else 0

    return rho, az, el, drho, daz, del_el


########################################################################################
# Site Topocentric (SEZ) Elements
########################################################################################


def _sez_transformation_matrix(lat: float, lon: float) -> np.ndarray:
    """Computes the transformation matrix from ECI to SEZ."""
    # Zenith component
    zvec = unit(
        np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    )

    # East component
    kvec = np.array([0, 0, 1])
    evec = unit(np.cross(kvec, zvec))

    # South component
    svec = unit(np.cross(evec, zvec))

    # Transformation matrix
    return np.vstack([svec, evec, zvec])


def rv2sez(
    reci: ArrayLike, veci: ArrayLike, lat: float, lon: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Converts position and velocity vectors into site topocentric (SEZ) coordinates.

    References:
        Vallado: 2022, p. 165-166

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        lat (float): Relative latitude of the second satellite w.r.t. the first
                     in radians
        lon (float): Relative longitude of the second satellite w.r.t the first
                     in radians

    Returns:
        tuple: (rsez, vsez, transmat)
            rsez (np.ndarray): SEZ position vector in km
            vsez (np.ndarray): SEZ velocity vector in km/s
            transmat (np.ndarray): Transformation matrix from ECI to SEZ
    """
    # Transformation matrix from ECI to SEZ
    transmat = _sez_transformation_matrix(lat, lon)

    # Transform position and velocity vectors
    rsez = transmat @ reci
    vsez = transmat @ veci

    return rsez, vsez, transmat


def sez2rv(
    rsez: ArrayLike, vsez: ArrayLike, lat: float, lon: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Converts site topocentric (SEZ) coordinates into ECI position and velocity
    vectors.

    References:
        Vallado: 2022, p. 165-166

    Args:
        rsez (array_like): SEZ position vector in km
        vsez (array_like): SEZ velocity vector in km/s
        lat (float): Relative latitude of the second satellite w.r.t. the first
                     in radians
        lon (float): Relative longitude of the second satellite w.r.t the first
                     in radians

    Returns:
        tuple: (reci, veci, transmat)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            transmat (np.ndarray): Transformation matrix from SEZ to ECI
    """
    # Transformation matrix from SEZ to ECI
    transmat = _sez_transformation_matrix(lat, lon).T

    # Transform back to ECI
    reci = transmat @ rsez
    veci = transmat @ vsez

    return reci, veci, transmat


########################################################################################
# Satellite Coordinate Systems
########################################################################################


def rv2rsw(
    reci: ArrayLike, veci: ArrayLike
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms position and velocity vectors into radial, tangential (in-track), and
    normal (cross-track) coordinates, i.e. RSW frame.

    Note: There are numerous nomenclatures for these systems. This is the RSW system of
    Vallado. The reverse values are found using the transmat transpose.

    References:
        Vallado: 2022, p. 166

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s

    Returns:
        tuple: (rrsw, vrsw, transmat)
            rrsw (np.ndarray): RSW position vector in km
            vrsw (np.ndarray): RSW velocity vector in km/s
            transmat (np.ndarray): Transformation matrix from ECI to RSW
    """
    # Radial component
    rvec = unit(reci)

    # Cross-track component
    wvec = unit(np.cross(reci, veci))

    # Along-track component
    svec = unit(np.cross(wvec, rvec))

    # Assemble transformation matrix from ECI to RSW frame
    transmat = np.vstack([rvec, svec, wvec])

    rrsw = np.dot(transmat, reci)
    vrsw = np.dot(transmat, veci)

    return rrsw, vrsw, transmat


def rv2ntw(
    reci: ArrayLike, veci: ArrayLike
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms position and velocity vectors into normal (in-radial), tangential
    (velocity), and normal (cross-track) coordinates.

    Note: Sometimes the first vector is called "along-radial". the tangential direction
    is always aligned with the velocity vector. This is the NTW system of Vallado.

    References:
        Vallado: 2022, p. 166

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s

    Returns:
        tuple: (rntw, vntw, transmat)
            rntw (np.ndarray): NTW position vector in km
            vntw (np.ndarray): NTW velocity vector in km/s
            transmat (np.ndarray): Transformation matrix from ECI to NTW
    """
    # In-velocity component
    tvec = unit(veci)

    # Cross-track component
    wvec = unit(np.cross(reci, veci))

    # Along-radial component
    nvec = unit(np.cross(tvec, wvec))

    # Assemble transformation matrix from ECI to NTW frame
    transmat = np.array([nvec, tvec, wvec])

    rntw = np.dot(transmat, reci)
    vntw = np.dot(transmat, veci)

    return rntw, vntw, transmat


def rv2pqw(reci: ArrayLike, veci: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
    """Transforms position and velocity vectors into perifocal (PQW) coordinates.

    References:
        Vallado: 2022, p. 166

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s

    Returns:
        tuple: (rpqw, vpqw)
            rpqw (np.ndarray): PQW position vector in km
            vpqw (np.ndarray): PQW velocity vector in km/s
    """
    # Covert to classical elements
    p, _, ecc, _, _, _, nu, *_ = rv2coe(reci, veci)

    # Return if true anomaly is undefined
    if np.isnan(nu):
        return np.zeros(3), np.zeros(3)

    # Form PQW position vector
    sinnu, cosnu = np.sin(nu), np.cos(nu)
    temp = p / (1 + ecc * cosnu)
    rpqw = temp * np.array([cosnu, sinnu, 0])

    # Form PQW velocity vector
    p = const.SMALL if np.abs(p) < const.SMALL else p
    vpqw = np.array(
        [-np.sqrt(const.MU / p) * sinnu, np.sqrt(const.MU / p) * (ecc + cosnu), 0]
    )

    return rpqw, vpqw


########################################################################################
# Geodetic Elements
########################################################################################


def ecef2ll(r: ArrayLike) -> Tuple[float, float, float, float]:
    """Converts an ECEF position vector into geodetic and geocentric latitude,
    longitude, and height above the ellipsoid using the Astronomical Almanac method.

    References:
        Vallado: 2022, p. 174, Algorithm 12

    Args:
        r (array_like): ECEF position vector in km

    Returns:
        tuple: (latgc, latgd, lon, hellp)
            latgc (float): Geocentric latitude in radians (-pi to pi)
            latgd (float): Geodetic latitude in radians (-pi to pi)
            lon (float): Longitude in radians (-2pi to 2pi)
            hellp (float): Height above the ellipsoid in km
    """
    # Compute magnitude of the position vector
    magr = np.linalg.norm(r)

    # Compute longitude (right ascension approximation)
    temp = np.sqrt(r[0] ** 2 + r[1] ** 2)
    if np.abs(temp) < const.SMALL:
        lon = np.sign(r[2]) * np.pi * 0.5
    else:
        lon = np.arctan2(r[1], r[0])

    # Adjust longitude to be within [-2π, 2π] range
    if np.abs(lon) >= np.pi:
        lon += const.TWOPI if lon < 0 else -const.TWOPI

    # Compute geodetic latitude as an initial approximation (declination)
    latgd = np.arcsin(r[2] / magr)

    # Iterate to refine geodetic latitude
    i, c = 0, 0
    olddelta = latgd + 10

    while np.abs(olddelta - latgd) >= const.SMALL and i < 10:
        olddelta = latgd
        sintemp = np.sin(latgd)
        c = const.RE / np.sqrt(1 - const.ECCEARTHSQRD * sintemp**2)
        latgd = np.arctan((r[2] + c * const.ECCEARTHSQRD * sintemp) / temp)
        i += 1

    # Calculate height above the ellipsoid
    if np.pi * 0.5 - np.abs(latgd) < np.radians(1):  # within 1 deg of poles
        hellp = (temp / np.cos(latgd)) - c
    else:
        s = c * (1 - const.ECCEARTHSQRD)
        hellp = r[2] / np.sin(latgd) - s

    # Compute geocentric latitude
    latgc = np.arcsin(r[2] / magr)

    return latgc, latgd, lon, hellp


def ecef2llb(r: ArrayLike) -> Tuple[float, float, float, float]:
    """Converts an ECEF position vector into geodetic and geocentric latitude,
    longitude, and height above the ellipsoid using the Borkowski method.

    References:
        Vallado: 2022, p. 175-176, Algorithm 13

    Args:
        r (array_like): ECEF position vector in km

    Returns:
        tuple: (latgc, latgd, lon, hellp)
            latgc (float): Geocentric latitude in radians (-pi to pi)
            latgd (float): Geodetic latitude in radians (-pi to pi)
            lon (float): Longitude in radians (-2pi to 2pi)
            hellp (float): Height above the ellipsoid in km
    """
    # Constants
    third = 1 / 3

    # Compute magnitude of the position vector
    magr = np.linalg.norm(r)

    # Earth semi-major and semi-minor axes
    a = const.RE
    b = np.sign(r[2]) * a * np.sqrt(1.0 - const.ECCEARTHSQRD)

    # Compute longitude (right ascension approximation)
    temp = np.sqrt(r[0] ** 2 + r[1] ** 2)
    if np.abs(temp) < const.SMALL:
        lon = np.sign(r[2]) * np.pi * 0.5
    else:
        lon = np.arctan2(r[1], r[0])

    # Adjust longitude to be within [-2π, 2π] range
    if np.abs(lon) >= np.pi:
        lon += const.TWOPI if lon < 0 else -const.TWOPI

    # Intermediate variables for polynomial solution
    atemp = 1 / (a * temp)
    e = (b * r[2] - a**2 + b**2) * atemp
    f = (b * r[2] + a**2 - b**2) * atemp
    p = 4 * third * (e * f + 1)
    q = 2 * (e**2 - f**2)
    d = p**3 + q**2

    # Solve the cubic equation based on the discriminant `d`
    if d > 0:
        nu = (np.sqrt(d) - q) ** third - (np.sqrt(d) + q) ** third
    else:
        sqrtp = np.sqrt(-p)
        nu = 2 * sqrtp * np.cos(third * np.arccos(q / (p * sqrtp)))

    # Intermediate variables for latitude and height computations
    g = 0.5 * (np.sqrt(e**2 + nu) + e)
    t = np.sqrt(g**2 + (f - nu * g) / (2 * g - e)) - g

    # Compute geodetic latitude and height above the ellipsoid
    latgd = np.arctan(a * (1 - t**2) / (2 * b * t))
    hellp = (temp - a * t) * np.cos(latgd) + (r[2] - b) * np.sin(latgd)

    # Compute geocentric latitude
    latgc = np.arcsin(r[2] / magr)

    return latgc, latgd, lon, hellp


########################################################################################
# Modified Equidistant Cylindrical (EQCM)
########################################################################################


def _compute_oes(r_rtn: np.ndarray, v_rtn: np.ndarray, ecc_tol: float):
    """Compute orbital elements"""
    # Calculate eccentricity vector and magnitude
    h_vec = np.cross(r_rtn, v_rtn)
    p = np.dot(h_vec, h_vec) / const.MU
    ecc_vec = np.cross(v_rtn, h_vec) / const.MU - r_rtn / np.linalg.norm(r_rtn)
    ecc = np.linalg.norm(ecc_vec)

    # Calculate semimajor axis and perigee unit vector
    a = p / (1 - ecc**2)
    perigee_unit = ecc_vec / ecc if ecc > ecc_tol else r_rtn / np.linalg.norm(r_rtn)
    lambda_perigee = np.arctan2(perigee_unit[1], perigee_unit[0])

    return p, ecc, a, perigee_unit, lambda_perigee


def _compute_rv_target(ecc, a, nu2, perigee_unit):
    """Compute the future position and velocity of the target in RTN2"""
    # Computes the future position and velocity of the target
    r2_tgt = a * (1 - ecc**2) / (1 + ecc * np.cos(nu2))
    p_vec = perigee_unit
    q_vec = np.cross([0, 0, 1], p_vec)

    r2_vec_tgt = r2_tgt * (np.cos(nu2) * p_vec + np.sin(nu2) * q_vec)
    v2_vec_tgt = np.sqrt(const.MU / (a * (1 - ecc**2))) * (
        -np.sin(nu2) * p_vec + (ecc + np.cos(nu2)) * q_vec
    )

    # Convert to RTN2 frame
    *_, rot_rtn1_to_rtn2 = rv2rsw(r2_vec_tgt, v2_vec_tgt)
    r_tgt_rtn2 = rot_rtn1_to_rtn2 @ r2_vec_tgt
    v_tgt_rtn2 = rot_rtn1_to_rtn2 @ v2_vec_tgt

    return r_tgt_rtn2, v_tgt_rtn2, r2_tgt


def _get_sez_rotation(phi, lambda_):
    """Computes the transformation matrix from RTN to SEZ frame."""
    sin_phi, cos_phi = np.sin(phi), np.cos(phi)
    sin_lambda, cos_lambda = np.sin(lambda_), np.cos(lambda_)

    return np.array(
        [
            [sin_phi * cos_lambda, sin_phi * sin_lambda, -cos_phi],
            [-sin_lambda, cos_lambda, 0],
            [cos_phi * cos_lambda, cos_phi * sin_lambda, sin_phi],
        ]
    )


def eci_to_eqcm_rtn(
    r_tgt_eci: ArrayLike,
    v_tgt_eci: ArrayLike,
    r_int_eci: ArrayLike,
    v_int_eci: ArrayLike,
    ecc_tol: float = 1e-7,
) -> tuple[np.ndarray, np.ndarray]:
    """Finds the relative position and velocity vectors in the Modified Equidistant
    Cylindrical (EQCM) frame given the ECI target and interceptor states with RTN (RSW)
    ordered components.

    References:
        Alfano: 2012

    Args:
        r_tgt_eci (array_like): ECI position vector of the target in km
        v_tgt_eci (array_like): ECI velocity vector of the target in km/s
        r_int_eci (array_like): ECI position vector of the interceptor in km
        v_int_eci (array_like): ECI velocity vector of the interceptor in km/s
        ecc_tol (float): Tolerance for eccentricity (defaults to 1e-7)

    Returns:
        tuple: (r_int_eqcm, v_int_eqcm)
            r_int_eqcm (np.ndarray): EQCM position vector of the interceptor in km
            v_int_eqcm (np.ndarray): EQCM velocity vector of the interceptor in km/s
    """
    # Compute rotation matrix from ECI to RTN1 frame for target
    r_tgt_rtn1, v_tgt_rtn1, rot_eci_to_rtn1 = rv2rsw(r_tgt_eci, v_tgt_eci)
    r_int_rtn1, v_int_rtn1 = rot_eci_to_rtn1 @ r_int_eci, rot_eci_to_rtn1 @ v_int_eci

    # Compute magnitudes
    mag_r_tgt, mag_r_int = np.linalg.norm(r_tgt_rtn1), np.linalg.norm(r_int_rtn1)

    # Compute lambda and phi rotation angles to go from target to interceptor
    # (lambda_tgt will be 0)
    sin_phi = r_int_rtn1[2] / mag_r_int
    phi = np.arcsin(sin_phi)
    lambda_ = np.arctan2(r_int_rtn1[1], r_int_rtn1[0])

    # Orbital elements of the target at present (nu1) and future (nu2) locations
    p_tgt, ecc_tgt, a_tgt, perigee_unit, lambda_perigee = _compute_oes(
        r_tgt_rtn1, v_tgt_rtn1, ecc_tol
    )
    nu1 = -lambda_perigee
    nu2 = lambda_ - lambda_perigee

    # Future position and velocity of target
    r_tgt_rtn2, v_tgt_rtn2, r2_tgt = _compute_rv_target(
        ecc_tgt, a_tgt, nu2, perigee_unit
    )

    # Convert interceptor components to SEZ frame
    rot_to_sez = _get_sez_rotation(phi, lambda_)
    r_int_sez = rot_to_sez @ r_int_rtn1
    v_int_sez = rot_to_sez @ v_int_rtn1

    # Position components in EQCM
    r_int_eqcm = np.zeros((3, 1))
    r_int_eqcm[0] = r_int_sez[2] - r_tgt_rtn2[0]

    ea0, _ = newtonnu(ecc_tgt, nu1)
    ea1, _ = newtonnu(ecc_tgt, nu2)

    # Fix quadrants for special cases
    if abs(ea1 - ea0) > np.pi:
        ea0 = ea0 + const.TWOPI if ea0 < 0 else ea0 - const.TWOPI

    # Calculate elliptic integrals of the second kind
    e0 = ellipeinc(ea0, ecc_tgt**2)
    e1 = ellipeinc(ea1, ecc_tgt**2)

    fixit = np.pi if e1 - e0 < 0 else 0
    r_int_eqcm[1] = a_tgt * (e1 - e0 + fixit)  # arclength value
    r_int_eqcm[2] = phi * r2_tgt

    # Velocity components in EQCM
    v_int_eqcm = np.zeros((3, 1))
    lambda_dot = v_int_sez[1] / (mag_r_int * np.cos(phi))
    v_int_eqcm[0] = v_int_sez[2] - v_tgt_rtn2[0]
    v_int_eqcm[1] = lambda_dot * r2_tgt - (v_tgt_rtn1[1] / mag_r_tgt) * mag_r_tgt
    v_int_eqcm[2] = (-v_int_sez[0] / mag_r_int) * r2_tgt

    return r_int_eqcm.ravel(), v_int_eqcm.ravel()


def eqcm_to_eci_rtn(
    r_tgt_eci: ArrayLike,
    v_tgt_eci: ArrayLike,
    r_int_eqcm: ArrayLike,
    v_int_eqcm: ArrayLike,
    ecc_tol: float = 1e-7,
) -> Tuple[np.ndarray, np.ndarray]:
    """Finds the interceptor position and velocity vectors in the ECI (RTN) frame given
    the target ECI state and interceptor Modified Equidistant Cylindrical (EQCM) state.

    References:
        Alfano: 2012

    Args:
        r_tgt_eci (array_like): ECI position vector of the target in km
        v_tgt_eci (array_like): ECI velocity vector of the target in km/s
        r_int_eqcm (array_like): EQCM position vector of the interceptor in km
        v_int_eqcm (array_like): EQCM velocity vector of the interceptor in km/s
        ecc_tol (float): Tolerance for eccentricity (defaults to 1e-7)

    Returns:
        tuple: (r_int_eci, v_int_eci)
            r_int_eci (np.ndarray): ECI position vector of the interceptor in km
            v_int_eci (np.ndarray): ECI velocity vector of the interceptor in km/s
    """
    # Convert target vectors to RTN1 frame
    r_tgt_rtn1, v_tgt_rtn1, rot_eci_to_rtn1 = rv2rsw(r_tgt_eci, v_tgt_eci)
    mag_r_tgt = np.linalg.norm(r_tgt_rtn1)

    # Compute orbital elements of the target at nu1
    p_tgt, ecc_tgt, a_tgt, perigee_unit, lambda_perigee = _compute_oes(
        r_tgt_rtn1, v_tgt_rtn1, ecc_tol
    )
    nu1 = -lambda_perigee

    # Compute nu2 and lambda from orbit arc
    arclength = r_int_eqcm[1]
    ea1, _ = newtonnu(ecc_tgt, nu1)

    # Tolerance and iteration definitions
    arclength_tol = 1e-3  # 1 m tolerance
    n_ea_its = 10

    # Compute eccentric anomaly correction
    if abs(float(arclength)) > arclength_tol:
        delta_ea = utils.inverse_elliptic2(arclength / a_tgt, ecc_tgt**2)
        e1 = ellipeinc(ea1, ecc_tgt**2)
        ea2 = ea1 + delta_ea

        # Iteratively refine eccentric anomaly correction
        for _ in range(n_ea_its):
            e2 = ellipeinc(ea2, ecc_tgt**2)
            arclength1 = a_tgt * (e2 - e1)
            if abs(arclength1 - arclength) < arclength_tol:
                break
            correction = arclength / (ea2 - ea1)
            ea2 -= (arclength1 - arclength) / correction

        _, nu2 = newtone(ecc_tgt, ea2)
    else:
        nu2 = nu1

    lambda_ = nu2 - nu1
    sin_lambda, cos_lambda = np.sin(lambda_), np.cos(lambda_)

    # Compute future target position and velocity at nu2
    r_tgt_rtn2, v_tgt_rtn2, r2_tgt = _compute_rv_target(
        ecc_tgt, a_tgt, nu2, perigee_unit
    )

    # Compute phi and interceptor position unit vector in RTN1
    phi = r_int_eqcm[2] / r2_tgt
    r_int_unit_rtn1 = np.array(
        [cos_lambda * np.cos(phi), sin_lambda * np.cos(phi), np.sin(phi)]
    )

    # Convert to SEZ frame
    rot_to_sez = _get_sez_rotation(phi, lambda_)
    r_int_unit_sez = rot_to_sez @ r_int_unit_rtn1

    # Compute proper scaling from Z component of interceptor
    r_int_sez = np.zeros(3)
    r_int_sez[2] = r_int_eqcm[0] + r_tgt_rtn2[0]
    mag_r_int = r_int_sez[2] / r_int_unit_sez[2]
    r_int_rtn1 = mag_r_int * r_int_unit_rtn1

    # Compute velocity components in RTN1
    lambda_dot = (v_int_eqcm[1] + (v_tgt_rtn1[1] / mag_r_tgt) * mag_r_tgt) / r2_tgt

    v_int_sez = np.array(
        [
            (-v_int_eqcm[2] / r2_tgt) * mag_r_int,
            lambda_dot * mag_r_int * np.cos(phi),
            v_int_eqcm[0] + v_tgt_rtn2[0],
        ]
    )

    v_int_rtn1 = rot_to_sez.T @ v_int_sez

    # Convert all back to original ECI frame
    r_int_eci = rot_eci_to_rtn1.T @ r_int_rtn1
    v_int_eci = rot_eci_to_rtn1.T @ v_int_rtn1

    return r_int_eci, v_int_eci


########################################################################################
# Miscellaneous
########################################################################################


def perifocal_transform(
    i: float, raan: float, w: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transform classical elements from equatorial (ECI) frame to perifocal frame.

    References:
        Vallado: 2022, p. 171, Equation 3-31

    Args:
        i (float): Inclination in radians
        raan (float): Right ascension of the ascending node in radians
        w (float): Argument of periapsis in radians

    Returns:
        tuple: (p_, q_, w_)
            p_ (np.ndarray): P vector in perifocal frame
            q_ (np.ndarray): Q vector in perifocal frame
            w_ (np.ndarray): W vector in perifocal frame
    """
    p_ = np.array(
        [
            np.cos(w) * np.cos(raan) - np.sin(w) * np.sin(raan) * np.cos(i),
            np.cos(w) * np.sin(raan) + np.sin(w) * np.cos(raan) * np.cos(i),
            np.sin(w) * np.sin(i),
        ]
    )

    q_ = np.array(
        [
            -np.sin(w) * np.cos(raan) - np.cos(w) * np.sin(raan) * np.cos(i),
            -np.sin(w) * np.sin(raan) + np.cos(w) * np.cos(raan) * np.cos(i),
            np.cos(w) * np.sin(i),
        ]
    )

    w_ = np.array([np.sin(raan) * np.sin(i), -np.cos(raan) * np.sin(i), np.cos(i)])

    return p_, q_, w_
