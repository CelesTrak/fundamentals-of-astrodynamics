# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 6 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import numpy as np

from ...constants import SMALL, MU, HALFPI, TWOPI, OBLIQUITYEARTH
from ...mathtime.vector import rot1, rot2, rot3, angle, unit
from ..time.frame_conversions import ecef2eci, eci2ecef
from .kepler import OrbitType, determine_orbit_type, newtonnu, newtonm
from .utils import site


###############################################################################
# Local Utility Functions
###############################################################################

def is_equatorial(inc):
    """Equatorial check for inclinations"""
    return inc < SMALL or abs(inc - np.pi) < SMALL


###############################################################################
# Spherical Elements
###############################################################################

def adbar2rv(rmag, vmag, rtasc, decl, fpav, az):
    """Conversion from spherical elements to position & velocity vectors

    This function transforms the orbital elements (rtasc, decl, fpav,
    azimuth, position and velocity magnitude) into ECI position and velocity
    vectors.

    References:
        Vallado: 2001, XX
        Chobotov       70

    Args:
        rmag (float): ECI position vector magnitude in km
        vmag: (float): ECI velocity vector magnitude in km/sec
        rtasc (float): Right ascension of satellite in radians
        decl (float): Declination of satellite in radians
        fpav (float): Satellite flight path angle from vertical in radians
        az (float): Satellite flight path azimuth in radians

    Returns:
        r (numpy array): ECI position vector
        v (numpy array): ECI velocity vector
    """
    # Form position vector
    r = np.array(
        [
            rmag * np.cos(decl) * np.cos(rtasc),
            rmag * np.cos(decl) * np.sin(rtasc),
            rmag * np.sin(decl)
        ]
    )

    # Form velocity vector
    v = np.array(
        [
            vmag * (np.cos(rtasc)
                    * (-np.cos(az) * np.sin(fpav) * np.sin(decl)
                       + np.cos(fpav) * np.cos(decl))
                    - np.sin(az) * np.sin(fpav) * np.sin(rtasc)),
            vmag * (np.sin(rtasc)
                    * (-np.cos(az) * np.sin(fpav) * np.sin(decl)
                       + np.cos(fpav) * np.cos(decl))
                    + np.sin(az) * np.sin(fpav) * np.cos(rtasc)),
            vmag * (np.cos(az) * np.cos(decl) * np.sin(fpav)
                    + np.cos(fpav) * np.sin(decl))
        ]
    )

    return r, v


def rv2adbar(r, v):
    """Conversion from position & velocity vectors to spherical elements

    This function transforms a position and velocity vector into the adbarv
    elements: rtasc, decl, fpav, azimuth, position and velocity magnitude.

    References:
        Vallado: 2001, xx
        Chobotov       70

    Args:
        r (array_like): ECI position vector
        v (array_like): ECI velocity vector

    Returns:
        rmag (float): ECI position vector magnitude in km
        vmag: (float): ECI velocity vector magnitude in km/sec
        rtasc (float): Right ascension of satellite in radians
        decl (float): Declination of satellite in radians
        fpav (float): Satellite flight path angle from vertical in radians
        az (float): Satellite flight path azimuth in radians
    """
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    rtemp = np.sqrt(r[0]**2 + r[1]**2)
    vtemp = np.sqrt(v[0]**2 + v[1]**2)

    # Right ascension of sateillite
    if rtemp < SMALL:
        rtasc = np.arctan2(v[1], v[0]) if vtemp > SMALL else 0
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


###############################################################################
# Classical Elements
###############################################################################

def coe2rv(p, ecc, incl, raan, nu=0, arglat=0, truelon=0, lonper=0):
    """Convert from classical elements to position & velocity vectors.

    This function finds the position and velocity vectors in geocentric
    equatorial (ECI) system given the classical orbit elements.

    References:
        Vallado: 2007, p. 126, Algorithm 10
        Chobotov       70

    Args:
        p (float): Semi-latus rectum of the orbit in km
        ecc (float): Eccentricity of the orbit
        incl (float): Inclination of the orbit in radians
        raan (float): Right ascension of the ascending node (RAAN) in radians
        nu (float): True anomaly in radians
        arglat (float, optional): Argument of latitude in radians
        truelon (float, optional): True longitude in radians
        lonper (float, optional): Longitude of periapsis in radians

    Returns:
        tuple: (r, v)
            r (numpy.ndarray): ECI position vector in km
            v (numpy.ndarray): ECI velocity vector in km/s
    """
    # Handle special cases for orbit type
    if ecc < SMALL:
        if is_equatorial(incl):
            # Circular equatorial
            argp, raan = 0.0, 0.0
            nu = truelon
        else:
            # Circular inclined
            argp = 0.0
            nu = arglat
    else:
        # Elliptical equatorial
        if is_equatorial(incl):
            argp = lonper
            raan = 0.0
        else:
            # Rectilinear orbits
            nu = arglat
            argp = 0.0

    # Compute position and velocity in the perifocal coordinate system
    cosnu = np.cos(nu)
    sinnu = np.sin(nu)
    r_pqw = np.array([cosnu, sinnu, 0.0]) * (p / (1 + ecc * cosnu))
    v_pqw = np.array([-sinnu, ecc + cosnu, 0.0]) * (np.sqrt(MU / p))

    # Transform from PQW to IJK (GEC)
    r = rot3(rot1(rot3(r_pqw, -argp), -incl), -raan)
    v = rot3(rot1(rot3(v_pqw, -argp), -incl), -raan)

    return r, v


def rv2coe(r, v):
    """Converts position and velocity vectors into classical orbital elements.

    Args:
        r (array_like): Position vector in km
        v (array_like): Velocity vector in km/s

    References:
        Vallado: 2007, p. 121, Algorithm 9

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
        return 2 * np.pi - ang

    (p, a, ecc, incl, raan, argp, nu,
     m, arglat, truelon, lonper) = (np.nan,) * 11
    orbit_type = None

    # Make sure position and velocity vectors are numpy arrays
    r = np.array(r)
    v = np.array(v)

    # Get magnitude of position and velocity vectors
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)

    # Get angular momentum
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    # Elements are undefined for negative angular momentum
    if h_mag < 0:
        return (
            p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper,
            orbit_type
        )

    # Define line of nodes vector
    n_vec = np.array([-h[1], h[0], 0])
    n_mag = np.linalg.norm(n_vec)

    # Get eccentricity vector
    e_vec = ((v_mag ** 2 - MU / r_mag) * r - np.dot(r, v) * v) / MU
    ecc = np.linalg.norm(e_vec)

    # find a, e, and p (semi-latus rectum)
    sme = (v_mag ** 2 / 2) - (MU / r_mag)
    if abs(sme) > SMALL:
        a = - MU / (2 * sme)
    else:
        a = np.inf

    # Semi-latus rectum
    p = h_mag ** 2 / MU

    # Find inclination
    incl = np.arccos(h[2] / h_mag)

    # Determine orbit type
    orbit_type = determine_orbit_type(ecc, incl, tol=SMALL)

    # Find right ascension of ascending node
    if n_mag > SMALL:
        raan = np.arccos(np.clip(n_vec[0] / n_mag, -1, 1))
        if n_vec[1] < 0:
            raan = adjust_angle(raan)

    # Find argument of periapsis
    if orbit_type is OrbitType.EPH_INCLINED:
        argp = angle(n_vec, e_vec)
        if e_vec[2] < 0:
            argp = adjust_angle(argp)

    # Find true anomaly at epoch
    if orbit_type in [OrbitType.EPH_INCLINED, OrbitType.EPH_EQUATORIAL]:
        nu = angle(e_vec, r)
        if np.dot(r, v) < 0:
            nu = adjust_angle(nu)

    # Find argument of latitude (inclined cases)
    if orbit_type in [OrbitType.CIR_INCLINED, OrbitType.EPH_INCLINED]:
        arglat = angle(n_vec, r)
        if r[2] < 0:
            arglat = adjust_angle(arglat)

    # Find longitude of periapsis
    if ecc > SMALL and orbit_type is OrbitType.EPH_EQUATORIAL:
        lonper = np.arccos(np.clip(e_vec[0] / ecc, -1, 1))
        if e_vec[1] < 0:
            lonper = adjust_angle(lonper)
        if incl > np.pi / 2:
            lonper = adjust_angle(lonper)

    # Find true longitude
    if r_mag > SMALL and orbit_type is OrbitType.CIR_EQUATORIAL:
        truelon = np.arccos(np.clip(r[0] / r_mag, -1, 1))
        if r[1] < 0:
            truelon = adjust_angle(truelon)
        if incl > np.pi / 2:
            truelon = adjust_angle(truelon)

    # Find mean anomaly for all orbits
    e, m = newtonnu(ecc, nu)

    return (
        p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper, orbit_type
    )


###############################################################################
# Equinoctial Elements
###############################################################################

def eq2rv(a, af, ag, chi, psi, meanlon, fr):
    """Convert from equinoctial elements to position & velocity vectors.

    This function finds the position and velocity vectors in geocentric
    equatorial (ijk) system given the equinoctial orbit elements.

    References:
        Vallado: 2013, p. 108

    Args:
        a (float): Semi-major axis in km
        af (float): Component of eccentricity vector (also called k)
        ag (float): Component of eccentricity vector (also called h)
        chi (float): Component of node vector in eqw (also called p)
        psi (float): Component of node vector in eqw (also called q)
        meanlon (float): Mean longitude in radians
        fr (int): Retrograde factor (+1 for prograde, -1 for retrograde)

    Returns:
        tuple:
            np.array: Position vector in km
            np.array: Velocity vector in km/s

    TODO:
        - Add vector option for conversion
    """
    # Initialize variables
    arglat, truelon, lonper = (0., ) * 3

    # Compute eccentricity
    ecc = np.sqrt(af**2 + ag**2)
    p = a * (1.0 - ecc**2)
    incl = (
        np.pi * ((1.0 - fr) * 0.5)
        + 2.0 * fr * np.arctan(np.sqrt(chi**2 + psi**2))
    )
    omega = np.arctan2(chi, psi)
    argp = np.arctan2(ag, af) - fr * omega

    if ecc < SMALL:
        # Circular orbits
        if is_equatorial(incl):
            # Circular equatorial
            truelon = omega
            omega = 0
        else:
            # Circular inclined
            arglat = argp
    else:
        # Elliptical equatorial
        if is_equatorial(incl):
            lonper = argp
            omega = 0

    # Mean anomaly
    m = meanlon - fr * omega - argp
    m = np.mod(m, 2.0 * np.pi)

    # Solve for eccentric anomaly and true anomaly
    e0, nu = newtonm(ecc, m)

    if ecc < SMALL:
        # Circular orbits
        if is_equatorial(incl):
            # Circular equatorial
            truelon = nu
        else:
            # Circular inclined
            arglat = nu - fr * omega

    # Convert back to position and velocity vectors
    return coe2rv(p, ecc, incl, omega, nu, arglat, truelon, lonper)


def rv2eq(r, v):
    """Convert from position & velocity vectors to equinoctial elements.

    References:
        Vallado: 2013, p. 108
        Chobotov:       30

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
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper, _ = (
        rv2coe(r, v)
    )

    # Determine retrograde factor
    fr = -1 if abs(incl - np.pi) < SMALL else 1

    if ecc < SMALL:
        # Circular orbits
        if is_equatorial(incl):
            # Circular Equatorial
            argp, omega = 0.0, 0.0
            nu, m = truelon, truelon
        else:
            # Circular inclined
            argp = 0.0
            nu, m = arglat
    else:
        # Elliptical equatorial
        if is_equatorial(incl):
            argp = lonper
            omega = 0.0

    # Calculate mean motion
    # TODO: put in separate utility function
    n = np.sqrt(MU / (a**3))

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
    meanlon = np.mod(meanlon + TWOPI, TWOPI)

    # Determine true longitude
    truelon = fr * omega + argp + nu
    truelon = np.mod(truelon + TWOPI, TWOPI)

    return a, n, af, ag, chi, psi, meanlon, truelon, fr


###############################################################################
# Topocentric Elements
###############################################################################

def tradec2rv(rho, trtasc, tdecl, drho, tdrtasc, tddecl, rseci, vseci):
    """Converts topocentric coordinates (range, right ascension, declination,
    and their rates) into geocentric equatorial (ECI) position and velocity
    vectors.

    References:
        Vallado: 2022, p. 254, Eqs. 4-1 to 4-2

    Args:
        rho (float): Satellite range from site in km
        trtasc (float): Topocentric right ascension in radians
        tdecl (float): Topocentric declination in radians
        drho (float): Range rate in km/s
        tdrtasc (float): Topocentric right ascension rate in rad/s
        tddecl (float): Topocentric declination rate in rad/s
        rseci (array_like): ECI site position vector in km
        vseci (np.array): ECI site velocity vector in km/s

    Returns:
        tuple: (reci, veci)
            reci (np.array): ECI position vector in km
            veci (np.array): ECI velocity vector in km/s
    """

    # Calculate topocentric slant range vectors
    rhov = np.array([
        rho * np.cos(tdecl) * np.cos(trtasc),
        rho * np.cos(tdecl) * np.sin(trtasc),
        rho * np.sin(tdecl)
    ])

    # Slant range rate vectors
    drhov = np.array([
        drho * np.cos(tdecl) * np.cos(trtasc) - rho * np.sin(tdecl) * np.cos(
            trtasc) * tddecl - rho * np.cos(tdecl) * np.sin(trtasc) * tdrtasc,
        drho * np.cos(tdecl) * np.sin(trtasc) - rho * np.sin(tdecl) * np.sin(
            trtasc) * tddecl + rho * np.cos(tdecl) * np.cos(trtasc) * tdrtasc,
        drho * np.sin(tdecl) + rho * np.cos(tdecl) * tddecl
    ])

    # ECI position and velocity vectors
    reci = rhov + rseci
    veci = drhov + vseci

    return reci, veci


def rv2tradec(reci, veci, rseci, vseci):
    """Converts geocentric equatorial (ECI) position and velocity vectors into
    range, topocentric right ascension, declination, and rates.

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
            trtasc (float): Topocentric right ascension in radians
            tdecl (float): Topocentric declination in radians
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
    if temp < SMALL:
        trtasc = np.arctan2(drhoveci[1], drhoveci[0])
    else:
        trtasc = np.arctan2(rhoveci[1], rhoveci[0])

    # Directly over the North Pole
    if temp < SMALL:
        tdecl = np.sign(rhoveci[2]) * np.pi / 2  # +- 90 deg
    else:
        tdecl = np.arcsin(rhoveci[2] / rho)

    if trtasc < 0.0:
        trtasc += 2.0 * np.pi

    # Calculate topocentric right ascension and declination rates
    temp1 = -rhoveci[1] ** 2 - rhoveci[0] ** 2
    drho = np.dot(rhoveci, drhoveci) / rho
    if abs(temp1) > SMALL:
        dtrtasc = (drhoveci[0] * rhoveci[1] - drhoveci[1] * rhoveci[0]) / temp1
    else:
        dtrtasc = 0.0

    if abs(temp) > SMALL:
        dtdecl = (drhoveci[2] - drho * np.sin(tdecl)) / temp
    else:
        dtdecl = 0.0

    return rho, trtasc, tdecl, drho, dtrtasc, dtdecl


###############################################################################
# Flight Elements
###############################################################################

def flt2rv(rmag, vmag, latgc, lon, fpa, az, ttt, jdut1, lod, xp, yp, ddpsi,
           ddeps, eqeterms=True):
    """Converts flight elements into ECI position and velocity vectors.

    References:
        Vallado: 2013, XX
        Escobal:      397
        Chobotov:      67

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
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (reci, veci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
    """
    # Form position vector
    recef = np.array([
        rmag * np.cos(latgc) * np.cos(lon),
        rmag * np.cos(latgc) * np.sin(lon),
        rmag * np.sin(latgc)
    ])

    # Convert r to ECI
    vecef = np.zeros(3)  # this is a dummy for now
    aecef = np.zeros(3)
    reci, veci, _ = ecef2eci(
        recef, vecef, aecef, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, eqeterms
    )

    # Calculate right ascension and declination
    if np.sqrt(reci[0]**2 + reci[1]**2) < SMALL:
        rtasc = np.arctan2(veci[1], veci[0])
    else:
        rtasc = np.arctan2(reci[1], reci[0])
    decl = np.arcsin(reci[2] / rmag)

    # Form velocity vector
    fpav = np.pi * 0.5 - fpa
    veci = vmag * np.array([
        # First element
        (-np.cos(rtasc) * np.sin(decl)
         * (np.cos(az) * np.cos(fpav)
            - np.sin(rtasc) * np.sin(az) * np.cos(fpav))
         + np.cos(rtasc) * np.sin(decl) * np.sin(fpav)),

        # Second element
        (-np.sin(rtasc) * np.sin(decl)
         * (np.cos(az) * np.cos(fpav)
            + np.cos(rtasc) * np.sin(az) * np.cos(fpav))
         + np.sin(rtasc) * np.cos(decl) * np.sin(fpav)),

        # Third element
        (np.sin(decl) * np.sin(fpav)
         + np.cos(decl) * np.cos(az) * np.cos(fpav))
    ])

    return reci, veci


def rv2flt(reci, veci, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, eqeterms):
    """Transforms a position and velocity vector to flight elements.

    References:
        Vallado: 2001, XX

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
    rmag = np.linalg.norm(reci)
    vmag = np.linalg.norm(veci)

    # Convert r to ECEF for lat/lon calculations
    aecef = np.zeros(3)
    recef, vecef, aecef = eci2ecef(
        reci, veci, aecef, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, eqeterms
    )

    # Calculate longitude
    if np.sqrt(recef[0]**2 + recef[1]**2) < SMALL:
        lon = np.arctan2(vecef[1], vecef[0])
    else:
        lon = np.arctan2(recef[1], recef[0])

    latgc = np.arcsin(recef[2] / rmag)

    # Calculate right ascension and declination
    if np.sqrt(reci[0]**2 + reci[1]**2) < SMALL:
        rtasc = np.arctan2(veci[1], veci[0])
    else:
        rtasc = np.arctan2(reci[1], reci[0])

    decl = np.arcsin(reci[2] / rmag)

    # Calculate flight path angle
    h = np.cross(reci, veci)
    hmag = np.linalg.norm(h)
    rdotv = np.dot(reci, veci)
    fpav = np.arctan2(hmag, rdotv)
    fpa = np.pi / 2 - fpav

    # Calculte azimuth
    hcrossr = np.cross(h, reci)
    az = np.arctan2(
        reci[0] * hcrossr[1] - reci[1] * hcrossr[0], hcrossr[2] * rmag
    )

    return lon, latgc, rtasc, decl, fpa, az, rmag, vmag


###############################################################################
# Ecliptic Elements
###############################################################################

def ell2rv(rr, ecllon, ecllat, drr, decllon, decllat):
    """Transforms ecliptic latitude and longitude to position and velocity
    vectors.

    References:
        Vallado: 2004, XX

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
    r = np.array([
        rr * np.cos(ecllat) * np.cos(ecllon),
        rr * np.cos(ecllat) * np.sin(ecllon),
        rr * np.sin(ecllat)
    ])

    # Calculate velocity vector in ecliptic coordinates
    v = np.array([
        # X component
        drr * np.cos(ecllat) * np.cos(ecllon)
        - rr * np.sin(ecllat) * np.cos(ecllon) * decllat
        - rr * np.cos(ecllat) * np.sin(ecllon) * decllon,
        # Y component
        drr * np.cos(ecllat) * np.sin(ecllon)
        - rr * np.sin(ecllat) * np.sin(ecllon) * decllat
        + rr * np.cos(ecllat) * np.cos(ecllon) * decllon,
        # Z component
        drr * np.sin(ecllat) + rr * np.cos(ecllat) * decllat
    ])

    # Rotate position and velocity vectors to the ECI frame
    reci = rot1(r, -OBLIQUITYEARTH)
    veci = rot1(v, -OBLIQUITYEARTH)

    return reci, veci


def rv2ell(reci, veci):
    """Transforms position and velocity vectors to ecliptic latitude and
    longitude.

    References:
        Vallado: 2004, XX

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
    r = rot1(reci, OBLIQUITYEARTH)
    v = rot1(veci, OBLIQUITYEARTH)

    # Calculate magnitudes
    rr = np.linalg.norm(r)
    temp = np.sqrt(r[0]**2 + r[1]**2)

    # Calculate ecliptic longitude
    if temp < SMALL:
        temp1 = np.sqrt(v[0]**2 + v[1]**2)
        ecllon = np.arctan2(v[1], v[0]) if abs(temp1) > SMALL else 0
    else:
        ecllon = np.arctan2(r[1], r[0])

    # Calculate ecliptic latitude
    ecllat = np.arcsin(r[2] / rr)

    # Calculate rates
    temp1 = -r[1]**2 - r[0]**2  # Different now
    drr = np.dot(r, v) / rr
    decllon = (v[0] * r[1] - v[1] * r[0]) / temp1 if abs(temp1) > SMALL else 0
    decllat = (v[2] - drr * np.sin(ecllat)) / temp if abs(temp) > SMALL else 0

    return rr, ecllon, ecllat, drr, decllon, decllat


###############################################################################
# Celestial Elements
###############################################################################

def radec2rv(rr, rtasc, decl, drr, drtasc, ddecl):
    """Transforms celestial (right ascension and declination) elements to
    position and velocity vectors.

    References:
        Vallado: 2001, p. 246-248, Algorithm 25

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
    r = np.array([
        rr * np.cos(decl) * np.cos(rtasc),
        rr * np.cos(decl) * np.sin(rtasc),
        rr * np.sin(decl)
    ])

    # Velocity vector
    v = np.array([
        # X component
        drr * np.cos(decl) * np.cos(rtasc)
        - rr * np.sin(decl) * np.cos(rtasc) * ddecl
        - rr * np.cos(decl) * np.sin(rtasc) * drtasc,
        # Y component
        drr * np.cos(decl) * np.sin(rtasc)
        - rr * np.sin(decl) * np.sin(rtasc) * ddecl
        + rr * np.cos(decl) * np.cos(rtasc) * drtasc,
        # Z component
        drr * np.sin(decl) + rr * np.cos(decl) * ddecl
    ])

    return r, v


def rv2radec(r, v):
    """Transforms position and velocity vectors to celestial (right
    ascension and declination) elements.

    References:
        Vallado: 2001, p. 246-248, Algorithm 25

    Args:
        r (array_like): Position vector in km
        v (array_like): Velocity vector in km/s

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
    rtasc = np.arctan2(v[1], v[0]) if temp < SMALL else np.arctan2(r[1], r[0])
    rtasc += TWOPI if rtasc < 0.0 else rtasc

    # Calculate declination
    decl = np.arcsin(r[2] / rr)

    # Calculate radius rate
    drr = np.dot(r, v) / rr
    temp1 = -r[1] * r[1] - r[0] * r[0]

    # Calculate right ascension rate
    drtasc = (v[0] * r[1] - v[1] * r[0]) / temp1 if abs(temp1) > SMALL else 0

    # Calculate declination rate
    ddecl = (v[2] - drr * np.sin(decl)) / temp if abs(temp) > SMALL else 0

    return rr, rtasc, decl, drr, drtasc, ddecl


###############################################################################
# Azimuth-Elevation Elements
###############################################################################

def raz2rvs(rho, az, el, drho, daz, del_el):
    """Converts range, azimuth, and elevation values with slant range and
    velocity vectors  for a satellite from a radar site in the topocentric
    horizon (SEZ) system.

    References:
        Vallado: 2001, p. 250-251, Eqs. 4-4 and 4-5

    Args:
        rho (float): Satellite range from site in km
        az (float): Azimuth in radians (0 to 2pi)
        el (float): Elevation in radians (-pi/2 to pi/2)
        drho (float): Range rate in km/s
        daz (float): Azimuth rate in rad/s
        del_el (float): Elevation rate in rad/s

    Returns:
        tuple:
            rhosez (np.ndarray): SEZ range vector in km
            drhosez (np.ndarray): SEZ velocity vector in km/s
    """
    # Initialize values
    sinel, cosel = np.sin(el), np.cos(el)
    sinaz, cosaz = np.sin(az), np.cos(az)

    # Form SEZ range vector
    rhosez = np.array([
        -rho * cosel * cosaz,
        rho * cosel * sinaz,
        rho * sinel
    ])

    # Form SEZ velocity vector
    drhosez = np.array([
        -drho * cosel * cosaz + rhosez[2] * del_el * cosaz + rhosez[1] * daz,
        drho * cosel * sinaz - rhosez[2] * del_el * sinaz - rhosez[0] * daz,
        drho * sinel + rho * del_el * cosel
    ])

    return rhosez, drhosez


def rvs2raz(rhosez, drhosez):
    """Transforms SEZ range and velocity vectors to range, azimuth, and
    elevation values and their rates.

    References:
        Vallado: 2001, p. 250-251, Eqs. 4-4 and 4-5

    Args:
        rhosez (array_like): SEZ range vector in km
        drhosez (array_like): SEZ velocity vector in km/s

    Returns:
        tuple:
            rho (float): Satellite range from site in km
            az (float): Azimuth in radians (0 to 2pi)
            el (float): Elevation in radians (-pi/2 to pi/2)
            drho (float): Range rate in km/s
            daz (float): Azimuth rate in rad/s
            del_el (float): Elevation rate in rad/s
    """
    # Calculate azimuth
    temp = np.sqrt(rhosez[0] ** 2 + rhosez[1] ** 2)
    if abs(rhosez[1]) < SMALL:
        if temp < SMALL:
            az = np.arctan2(drhosez[1], -drhosez[0])
        else:
            az = np.pi if rhosez[0] > 0.0 else 0.0
    else:
        az = np.arctan2(rhosez[1], -rhosez[0])

    # Calculate elevation
    el = (
        np.sign(rhosez[2]) * SMALL if temp < SMALL
        else np.arcsin(rhosez[2] / np.linalg.norm(rhosez))
    )

    # Calculate range
    rho = np.linalg.norm(rhosez)

    # Range rate
    drho = np.dot(rhosez, drhosez) / rho

    # Azimuth rate
    daz = (
        (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) / (temp ** 2)
        if abs(temp ** 2) > SMALL else 0.0
    )

    # Elevation rate
    del_el = (
        (drhosez[2] - drho * np.sin(el)) / temp if abs(temp) > SMALL else 0.0
    )

    return rho, az, el, drho, daz, del_el


def razel2rv(rho, az, el, drho, daz, del_el, latgd, lon, alt, ttt, jdut1, lod,
             xp, yp, ddpsi, ddeps, eqeterms):
    """Transforms range, azimuth, elevation, and their rates to the geocentric
    equatorial (ECI) position and velocity vectors.

    References:
        Vallado: 2001, p. 250-255, Algorithm 27

    Args:
        rho (float): Satellite range from site [km]
        az (float): Azimuth [rad]
        el (float): Elevation [rad]
        drho (float): Range rate [km/s]
        daz (float): Azimuth rate [rad/s]
        del_el (float): Elevation rate [rad/s]
        latgd (float): Geodetic latitude of site in radians
        lon (float): Longitude of site in radians
        alt (float): Altitude of site in km
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (reci, veci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
    """
    # Find SEZ range and velocity vectors
    rhosez, drhosez = raz2rvs(rho, az, el, drho, daz, del_el)

    # Perform SEZ to ECEF transformation
    rhoecef = rot3(rot2(rhosez, latgd - HALFPI), -lon).T
    drhoecef = rot3(rot2(drhosez, latgd - HALFPI), -lon).T

    # Find ECEF range and velocity vectors
    rs, vs = site(latgd, lon, alt)
    recef = rhoecef + rs
    vecef = drhoecef

    # Convert ECEF to ECI
    a = np.array([0, 0, 0])
    reci, veci, aeci = ecef2eci(
        recef, vecef, a, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, eqeterms
    )

    return reci, veci


def rv2razel(reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp, ddpsi,
             ddeps, eqeterms=True):
    """Transforms ECI position and velocity vectors to range, azimuth,
    elevation, and their rates.

    The value of `SMALL` can affect the rate term calculations. the solution
    uses the velocity vector to find the singular cases. also, the elevation
    and azimuth rate terms are not observable unless the acceleration vector is
    available.

    References:
        Vallado: 2007, p. 268-269, Algorithm 27

    Args:
        reci (array): ECI position vector in km
        veci (array): ECI velocity vector in km/s
        latgd (float): Geodetic latitude of site in radians
        lon (float): Longitude of site in radians
        alt (float): Altitude of site in km
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        eqeterms (bool, optional): Add terms for ast calculation (default True)

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
    rsecef, vsecef = site(latgd, lon, alt)

    # Convert ECI to ECEF
    a = np.array([0, 0, 0])
    recef, vecef, aecef = eci2ecef(
        reci, veci, a, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, eqeterms
    )

    # Find ECEF range vector from site to satellite
    rhoecef = recef - rsecef
    drhoecef = vecef
    rho = np.linalg.norm(rhoecef)

    # Convert to SEZ for calculations
    tempvec = rot3(rhoecef, lon)
    rhosez = rot2(tempvec, HALFPI - latgd)

    tempvec = rot3(drhoecef, lon)
    drhosez = rot2(tempvec, HALFPI - latgd)

    # Calculate azimuth and elevation
    temp = np.sqrt(rhosez[0]**2 + rhosez[1]**2)
    if temp < SMALL:
        el = np.sign(rhosez[2]) * HALFPI
    else:
        magrhosez = np.linalg.norm(rhosez)
        el = np.arcsin(rhosez[2] / magrhosez)

    if temp < SMALL:
        az = np.arctan2(drhosez[1], -drhosez[0])
    else:
        az = np.arctan2(rhosez[1] / temp, -rhosez[0] / temp)

    # Calculate range, azimuth, and elevation rates
    drho = np.dot(rhosez, drhosez) / rho
    daz = (
        (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) / (temp * temp)
        if temp > SMALL else 0.0
    )
    del_el = (
        (drhosez[2] - drho * np.sin(el)) / temp if abs(temp) > SMALL else 0.0
    )

    return rho, az, el, drho, daz, del_el


###############################################################################
# Satellite Coordinate Systems
###############################################################################

def rv2rsw(reci, veci):
    """Transforms position and velocity vectors into radial, tangential
    (in-track), and normal (cross-track) coordinates, i.e. RSW frame.

    Note: There are numerous nomenclatures for these systems. This is the RSW
    system of Vallado. Rhe reverse values are found using the transmat
    transpose.

    References:
        Vallado: 2007, p. 172

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


def rv2ntw(reci, veci):
    """Transforms position and velocity vectors into normal (in-radial),
    tangential (velocity), and normal (cross-track) coordinates.

    Note that sometimes the first vector is called "along-radial". the
    tangential direction is always aligned with the velocity vector. this is
    the NTW system of Vallado.

    References:
        Vallado: 2007, p. 172

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
