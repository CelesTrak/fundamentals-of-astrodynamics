# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import logging
from typing import Tuple

import numpy as np

from ... import constants as const


# Set up logging
logger = logging.getLogger(__name__)


def pathm(llat: float, llon: float, range_: float, az: float) -> Tuple[float, float]:
    """Determines the end position (latitude and longitude) for a given range and
    azimuth from a given starting point.

    References:
        Vallado: 2022, p. 872, Eq. 11-6 and 11-7

    Args:
        llat (float): Start geocentric latitude in radians (-pi/2 to pi/2)
        llon (float): Start longitude in radians (0 to 2pi).
        range_(float): Range between points in Earth radii
        az (float): Azimuth in radians (0 to 2pi)

    Returns:
        tuple: (tlat, tlon)
            tlat (float): End geocentric latitude in radians (-pi/2 to pi/2)
            tlon (float): End longitude in radians (0 to 2pi)
    """
    # Normalize inputs
    az = az % const.TWOPI
    llon = (llon + const.TWOPI) % const.TWOPI
    range_ = range_ % const.TWOPI

    # Find geocentric latitude
    tlat = np.arcsin(
        np.sin(llat) * np.cos(range_) + np.cos(llat) * np.sin(range_) * np.cos(az)
    )

    # Find delta n, the angle between the points
    deltan = 0
    if abs(np.cos(tlat)) > const.SMALL and abs(np.cos(llat)) > const.SMALL:
        sindn = np.sin(az) * np.sin(range_) / np.cos(tlat)
        cosdn = (np.cos(range_) - np.sin(tlat) * np.sin(llat)) / (
            np.cos(tlat) * np.cos(llat)
        )
        deltan = np.arctan2(sindn, cosdn)
    else:
        # Case where launch is within a small distance of a pole
        if abs(np.cos(llat)) <= const.SMALL:
            deltan = az + np.pi if np.pi < range_ < const.TWOPI else az

        # Case where end point is within a small distance of a pole
        elif abs(np.cos(tlat)) <= const.SMALL:
            deltan = 0

    # Compute end longitude
    tlon = (llon + deltan) % const.TWOPI

    return tlat, tlon


def rngaz(
    llat: float, llon: float, tlat: float, tlon: float, tof: float = 0.0
) -> Tuple[float, float]:
    """Calculates the range and azimuth between two specified ground points
    on a spherical Earth.

    References:
        Vallado: 2022, p. 872, Eq. 11-3 to 11-5

    Args:
        llat (float): Start geocentric latitude in radians (-pi/2 to pi/2)
        llon (float): Start longitude in radians (0 to 2pi)
        tlat (float): End geocentric latitude in radians (-pi/2 to pi/2)
        tlon (float): End longitude in radians (0 to 2pi)
        tof (float): Time of flight if applicable, in seconds (default is 0)

    Returns:
        tuple: (range_, az)
            range_ (float): Range between points in km
            az (float): Azimuth in radians (0 to 2pi)
    """
    # Calculate the spherical range
    range_ = np.arccos(
        np.sin(llat) * np.sin(tlat)
        + np.cos(llat) * np.cos(tlat) * np.cos(tlon - llon + const.EARTHROT * tof)
    )

    # Check for special cases where range is 0 or half the Earth
    if abs(np.sin(range_) * np.cos(llat)) < const.SMALL:
        az = np.pi if abs(range_ - np.pi) < const.SMALL else 0
    else:
        az = np.arccos(
            (np.sin(tlat) - np.cos(range_) * np.sin(llat))
            / (np.sin(range_) * np.cos(llat))
        )

    # Adjust azimuth if it is greater than pi
    if np.sin(tlon - llon + const.EARTHROT * tof) < 0:
        az = const.TWOPI - az

    return range_ * const.RE, az


def satfov(
    az: float,
    slatgd: float,
    slon: float,
    salt: float,
    tfov: float,
    etactr: float,
    tol_fov: float = 1e-5,
) -> Tuple[float, float, float, float]:
    """Finds parameters relating to a satellite's field of view (FOV).

    References:
        Vallado: 2022, p. 874-876, Eq. 11-8 to 11-13

    Args:
        az (float): Azimuth in radians (0.0 to 2pi)
        slatgd (float): Geodetic latitude of the satellite in radians (-pi/2 to pi/2)
        slon (float): Longitude of the satellite in radians
        salt (float): Altitude of the satellite in km
        tfov (float): Total field of view in radians
        etactr (float): Center where the sensor is pointed in radians
        tol_fov (float): Tolerance for the center of the FOV (default is 1e-5)

    Returns:
        tuple: (rhomin, rhomax, lat, lon)
            rhomin (float): Minimum slant range in km
            rhomax (float): Maximum slant range in km
            lat (float): Latitude of the center of the FOV in radians
            lon (float): Longitude of the center of the FOV in radians
    """
    # Satellite parameters and limiting cases
    r = const.RE + salt

    # Maximum field of view
    fovmin = etactr + tfov * 0.5
    gamma = np.pi - np.arcsin(r * np.sin(fovmin) / const.RE)  # use larger angle
    rhomax = const.RE * np.cos(gamma) + r * np.cos(fovmin)

    # Slant range
    lambda_ = np.arcsin(rhomax * np.sin(fovmin) / const.RE)
    rhomin = lambda_ * const.RE

    # Location of center of FOV
    if abs(etactr) > tol_fov:
        lat, lon = pathm(slatgd, slon, lambda_, az)
    else:
        lat, lon = slatgd, slon

    return rhomin, rhomax, lat, lon
