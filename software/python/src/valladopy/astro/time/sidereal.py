# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 7 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import numpy as np
from typing import Tuple

from ... import constants as const
from ...mathtime.julian_date import jday


def gstime(jdut1: float) -> float:
    """Calculates the Greenwich Sidereal Time (IAU-82).

    References:
        Vallado: 2022, p. 189, Eq. 3-48

    Args:
        jdut1 (float): Julian date of UT1 (days from 4713 BC)

    Returns:
        float: Greenwich Sidereal Time in radians (0 to 2pi)
    """
    # Julian centuries from the J2000 epoch
    tut1 = (jdut1 - const.J2000) / const.CENT2DAY

    # Calculate Greenwich Sidereal Time in seconds
    gst = (
        -6.2e-6 * tut1**3
        + 0.093104 * tut1**2
        + (876600 * const.HR2SEC + 8640184.812866) * tut1
        + 67310.54841
    )

    # Convert to radians
    return np.remainder(gst * const.EARTHROT_APPROX, const.TWOPI)


def gstime0(year: int) -> float:
    """Calculates the Greenwich Sidereal Time at the beginning (0 hr UT1 on January 1)
    of the given year.

    References:
        Vallado: 2022, p. 189-190, Eq. 3-48

    Args:
        year (int): Year (e.g., 1998, 1999, etc.)

    Returns:
        float: Greenwich Sidereal Time in radians (0 to 2pi)
    """
    # Calculate Julian Date at 0 hr UT1 on January 1
    jd, _ = jday(year, month=1, day=1)

    return gstime(jd)


def lstime(lon: float, jdut1: float) -> Tuple[float, float]:
    """Calculates the local sidereal time (LST) and Greenwich sidereal time (GST) at a
    given location (GST from IAU-82).

    References:
        Vallado: 2022, p. 190, Algorithm 15

    Args:
        lon (float): Longitude of the site in radians (-2pi to 2pi) (West is negative)
        jdut1 (float): Julian date of UT1 (days from 4713 BC)

    Returns:
        tuple: (lst, gst)
            lst (float): Local sidereal time (LST) in radians (0 to 2pi)
            gst (float): Greenwich sidereal time (GST) in radians (0 to 2pi)
    """
    # Calculate GST
    gst = gstime(jdut1)

    # Calculate LST
    lst = np.mod(lon + gst, const.TWOPI)

    return lst, gst


def sidereal(
    jdut1: float,
    deltapsi: float,
    meaneps: float,
    omega: float,
    lod: float,
    use_iau80: bool = True,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculates the transformation matrix that accounts for the effects of
    sidereal time.

    References:
        Vallado: 2022, p. 224-225

    Args:
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        deltapsi (float): Nutation angle in radians
        meaneps (float): Mean obliquity of the ecliptic in radians
        omega (float): Longitude of ascending node of the moon in radians
        lod (float): Length of day in seconds
        use_iau80 (bool, optional): Use IAU-80 theory (default True)
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (st, stdot)
            st (np.ndarray): Transformation matrix for PEF to TOD
            stdot (np.ndarray): Transformation rate matrix
    """
    # Calculate apparent GMST
    if use_iau80:
        # Find GMST
        gmst = gstime(jdut1)

        # Find mean apparent sidereal time
        if jdut1 > 2450449.5 and eqeterms:
            ast = (
                gmst
                + deltapsi * np.cos(meaneps)
                + 0.00264 * const.ARCSEC2RAD * np.sin(omega)
                + 0.000063 * const.ARCSEC2RAD * np.sin(2 * omega)
            )
        else:
            ast = gmst + deltapsi * np.cos(meaneps)

        ast = np.remainder(ast, const.TWOPI)
    else:
        tut1d = jdut1 - const.J2000
        era = const.TWOPI * (0.7790572732640 + 1.00273781191135448 * tut1d)
        ast = np.remainder(era, const.TWOPI)

    # Transformation matrix for PEF to TOD
    st = np.array(
        [[np.cos(ast), -np.sin(ast), 0], [np.sin(ast), np.cos(ast), 0], [0, 0, 1]]
    )

    # Sidereal time rate matrix
    omegaearth = const.EARTHROT * (1 - lod / const.DAY2SEC)
    stdot = np.array(
        [
            [-omegaearth * np.sin(ast), -omegaearth * np.cos(ast), 0],
            [omegaearth * np.cos(ast), -omegaearth * np.sin(ast), 0],
            [0, 0, 0],
        ]
    )

    return st, stdot
