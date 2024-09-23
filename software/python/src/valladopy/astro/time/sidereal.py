# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 7 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------


import numpy as np

from ... import constants as const


def gstime(jdut1):
    """Calculates the Greenwich Sidereal Time (IAU-82).

    References:
        Vallado: 2022, p. 189, Eq. 3-48

    Args:
        jdut1 (float): Julian date of UT1 (days from 4713 BC)

    Returns:
        float: Greenwich Sidereal Time (0 to 2pi radians)
    """
    # Julian centuries from the J2000.0 epoch
    tut1 = (jdut1 - const.J2000) / const.CENT2DAY

    # Calculate Greenwich Sidereal Time in seconds
    gst = (
        -6.2e-6 * tut1**3
        + 0.093104 * tut1**2
        + (876600.0 * 3600 + 8640184.812866) * tut1
        + 67310.54841
    )

    # Convert to radians
    return np.remainder(gst * const.EARTHROT_APPROX, const.TWOPI)


def sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms=True):
    """Calculates the transformation matrix that accounts for the effects of
    sidereal time.

    References:
        Vallado: 2013, p. 223-224

    Args:
        jdut1 (float): Julian date of UT1
        deltapsi (float): Nutation angle in radians
        meaneps (float): Mean obliquity of the ecliptic in radians
        omega (float): Longitude of ascending node of the moon in radians
        lod (float): Length of day in seconds
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple:
            st (np.ndarray): Transformation matrix for PEF to TOD
            stdot (np.ndarray): Transformation rate matrix
    """
    # Find GMST
    gmst = gstime(jdut1)

    # Find mean apparent sidereal time
    if jdut1 > 2450449.5 and eqeterms > 0:
        ast = (
            gmst
            + deltapsi * np.cos(meaneps)
            + 0.00264 * const.ARCSEC2RAD * np.sin(omega)
            + 0.000063 * const.ARCSEC2RAD * np.sin(2.0 * omega)
        )
    else:
        ast = gmst + deltapsi * np.cos(meaneps)

    ast = np.remainder(ast, const.TWOPI)
    omegaearth = const.EARTHROT * (1.0 - lod / const.DAY2SEC)

    # Transformation matrix for PEF to TOD
    st = np.array(
        [
            [np.cos(ast), -np.sin(ast), 0.0],
            [np.sin(ast), np.cos(ast), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    # Sidereal time rate matrix
    stdot = np.array(
        [
            [-omegaearth * np.sin(ast), -omegaearth * np.cos(ast), 0.0],
            [omegaearth * np.cos(ast), -omegaearth * np.sin(ast), 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    return st, stdot
