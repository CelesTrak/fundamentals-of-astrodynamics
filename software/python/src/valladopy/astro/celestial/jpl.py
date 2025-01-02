# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 22 January 2018
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

from enum import Enum
from typing import Any, Tuple, Dict

import numpy as np
from scipy.interpolate import CubicSpline

from ... import constants as const
from ...mathtime.julian_date import jday


class JPLInterp(Enum):
    NONE = "n"
    LINEAR = "l"
    SPLINE = "s"


def init_jplde(filepath: str) -> Tuple[Dict[str, np.ndarray], float, float]:
    """Initializes the JPL planetary ephemeris data by loading the sun and moon
    positions.

    Args:
        filepath (str): Path to the input text file containing ephemeris data

    Returns:
        tuple: (jpldearr, jdjpldestart, jdjpldestart_frac)
            jpldearr (dict[str, np.ndarray]): Dictionary of JPL DE data records
            jdjpldestart (float): Julian date of the start of the JPL DE data
            jdjpldestart_frac (float): Fractional part of the Julian date at the start
    """
    # Load the input file data
    file_data = np.loadtxt(filepath)

    # Initialize the JPL DE data dictionary
    jpldearr = {
        "year": file_data[:, 0].astype(int),
        "month": file_data[:, 1].astype(int),
        "day": file_data[:, 2].astype(int),
        "rsun1": file_data[:, 3],
        "rsun2": file_data[:, 4],
        "rsun3": file_data[:, 5],
        "rsmag": file_data[:, 6],
        "rmoon1": file_data[:, 8],
        "rmoon2": file_data[:, 9],
        "rmoon3": file_data[:, 10],
    }

    # Calculate Modified Julian Date (MJD)
    jd, jd_frac = jday(jpldearr["year"], jpldearr["month"], jpldearr["day"])
    jpldearr["mjd"] = jd + jd_frac - const.JD_TO_MJD_OFFSET

    # Find the start epoch date
    jdjpldestart, jdjpldestart_frac = jday(
        jpldearr["year"][0], jpldearr["month"][0], jpldearr["day"][0]
    )

    return jpldearr, jdjpldestart, jdjpldestart_frac


def find_jplde_param(
    jdtdb: float,
    jdtdb_f: float,
    jpldearr: dict[str, Any],
    jdjpldestart: float,
    interp: JPLInterp = JPLInterp.NONE,
) -> Tuple[np.ndarray, np.ndarray]:
    """Finds the JPL DE parameters for a given time using interpolation.

    References:
        Vallado, 2022

    Args:
        jdtdb (float): Epoch Julian date (days from 4713 BC)
        jdtdb_f (float): Fractional part of the epoch Julian date
        jpldearr (dict[str, Any]): Dictionary of JPL DE data records
        jdjpldestart (float): Julian date of the start of the JPL DE data
        interp (JPLInterp, optional): Interpolation method to use
                                      (default is JPLInterp.NONE)

    Returns:
        tuple: (rsun, rmoon)
            rsun (np.ndarray): ECI sun position vector in km
            rmoon (np.ndarray): ECI moon position vector in km
    """
    # Compute whole-day Julian date and minutes from midnight
    jdb = np.floor(jdtdb + jdtdb_f) + 0.5
    mfme = (jdtdb + jdtdb_f - jdb) * const.DAY2MIN
    if mfme < 0.0:
        mfme += const.DAY2MIN

    # Determine record index
    jdjpldestarto = np.floor(jdtdb + jdtdb_f - jdjpldestart)
    recnum = int(jdjpldestarto) - 1

    # Default values if out of bounds
    if not (0 <= recnum <= len(jpldearr["rsun1"]) - 1):
        return np.zeros(3), np.zeros(3)

    # Non-interpolated values
    rsun = np.array(
        [
            jpldearr["rsun1"][recnum],
            jpldearr["rsun2"][recnum],
            jpldearr["rsun3"][recnum],
        ]
    )
    rmoon = np.array(
        [
            jpldearr["rmoon1"][recnum],
            jpldearr["rmoon2"][recnum],
            jpldearr["rmoon3"][recnum],
        ]
    )

    # Fix for interpolation (get back to minutes)
    fixf = mfme / const.DAY2MIN

    # Linear interpolation
    if interp == JPLInterp.LINEAR:
        rsun += fixf * (
            np.array(
                [
                    jpldearr["rsun1"][recnum + 1],
                    jpldearr["rsun2"][recnum + 1],
                    jpldearr["rsun3"][recnum + 1],
                ]
            )
            - rsun
        )
        rmoon += fixf * (
            np.array(
                [
                    jpldearr["rmoon1"][recnum + 1],
                    jpldearr["rmoon2"][recnum + 1],
                    jpldearr["rmoon3"][recnum + 1],
                ]
            )
            - rmoon
        )

    # Cubic spline interpolation
    elif interp == JPLInterp.SPLINE:
        idx1, idx2 = recnum - 1, recnum + 3
        mjds = jpldearr["mjd"][idx1:idx2]

        # Interpolate each component of rsun and rmoon separately
        rsun[0] = CubicSpline(mjds, jpldearr["rsun1"][idx1:idx2])(mjds[1] + fixf)
        rsun[1] = CubicSpline(mjds, jpldearr["rsun2"][idx1:idx2])(mjds[1] + fixf)
        rsun[2] = CubicSpline(mjds, jpldearr["rsun3"][idx1:idx2])(mjds[1] + fixf)
        rmoon[0] = CubicSpline(mjds, jpldearr["rmoon1"][idx1:idx2])(mjds[1] + fixf)
        rmoon[1] = CubicSpline(mjds, jpldearr["rmoon2"][idx1:idx2])(mjds[1] + fixf)
        rmoon[2] = CubicSpline(mjds, jpldearr["rmoon3"][idx1:idx2])(mjds[1] + fixf)

    return rsun, rmoon
