# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 22 January 2018
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from enum import Enum
from typing import Any, Tuple, Dict

import numpy as np
from scipy.interpolate import CubicSpline

from ... import constants as const
from ...mathtime.julian_date import jday


class JPLInterp(Enum):
    LINEAR = "l"
    SPLINE = "s"


def read_jplde(
    filepath: str, include_hr: bool = True
) -> Tuple[Dict[str, np.ndarray], float, float]:
    """Initializes the JPL planetary ephemeris data by loading the sun and moon
    positions.

    Args:
        filepath (str): Path to the input text file containing ephemeris data
        include_hr (bool, optional): Set to True if the file includes the hour column
                                     (defaults to True)

    Returns:
        tuple: (jpldearr, jdjpldestart, jdjpldestart_frac)
            jpldearr (dict[str, np.ndarray]): Dictionary of JPL DE data records
            jdjpldestart (float): Julian date of the start of the JPL DE data
            jdjpldestart_frac (float): Fractional part of the Julian date at the start
    """
    # Load the input file data
    file_data = np.loadtxt(filepath)

    # Adjust indices based on the presence of the hour column
    offset = 1 if include_hr else 0

    # Extract common fields
    year = file_data[:, 0].astype(int)
    month = file_data[:, 1].astype(int)
    day = file_data[:, 2].astype(int)
    hr = file_data[:, 3] if include_hr else np.zeros_like(year)

    # Construct the JPL DE data dictionary
    jpldearr = {
        "year": year,
        "month": month,
        "day": day,
        "hour": hr,
        "rsun1": file_data[:, 3 + offset],
        "rsun2": file_data[:, 4 + offset],
        "rsun3": file_data[:, 5 + offset],
        "rsmag": file_data[:, 6 + offset],
        "rmoon1": file_data[:, 8 + offset],
        "rmoon2": file_data[:, 9 + offset],
        "rmoon3": file_data[:, 10 + offset],
        "mjd": np.zeros(len(year)),
    }

    # Calculate Modified Julian Date (MJD)
    for i in range(len(year)):
        jd, jd_frac = jday(
            jpldearr["year"][i],
            jpldearr["month"][i],
            jpldearr["day"][i],
            jpldearr["hour"][i],
        )
        jpldearr["mjd"][i] = jd + jd_frac - const.JD_TO_MJD_OFFSET

    # Find the start epoch date
    jdjpldestart, jdjpldestart_frac = jday(
        jpldearr["year"][0],
        jpldearr["month"][0],
        jpldearr["day"][0],
        jpldearr["hour"][0],
    )

    return jpldearr, jdjpldestart, jdjpldestart_frac


def find_jplde_param(
    jdtdb: float,
    jdtdb_f: float,
    jpldearr: dict[str, Any],
    interp: JPLInterp | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Finds the JPL DE parameters for a given time using interpolation.

    References:
        Vallado, 2022

    Args:
        jdtdb (float): Epoch Julian date (days from 4713 BC)
        jdtdb_f (float): Fractional part of the epoch Julian date
        jpldearr (dict[str, Any]): Dictionary of JPL DE data records
        interp (JPLInterp, optional): Interpolation method to use (default is None)

    Returns:
        tuple: (rsun, rmoon)
            rsun (np.ndarray): ECI sun position vector in km
            rmoon (np.ndarray): ECI moon position vector in km

    TODO:
        - Fix for the 1/2 day offset case
    """
    # Compute whole-day Julian date and minutes from midnight
    jdb = np.floor(jdtdb + jdtdb_f) + 0.5
    mfme = (jdtdb + jdtdb_f - jdb) * const.DAY2MIN
    if mfme < 0:
        mfme += const.DAY2MIN

    # Determine record index
    jdjpldestarto = np.floor(
        jdtdb + jdtdb_f - jpldearr["mjd"][0] - const.JD_TO_MJD_OFFSET
    )
    if np.any(jpldearr["hour"]):
        # TODO: this works better when the date is closer to the hour, but not so great
        #       when it is close to the 1/2 day; Matlab and C# use a +1 additional
        #       offset (i.e. +2 instead of +1) but this doesn't seem quite right either
        recnum = int(jdjpldestarto) * 2 + 1  # 12-hr data
        if mfme > 720:
            mfme -= 720
    else:
        recnum = int(jdjpldestarto)  # 1-day data

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


def sunmoon(
    jdtdb: float,
    jdtdb_f: float,
    jpldearr: dict[str, Any],
    interp: JPLInterp | None = None,
) -> Tuple[np.ndarray, float, float, np.ndarray, float, float]:
    """Calculates the geocentric equatorial position vectors of the Sun and Moon.

    Args:
        jdtdb (float): Epoch Julian date (days from 4713 BC)
        jdtdb_f (float): Fractional part of the epoch Julian date
        jpldearr (dict[str, Any]): Dictionary of JPL DE data records
        interp (JPLInterp, optional): Interpolation method to use (default is None)

    Returns:
        tuple: (rsun, rtascs, decls, rmoon, rtascm, declm)
            rsun (np.ndarray): ECI sun position vector in km
            rtascs (float): Sun right ascension in radians
            decls (float): Sun declination in radians
            rmoon (np.ndarray): ECI moon position vector in km
            rtascm (float): Moon right ascension in radians
            declm (float): Moon declination in radians
    """
    # Get Sun and Moon position vectors and magnitudes
    rsun, rmoon = find_jplde_param(jdtdb, jdtdb_f, jpldearr, interp)

    # Sun right ascension
    temp_sun = np.hypot(rsun[0], rsun[1])
    rtascs = 0 if temp_sun < const.SMALL else np.arctan2(rsun[1], rsun[0])

    # Sun declination
    decls = np.arcsin(rsun[2] / np.linalg.norm(rsun))

    # Moon right ascension
    temp_moon = np.hypot(rmoon[0], rmoon[1])
    rtascm = 0 if temp_moon < const.SMALL else np.arctan2(rmoon[1], rmoon[0])

    # Moon declination
    declm = np.arcsin(rmoon[2] / np.linalg.norm(rmoon))

    return rsun, rtascs, decls, rmoon, rtascm, declm
