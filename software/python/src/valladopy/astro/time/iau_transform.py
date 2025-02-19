# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 16 July 2004
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from enum import Enum
from typing import Tuple

import numpy as np
from scipy.interpolate import CubicSpline

from .data import IAU06pnOldArray, IAU06Array, IAU06xysArray, EOPArray
from .utils import FundArgs, fundarg, precess
from ... import constants as const
from ...mathtime.julian_date import jday
from ...mathtime.vector import rot1mat, rot2mat, rot3mat


class InterpolationMode(Enum):
    LINEAR = "linear"
    SPLINE = "spline"


def iau06era(jdut1: float) -> np.ndarray:
    """Calculates the transformation matrix that accounts for sidereal time via the
    Earth Rotation Angle (ERA).

    References:
        Vallado, 2022, p. 214

    Args:
        jdut1 (float): Julian date of UT1 (days)

    Returns:
        np.ndarray: 3x3 transformation matrix for PEF to IRE
    """
    # Julian centuries of UT1 (in days from J2000 epoch)
    tut1d = jdut1 - const.J2000

    # Earth rotation angle (ERA) in radians
    era = const.TWOPI * (0.779057273264 + 1.00273781191135448 * tut1d)
    era = np.mod(era, const.TWOPI)

    # Transformation matrix from PEF to IRE
    return np.array(
        [[np.cos(era), -np.sin(era), 0], [np.sin(era), np.cos(era), 0], [0, 0, 1]]
    )


def iau06gst(
    jdut1: float, ttt: float, deltapsi: float, fundargs: FundArgs, iau06arr: IAU06Array
) -> Tuple[float, np.ndarray]:
    """Calculates the IAU 2006 Greenwich Sidereal Time (GST) and transformation matrix.

    References:
        Vallado, 2022, p. 217

    Args:
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        ttt (float): Julian centuries of TT
        deltapsi (float): Change in longitude in radians
        fundargs (FundArgs): Delaunay and planetary arguments
        iau06arr (IAU06Array): IAU 2006 data

    Returns:
        tuple: (gst, st)
            gst (float): Greenwich Sidereal Time in radians (0 to 2pi)
            st (np.ndarray): 3x3 transformation matrix
    """
    # Mean obliquity of the ecliptic
    epsa = (
        84381.406
        - 46.836769 * ttt
        - 0.0001831 * ttt**2
        + 0.0020034 * ttt**3
        - 0.000000576 * ttt**4
        - 0.0000000434 * ttt**5
    )  # arcseconds
    epsa = np.mod(np.radians(epsa / const.DEG2ARCSEC), const.TWOPI)

    # Evaluate the EE complementary terms
    gstsum0, gstsum1 = 0, 0
    n_elem = len(iau06arr.agsti) - 1
    for i in range(n_elem):
        tempval = (
            iau06arr.agsti[i, 0] * fundargs.l
            + iau06arr.agsti[i, 1] * fundargs.l1
            + iau06arr.agsti[i, 2] * fundargs.f
            + iau06arr.agsti[i, 3] * fundargs.d
            + iau06arr.agsti[i, 4] * fundargs.omega
            + iau06arr.agsti[i, 5] * fundargs.lonmer
            + iau06arr.agsti[i, 6] * fundargs.lonven
            + iau06arr.agsti[i, 7] * fundargs.lonear
            + iau06arr.agsti[i, 8] * fundargs.lonmar
            + iau06arr.agsti[i, 9] * fundargs.lonjup
            + iau06arr.agsti[i, 10] * fundargs.lonsat
            + iau06arr.agsti[i, 11] * fundargs.lonurn
            + iau06arr.agsti[i, 12] * fundargs.lonnep
            + iau06arr.agsti[i, 13] * fundargs.precrate
        )
        gstsum0 += iau06arr.agst[i, 0] * np.sin(tempval) + iau06arr.agst[i, 1] * np.cos(
            tempval
        )

    tempval = (
        iau06arr.agsti[n_elem, 0] * fundargs.l
        + iau06arr.agsti[n_elem, 1] * fundargs.l1
        + iau06arr.agsti[n_elem, 2] * fundargs.f
        + iau06arr.agsti[n_elem, 3] * fundargs.d
        + iau06arr.agsti[n_elem, 4] * fundargs.omega
        + iau06arr.agsti[n_elem, 5] * fundargs.lonmer
        + iau06arr.agsti[n_elem, 6] * fundargs.lonven
        + iau06arr.agsti[n_elem, 7] * fundargs.lonear
        + iau06arr.agsti[n_elem, 8] * fundargs.lonmar
        + iau06arr.agsti[n_elem, 9] * fundargs.lonjup
        + iau06arr.agsti[n_elem, 10] * fundargs.lonsat
        + iau06arr.agsti[n_elem, 11] * fundargs.lonurn
        + iau06arr.agsti[n_elem, 12] * fundargs.lonnep
        + iau06arr.agsti[n_elem, 13] * fundargs.precrate
    )
    gstsum1 += iau06arr.agst[n_elem, 0] * ttt * np.sin(tempval) + iau06arr.agst[
        n_elem, 1
    ] * ttt * np.cos(tempval)
    eect2000 = gstsum0 + gstsum1 * ttt

    # Equation of the equinoxes
    ee2000 = deltapsi * np.cos(epsa) + eect2000

    # Earth rotation angle (ERA)
    tut1d = jdut1 - const.J2000  # days from the Jan 1, 2000 12h epoch (UT1)
    era = const.TWOPI * (0.779057273264 + 1.00273781191135448 * tut1d)
    era = np.mod(era, const.TWOPI)

    # Greenwich Mean Sidereal Time (GMST), IAU 2000
    gmst2000 = era + (
        (
            0.014506
            + 4612.156534 * ttt
            + 1.3915817 * ttt**2
            - 0.00000044 * ttt**3
            + 0.000029956 * ttt**4
            + 0.0000000368 * ttt**5
        )
        * const.ARCSEC2RAD
    )

    # Greenwich Sidereal Time (GST)
    gst = gmst2000 + ee2000

    # Transformation matrix
    st = np.array(
        [[np.cos(gst), -np.sin(gst), 0], [np.sin(gst), np.cos(gst), 0], [0, 0, 1]]
    )

    return gst, st


########################################################################################
# IAU 2006 Precession-Nutation Theories (IAU2006/2000A and IAU2006/2000B)
########################################################################################


def _build_transformation_matrices(
    ttt: float, deltaeps: float, deltapsi: float, use_extended_prec: bool
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constructs nutation, precession, and combined precession-nutation matrices.

    Args:
        ttt (float): Julian centuries of TT
        deltaeps (float): Nutation in obliquity in radians
        deltapsi (float): Nutation in longitude in radians
        use_extended_prec (bool): Whether to include extended precession terms

    Returns:
        tuple:
            nut (np.ndarray): Nutation matrix (mean to true transformation)
            prec (np.ndarray): Precession matrix (J2000 to date transformation)
            pnb (np.ndarray): Combined precession-nutation matrix (ICRS to GCRF)
    """
    # Get precession angles
    _, psia, wa, ea, xa = precess(ttt, opt="06")

    # Obliquity of the ecliptic
    oblo = 84381.406 * const.ARCSEC2RAD

    # Nutation matrix
    a1 = rot1mat(ea + deltaeps)
    a2 = rot3mat(deltapsi)
    a3 = rot1mat(-ea)
    nut = a3 @ a2 @ a1

    # Precession matrix
    a4 = rot3mat(-xa)
    a5 = rot1mat(wa)
    a6 = rot3mat(psia)
    a7 = rot1mat(-oblo)

    # ICRS to J2000
    a8 = rot1mat(-0.0068192 * const.ARCSEC2RAD)
    a9 = rot2mat(0.041775 * np.sin(oblo) * const.ARCSEC2RAD)
    a10 = rot3mat(0.0146 * const.ARCSEC2RAD)

    # Precession and combined matrices
    if use_extended_prec:
        prec = a10 @ a9 @ a8 @ a7 @ a6 @ a5 @ a4
        pnb = prec @ nut
    else:
        prec = a7 @ a6 @ a5 @ a4
        pnb = a10 @ a9 @ a8 @ prec @ nut

    return nut, prec, pnb


def iau06pna(
    ttt: float, iau06arr: IAU06pnOldArray
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, FundArgs]:
    """Calculates the transformation matrix that accounts for the effects of
    precession-nutation using the IAU2006 precession theory and the IAU2000A nutation
    model.

    References:
        Vallado, 2022, p. 214-216

    Args:
        ttt (float): Julian centuries of TT
        iau06arr (IAU06pnOldArray): IAU 2006 data (old nutation coefficients)

    Returns:
        tuple: (deltapsi, pnb, prec, nut, fundargs)
            deltapsi (float): Change in longitude in radians
            pnb (np.ndarray): Combined precession-nutation matrix
            prec (np.ndarray): Precession transformation matrix (MOD to J2000)
            nut (np.ndarray): Nutation transformation matrix (IRE to GCRF)
            fundargs (FundArgs): Delaunay and planetary arguments
    """
    # Obtain data for calculations from the IAU 2006 nutation theory
    fundargs = fundarg(ttt, opt="06")

    # Compute luni-solar nutation
    pnsum, ensum = 0, 0
    for i in range(len(iau06arr.apni) - 1, -1, -1):
        tempval = (
            iau06arr.apni[i, 0] * fundargs.l
            + iau06arr.apni[i, 1] * fundargs.l1
            + iau06arr.apni[i, 2] * fundargs.f
            + iau06arr.apni[i, 3] * fundargs.d
            + iau06arr.apni[i, 4] * fundargs.omega
        )
        tempval = np.mod(tempval, const.TWOPI)
        pnsum += (iau06arr.apn[i, 0] + iau06arr.apn[i, 1] * ttt) * np.sin(
            tempval
        ) + iau06arr.apn[i, 4] * np.cos(tempval)
        ensum += (iau06arr.apn[i, 2] + iau06arr.apn[i, 3] * ttt) * np.cos(
            tempval
        ) + iau06arr.apn[i, 6] * np.sin(tempval)

    # Compute planetary nutation
    pplnsum, eplnsum = 0, 0
    for i in range(len(iau06arr.appli)):
        tempval = (
            iau06arr.appli[i, 0] * fundargs.l
            + iau06arr.appli[i, 1] * fundargs.l1
            + iau06arr.appli[i, 2] * fundargs.f
            + iau06arr.appli[i, 3] * fundargs.d
            + iau06arr.appli[i, 4] * fundargs.omega
            + iau06arr.appli[i, 5] * fundargs.lonmer
            + iau06arr.appli[i, 6] * fundargs.lonven
            + iau06arr.appli[i, 7] * fundargs.lonear
            + iau06arr.appli[i, 8] * fundargs.lonmar
            + iau06arr.appli[i, 9] * fundargs.lonjup
            + iau06arr.appli[i, 10] * fundargs.lonsat
            + iau06arr.appli[i, 11] * fundargs.lonurn
            + iau06arr.appli[i, 12] * fundargs.lonnep
            + iau06arr.appli[i, 13] * fundargs.precrate
        )
        pplnsum += iau06arr.appl[i, 0] * np.sin(tempval) + iau06arr.appl[i, 1] * np.cos(
            tempval
        )
        eplnsum += iau06arr.appl[i, 2] * np.sin(tempval) + iau06arr.appl[i, 3] * np.cos(
            tempval
        )

    # Combine nutation components
    deltapsi = pnsum + pplnsum
    deltaeps = ensum + eplnsum

    # Apply IAU 2006 corrections
    j2d = -2.7774e-6 * ttt * const.ARCSEC2RAD
    deltapsi += deltapsi * (0.4697e-6 + j2d)
    deltaeps += deltaeps * j2d

    # Build transformation matrices
    nut, prec, pnb = _build_transformation_matrices(ttt, deltaeps, deltapsi, False)

    return deltapsi, pnb, prec, nut, fundargs


def iau06pnb(
    ttt: float, iau06arr: IAU06pnOldArray
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, FundArgs]:
    """Calculates the transformation matrix that accounts for the effects of
    precession-nutation using the IAU2006 precession theory and a simplified nutation
    model based on IAU2000B.

    References:
        Vallado, 2022, p. 214-216

    Args:
        ttt (float): Julian centuries of TT
        iau06arr (IAU06pnOldArray): IAU 2006 data (old nutation coefficients)

    Returns:
        tuple:
            deltapsi (float): Change in longitude in radians
            pnb (np.ndarray): Combined precession-nutation matrix
            prec (np.ndarray): Precession transformation matrix (MOD to J2000)
            nut (np.ndarray): Nutation transformation matrix (IRE to GCRF)
            fundargs (FundArgs): Delaunay and planetary arguments
    """
    # Definitions
    iau2000b_terms = 77

    # Obtain data for calculations from the 2000b theory
    fundargs = fundarg(ttt, opt="02")

    # Compute luni-solar nutation
    pnsum, ensum = 0, 0
    for i in range(iau2000b_terms - 1, -1, -1):
        tempval = (
            iau06arr.apni[i, 0] * fundargs.l
            + iau06arr.apni[i, 1] * fundargs.l1
            + iau06arr.apni[i, 2] * fundargs.f
            + iau06arr.apni[i, 3] * fundargs.d
            + iau06arr.apni[i, 4] * fundargs.omega
        )
        pnsum += (iau06arr.apn[i, 0] + iau06arr.apn[i, 1] * ttt) * np.sin(tempval) + (
            iau06arr.apn[i, 4] + iau06arr.apn[i, 5] * ttt
        ) * np.cos(tempval)
        ensum += (iau06arr.apn[i, 2] + iau06arr.apn[i, 3] * ttt) * np.cos(tempval) + (
            iau06arr.apn[i, 6] + iau06arr.apn[i, 7] * ttt
        ) * np.sin(tempval)

    # Planetary nutation constants
    pplnsum = -0.000135 * const.ARCSEC2RAD
    eplnsum = 0.000388 * const.ARCSEC2RAD

    # Combine nutation components
    deltapsi = pnsum + pplnsum
    deltaeps = ensum + eplnsum

    # Build transformation matrices
    nut, prec, pnb = _build_transformation_matrices(ttt, deltaeps, deltapsi, True)

    return deltapsi, pnb, prec, nut, fundargs


########################################################################################
# IAU 2006 XYS Parameters
########################################################################################


def _get_mfme_recnum(jd: float, jdf: float, mjd0: float) -> Tuple[float, int]:
    # Calculate the Julian day at 0000 hr
    jdb = np.floor(jd + jdf) + 0.5
    mfme = (jd + jdf - jdb) * const.DAY2MIN
    if mfme < 0:
        mfme += const.DAY2MIN

    # Find the record number corresponding to the desired day
    recnum = int(np.floor(jd + jdf - mjd0 - const.JD_TO_MJD_OFFSET))

    return mfme, recnum


def findxysparam(
    jdtt: float,
    jdttf: float,
    iau06xysarr: IAU06xysArray,
    interp: InterpolationMode | None = None,
) -> Tuple[float, float, float]:
    """Finds the X, Y, S parameters for a given time with optional interpolation.

    References:
        Vallado, 2013

    Args:
        jdtt (float): Epoch Julian day (days from 4713 BC)
        jdttf (float): Epoch Julian day fraction (day fraction from jdutc)
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        interp (InterpolationMode, optional): Interpolation mode (default: None)

    Returns:
        tuple: (x, y, s)
            x (float): X component of CIO in radians
            y (float): Y component of CIO in radians
            s (float): S component in radians
    """
    # Find the record number corresponding to the desired day
    mfme, recnum = _get_mfme_recnum(jdtt, jdttf, iau06xysarr.mjd_tt[0])

    # Check for out-of-bound values
    if 0 <= recnum <= len(iau06xysarr.x) - 1:
        if interp == InterpolationMode.LINEAR:
            # Perform linear interpolation
            target_time = iau06xysarr.mjd_tt[recnum] + mfme / const.DAY2MIN
            x = np.interp(
                target_time,
                iau06xysarr.mjd_tt[recnum : recnum + 2],
                iau06xysarr.x[recnum : recnum + 2],
            )
            y = np.interp(
                target_time,
                iau06xysarr.mjd_tt[recnum : recnum + 2],
                iau06xysarr.y[recnum : recnum + 2],
            )
            s = np.interp(
                target_time,
                iau06xysarr.mjd_tt[recnum : recnum + 2],
                iau06xysarr.s[recnum : recnum + 2],
            )
        elif interp == InterpolationMode.SPLINE:
            # Perform cubic spline interpolation
            start_idx = max(0, recnum - 1)
            end_idx = min(len(iau06xysarr.x), recnum + 3)
            cs_x = CubicSpline(
                iau06xysarr.mjd_tt[start_idx:end_idx], iau06xysarr.x[start_idx:end_idx]
            )
            cs_y = CubicSpline(
                iau06xysarr.mjd_tt[start_idx:end_idx], iau06xysarr.y[start_idx:end_idx]
            )
            cs_s = CubicSpline(
                iau06xysarr.mjd_tt[start_idx:end_idx], iau06xysarr.s[start_idx:end_idx]
            )
            target_time = iau06xysarr.mjd_tt[recnum] + mfme / const.DAY2MIN
            x = cs_x(target_time).item()
            y = cs_y(target_time).item()
            s = cs_s(target_time).item()
        else:
            # No interpolation
            x = iau06xysarr.x[recnum]
            y = iau06xysarr.y[recnum]
            s = iau06xysarr.s[recnum]
    else:
        # Default values for out-of-bound requests
        x, y, s = 0, 0, 0

    return x, y, s


def create_xys(
    iau06arr: IAU06Array,
    directory: str | None = None,
    filename: str = "xysdata.dat",
    yr_span: int = 1,
    dt_day: int = 1,
    ymdhms: Tuple[int, int, int, int, int, float] = (1957, 1, 1, 0, 0, 0.0),
) -> np.ndarray:
    """Generate the XYS data array and optionally save to a file.

    This function precalculates the XYS parameters and optionally stores in a data file
    for efficient access in the future.

    Args:
        iau06arr (IAU06Array): IAU 2006 data
        directory (str, optional): Directory to save the output file (default: None)
        filename (str, optional): Output filename (default: 'xysdata.dat')
        yr_span (int, optional): Number of years to generate data (default: 142 years)
        dt_day (int, optional): Time step in days (default: 1 day)
        ymdhms (tuple, optional): Initial date in (yr, mo, day, hr, min, sec) format
                                  (default: (1957, 1, 1, 0, 0, 0.0))

    Returns:
        np.ndarray: Array of XYS data with columns [jdtt, jdftt, x, y, s]

    Notes:
        - This is pretty slow due to `iau06xys_series` being called for each day and
          could use some optimization.
        - MATLAB and C# versions hardcode the year start and duration; the default
          duration is set to 1 year here instead of the 142 years used in those

    TODO:
        - Look into using `jit` for performance improvements of downstream functions
    """
    # Initialize the starting Julian date
    jdtt, jdftt = jday(*ymdhms)

    # Calculate the number of rows for the array
    num_rows = yr_span * int(const.YR2DAY) // dt_day + 1

    # Pre-initialize the data array
    xys_data = np.zeros((num_rows, 5))

    # Generate the data
    for i in range(num_rows):
        ttt = (jdtt + jdftt - const.J2000) / const.CENT2DAY
        x, y, s = iau06xys_series(ttt, iau06arr)

        # Store data in the array
        xys_data[i] = [jdtt, jdftt, x, y, s]

        # Increment the Julian date
        jdtt += dt_day

    # Optionally save to a file
    if directory:
        np.savetxt(
            f"{directory}/{filename}",
            xys_data,
            fmt=["%15.6f", "%13.11f", "%15.12f", "%15.12f", "%15.12f"],
            header="jdtt jdftt x y s",
            comments="",
        )

    return xys_data


def iau06xys_series(ttt: float, iau06arr: IAU06Array) -> Tuple[float, float, float]:
    """Calculates the XYS parameters for the IAU2006 CIO theory.

    This is the series implementation of the XYS parameters, which are used to compute
    the Celestial Intermediate Origin (CIO) locator.

    References:
        Vallado, 2022, p. 214-216

    Args:
        ttt (float): Julian centuries of TT
        iau06arr (IAU06Array): IAU 2006 data

    Returns:
        tuple: (x, y, s)
            x (float): Coordinate of CIP in radians
            y (float): Coordinate of CIP in radians
            s (float): Coordinate in radians
    """
    # Fundamental arguments from the IAU 2006 nutation theory
    fundargs = fundarg(ttt, opt="06")

    # Powers of TTT
    ttt2, ttt3, ttt4, ttt5 = ttt**2, ttt**3, ttt**4, ttt**5

    # Limits for the x, y, and s series. These numbers correspond to the ranges of
    # terms used in the calculations for each group:
    # - Group 1: Main series (1306 terms for x, 962 for y, etc.)
    # - Group 2: Secondary contributions (253 terms for x, 277 for y, etc.)
    # - Group 3: Smaller corrections (36 terms for x, 30 for y, etc.)
    # - Group 4: Even smaller corrections (4 terms for x, 5 for y, etc.)
    # - Group 5: Minimal corrections (1 term each for both x and y, 0 for s)

    # Compute X
    limits_x = [1306, 253, 36, 4, 1]  # total sum = 1600 (axs0 and a0xi length)
    x_sums = [0] * len(limits_x)

    # Loop over each group
    for group, limit in enumerate(limits_x):
        start_index = sum(limits_x[:group])
        for i in range(limit):
            idx = start_index + i
            tempval = (
                iau06arr.ax0i[idx, 0] * fundargs.l
                + iau06arr.ax0i[idx, 1] * fundargs.l1
                + iau06arr.ax0i[idx, 2] * fundargs.f
                + iau06arr.ax0i[idx, 3] * fundargs.d
                + iau06arr.ax0i[idx, 4] * fundargs.omega
                + iau06arr.ax0i[idx, 5] * fundargs.lonmer
                + iau06arr.ax0i[idx, 6] * fundargs.lonven
                + iau06arr.ax0i[idx, 7] * fundargs.lonear
                + iau06arr.ax0i[idx, 8] * fundargs.lonmar
                + iau06arr.ax0i[idx, 9] * fundargs.lonjup
                + iau06arr.ax0i[idx, 10] * fundargs.lonsat
                + iau06arr.ax0i[idx, 11] * fundargs.lonurn
                + iau06arr.ax0i[idx, 12] * fundargs.lonnep
                + iau06arr.ax0i[idx, 13] * fundargs.precrate
            )
            x_sums[group] += iau06arr.ax0[idx, 0] * np.sin(tempval) + iau06arr.ax0[
                idx, 1
            ] * np.cos(tempval)

    # Final value for x
    x = (
        -0.016617
        + 2004.191898 * ttt
        - 0.4297829 * ttt2
        - 0.19861834 * ttt3
        - 0.000007578 * ttt4
        + 0.0000059285 * ttt5
    )
    x = x * const.ARCSEC2RAD + sum(x_sums * np.array([1, ttt, ttt2, ttt3, ttt4]))

    # Compute Y
    limits_y = [962, 277, 30, 5, 1]  # total sum = 1275 (ays0 and a0yi length)
    y_sums = [0] * len(limits_y)

    # Loop over each group
    for group, limit in enumerate(limits_y):
        start_index = sum(limits_y[:group])
        for i in range(limit):
            idx = start_index + i
            tempval = (
                iau06arr.ay0i[idx, 0] * fundargs.l
                + iau06arr.ay0i[idx, 1] * fundargs.l1
                + iau06arr.ay0i[idx, 2] * fundargs.f
                + iau06arr.ay0i[idx, 3] * fundargs.d
                + iau06arr.ay0i[idx, 4] * fundargs.omega
                + iau06arr.ay0i[idx, 5] * fundargs.lonmer
                + iau06arr.ay0i[idx, 6] * fundargs.lonven
                + iau06arr.ay0i[idx, 7] * fundargs.lonear
                + iau06arr.ay0i[idx, 8] * fundargs.lonmar
                + iau06arr.ay0i[idx, 9] * fundargs.lonjup
                + iau06arr.ay0i[idx, 10] * fundargs.lonsat
                + iau06arr.ay0i[idx, 11] * fundargs.lonurn
                + iau06arr.ay0i[idx, 12] * fundargs.lonnep
                + iau06arr.ay0i[idx, 13] * fundargs.precrate
            )
            y_sums[group] += iau06arr.ay0[idx, 0] * np.sin(tempval) + iau06arr.ay0[
                idx, 1
            ] * np.cos(tempval)

    # Final value for y
    y = (
        -0.006951
        - 0.025896 * ttt
        - 22.4072747 * ttt2
        + 0.00190059 * ttt3
        + 0.001112526 * ttt4
        + 0.0000001358 * ttt5
    )
    y = y * const.ARCSEC2RAD + sum(y_sums * np.array([1, ttt, ttt2, ttt3, ttt4]))

    # Compute S
    limits_s = [33, 3, 25, 4, 1]  # total sum = 66 (ass0 and a0si length)
    s_sums = [0] * len(limits_s)

    # Loop over each group
    for group, limit in enumerate(limits_s):
        start_index = sum(limits_s[:group])
        for i in range(limit):
            idx = start_index + i
            tempval = (
                iau06arr.as0i[idx, 0] * fundargs.l
                + iau06arr.as0i[idx, 1] * fundargs.l1
                + iau06arr.as0i[idx, 2] * fundargs.f
                + iau06arr.as0i[idx, 3] * fundargs.d
                + iau06arr.as0i[idx, 4] * fundargs.omega
                + iau06arr.as0i[idx, 5] * fundargs.lonmer
                + iau06arr.as0i[idx, 6] * fundargs.lonven
                + iau06arr.as0i[idx, 7] * fundargs.lonear
                + iau06arr.as0i[idx, 8] * fundargs.lonmar
                + iau06arr.as0i[idx, 9] * fundargs.lonjup
                + iau06arr.as0i[idx, 10] * fundargs.lonsat
                + iau06arr.as0i[idx, 11] * fundargs.lonurn
                + iau06arr.as0i[idx, 12] * fundargs.lonnep
                + iau06arr.as0i[idx, 13] * fundargs.precrate
            )
            s_sums[group] += iau06arr.as0[idx, 0] * np.sin(tempval) + iau06arr.as0[
                idx, 1
            ] * np.cos(tempval)

    # Final value for s
    s = (
        0.000094
        + 0.00380865 * ttt
        - 0.00012268 * ttt2
        - 0.07257411 * ttt3
        + 0.00002798 * ttt4
        + 0.00001562 * ttt5
    )
    s = (
        -x * y * 0.5
        + s * const.ARCSEC2RAD
        + sum(s_sums * np.array([1, ttt, ttt2, ttt3, ttt4]))
    )

    return x, y, s


def iau06xys(
    ttt: float,
    iau06arr: IAU06Array,
    ddx: float = 0.0,
    ddy: float = 0.0,
    iau06xysarr: IAU06xysArray | None = None,
    use_full_series: bool = True,
) -> Tuple[float, float, float, np.ndarray]:
    """Calculates the transformation matrix that accounts for the effects of
    precession-nutation using the IAU2006 theory.

    References:
        Vallado, 2022, pp. 214, 221

    Args:
        ttt (float): Julian centuries of TT
        iau06arr (IAU06Array): IAU 2006 data
        ddx (float, optional): EOP correction for x in radians (default is 0)
        ddy (float, optional): EOP correction for y in radians (default is 0)
        iau06xysarr (IAU06xysArray, optional): IAU 2006 XYS data (default is None)
        use_full_series (bool, optional): Whether to use the full series implementation
                                          for XYS parameters (default is True)

    Returns:
        tuple: (x, y, s, pn)
            x (float): Coordinate of CIP in radians
            y (float): Coordinate of CIP in radians
            s (float): Coordinate in radians
            pn (np.ndarray): Transformation matrix for TIRS-GCRF
    """
    # Calculate X, Y, and S components
    if use_full_series:
        # Use the full series implementation
        x, y, s = iau06xys_series(ttt, iau06arr)
    else:
        # Check that the XYS array is provided
        if not iau06xysarr:
            raise ValueError("IAU 2006 XYS array must be provided for interpolation!")

        # Find the X, Y, and S components using spline interpolation
        # TODO: allow for user to specify interpolation method?
        jdtt = ttt * const.CENT2DAY + const.J2000
        x, y, s = findxysparam(jdtt, 0, iau06xysarr, InterpolationMode.SPLINE)

    # Apply any corrections for x and y
    x += ddx
    y += ddy

    # Calculate the 'a' parameter based on x and y
    a = 0.5 + 0.125 * (x**2 + y**2)

    # Build nutation matrices
    nut1 = np.array(
        [
            [1 - a * x**2, -a * x * y, x],
            [-a * x * y, 1 - a * y**2, y],
            [-x, -y, 1 - a * (x**2 + y**2)],
        ]
    )
    nut2 = np.array([[np.cos(s), np.sin(s), 0], [-np.sin(s), np.cos(s), 0], [0, 0, 1]])

    # Combine to form the final transformation matrix
    pn = np.dot(nut1, nut2)

    return x, y, s, pn


########################################################################################
# Earth Orientation Parameters (EOP)
########################################################################################


def findeopparam(
    jd: float, jdf: float, eoparr: EOPArray, interp: InterpolationMode | None = None
) -> Tuple[float, float, float, float, float, float, float, float, float]:
    """Finds the EOP parameters for a given time with optional interpolation.

    References:
        Vallado, 2013

    Args:
        jd (float): Epoch Julian day (days from 4713 BC)
        jdf (float): Epoch Julian day fraction (day fraction from jdutc)
        eoparr (EOPArray): EOP data
        interp (InterpolationMode, optional): Interpolation mode (default: None)

    Returns:
        tuple: (dut1, dat, lod, xp, yp, ddpsi, ddeps, dx, dy)
            dut1 (float): Julian date of UT1 (days from 4713 BC)
            dat (int): TAI - UTC in seconds
            lod (float): Length of day in seconds
            xp (float): Polar motion coefficient in radians
            yp (float): Polar motion coefficient in radians
            ddpsi (float): Delta psi (nutation in longitude) correction in radians
            ddeps (float): Delta epsilon (nutation in obliquity) correction in radians
            dx (float): Celestial pole (CIP) x offset in radians
            dy (float): Celestial pole (CIP) y offset in radians
    """
    # Find the record number corresponding to the desired day
    mfme, recnum = _get_mfme_recnum(jd, jdf, eoparr.mjd[0])

    # Ensure recnum is within valid bounds
    if 0 <= recnum < len(eoparr.mjd) - 1:
        mjd = eoparr.mjd
        params = np.vstack(
            [
                eoparr.dut1,
                eoparr.dat,
                eoparr.lod,
                eoparr.xp,
                eoparr.yp,
                eoparr.ddpsi,
                eoparr.ddeps,
                eoparr.dx,
                eoparr.dy,
            ]
        ).T

        if interp == InterpolationMode.LINEAR:
            # Linear interpolation
            fixf = mfme / const.DAY2MIN
            weights = [1 - fixf, fixf]
            interp_params = params[recnum : recnum + 2].T @ weights
        elif interp == InterpolationMode.SPLINE:
            # Cubic spline interpolation
            cs = CubicSpline(
                mjd[recnum - 1 : recnum + 3], params[recnum - 1 : recnum + 3], axis=0
            )
            interp_params = cs(mjd[recnum] + mfme / const.DAY2MIN)
        else:
            # No interpolation
            interp_params = params[recnum]
    else:
        # Default values for out-of-bounds requests
        interp_params = np.zeros(9)

    # Convert units for certain parameters
    interp_params[1] = int(interp_params[1])
    interp_params[3:9] *= const.ARCSEC2RAD

    return tuple(interp_params)
