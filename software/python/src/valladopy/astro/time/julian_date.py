# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------


import numpy as np

from ... import constants as const


def jday(year, month, day, hour, minute, second):
    """Calculate the Julian Date (JD) and fractional day

    This function finds the Julian Date given the year, month, day, and time.

    References:
        Vallado: 2007, p. 189, Algorithm 14, Example 3-14

    Args:
        year (int): Year (e.g., 2024)
        month (int): Month (1 to 12)
        day (int): Day (1 to 31)
        hour (int): Universal Time hour (0 to 23)
        minute (int): Universal Time minute (0 to 59)
        second (float): Universal Time second (0.0 to 59.999)

    Returns:
        jd (float): Julian Date
        jd_frac (float): Fractional part of the Julian Date
    """
    # Calculate Julian Date
    jd = (
        367.0 * year
        - np.floor((7 * (year + np.floor((month + 9) / 12.0))) * 0.25)
        + np.floor(275 * month / 9.0)
        + day
        + 1721013.5
    )  # Use - 678987.0 to go to MJD directly

    # Calculate fractional part of the day
    jd_frac = (second + minute * const.MIN2SEC + hour * const.HR2SEC) / const.DAY2SEC

    # Adjust if jd_frac > 1
    if jd_frac > 1.0:
        jd += np.floor(jd_frac)
        jd_frac -= np.floor(jd_frac)

    return jd, jd_frac
