# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import calendar
import logging
from enum import Enum
from typing import Tuple

import numpy as np

from .calendar import days2mdh
from .utils import hms2sec, sec2hms
from .. import constants as const


# Initialize the logger
logger = logging.getLogger(__name__)


# Constants
JULIAN_DATE_REFERENCE_YEAR = 1900  # reference year for the Julian Date
JULIAN_DATE_1900 = 2415019.5  # Julian date for January 0, 1900
JULIAN_DATE_EPOCH_OFFSET = 4716  # Julian date offset from the Gregorian calendar
DST_RULE_CHANGE_YEAR = 2007  # year when DST rules changed in the US
DST_CHANGE_UTC_HOUR = -10  # hour when DST changes in UTC (02:00)


class CalendarType(Enum):
    JULIAN = "j"
    GREGORIAN = "g"


def jd2sse(julian_date: float) -> float:
    """Converts Julian Date to seconds since epoch.

    Args:
        julian_date (float): The Julian Date (days from 4713 BC)

    Returns:
        float: Seconds since epoch (1 Jan 2000 00:00:00 UTC)
    """
    return (julian_date - const.J2000_UTC) * const.DAY2SEC


def jday(
    year: int, month: int, day: int, hour: int = 0, minute: int = 0, second: float = 0.0
) -> Tuple[float, float]:
    """Calculate the Julian Date (JD) and fractional day.

    This function finds the Julian Date given the year, month, day, and time.

    References:
        Vallado: 2022, p. 184-185, Algorithm 14

    Args:
        year (int): Year (e.g., 2024)
        month (int): Month (1 to 12)
        day (int): Day (1 to 31)
        hour (int, optional): Universal Time hour (0 to 23) (defaults to 0)
        minute (int, optional): Universal Time minute (0 to 59) (defaults to 0)
        second (float, optional): Universal Time second (0.0 to 59.999)
                                  (defaults to 0.0)

    Returns:
        tuple: (jd, jd_frac)
            jd (float): Julian Date
            jd_frac (float): Fractional part of the Julian Date
    """
    # Calculate Julian Date
    jd = (
        367 * year
        - np.floor((7 * (year + np.floor((month + 9) / 12))) * 0.25)
        + np.floor(275 * month / 9)
        + day
        + 1721013.5
    )  # use - 678987.0 to go to MJD directly

    # Calculate fractional part of the day
    jd_frac = (second + minute * const.MIN2SEC + hour * const.HR2SEC) / const.DAY2SEC

    # Adjust if jd_frac > 1
    if jd_frac > 1:
        jd += np.floor(jd_frac)
        jd_frac -= np.floor(jd_frac)

    return jd, jd_frac


def jday_all(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    calendar_type: CalendarType = CalendarType.GREGORIAN,
) -> Tuple[float, float]:
    """Calculates the Julian Date and fractional Julian Day.

    The Julian Date is defined by each elapsed day since noon, Jan 1, 4713 BC.

    References:
        Vallado: 2022, p. 185

    Args:
        year (int): Year (e.g., 1900 .. 2100)
        month (int): Month (1 .. 12)
        day (int): Day (1 .. 31)
        hour (int): Universal Time hour (0 .. 23)
        minute (int): Universal Time minute (0 .. 59)
        second (float): Universal Time second (0.0 .. 59.999)
        calendar_type (CalendarType, optional): Calendar type (Julian or Gregorian)
                                                (Defaults to CalendarType.GREGORIAN)

    Returns:
        tuple: (jd, jdfrac)
            jd (float): Julian Date
            jdfrac (float): Fractional part of the Julian Date
    """
    if month <= 2:
        year -= 1
        month += 12

    # Determine B based on the calendar type
    if calendar_type == CalendarType.JULIAN:
        b = 0
    elif calendar_type == CalendarType.GREGORIAN:
        b = 2 - (year // 100) + (year // 400)
    else:
        raise ValueError(
            "Invalid calendar type. Must be either CalendarType.JULIAN or"
            "CalendarType.GREGORIAN"
        )

    # Compute Julian Date
    jd = (
        int(const.YR2DAY * (year + JULIAN_DATE_EPOCH_OFFSET))
        + int(30.6001 * (month + 1))
        + day
        + b
        - 1524.5
    )

    # Compute fractional day
    jdfrac = (hour * const.HR2SEC + minute * const.MIN2SEC + second) / const.DAY2SEC

    # Normalize jdfrac if it exceeds 1.0
    if jdfrac >= 1:
        jd += int(jdfrac)
        jdfrac %= 1

    return jd, jdfrac


def invjday(jd: float, jdfrac: float) -> Tuple[int, int, int, int, int, float]:
    """Converts Julian Date and fractional day to calendar date and time.

    References:
        Vallado: 2022, p. 203-204, Algorithm 22

    Args:
        jd (float): Julian Date (days from 4713 BC)
        jdfrac (float): Fractional part of the Julian Date

    Returns:
        tuple: (year, month, day, hour, minute, second)
            year (int): Year
            month (int): Month
            day (int): Day
            hour (int): Hour
            minute (int): Minute
            second (float): Second

    Notes:
        - This assumes the Gregorian calendar type.
    """
    # Normalize jdfrac if it spans multiple days
    if abs(jdfrac) >= 1:
        jd += int(jdfrac)
        jdfrac %= 1

    # Adjust for fraction of a day in the Julian Date
    dt = jd - int(jd) - 0.5
    if abs(dt) > 1e-8:
        jd -= dt
        jdfrac += dt

    # Compute year and day of year
    temp = jd - JULIAN_DATE_1900
    tu = temp / const.YR2DAY
    year = JULIAN_DATE_REFERENCE_YEAR + int(tu)
    leap_years = (year - (JULIAN_DATE_REFERENCE_YEAR + 1)) // 4
    days = int(
        temp
        - ((year - JULIAN_DATE_REFERENCE_YEAR) * np.floor(const.YR2DAY) + leap_years)
    )

    # Handle start-of-year edge case
    if days + jdfrac < 1:
        year -= 1
        leap_years = (year - (JULIAN_DATE_REFERENCE_YEAR + 1)) // 4
        days = int(
            temp
            - (
                (year - JULIAN_DATE_REFERENCE_YEAR) * np.floor(const.YR2DAY)
                + leap_years
            )
        )

    # Convert days of year + fractional day to calendar date and time
    month, day, hour, minute, second = days2mdh(year, days + jdfrac)

    return year, month, day, hour, minute, second


def day_of_week(jd: float) -> int:
    """Finds the day of the week for a given Julian date.

    References:
        Vallado: 2022, p. 183, Eq. 3-42

    Args:
        jd (float): Julian date (days from 4713 BC)

    Returns:
        int: Day of the week (1 for Sunday, 2 for Monday, etc.)
    """
    # Ensure the Julian date corresponds to 0.0 hours of the day
    jd = int(jd + 0.5)

    # Calculate the day of the week (1 = Sunday, ..., 7 = Saturday)
    return ((jd + 1) % 7) + 1


def find_dst_date(year: int, month: int, target_week: int, weekday: int) -> int | None:
    """Find the nth occurrence or the last occurrence of a specific weekday in a given
    month.

    Args:
        year (int): Year of interest
        month (int): Month of interest
        target_week (int): Week index (1-based for nth week, -1 for last week)
        weekday (int): Day of the week (0=Monday, 6=Sunday)

    Returns:
        int: Day of the month for the specified weekday and week (None if invalid)
    """
    weeks = calendar.monthcalendar(year, month)
    if target_week == -1:  # last occurrence
        for week in reversed(weeks):
            if week[weekday] != 0:  # non-zero means the day exists
                return week[weekday]
    else:  # nth occurrence
        count = 0
        for week in weeks:
            if week[weekday] != 0:  # non-zero means the day exists
                count += 1
                if count == target_week:
                    return week[weekday]
    return None


def daylight_savings(year: int, lon: float) -> Tuple[int, int, float, float]:
    """Find the start and stop dates for daylight savings time (DST) in a given year.

    References:
        Vallado: 2022, p. 183

    This function uses U.S.-specific DST rules:
        - Before 2007: DST starts on the first Sunday of April and ends on the last
          Sunday of October.
        - From 2007 onward: DST starts on the second Sunday of March and ends on the
          first Sunday of November.

    Args:
        year (int): The year (1900 .. 2100)
        lon (float): Longitude of the site (WEST is negative) in radians

    Returns:
        tuple: (startday, stopday, jdstartdst, jdstopdst)
            startday (int): Day in March when DST starts
            stopday (int): Day in November when DST ends
            jdstartdst (float): Julian date of DST start
            jdstopdst (float): Julian date of DST end
    """
    # Check if the longitude corresponds to a U.S. time zone
    lon_deg = np.degrees(lon)
    if lon_deg < -125 or lon_deg > -66:
        logger.warning(
            "Longitude does not correspond to a U.S. time zone. DST dates may not align"
            "with local transitions."
        )

    # Determine the site zone (0.0 gives Greenwich time)
    zone = int(lon_deg / const.DEG2HR)
    if zone > 0:
        zone -= const.DAY2HR

    # Find the start and stop days for DST
    if year < DST_RULE_CHANGE_YEAR:
        # Pre-2007: DST starts 1st Sunday in April, ends last Sunday in October
        startday = find_dst_date(year, 4, 1, calendar.SUNDAY)
        stopday = find_dst_date(year, 10, -1, calendar.SUNDAY)
    else:
        # 2007 and later: DST starts 2nd Sunday in March, ends 1st Sunday in November
        startday = find_dst_date(year, 3, 2, calendar.SUNDAY)
        stopday = find_dst_date(year, 11, 1, calendar.SUNDAY)

    # Check for invalid dates
    if startday is None or stopday is None:
        logger.warning("Invalid DST dates. Using default values.")
        startday = 1 if year < DST_RULE_CHANGE_YEAR else 8
        stopday = 31

    # Julian date for DST start
    jdstartdst, jdstartdst_frac = jday(
        year, 3 if year >= DST_RULE_CHANGE_YEAR else 4, startday, 12, 0, 0
    )
    jdstartdst += jdstartdst_frac
    jdstartdst += (DST_CHANGE_UTC_HOUR - zone) / const.DAY2HR

    # Julian date for DST end
    jdstopdst, jdstopdst_frac = jday(
        year, 11 if year >= DST_RULE_CHANGE_YEAR else 10, stopday, 12, 0, 0
    )
    jdstopdst += jdstopdst_frac
    jdstopdst += (DST_CHANGE_UTC_HOUR - zone) / const.DAY2HR

    return startday, stopday, jdstartdst, jdstopdst


def convtime(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    timezone: int,
    dut1: float,
    dat: float,
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
    float,
    float,
    float,
    float,
]:
    """Convert UTC or UT1 to various time systems.

    This function finds the time parameters and Julian Century values for inputs of UTC
    or UT1. Since calucations are in UTC, you must include timezone if you enter a local
    time; otherwise it should be zero.

    References:
        Vallado: 2022, p. 196, Algorithm 16

    Args:
        year (int): Year (1900 .. 2100)
        month (int): Month (1 .. 12)
        day (int): Day (1 .. 31)
        hour (int): Universal Time hour (0 .. 23)
        minute (int): Universal Time minute (0 .. 59)
        second (float): Universal Time second (0.0 .. 59.999)
        timezone (int): Offset to UTC from local time (0 .. 23 hr)
        dut1 (float): Delta of UT1 - UTC in seconds
        dat (float): Delta of TAI - UTC in seconds

    Returns:
        tuple: (ut1, tut1, jdut1, jdut1frac, utc, tai, gps, tt, ttt, jdtt, jdttfrac,
                tdb, ttdb, jdtdb, jdtdbfrac)
            ut1 (float): UT1 in seconds (from start of day)
            tut1 (float): Julian centuries since J2000 for UT1
            jdut1 (float): Julian Date for UT1 (days only, days from 4713 BC)
            jdut1frac (float): Julian Date for UT1 (fractional part)
            utc (float): Coordinated Universal Time in seconds
            tai (float): International Atomic Time in seconds
            gps (float): GPS Time in seconds
            tt (float): Terrestrial Time in seconds (from start of day)
            ttt (float): Julian centuries since J2000 for TT
            jdtt (float): Julian Date for TT (days only, days from 4713 BC)
            jdttfrac (float): Julian Date for TT (fractional part)
            tdb (float): Terrestrial Barycentric Time in seconds (from start of day)
            ttdb (float): Julian centuries since J2000 for TDB
            jdtdb (float): Julian Date for TDB (days only, days from 4713 BC)
            jdtdbfrac (float): Julian Date for TDB (fractional part)
    """
    # UTC in seconds from the start of the day
    local_hour = timezone + hour
    utc = hms2sec(local_hour, minute, second)

    # UT1
    ut1 = utc + dut1
    hrtemp, mintemp, sectemp = sec2hms(ut1)
    jdut1, jdut1frac = jday(year, month, day, hrtemp, mintemp, sectemp)
    tut1 = (jdut1 + jdut1frac - const.J2000) / const.CENT2DAY

    # TAI
    tai = utc + dat

    # GPS
    gps = tai + 19  # seconds

    # TT
    tt = tai + 32.184  # seconds
    hrtemp, mintemp, sectemp = sec2hms(tt)
    jdtt, jdttfrac = jday(year, month, day, hrtemp, mintemp, sectemp)
    ttt = (jdtt + jdttfrac - const.J2000) / const.CENT2DAY

    # TDB
    tdb = (
        tt
        + 0.001657 * np.sin(628.3076 * ttt + 6.2401)
        + 0.000022 * np.sin(575.3385 * ttt + 4.297)
        + 0.000014 * np.sin(1256.6152 * ttt + 6.1969)
        + 0.000005 * np.sin(606.9777 * ttt + 4.0212)
        + 0.000005 * np.sin(52.9691 * ttt + 0.4444)
        + 0.000002 * np.sin(21.3299 * ttt + 5.5431)
        + 0.00001 * ttt * np.sin(628.3076 * ttt + 4.249)
    )
    hrtemp, mintemp, sectemp = sec2hms(tdb)
    jdtdb, jdtdbfrac = jday(year, month, day, hrtemp, mintemp, sectemp)
    ttdb = (jdtdb + jdtdbfrac - const.J2000) / const.CENT2DAY

    return (
        ut1,
        tut1,
        jdut1,
        jdut1frac,
        utc,
        tai,
        gps,
        tt,
        ttt,
        jdtt,
        jdttfrac,
        tdb,
        ttdb,
        jdtdb,
        jdtdbfrac,
    )
