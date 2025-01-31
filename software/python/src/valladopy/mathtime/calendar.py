# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import calendar
from datetime import datetime, timedelta
from typing import Dict, List, Tuple

from .. import constants as const


# Constants
SEC2MICROSEC = 1e6


def initialize_time(year: int) -> Dict[str, List]:
    """Initializes time-related information for a given year.

    Args:
        year (int): The year for which the data is initialized
                    (to account for leap years)

    Returns:
        dict[str, List]: A dictionary with three keys:
              - 'lmonth': List of days in each month for the given year
              - 'monthtitle': List of abbreviated month names
              - 'daytitle': List of abbreviated day names (starting with Sunday)
    """
    return {
        "lmonth": [calendar.monthrange(year, i)[1] for i in range(1, 13)],
        "monthtitle": [calendar.month_abbr[i].lower() for i in range(1, 13)],
        "daytitle": [calendar.day_abbr[(i + 6) % 7].lower() for i in range(7)],
    }


def get_int_month(month_str: str):
    """Converts a 3-character month abbreviation to its integer equivalent.

    Args:
        month_str (str): Abbreviation of the month (e.g., 'jan', 'feb')

    Returns:
        int: Integer representation of the month (1 for January, 2 for February, etc.)

    Raises:
        ValueError: If the input string is not a valid month abbreviation.
    """
    try:
        return datetime.strptime(month_str[:3].capitalize(), "%b").month
    except ValueError as e:
        raise ValueError(f"Invalid month abbreviation: {month_str}") from e


def get_int_day(day_str: str) -> int:
    """Converts a 3-character day abbreviation to its integer equivalent.

    Args:
        day_str (str): Abbreviation of the day (e.g., 'sun', 'mon')

    Returns:
        int: Integer representation of the day (1 for Sunday, 2 for Monday, etc.)

    Raises:
        ValueError: If the day abbreviation is invalid
    """
    # Shift days to start with Sunday
    shifted_days = [calendar.day_abbr[(i + 6) % 7].lower() for i in range(7)]
    try:
        return shifted_days.index(day_str.lower()) + 1
    except ValueError as e:
        raise ValueError(f"Invalid day abbreviation: {day_str}") from e


def days2mdh(year: int, days: float) -> Tuple[int, int, int, int, float]:
    """Converts day of the year to month, day, hour, minute, and second.

    Args:
        year (int): Year (e.g., 1900 .. 2100)
        days (float): Day of the year (1.0 .. 366.0)

    Returns:
        tuple: (month, day, hour, minute, second)
            month (int): Integer month (1 .. 12)
            day (int): Integer day (1 .. 31)
            hour (int): Integer hour (0 .. 23)
            minute (int): Integer minute (0 .. 59)
            second (float): Seconds (0.0 .. 59.999999)
    """
    # Get the start of the year
    start_of_year = datetime(year, 1, 1)

    # Calculate the full date and time by adding the fractional days
    delta = timedelta(days=days - 1)  # subtract 1 since days=1.0 means Jan 1
    full_date = start_of_year + delta

    # Extract the components
    month = full_date.month
    day = full_date.day
    hour = full_date.hour
    minute = full_date.minute
    second = full_date.second + full_date.microsecond / SEC2MICROSEC

    return month, day, hour, minute, second


def find_days(
    year: int, month: int, day: int, hour: int, minute: int, second: float
) -> float:
    """Finds the fractional days through a year.

    References:
        Vallado: 2022, p. 202

    Args:
        year (int): Year (1900 .. 2100)
        month (int): Month (1 .. 12)
        day (int): Day (1 .. 28, 29, 30, 31)
        hour (int): Hour (0 .. 23)
        minute (int): Minute (0 .. 59)
        second (float): Second (0.0 .. 59.999)

    Returns:
        float: Day of the year plus fraction of the day
    """
    # Start of the year
    start_of_year = datetime(year, 1, 1)

    # Current date and time
    current_date = datetime(
        year, month, day, hour, minute, int(second), int((second % 1) * SEC2MICROSEC)
    )

    # Compute the difference in days
    delta = current_date - start_of_year

    # Fractional days
    return (
        delta.days
        + delta.seconds / const.DAY2SEC
        + delta.microseconds / const.DAY2SEC / SEC2MICROSEC
    ) + 1  # add 1 since days=1.0 means Jan 1
