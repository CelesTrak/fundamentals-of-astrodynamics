# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from typing import Tuple

import numpy as np

from .. import constants as const


def safe_sqrt(value: float, context: str = "") -> float:
    """Safe square root function that checks for negative values.

    Args:
        value (float): The value to take the square root of
        context (str, optional): A description of the current step or variable
                                 being used

    Returns:
        float: The square root of the value

    Raises:
        ValueError: If the value is negative
    """
    if value < 0:
        error_message = f"Cannot take square root of negative value: ({value})"
        if context:
            error_message += f"\nContext: {context}"
        raise ValueError(error_message)
    return np.sqrt(value)


def hms2sec(hours: int, minutes: int, seconds: float) -> float:
    """Convert hours, minutes, and seconds to seconds.

    Args:
        hours (int): The number of hours
        minutes (int): The number of minutes
        seconds (float): The number of seconds

    Returns:
        float: The total number of seconds
    """
    return hours * const.HR2SEC + minutes * const.MIN2SEC + seconds


def sec2hms(seconds: float) -> Tuple[int, int, float]:
    """Convert seconds to hours, minutes, and seconds.

    Args:
        seconds (float): The total number of seconds

    Returns:
        tuple: (hours, minutes, seconds)
            hours (int): The number of hours
            minutes (int): The number of minutes
            seconds (float): The number of seconds
    """
    # Get the hours and the fraction of hours
    total_hours = seconds / const.HR2SEC
    hours = int(total_hours)
    hours_fraction = total_hours - hours

    # Get the minutes and seconds
    minutes = int(hours_fraction * const.MIN2SEC)
    secs = (hours_fraction - minutes / const.MIN2SEC) * const.HR2SEC

    # Adjust seconds to avoid floating point errors
    secs = round(secs) if abs(secs - round(secs)) < const.SMALL else secs

    return hours, minutes, secs


def hms2ut(hours: int, minutes: int, seconds: float) -> float:
    """Converts hours, minutes, and seconds into universal time.

    Args:
        hours (int): Hours (0 .. 23)
        minutes (int): Minutes (0 .. 59)
        seconds (float): Seconds (0.0 .. 59.999)

    Returns:
        float: Universal time in hrmin.sec format
    """
    return hours * 100 + minutes + seconds * 0.01


def ut2hms(ut: float) -> Tuple[int, int, float]:
    """Converts universal time (hhmm.sec format) into hours, minutes, and seconds.

    Args:
        ut (float): Universal time in hrmin.sec format

    Returns:
        tuple: (hours, minutes, seconds)
            hours (int): The number of hours (0 .. 23)
            minutes (int): The number of minutes (0 .. 59)
            seconds (float): The number of seconds (0.0 .. 59.999)
    """
    hr = int(np.floor(ut * 0.01))
    minute = int(np.floor(ut - hr * 100))
    second = (ut - hr * 100 - minute) * 100
    return hr, minute, second


def hms2rad(hours: int, minutes: int, seconds: float) -> float:
    """Convert hours, minutes, and seconds to radians.

    References:
        Vallado: 2022, p. 199, Algorithm 19

    Args:
        hours (int): The number of hours
        minutes (int): The number of minutes
        seconds (float): The number of seconds

    Returns:
        float: The total number of radians
    """
    return (hours + minutes / const.MIN2SEC + seconds / const.HR2SEC) * const.HR2RAD


def rad2hms(radians: float) -> Tuple[int, int, float]:
    """Convert radians to hours, minutes, and seconds.

    References:
        Vallado: 2022, p. 199-200, Algorithm 20

    Args:
        radians (float): The total number of radians

    Returns:
        tuple: (hours, minutes, seconds)
            hours (int): The number of hours
            minutes (int): The number of minutes
            seconds (float): The number of seconds
    """
    # Get the total seconds from the radians
    total_seconds = radians / const.HR2RAD * const.HR2SEC

    # Convert the total seconds to hours, minutes, and seconds
    return sec2hms(total_seconds)


def dms2rad(degrees: int, minutes: int, seconds: float) -> float:
    """Convert degrees, minutes, and seconds to radians.

    References:
        Vallado: 2022, p. 198, Algorithm 17

    Args:
        degrees (int): The number of degrees
        minutes (int): The number of minutes
        seconds (float): The number of seconds

    Returns:
        float: The total number of radians
    """
    return float(
        np.radians(degrees + minutes / const.DEG2MIN + seconds / const.DEG2ARCSEC)
    )


def rad2dms(radians: float) -> Tuple[int, int, float]:
    """Convert radians to degrees, minutes, and seconds.

    References:
        Vallado: 2022, p. 199, Algorithm 18

    Args:
        radians (float): The total number of radians

    Returns:
        tuple: (degrees, minutes, seconds)
            degrees (int): The number of degrees
            minutes (int): The number of minutes
            seconds (float): The number of seconds
    """
    # Convert radians to total degrees
    total_degrees = np.degrees(radians)
    degrees = int(total_degrees)
    degrees_fraction = total_degrees - degrees

    # Get minutes and seconds
    minutes = int(degrees_fraction * const.DEG2MIN)
    secs = (degrees_fraction - minutes / const.DEG2MIN) * const.DEG2ARCSEC

    # Adjust seconds to avoid floating-point errors
    secs = round(secs) if abs(secs - round(secs)) < const.SMALL else secs

    return degrees, minutes, secs
