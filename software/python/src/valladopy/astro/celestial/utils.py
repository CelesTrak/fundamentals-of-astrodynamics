# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 31 Oct 2003
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import logging
from enum import Enum
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const


# Set up logging
logger = logging.getLogger(__name__)


class EarthModel(Enum):
    SPHERICAL = "s"
    ELLIPSOIDAL = "e"


def in_sight(
    r1: ArrayLike, r2: ArrayLike, earth_model: EarthModel = EarthModel.ELLIPSOIDAL
) -> bool:
    """Determines if there is line-of-sight (LOS) between two satellites, considering
    the Earth's shape.

    References:
        Vallado: 2022, p. 312-315, Algorithm 35

    Args:
        r1 (array_like): Position vector of the first satellite in km
        r2 (array_like): Position vector of the second satellite in km
        earth_model (EarthModel, optional): Earth model to use (default is ELLIPSOIDAL)

    Returns:
        bool: True if there is line-of-sight, False otherwise
    """
    # Magnitudes
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)

    # Scale z-components for ellipsoidal Earth
    temp = (
        1 / np.sqrt(1 - const.ECCEARTHSQRD)
        if earth_model == EarthModel.ELLIPSOIDAL
        else 1
    )
    tr1 = np.array([r1[0], r1[1], r1[2] * temp])
    tr2 = np.array([r2[0], r2[1], r2[2] * temp])

    # Compute magnitudes and dot product
    asqrd = magr1**2
    bsqrd = magr2**2
    adotb = np.dot(tr1, tr2)

    # Compute minimum parametric value
    if abs(asqrd + bsqrd - 2 * adotb) < 1e-4:
        tmin = 0
    else:
        tmin = (asqrd - adotb) / (asqrd + bsqrd - 2 * adotb)
    logger.debug(f"Minimum parametric value (tmin): {tmin}")

    # Check line-of-sight (LOS)
    if tmin < 0 or tmin > 1:
        return True
    else:
        distsqrd = ((1 - tmin) * asqrd + adotb * tmin) / const.RE**2
        return True if distsqrd > 1 else False


def sun_ecliptic_parameters(t: float) -> Tuple[float, float, float]:
    """Compute the mean longitude, mean anomaly, and ecliptic longitude of the Sun.

    Args:
        t (float): Time since J2000 in Julian centuries (e.g. 'tut1' or 'ttdb')

    Returns:
        tuple: (mean_lon, mean_anomaly, ecliptic_lon)
            mean_lon (float): Mean longitude of the Sun in radians
            mean_anomaly (float): Mean anomaly of the Sun in radians
            ecliptic_lon (float): Ecliptic longitude of the Sun in radians
    """
    mean_lon = np.radians(280.4606184 + 36000.77005361 * t) % const.TWOPI
    mean_anomaly = np.radians(357.5277233 + 35999.05034 * t) % const.TWOPI
    ecliptic_lon = (
        np.radians(
            np.degrees(mean_lon)
            + 1.914666471 * np.sin(mean_anomaly)
            + 0.019994643 * np.sin(2 * mean_anomaly)
        )
        % const.TWOPI
    )

    return float(mean_lon), float(mean_anomaly), ecliptic_lon


def obliquity_ecliptic(t: float) -> float:
    """Compute the obliquity of the ecliptic.

    Args:
        t (float): Time since J2000 in Julian centuries (e.g. 'tut1' or 'ttdb')

    Returns:
        float: Obliquity of the ecliptic in radians
    """
    return float(np.radians(np.degrees(const.OBLIQUITYEARTH) - 0.0130042 * t))
