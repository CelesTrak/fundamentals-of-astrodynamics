# --------------------------------------------------------------------------------------
# Authors: David Vallado, Jeff Beck
# Date: 28 June 2005
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from dataclasses import dataclass

import numpy as np

from ... import constants as const
from ..time.sidereal import gstime


@dataclass
class SGP4InitOutput:
    use_deep_space: bool = False
    ainv: float = 0.0  # may not be used
    ao: float = 0.0
    con41: float = 0.0
    con42: float = 0.0
    cosio: float = 0.0
    cosio2: float = 0.0
    eccsq: float = 0.0
    omeosq: float = 0.0
    posq: float = 0.0
    rp: float = 0.0
    rteosq: float = 0.0
    sinio: float = 0.0
    gsto: float = 0.0
    no_unkozai: float = 0.0


def initl(
    xke: float,
    j2: float,
    epoch: float,
    ecco: float,
    inclo: float,
    no_kozai: float,
    use_afspc_mode: bool = True,
) -> SGP4InitOutput:
    """Initialize parameters for the SPG4 propagator.

    References:
        - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
        - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
        - Hoots, Schumacher, and Glover, 2004
        - Vallado, Crawford, Hujsak, Kelso, 2006

    Args:
        xke (float): Earth gravitational constant in radians/minute  # TODO: check units
        j2 (float): Earth second zonal harmonic
        epoch (float): Epoch time in days from Jan 0, 1950, 0 hr
        ecco (float): Eccentricity
        inclo (float): Inclination of the satellite
        no_kozai (float): Mean motion of the satellite
        use_afspc_mode (bool): Flag to use AFSPC mode (default = True)

    Returns:
        SGP4InitOutput: Dataclass encapsulating the initialized values
    """
    # Initialize output dataclass
    sgp4init_output = SGP4InitOutput()

    # Constants
    x2o3 = 2 / 3

    # Calculate auxiliary epoch quantities
    sgp4init_output.eccsq = ecco**2
    sgp4init_output.omeosq = 1.0 - sgp4init_output.eccsq
    sgp4init_output.rteosq = np.sqrt(sgp4init_output.omeosq)
    sgp4init_output.cosio = np.cos(inclo)
    sgp4init_output.cosio2 = sgp4init_output.cosio**2

    # Un-Kozai the mean motion
    ak = (xke / no_kozai) ** x2o3
    d1 = (
        0.75
        * j2
        * (3 * sgp4init_output.cosio2 - 1)
        / (sgp4init_output.rteosq * sgp4init_output.omeosq)
    )
    delta = d1 / (ak**2)
    adel = ak * (1 - delta**2 - delta * (1 / 3 + 134 * delta**2 / 81))
    delta = d1 / (adel**2)
    sgp4init_output.no_unkozai = no_kozai / (1 + delta)

    # Calculate other terms
    sgp4init_output.ao = (xke / sgp4init_output.no_unkozai) ** x2o3
    sgp4init_output.sinio = np.sin(inclo)
    po = sgp4init_output.ao * sgp4init_output.omeosq
    sgp4init_output.con42 = 1 - 5 * sgp4init_output.cosio2
    sgp4init_output.con41 = (
        -sgp4init_output.con42 - sgp4init_output.cosio2 - sgp4init_output.cosio2
    )
    sgp4init_output.ainv = 1 / sgp4init_output.ao
    sgp4init_output.posq = po**2
    sgp4init_output.rp = sgp4init_output.ao * (1 - ecco)

    # Calculate Greenwich Sidereal Time
    if use_afspc_mode:
        sgp4init_output.gsto = gstime(epoch + 2433281.5) % const.TWOPI
    else:
        # SGP4 fix - use old way of finding GST
        # Count integer number of days from 0 Jan 1970
        ts70 = epoch - 7305
        ids70 = np.floor(ts70 + 1e-8)
        tfrac = ts70 - ids70
        c1 = 1.72027916940703639e-2
        thgr70 = 1.7321343856509374
        fk5r = 5.07551419432269442e-15
        c1p2p = c1 + const.TWOPI
        sgp4init_output.gsto = np.mod(
            thgr70 + c1 * ids70 + c1p2p * tfrac + ts70 * ts70 * fk5r, const.TWOPI
        )

    return sgp4init_output
