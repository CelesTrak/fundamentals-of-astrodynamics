# --------------------------------------------------------------------------------------
# Authors: David Vallado, Jeff Beck
# Date: 28 June 2005
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from dataclasses import dataclass
from enum import Enum
from typing import Tuple

import numpy as np
from pydantic import BaseModel

from ... import constants as const
from ...mathtime.julian_date import jday
from ...mathtime.calendar import days_to_mdh
from ..time.sidereal import gstime
from .utils import WGSModel, getgravc


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


class TypeRun(Enum):
    """Character for mode of SGP4 Execution."""

    # fmt: off
    Catalog = "c"       # +/- 1 day from epoch, 20 min steps
    Verification = "v"  # start/stop/timestep from TLE input (Line 2)
    FromJD = "j"        # start/stop/timestep provided from start and stop Julian dates
    Manual = "m"        # custom start/stop/timestep provided


class Classification(Enum):
    """Classification of the satellite."""

    # fmt: off
    Unclassified = "U"
    Classified = "C"


class SatRec(BaseModel):
    radiusearthkm: float = const.RE
    xke: float = 0.0743669161
    mo: float = 0.0
    mdot: float = 0.0
    argpo: float = 0.0
    argpdot: float = 0.0
    nodeo: float = 0.0
    nodedot: float = 0.0
    nodecf: float = 0.0
    cc1: float = 0.0
    cc4: float = 0.0
    t2cof: float = 0.0
    omgcof: float = 0.0
    xmcof: float = 0.0
    eta: float = 0.0
    delmo: float = 0.0
    d2: float = 0.0
    d3: float = 0.0
    d4: float = 0.0
    t3cof: float = 0.0
    t4cof: float = 0.0
    t5cof: float = 0.0
    no: float = 0.0
    ecco: float = 0.0
    inclo: float = 0.0
    isimp: int = 0
    bstar: float = 0.0
    j3oj2: float = 0.0
    gsto: float = 0.0
    xfact: float = 0.0
    xlamo: float = 0.0
    atime: float = 0.0
    error: int = 0
    t: float = 0.0
    aycof: float = 0.0
    xlcof: float = 0.0
    con41: float = 0.0
    x1mth2: float = 0.0
    x7thm1: float = 0.0
    satnum: int = 0
    intldesg: str = ""
    epochyr: int = 0
    epochdays: float = 0.0
    ndot: float = 0.0
    nddot: float = 0.0
    elnum: int = 0
    revnum: int = 0
    no_kozai: float = 0.0
    jdsatepoch: float = 0.0
    jdsatepochf: float = 0.0
    classification: Classification = Classification.Unclassified


class SGP4:
    def __init__(
        self, wgs_model: WGSModel = WGSModel.WGS_84, use_afspc_mode: bool = True
    ):
        self.wgs_model = wgs_model
        self.use_afspc_mode = use_afspc_mode
        self.grav_const = getgravc(wgs_model)
        self.satrec = SatRec()
        self.use_deep_space = False

        # TLE attributes
        self.jdstart_full = None
        self.jdstop_full = None

    @staticmethod
    def preprocess_tle(tle_line1: str, tle_line2: str) -> Tuple[str, str]:
        # Fix line 1 issues
        tle_line1 = list(tle_line1)
        for j in range(10, 16):
            if tle_line1[j] == " ":
                tle_line1[j] = "_"
        if tle_line1[44] == " ":
            tle_line1[43] = tle_line1[44]
        tle_line1[44] = "."
        if tle_line1[7] == " ":
            tle_line1[7] = "U"
        if tle_line1[9] == " ":
            tle_line1[9] = "."
        for j in range(45, 50):
            if tle_line1[j] == " ":
                tle_line1[j] = "0"
        if tle_line1[51] == " ":
            tle_line1[51] = "0"
        if tle_line1[53] != " ":
            tle_line1[52] = tle_line1[53]
        tle_line1[53] = "."
        if tle_line1[62] == " ":
            tle_line1[62] = "0"
        if len(tle_line1) < 68 or tle_line1[67] == " ":
            tle_line1[67] = "0"

        # Fix line 2 issues
        tle_line2 = list(tle_line2)
        tle_line2[25] = "."
        for j in range(26, 33):
            if tle_line2[j] == " ":
                tle_line2[j] = "0"

        # Convert back to strings
        return "".join(tle_line1), "".join(tle_line2)

    def set_jd_from_from_ymdhms(
        self,
        start_ymdhms: Tuple[int, int, int, int, int, float],
        stop_ymdhms: Tuple[int, int, int, int, int, float],
    ):
        jdstart, jdstartf = jday(*start_ymdhms)
        jdstop, jdstopf = jday(*stop_ymdhms)
        self.jdstart_full = jdstart + jdstartf
        self.jdstop_full = jdstop + jdstopf

    def set_jd_from_yr_doy(
        self, start_yr: int, start_doy: float, stop_yr: int, stop_doy: float
    ):
        start_mdhms = days_to_mdh(start_yr, start_doy)
        stop_mdhms = days_to_mdh(stop_yr, stop_doy)

        self.set_jd_from_from_ymdhms((start_yr, *start_mdhms), (stop_yr, *stop_mdhms))

    def twoline2rv(
        self,
        tle_line1: str,
        tle_line2: str,
        typerun: TypeRun = TypeRun.Catalog,
        start: float | None = None,
        stop: float | None = None,
        step: float | None = None,
    ):
        """Parse TLE lines and populate SGP4 variables.

        This function converts the two line element (TLE) set character string data to
        variables and initializes the sgp4 variables. several intermediate variables
        and quantities are determined. The Verification mode permits quick checks of any
        changes to the underlying technical theory and works using a
        modified tle file in which the start, stop, and delta time values are
        included at the end of the second line of data. The Catalog mode simply
        propagates from -1440 to 1440 min from epoch and is useful when performing
        entire catalog runs.

        If using the FromJD mode, the start and stop Julian dates must be set before
        calling this function (see `set_jd_from_from_ymdhms` or `set_jd_from_yr_doy`).

        References:
            - NORAD Spacetrack Report #3
            - Vallado, Crawford, Hujsak, Kelso 2006

        Args:
            tle_line1 (str): First line of the TLE set
            tle_line2 (str): Second line of the TLE set
            typerun (TypeRun): Mode of execution (default = TypeRun.Catalog)
            start (float, optional): Start time in minutes from epoch (default = None)
            stop (float, optional): Stop time in minutes from epoch (default = None)
            step (float, optional): Time step in minutes (default = None)

        Returns:
            tuple (startmfe, stopmfe, deltamin)
                startmfe (float): Start time in minutes from epoch
                stopmfe (float): Stop time in minutes from epoch
                deltamin (float): Time step in minutes
        """
        # Constants
        xpdotp = const.DAY2MIN / const.TWOPI  # rev/day / rad/min

        # Preprocess TLE lines
        tle_line1, tle_line2 = self.preprocess_tle(tle_line1, tle_line2)

        # Parse the first line
        self.satrec.satnum = int(tle_line1[2:7])
        self.satrec.classification = Classification(tle_line1[7])
        self.satrec.intldesg = tle_line1[9:17].strip()
        self.satrec.epochyr = int(tle_line1[18:20])
        self.satrec.epochdays = float(tle_line1[20:32])
        self.satrec.ndot = float(tle_line1[33:43])
        self.satrec.nddot = float(tle_line1[44:50]) * 10 ** int(tle_line1[50:52])
        self.satrec.bstar = float(tle_line1[53:59]) * 10 ** int(tle_line1[59:61])
        self.satrec.elnum = int(tle_line1[64:68])

        # Parse the second line
        self.satrec.inclo = np.radians(float(tle_line2[8:16].strip()))
        self.satrec.nodeo = np.radians(float(tle_line2[17:25].strip()))
        self.satrec.ecco = float(f"0.{tle_line2[26:33].strip()}")
        self.satrec.argpo = np.radians(float(tle_line2[34:42].strip()))
        self.satrec.mo = np.radians(float(tle_line2[43:51].strip()))
        self.satrec.no_kozai = float(tle_line2[52:63].strip()) / xpdotp
        self.satrec.revnum = int(tle_line2[63:68].strip())

        # Convert epoch year to full year
        year = self.satrec.epochyr + 2000 if self.satrec.epochyr < 57 else 1900

        # Adjust ndot and nddot units
        self.satrec.ndot /= xpdotp * const.DAY2MIN  # rad/min^2
        self.satrec.nddot /= xpdotp * const.DAY2MIN**2  # rad/min^3

        # Compute Julian date of the epoch
        mdhms = days_to_mdh(year, self.satrec.epochdays)
        self.satrec.jdsatepoch, self.satrec.jdsatepochf = jday(year, *mdhms)

        # Default values for start, stop, and step
        startmfe, stopmfe, deltamin = 0, const.DAY2MIN, 1

        # Set start, stop, and step based on the type of run
        # Complete catalog evaluation
        if typerun == TypeRun.Catalog:
            startmfe, stopmfe, deltamin = -const.DAY2MIN, const.DAY2MIN, 20

        # Verification - use TLE start/stop/step values
        elif typerun == TypeRun.Verification:
            try:
                startmfe = float(tle_line2[69:81].strip())
                stopmfe = float(tle_line2[82:96].strip())
                deltamin = float(tle_line2[96:105].strip())
            except ValueError:
                raise ValueError("Input TLE does not support verification mode.")

        # From Julian dates (these must be set before calling this function)
        elif typerun == TypeRun.FromJD:
            if any(value is None for value in (self.jdstart_full, self.jdstop_full)):
                raise ValueError(
                    "FromJD mode requires start and stop Julian dates. "
                    "Please set them prior to calling this function!"
                )
            jdsatepoch = self.satrec.jdsatepoch + self.satrec.jdsatepochf
            startmfe = (self.jdstart_full - jdsatepoch) * const.DAY2MIN
            stopmfe = (self.jdstop_full - jdsatepoch) * const.DAY2MIN
            deltamin = step or deltamin

        # Manual mode - use provided start/stop/step values
        elif typerun == TypeRun.Manual:
            startmfe = start or startmfe
            stopmfe = stop or stopmfe
            deltamin = step or deltamin

        # Invalid mode
        else:
            raise ValueError(f"Invalid mode: {typerun}")

        # Initialize SGP4
        self.sgp4init()

        return startmfe, stopmfe, deltamin

    def sgp4init(self):
        """
        Initialize the satellite parameters.
        """
        pass

    def propagate(self, t: float):
        """
        Perform the propagation for the satellite at time t.
        """
        pass


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
