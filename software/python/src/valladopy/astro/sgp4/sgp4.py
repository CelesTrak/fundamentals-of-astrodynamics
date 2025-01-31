# --------------------------------------------------------------------------------------
# Authors: David Vallado, Jeff Beck
# Date: 28 June 2005
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import logging
import math
from dataclasses import dataclass
from enum import Enum
from typing import Tuple

import numpy as np

from ... import constants as const
from ...mathtime.julian_date import jday
from ...mathtime.calendar import days2mdh
from ..time.sidereal import gstime
from .deep_space import DeepSpace
from .utils import SatRec, WGSModel, getgravc


# Set up logging
logger = logging.getLogger(__name__)

# Constants
JD_EPOCH_1950 = 2433281.5


@dataclass
class SGP4InitOutput:
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
    # fmt: off
    """Character for mode of SGP4 Execution."""
    Catalog = "c"       # +/- 1 day from epoch, 20 min steps
    Verification = "v"  # start/stop/timestep from TLE input (Line 2)
    FromJD = "j"        # start/stop/timestep provided from start and stop Julian dates
    Manual = "m"        # custom start/stop/timestep provided


class Classification(Enum):
    Unclassified = "U"
    Classified = "C"


class PropagationError(Enum):
    INVALID_ELEMENTS = 1
    NEGATIVE_MEAN_MOTION = 2
    ECCENTRICITY_OUT_OF_RANGE = 3
    NEGATIVE_SEMILATUS_RECTUM = 4
    ORBITAL_DECAY = 5


class SGP4:
    """Class for running the SGP4 propagator.

    This class provides the necessary functions to propagate a satellite's orbit using
    the Simplified General Perturbations 4 (SGP4) model.

    References:
        - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
        - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
        - Hoots, Schumacher, and Glover, 2004
        - Vallado, Crawford, Hujsak, Kelso, 2006

    Args:
        wgs_model (WGSModel): The WGS model to use (default = WGSModel.WGS_84)
        use_afspc_mode (bool): Flag to use AFSPC mode for GST calculation
                               (default = True)

    Attributes:
        wgs_model (WGSModel): The WGS model to use
        use_afspc_mode (bool): Flag to use AFSPC mode for GST calculation
        grav_const (GravitationalConstants): Gravitational constants for the Earth
        satrec (SatRec): Dataclass for satellite elements
        use_deep_space (bool): Flag to use deep space model
        x2o3 (float): 2/3 constant for deep space model
        jdstart_full (float): Start Julian date for TLE start time calculation
        jdstop_full (float): Stop Julian date for TLE stop time calculation
        sgp4init_out (SGP4InitOutput): Output dataclass for SGP4 initialization
        ds (DeepSpace): Deep space object for deep space model

    TODO:
        - This class could be further refactored/cleaned up for better readability.
    """

    def __init__(
        self, wgs_model: WGSModel = WGSModel.WGS_84, use_afspc_mode: bool = True
    ):
        self.wgs_model = wgs_model
        self.use_afspc_mode = use_afspc_mode
        self.grav_const = getgravc(wgs_model)
        self.satrec = SatRec()
        self.use_deep_space = False
        self.x2o3 = 2 / 3

        # TLE attributes
        self.jdstart_full = None
        self.jdstop_full = None

        # SGP4 attributes
        self.sgp4init_out = None

        # Deep space variables
        self.ds = None

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
        start_mdhms = days2mdh(start_yr, start_doy)
        stop_mdhms = days2mdh(stop_yr, stop_doy)

        self.set_jd_from_from_ymdhms((start_yr, *start_mdhms), (stop_yr, *stop_mdhms))

    def twoline2rv(
        self,
        tle_line1: str,
        tle_line2: str,
        typerun: TypeRun = TypeRun.Catalog,
        start: float | None = None,
        stop: float | None = None,
        step: float | None = None,
    ) -> Tuple[float, float, float, np.ndarray, np.ndarray]:
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

        Args:
            tle_line1 (str): First line of the TLE set
            tle_line2 (str): Second line of the TLE set
            typerun (TypeRun): Mode of execution (default = TypeRun.Catalog)
            start (float, optional): Start time in minutes from epoch (default = None)
            stop (float, optional): Stop time in minutes from epoch (default = None)
            step (float, optional): Time step in minutes (default = None)

        Returns:
            tuple (r_init, v_init, startmfe, stopmfe, deltamin)
                startmfe (float): Start time in minutes from epoch
                stopmfe (float): Stop time in minutes from epoch
                deltamin (float): Time step in minutes
                r_init (np.ndarray): Initial position vector in TEME frame in km
                v_init (np.ndarray): Initial velocity vector in TEME frame in km/s
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
        mdhms = days2mdh(year, self.satrec.epochdays)
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
        epoch = self.satrec.jdsatepoch + self.satrec.jdsatepochf - JD_EPOCH_1950
        r_init, v_init = self.sgp4init(epoch)

        return startmfe, stopmfe, deltamin, r_init, v_init

    def initl(self, epoch: float):
        """Initialize parameters for the SPG4 propagator.

        Args:
            epoch (float): Epoch time in days from Jan 0, 1950, 0 hr

        Returns:
            None
        """
        # Initialize output dataclass
        out = SGP4InitOutput()

        # Calculate auxiliary epoch quantities
        out.eccsq = self.satrec.ecco**2
        out.omeosq = 1 - out.eccsq
        out.rteosq = np.sqrt(out.omeosq)
        out.cosio = np.cos(self.satrec.inclo)
        out.cosio2 = out.cosio**2

        # Un-Kozai the mean motion
        ak = (self.grav_const.xke / self.satrec.no_kozai) ** self.x2o3
        d1 = (
            0.75 * self.grav_const.j2 * (3 * out.cosio2 - 1) / (out.rteosq * out.omeosq)
        )
        delta = d1 / (ak**2)
        adel = ak * (1 - delta**2 - delta * (1 / 3 + 134 * delta**2 / 81))
        delta = d1 / (adel**2)
        out.no_unkozai = self.satrec.no_kozai / (1 + delta)

        # Calculate other terms
        out.ao = (self.grav_const.xke / out.no_unkozai) ** self.x2o3
        out.sinio = np.sin(self.satrec.inclo)
        po = out.ao * out.omeosq
        out.con42 = 1 - 5 * out.cosio2
        out.con41 = -out.con42 - out.cosio2 - out.cosio2
        out.ainv = 1 / out.ao
        out.posq = po**2
        out.rp = out.ao * (1 - self.satrec.ecco)

        # Calculate Greenwich Sidereal Time
        if self.use_afspc_mode:
            out.gsto = gstime(epoch + JD_EPOCH_1950) % const.TWOPI
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
            out.gsto = np.mod(
                thgr70 + c1 * ids70 + c1p2p * tfrac + ts70 * ts70 * fk5r, const.TWOPI
            )

        self.sgp4init_out = out

    def _adjust_perigee(self, ss, qzms2t):
        """Adjusts sfour and qzms24 for perigees below 156 km."""
        perigee = (self.sgp4init_out.rp - 1) * self.grav_const.radiusearthkm

        # Adjust sfour and qzms24 for perigees below 156 km
        sfour, qzms24 = ss, qzms2t
        if perigee < 156:
            sfour = 20 if perigee < 98 else perigee - 78
            qzms24 = ((120 - sfour) / self.grav_const.radiusearthkm) ** 4
            sfour = sfour / self.grav_const.radiusearthkm + 1

        return sfour, qzms24

    def _compute_perturbation_constants(self, coef, coef1, etasq, eeta, psisq, tsi):
        """Computes perturbation constants for SGP4."""
        cc2 = (
            coef1
            * self.satrec.no
            * (
                self.sgp4init_out.ao * (1 + 1.5 * etasq + eeta * (4 + etasq))
                + 0.375
                * self.grav_const.j2
                * tsi
                / psisq
                * self.sgp4init_out.con41
                * (8 + 3 * etasq * (8 + etasq))
            )
        )
        self.satrec.cc1 = self.satrec.bstar * cc2

        # Compute CC3
        cc3 = 0
        if self.satrec.ecco > 1e-4:
            cc3 = (
                -2
                * coef
                * tsi
                * self.grav_const.j3oj2
                * self.satrec.no
                * self.sgp4init_out.sinio
                / self.satrec.ecco
            )

        # Additional short-periodics parameters
        self.satrec.x1mth2 = 1 - self.sgp4init_out.cosio2
        self.satrec.cc4 = (
            2
            * self.satrec.no
            * coef1
            * self.sgp4init_out.ao
            * self.sgp4init_out.omeosq
            * (
                self.satrec.eta * (2 + 0.5 * etasq)
                + self.satrec.ecco * (0.5 + 2 * etasq)
                - self.grav_const.j2
                * tsi
                / (self.sgp4init_out.ao * psisq)
                * (
                    -3
                    * self.sgp4init_out.con41
                    * (1 - 2 * eeta + etasq * (1.5 - 0.5 * eeta))
                    + 0.75
                    * self.satrec.x1mth2
                    * (2 * etasq - eeta * (1 + etasq))
                    * np.cos(2 * self.satrec.argpo)
                )
            )
        )
        self.satrec.cc5 = (
            2
            * coef1
            * self.sgp4init_out.ao
            * self.sgp4init_out.omeosq
            * (1 + 2.75 * (etasq + eeta) + eeta * etasq)
        )

        return cc3

    def _update_motion_rates(self, pinvsq):
        # Compute higher-order terms
        cosio4 = self.sgp4init_out.cosio2**2
        temp1 = 1.5 * self.grav_const.j2 * pinvsq * self.satrec.no
        temp2 = 0.5 * temp1 * self.grav_const.j2 * pinvsq
        temp3 = -0.46875 * self.grav_const.j4 * pinvsq**2 * self.satrec.no

        # Update motion rates
        self.satrec.mdot = (
            self.satrec.no
            + 0.5 * temp1 * self.sgp4init_out.rteosq * self.sgp4init_out.con41
            + 0.0625
            * temp2
            * self.sgp4init_out.rteosq
            * (13 - 78 * self.sgp4init_out.cosio2 + 137 * cosio4)
        )
        self.satrec.argpdot = (
            -0.5 * temp1 * self.sgp4init_out.con42
            + 0.0625 * temp2 * (7 - 114 * self.sgp4init_out.cosio2 + 395 * cosio4)
            + temp3 * (3 - 36 * self.sgp4init_out.cosio2 + 49 * cosio4)
        )

        xhdot1 = -temp1 * self.sgp4init_out.cosio
        self.satrec.nodedot = (
            xhdot1
            + (
                0.5 * temp2 * (4 - 19 * self.sgp4init_out.cosio2)
                + 2 * temp3 * (3 - 7 * self.sgp4init_out.cosio2)
            )
            * self.sgp4init_out.cosio
        )
        xpidot = self.satrec.argpdot + self.satrec.nodedot

        return xhdot1, xpidot

    def _initialize_deep_space(self, epoch, xpidot):
        # Initialize deep space object
        self.ds = DeepSpace(
            epoch,
            self.satrec.ecco,
            self.satrec.inclo,
            self.satrec.nodeo,
            self.satrec.argpo,
            self.satrec.no,
            self.satrec.mo,
            self.use_afspc_mode,
        )

        # Set attributes and define variables
        self.use_deep_space = True
        self.satrec.isimp = True
        tc = argpm = nodem = mm = 0
        inclm = self.satrec.inclo

        # Call dscom function to compute deep-space common variables
        self.ds.dscom(tc)

        # Call dpper function to adjust for perturbations, if necessary
        if not self.satrec.init:
            self.ds.dpper(self.satrec.t)
            self.satrec.ecco = self.ds.ep
            self.satrec.inclo = self.ds.inclp
            self.satrec.nodeo = self.ds.nodep
            self.satrec.argpo = self.ds.argpp
            self.satrec.mo = self.ds.mp

        # Call dsinit function
        self.ds.dsinit(
            self.satrec,
            self.grav_const.xke,
            tc,
            self.sgp4init_out.gsto,
            xpidot,
            self.sgp4init_out.eccsq,
            inclm,
            nodem,
            argpm,
            mm,
        )

    def _initialize_non_deep_space(self, tsi, sfour):
        cc1sq = self.satrec.cc1**2
        self.satrec.d2 = 4 * self.sgp4init_out.ao * tsi * cc1sq
        temp = self.satrec.d2 * tsi * self.satrec.cc1 / 3
        self.satrec.d3 = (17 * self.sgp4init_out.ao + sfour) * temp
        self.satrec.d4 = (
            0.5
            * temp
            * self.sgp4init_out.ao
            * tsi
            * (221 * self.sgp4init_out.ao + 31 * sfour)
            * self.satrec.cc1
        )
        self.satrec.t3cof = self.satrec.d2 + 2 * cc1sq
        self.satrec.t4cof = 0.25 * (
            3 * self.satrec.d3 + self.satrec.cc1 * (12 * self.satrec.d2 + 10 * cc1sq)
        )
        self.satrec.t5cof = 0.2 * (
            3 * self.satrec.d4
            + 12 * self.satrec.cc1 * self.satrec.d3
            + 6 * self.satrec.d2**2
            + 15 * cc1sq * (2 * self.satrec.d2 + cc1sq)
        )

    def sgp4init(
        self, epoch: float, tol: float = const.SMALL
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Initializes variables for SGP4.

        Args:
            epoch (float): Epoch time in days from Jan 0, 1950 0 hr
            tol (float, optional): Tolerance for small values (default = const.SMALL)

        Returns:
            tuple (r_init, v_init)
                r_init (np.ndarray): Initial position vector in TEME frame in km
                v_init (np.ndarray): Initial velocity vector in TEME frame in km/s

        TODO:
            - Define magic numbers
        """
        # Earth constants
        ss = 78 / self.grav_const.radiusearthkm + 1
        qzms2t = ((120 - 78) / self.grav_const.radiusearthkm) ** 4

        # Initialize SGP4 variables
        self.initl(epoch)

        # Calculate derived orbital parameters
        self.satrec.no = self.sgp4init_out.no_unkozai
        self.satrec.a = (self.satrec.no * self.grav_const.tumin) ** (-2 / 3)
        self.satrec.alta = self.satrec.a * (1 + self.satrec.ecco) - 1
        self.satrec.altp = self.satrec.a * (1 - self.satrec.ecco) - 1

        # Ensure valid orbital elements and positive mean motion
        if not (self.sgp4init_out.omeosq >= 0 and self.satrec.no >= 0):
            r_init, v_init = self.propagate(0)
            return r_init, v_init

        # Determine if perigee is less than 220 km
        if self.sgp4init_out.rp < (220 / self.grav_const.radiusearthkm + 1):
            self.satrec.isimp = True

        # Adjust constants for perigee below 156 km
        sfour, qzms24 = self._adjust_perigee(ss, qzms2t)

        # Definitions
        pinvsq = 1 / self.sgp4init_out.posq
        tsi = 1 / (self.sgp4init_out.ao - sfour)
        self.satrec.eta = self.sgp4init_out.ao * self.satrec.ecco * tsi
        etasq = self.satrec.eta**2
        eeta = self.satrec.ecco * self.satrec.eta
        psisq = abs(1 - etasq)
        coef = qzms24 * tsi**4
        coef1 = coef / psisq**3.5

        # Compute perturbation constants
        cc3 = self._compute_perturbation_constants(coef, coef1, etasq, eeta, psisq, tsi)

        # Compute higher-order terms
        xhdot1, xpidot = self._update_motion_rates(pinvsq)

        # Update other coefficients
        self.satrec.omgcof = self.satrec.bstar * cc3 * np.cos(self.satrec.argpo)
        if self.satrec.ecco > 1e-4:
            self.satrec.xmcof = -self.x2o3 * coef * self.satrec.bstar / eeta
        self.satrec.nodecf = 3.5 * self.sgp4init_out.omeosq * xhdot1 * self.satrec.cc1
        self.satrec.t2cof = 1.5 * self.satrec.cc1

        # Handle divide-by-zero for xinco = 180 degrees
        den = tol
        if abs(self.sgp4init_out.cosio + 1) > tol:
            den = 1 + self.sgp4init_out.cosio

        self.satrec.xlcof = (
            -0.25
            * self.grav_const.j3oj2
            * self.sgp4init_out.sinio
            * (3 + 5 * self.sgp4init_out.cosio)
            / den
        )

        self.satrec.aycof = -0.5 * self.grav_const.j3oj2 * self.sgp4init_out.sinio
        self.satrec.delmo = (1 + self.satrec.eta * np.cos(self.satrec.mo)) ** 3
        self.satrec.sinmao = np.sin(self.satrec.mo)
        self.satrec.x7thm1 = 7 * self.sgp4init_out.cosio2 - 1

        # Determine initialization type
        if (const.TWOPI / self.satrec.no) >= 225:
            # Deep space initialization
            self._initialize_deep_space(epoch, xpidot)
        elif not self.satrec.isimp:
            # Non-deep space initialization
            self._initialize_non_deep_space(tsi, sfour)
        else:
            logger.warning(
                "Neither deep space nor non-deep space coefficients were initialized."
            )

        # Propagate to zero epoch
        r_init, v_init = self.propagate(0)

        return r_init, v_init

    def _apply_secular_gravity_drag(self, t):
        """Apply updates for secular gravity and atmospheric drag."""
        # Initialize secular terms
        t2 = t**2
        xmdf = self.satrec.mo + self.satrec.mdot * t
        argpdf = self.satrec.argpo + self.satrec.argpdot * t
        nodedf = self.satrec.nodeo + self.satrec.nodedot * t
        argpm, mm, nodem = argpdf, xmdf, nodedf + self.satrec.nodecf * t2
        nm, em, inclm = self.satrec.no, self.satrec.ecco, self.satrec.inclo
        tempa = 1 - self.satrec.cc1 * t
        tempe = self.satrec.bstar * self.satrec.cc4 * t
        templ = self.satrec.t2cof * t2

        # Check if non-deep space application
        if not self.satrec.isimp:
            t3, t4 = t**3, t**4
            d2, d3, d4 = self.satrec.d2, self.satrec.d3, self.satrec.d4
            delomg = self.satrec.omgcof * t
            delm = self.satrec.xmcof * (
                (1 + self.satrec.eta * np.cos(xmdf)) ** 3 - self.satrec.delmo
            )
            temp = delomg + delm
            mm, argpm = xmdf + temp, argpdf - temp
            tempa = tempa - d2 * t2 - d3 * t3 - d4 * t4
            tempe = tempe + self.satrec.bstar * self.satrec.cc5 * (
                np.sin(mm) - self.satrec.sinmao
            )
            templ = (
                templ
                + self.satrec.t3cof * t3
                + t4 * (self.satrec.t4cof + t * self.satrec.t5cof)
            )

        # Check if deep space application
        if self.use_deep_space:
            # Update orbital elements
            self.ds.set_mean_elems(em=em, inclm=inclm, nm=nm)

            # Add deep space contributions to mean elements for perturbing 3rd body
            self.ds.dspace(self.satrec, t, self.sgp4init_out.gsto)

            # Unpack updated values directly from `dsinit_out`
            em, inclm, nodem, argpm, nm, mm = self.ds.get_mean_elems()

        # Check if mean motion is less than zero
        if nm <= 0:
            logger.error("Mean motion is less than zero!")
            self.satrec.error = PropagationError.NEGATIVE_MEAN_MOTION

        # Compute some orbital elements
        am = (self.grav_const.xke / nm) ** self.x2o3 * tempa**2
        nm = self.grav_const.xke / am**1.5
        em -= tempe

        # Check if elements are within valid ranges
        if (em >= 1) or (em < -0.001) or (am < 0.95):
            logger.error("Mean elements are out of range!")
            self.satrec.error = PropagationError.INVALID_ELEMENTS

        # Update mean elements
        em = max(em, 1e-6)
        mm += self.satrec.no * templ
        xlm = mm + argpm + nodem
        nodem = math.remainder(nodem, const.TWOPI)
        argpm = math.remainder(argpm, const.TWOPI)
        xlm = math.remainder(xlm, const.TWOPI)
        mm = math.remainder(xlm - argpm - nodem, const.TWOPI)

        return nm, em, inclm, nodem, argpm, mm, am, templ

    def _compute_ds_periodics(self, t, ep, xincp, nodep, argpp, mp, tol):
        """Compute deep space long period periodic contributions"""
        # Set input values
        self.ds.set_attributes(ep=ep, inclp=xincp, nodep=nodep, argpp=argpp, mp=mp)

        # Add deep space periodics
        self.ds.dpper(t)

        # Retrieve updated values
        ep, xincp, nodep, argpp, mp = self.ds.get_attributes(
            "ep", "inclp", "nodep", "argpp", "mp"
        )

        # Handle inclination and node adjustments
        if xincp < 0:
            xincp = -xincp
            nodep += np.pi
            argpp -= np.pi

        # Check eccentricity range
        if (ep < 0) or (ep > 1):
            logger.error(f"Perturbed eccentricity is out of range!: {ep: .2f}")
            self.satrec.error = PropagationError.ECCENTRICITY_OUT_OF_RANGE

        # Compute long period periodics
        sinip, cosip = np.sin(xincp), np.cos(xincp)
        self.satrec.aycof = -0.5 * self.grav_const.j3oj2 * sinip
        den = 1 + cosip if abs(cosip + 1) > tol else const.SMALL
        self.satrec.xlcof = (
            -0.25 * self.grav_const.j3oj2 * sinip * (3 + 5 * cosip) / den
        )

        return ep, xincp, nodep, argpp, mp, sinip, cosip

    def _compute_long_periodics(self, t, em, inclm, nodem, argpm, mm, am, tol):
        """Compute lunar-solar and long period periodics."""
        # Initialize periodics
        ep, xincp, argpp, nodep, mp = em, inclm, argpm, nodem, mm
        sinip, cosip = np.sin(inclm), np.cos(inclm)

        # Update for deep space periodics
        if self.use_deep_space:
            ep, xincp, nodep, argpp, mp, sinip, cosip = self._compute_ds_periodics(
                t, ep, xincp, nodep, argpp, mp, tol
            )

        # Additional periodics
        axnl = ep * np.cos(argpp)
        temp = 1 / (am * (1 - ep**2))
        aynl = ep * np.sin(argpp) + temp * self.satrec.aycof
        xl = mp + argpp + nodep + temp * self.satrec.xlcof * axnl

        return axnl, aynl, xl, xincp, nodep, sinip, cosip

    @staticmethod
    def _solve_keplers_equation(xl, nodep, axnl, aynl, n_iter, tol, lim=0.95):
        """Solve Kepler's equation"""
        u = math.remainder(xl - nodep, const.TWOPI)
        eo1, tem5, ktr = u, np.inf, 1
        sineo1, coseo1 = np.sin(eo1), np.cos(eo1)

        while (abs(tem5) >= tol) and (ktr <= n_iter):
            # Compute correction
            tem5 = 1 - coseo1 * axnl - sineo1 * aynl
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5

            # Limit correction to avoid divergence
            if abs(tem5) >= lim:
                tem5 = lim if tem5 > 0 else -lim

            eo1 += tem5
            sineo1, coseo1 = np.sin(eo1), np.cos(eo1)
            ktr += 1

        return sineo1, coseo1

    def _compute_short_periodics(
        self, am, nm, pl, el2, axnl, aynl, coseo1, sineo1, sinip, cosip, nodep, xincp
    ):
        """Compute short period periodics."""
        ecose = axnl * coseo1 + aynl * sineo1
        esine = axnl * sineo1 - aynl * coseo1
        rl = am * (1 - ecose)
        rdotl = np.sqrt(am) * esine / rl
        rvdotl = np.sqrt(pl) / rl
        betal = np.sqrt(1 - el2)
        temp = esine / (1 + betal)
        sinu = am / rl * (sineo1 - aynl - axnl * temp)
        cosu = am / rl * (coseo1 - axnl + aynl * temp)
        su = np.arctan2(sinu, cosu)
        sin2u, cos2u = 2 * sinu * cosu, 1 - 2 * sinu**2
        temp = 1 / pl
        temp1 = 0.5 * self.grav_const.j2 * temp
        temp2 = temp1 * temp

        # Deep space for short period periodics
        if self.use_deep_space:
            cosisq = cosip**2
            self.sgp4init_out.con41 = 3 * cosisq - 1
            self.satrec.x1mth2 = 1 - cosisq
            self.satrec.x7thm1 = 7 * cosisq - 1

        mrt = (
            rl * (1 - 1.5 * temp2 * betal * self.sgp4init_out.con41)
            + 0.5 * temp1 * self.satrec.x1mth2 * cos2u
        )
        su -= 0.25 * temp2 * self.satrec.x7thm1 * sin2u
        xnode = nodep + 1.5 * temp2 * cosip * sin2u
        xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u
        mvt = rdotl - nm * temp1 * self.satrec.x1mth2 * sin2u / self.grav_const.xke
        rvdot = (
            rvdotl
            + nm
            * temp1
            * (self.satrec.x1mth2 * cos2u + 1.5 * self.sgp4init_out.con41)
            / self.grav_const.xke
        )

        return mrt, su, xnode, xinc, mvt, rvdot

    def propagate(
        self, t: float, n_iter: int = 10, tol: float = const.SMALL
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Simplified General Perturbations 4 (SGP4) model from Space Command.

        This is an updated and combined version of SGP4 and SDP4 algorithms, which
        were originally published separately in SpaceTrack Report #3. This version
        follows the methodology from the AIAA 2006 paper describing the history and
        development of the algorithm.

        Note: `sgp4init` must be called prior to running this function!

        Args:
            t (float): Time since epoch in minutes
            n_iter (int, optional): Number of iterations for solving Kepler's equation
                                    (default = 10)
            tol (float, optional): Tolerance for small values (default = const.SMALL)

        Returns:
            tuple: (r, v)
                r (np.ndarray): Position vector in TEME frame in km
                v (np.ndarray): Velocity vector in TEME frame in km/s

        Return codes for `satrec.error` (non-zero indicates an error)
            1 - Mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
            2 - Mean motion < 0.0
            3 - Pert elements, ecc < 0.0  or  ecc > 1.0
            4 - Semi-latus rectum < 0.0
            5 - Satellite has decayed

        TODO:
            - Decide whether to return early if `satrec.error` is non-zero
        """
        # Initialize position and velocity vectors
        r, v = np.zeros(3), np.zeros(3)

        # Compute vkmpersec
        vkmpersec = self.grav_const.radiusearthkm * self.grav_const.xke / const.MIN2SEC

        # Clear error flag and compute time-based quantities
        self.satrec.t, self.satrec.error = t, None

        # Update for secular gravity and atmospheric drag
        nm, em, inclm, nodem, argpm, mm, am, templ = self._apply_secular_gravity_drag(t)

        # Add lunar-solar periodics
        axnl, aynl, xl, xincp, nodep, sinip, cosip = self._compute_long_periodics(
            t, em, inclm, nodem, argpm, mm, am, tol
        )

        # Solve Kepler's equation
        sineo1, coseo1 = self._solve_keplers_equation(
            xl, nodep, axnl, aynl, n_iter, tol
        )

        # Short period preliminary quantities
        el2 = axnl**2 + aynl**2
        pl = am * (1 - el2)

        # Check semi-latus rectum and return if negative
        if pl < 0:
            logger.error("Semi-latus rectum is negative!")
            self.satrec.error = PropagationError.NEGATIVE_SEMILATUS_RECTUM
            return r, v

        # Compute short period periodics
        mrt, su, xnode, xinc, mvt, rvdot = self._compute_short_periodics(
            am, nm, pl, el2, axnl, aynl, coseo1, sineo1, sinip, cosip, nodep, xincp
        )

        # Orientation vectors
        xmx, xmy = -np.sin(xnode) * np.cos(xinc), np.cos(xnode) * np.cos(xinc)
        ux = xmx * np.sin(su) + np.cos(xnode) * np.cos(su)
        uy = xmy * np.sin(su) + np.sin(xnode) * np.cos(su)
        uz = np.sin(xinc) * np.sin(su)
        vx = xmx * np.cos(su) - np.cos(xnode) * np.sin(su)
        vy = xmy * np.cos(su) - np.sin(xnode) * np.sin(su)
        vz = np.sin(xinc) * np.cos(su)

        # Position and velocity vectors
        r = np.array([ux, uy, uz]) * mrt * self.grav_const.radiusearthkm
        v = (np.array([ux, uy, uz]) * mvt + np.array([vx, vy, vz]) * rvdot) * vkmpersec

        # Check for decay condition
        if mrt < 1:
            self.satrec.error = PropagationError.ORBITAL_DECAY

        return r, v
