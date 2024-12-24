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

from ... import constants as const
from ...mathtime.julian_date import jday
from ...mathtime.calendar import days_to_mdh
from ..time.sidereal import gstime
from .deep_space import DeepSpace
from .utils import SatRec, WGSModel, getgravc


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
    # fmt: off
    """Classification of the satellite."""
    Unclassified = "U"
    Classified = "C"


class SGP4:
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
        epoch = self.satrec.jdsatepoch + self.satrec.jdsatepochf - JD_EPOCH_1950
        self.sgp4init(epoch)

        return startmfe, stopmfe, deltamin

    def initl(self, epoch: float):
        """Initialize parameters for the SPG4 propagator.

        References:
            - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
            - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
            - Hoots, Schumacher, and Glover, 2004
            - Vallado, Crawford, Hujsak, Kelso, 2006

        Args:
            epoch (float): Epoch time in days from Jan 0, 1950, 0 hr

        Returns:
            None
        """
        # Initialize output dataclass
        out = SGP4InitOutput()

        # Calculate auxiliary epoch quantities
        out.eccsq = self.satrec.ecco**2
        out.omeosq = 1.0 - out.eccsq
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

    def sgp4init(self, epoch: float, tol: float = const.SMALL):
        """Initializes variables for SGP4.

        References:
            - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
            - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
            - Hoots, Schumacher, and Glover, 2004
            - Vallado, Crawford, Hujsak, Kelso, 2006

        Args:
            epoch (float): Epoch time in days from Jan 0, 1950 0 hr
            tol (float, optional): Tolerance for small values (default = const.SMALL)

        Returns:
            None (updates self.satrec)

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
            self.propagate(0)
            return

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
            # Handle unexpected cases
            # TODO: Check if it really makes sense to raise an error here
            raise ValueError(
                "Initialization skipped: Satellite is neither deep-space nor"
                "non-deep-space. Check input parameters for inconsistencies."
            )

        # Propagate to zero epoch
        self.propagate(0)

    def propagate(
        self, tsince: float, n_iter: int = 10, tol: float = const.SMALL
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Simplified General Perturbations 4 (SGP4) model.

        Args:
            tsince (float): Time since epoch in minutes.
            n_iter (int, optional): Number of iterations for solving Kepler's equation
                                    (default = 10)
            tol (float, optional): Tolerance for small values (default = const.SMALL)

        Returns:
            tuple: (r, v)
                r (np.ndarray): ECI position vector in km.
                v (np.ndarray): ECI Velocity vector in km/s.
        """
        # Initialize position and velocity vectors
        r, v = np.zeros(3), np.zeros(3)

        # Compute vkmpersec
        vkmpersec = self.grav_const.radiusearthkm * self.grav_const.xke / 60

        # Clear sgp4 error flag
        self.satrec.t = tsince
        self.satrec.error = 0

        # Update for secular gravity and atmospheric drag
        xmdf = self.satrec.mo + self.satrec.mdot * self.satrec.t
        argpdf = self.satrec.argpo + self.satrec.argpdot * self.satrec.t
        nodedf = self.satrec.nodeo + self.satrec.nodedot * self.satrec.t
        argpm, mm = argpdf, xmdf
        t2 = self.satrec.t**2
        nodem = nodedf + self.satrec.nodecf * t2
        tempa = 1 - self.satrec.cc1 * self.satrec.t
        tempe = self.satrec.bstar * self.satrec.cc4 * self.satrec.t
        templ = self.satrec.t2cof * t2

        if not self.satrec.isimp:
            delomg = self.satrec.omgcof * self.satrec.t
            delm = self.satrec.xmcof * (
                (1 + self.satrec.eta * np.cos(xmdf)) ** 3 - self.satrec.delmo
            )
            temp = delomg + delm
            mm = xmdf + temp
            argpm = argpdf - temp
            t3 = t2 * self.satrec.t
            t4 = t3 * self.satrec.t
            tempa = (
                tempa - self.satrec.d2 * t2 - self.satrec.d3 * t3 - self.satrec.d4 * t4
            )
            tempe = tempe + self.satrec.bstar * self.satrec.cc5 * (
                np.sin(mm) - self.satrec.sinmao
            )
            templ = (
                templ
                + self.satrec.t3cof * t3
                + t4 * (self.satrec.t4cof + self.satrec.t * self.satrec.t5cof)
            )

        nm, em, inclm = self.satrec.no, self.satrec.ecco, self.satrec.inclo
        if self.use_deep_space:
            self.ds.dsinit_out.em = em
            self.ds.dsinit_out.inclm = inclm
            self.ds.dsinit_out.nm = nm
            self.ds.dspace(self.satrec, self.satrec.t, self.sgp4init_out.gsto)
            out = self.ds.dsinit_out
            em, inclm, nodem, argpm, nm, mm = (
                out.em,
                out.inclm,
                out.nodem,
                out.argpm,
                out.nm,
                out.mm,
            )

        if nm <= 0:
            # Return early with error
            self.satrec.error = 2
            return np.zeros(3), np.zeros(3)

        am = (self.grav_const.xke / nm) ** self.x2o3 * tempa * tempa
        nm = self.grav_const.xke / am**1.5
        em -= tempe

        if (em >= 1) or (em < -0.001) or (am < 0.95):
            # Return early with error
            self.satrec.error = 1
            return r, v

        em = max(em, 1e-6)
        mm = mm + self.satrec.no * templ
        xlm = mm + argpm + nodem
        nodem = np.remainder(nodem, const.TWOPI)
        argpm = np.remainder(argpm, const.TWOPI)
        xlm = np.remainder(xlm, const.TWOPI)
        mm = np.remainder(xlm - argpm - nodem, const.TWOPI)

        # Compute extra mean quantities
        sinim, cosim = np.sin(inclm), np.cos(inclm)

        # Add lunar-solar periodics
        ep, xincp, argpp, nodep, mp = em, inclm, argpm, nodem, mm
        sinip, cosip = sinim, cosim

        if self.use_deep_space:
            # Add lunar-solar periodics
            self.ds.ep = ep
            self.ds.inclp = xincp
            self.ds.nodep = nodep
            self.ds.argpp = argpp
            self.ds.mp = mp

            self.ds.dpper(self.satrec.t)

            ep, xincp, nodep, argpp, mp = (
                self.ds.ep,
                self.ds.inclp,
                self.ds.nodep,
                self.ds.argpp,
                self.ds.mp,
            )

            if xincp < 0:
                xincp = -xincp
                nodep += np.pi
                argpp -= np.pi

            if (ep < 0) or (ep > 1):
                # Return early with error
                self.satrec.error = 3
                return r, v

            # Long period periodics
            sinip, cosip = np.sin(xincp), np.cos(xincp)
            self.satrec.aycof = -0.5 * self.grav_const.j3oj2 * sinip
            den = 1 + cosip if abs(cosip + 1) > tol else const.SMALL
            self.satrec.xlcof = (
                -0.25 * self.grav_const.j3oj2 * sinip * (3 + 5 * cosip) / den
            )

        axnl = ep * np.cos(argpp)
        temp = 1 / (am * (1 - ep**2))
        aynl = ep * np.sin(argpp) + temp * self.satrec.aycof
        xl = mp + argpp + nodep + temp * self.satrec.xlcof * axnl

        # Solve Kepler's equation
        u = np.remainder(xl - nodep, const.TWOPI)
        eo1, tem5, ktr = u, np.inf, 1
        sineo1, coseo1 = np.sin(eo1), np.cos(eo1)

        while (abs(tem5) >= tol) and (ktr <= n_iter):
            # Upate sine and cosine values for eo1
            # TODO: this should be done at the end of the loop instead?
            sineo1, coseo1 = np.sin(eo1), np.cos(eo1)

            # Compute correction
            tem5 = 1 - coseo1 * axnl - sineo1 * aynl
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5

            # Limit correction to avoid divergence
            lim = 0.95
            if abs(tem5) >= lim:
                tem5 = lim if tem5 > 0 else -lim

            eo1 += tem5
            ktr += 1

        # Short period preliminary quantities
        ecose = axnl * coseo1 + aynl * sineo1
        esine = axnl * sineo1 - aynl * coseo1
        el2 = axnl * axnl + aynl * aynl
        pl = am * (1 - el2)

        if pl < 0:
            # Return early with error
            self.satrec.error = 4
            return r, v

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

        # Update for short period periodics
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

        # Orientation vectors
        sinsu, cossu = np.sin(su), np.cos(su)
        snod, cnod = np.sin(xnode), np.cos(xnode)
        sini, cosi = np.sin(xinc), np.cos(xinc)
        xmx, xmy = -snod * cosi, cnod * cosi
        ux = xmx * sinsu + cnod * cossu
        uy = xmy * sinsu + snod * cossu
        uz = sini * sinsu
        vx = xmx * cossu - cnod * sinsu
        vy = xmy * cossu - snod * sinsu
        vz = sini * cossu

        # Position and velocity vectors
        r = np.array([ux, uy, uz]) * mrt * self.grav_const.radiusearthkm
        v = (np.array([ux, uy, uz]) * mvt + np.array([vx, vy, vz]) * rvdot) * vkmpersec

        # Check for decay condition
        if mrt < 1:
            self.satrec.error = 6

        return r, v
