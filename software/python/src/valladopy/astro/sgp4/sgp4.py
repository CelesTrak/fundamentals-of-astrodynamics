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

    def initl(self, epoch: float) -> SGP4InitOutput:
        """Initialize parameters for the SPG4 propagator.

        References:
            - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
            - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
            - Hoots, Schumacher, and Glover, 2004
            - Vallado, Crawford, Hujsak, Kelso, 2006

        Args:
            epoch (float): Epoch time in days from Jan 0, 1950, 0 hr

        Returns:
            SGP4InitOutput: Dataclass encapsulating the initialized values
        """
        # Initialize output dataclass
        sgp4init_output = SGP4InitOutput()

        # Calculate auxiliary epoch quantities
        sgp4init_output.eccsq = self.satrec.ecco**2
        sgp4init_output.omeosq = 1.0 - sgp4init_output.eccsq
        sgp4init_output.rteosq = np.sqrt(sgp4init_output.omeosq)
        sgp4init_output.cosio = np.cos(self.satrec.inclo)
        sgp4init_output.cosio2 = sgp4init_output.cosio**2

        # Un-Kozai the mean motion
        ak = (self.grav_const.xke / self.satrec.no_kozai) ** self.x2o3
        d1 = (
            0.75
            * self.grav_const.j2
            * (3 * sgp4init_output.cosio2 - 1)
            / (sgp4init_output.rteosq * sgp4init_output.omeosq)
        )
        delta = d1 / (ak**2)
        adel = ak * (1 - delta**2 - delta * (1 / 3 + 134 * delta**2 / 81))
        delta = d1 / (adel**2)
        sgp4init_output.no_unkozai = self.satrec.no_kozai / (1 + delta)

        # Calculate other terms
        sgp4init_output.ao = (
            self.grav_const.xke / sgp4init_output.no_unkozai
        ) ** self.x2o3
        sgp4init_output.sinio = np.sin(self.satrec.inclo)
        po = sgp4init_output.ao * sgp4init_output.omeosq
        sgp4init_output.con42 = 1 - 5 * sgp4init_output.cosio2
        sgp4init_output.con41 = (
            -sgp4init_output.con42 - sgp4init_output.cosio2 - sgp4init_output.cosio2
        )
        sgp4init_output.ainv = 1 / sgp4init_output.ao
        sgp4init_output.posq = po**2
        sgp4init_output.rp = sgp4init_output.ao * (1 - self.satrec.ecco)

        # Calculate Greenwich Sidereal Time
        if self.use_afspc_mode:
            sgp4init_output.gsto = gstime(epoch + JD_EPOCH_1950) % const.TWOPI
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

    def sgp4init(self, epoch: float, tol: float = const.SMALL):
        """Initializes variables for SGP4.

        Args:
            epoch (float): Epoch time in days from Jan 0, 1950 0 hr.
            tol (float, optional): Tolerance for small values (default = 1e-10)
        """
        # Earth constants
        ss = 78 / self.grav_const.radiusearthkm + 1
        qzms2t = ((120 - 78) / self.grav_const.radiusearthkm) ** 4

        # Initialize epoch-dependent values
        self.satrec.t = 0
        self.sgp4init_out = self.initl(epoch)
        self.satrec.no = self.sgp4init_out.no_unkozai

        # Calculate derived orbital parameters
        self.satrec.a = (self.satrec.no * self.grav_const.tumin) ** (-2 / 3)
        self.satrec.alta = self.satrec.a * (1 + self.satrec.ecco) - 1
        self.satrec.altp = self.satrec.a * (1 - self.satrec.ecco) - 1

        if self.sgp4init_out.omeosq >= 0 and self.satrec.no >= 0:
            self.satrec.isimp = 0
            if self.sgp4init_out.rp < (220 / self.grav_const.radiusearthkm + 1):
                self.satrec.isimp = 1
            sfour = ss
            qzms24 = qzms2t
            perigee = (self.sgp4init_out.rp - 1) * self.grav_const.radiusearthkm

            # Adjust sfour and qzms24 for perigees below 156 km
            if perigee < 156:
                sfour = 20 if perigee < 98 else perigee - 78
                qzms24 = ((120 - sfour) / self.grav_const.radiusearthkm) ** 4
                sfour /= self.grav_const.radiusearthkm + 1

            pinvsq = 1 / self.sgp4init_out.posq
            tsi = 1 / (self.sgp4init_out.ao - sfour)
            self.satrec.eta = self.sgp4init_out.ao * self.satrec.ecco * tsi
            etasq = self.satrec.eta**2
            eeta = self.satrec.ecco * self.satrec.eta
            psisq = abs(1 - etasq)
            coef = qzms24 * tsi**4
            coef1 = coef / psisq**3.5

            # Compute CC1 and related parameters
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

            # Update other coefficients
            self.satrec.omgcof = self.satrec.bstar * cc3 * np.cos(self.satrec.argpo)
            self.satrec.xmcof = 0
            if self.satrec.ecco > 1e-4:
                self.satrec.xmcof = -self.x2o3 * coef * self.satrec.bstar / eeta
            self.satrec.nodecf = (
                3.5 * self.sgp4init_out.omeosq * xhdot1 * self.satrec.cc1
            )
            self.satrec.t2cof = 1.5 * self.satrec.cc1

            # Handle divide-by-zero for xinco = 180 degrees
            if abs(self.sgp4init_out.cosio + 1) > tol:
                den = 1.0 + self.sgp4init_out.cosio
            else:
                den = tol
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

            # Deep space initialization
            if (const.TWOPI / self.satrec.no) >= 225:

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
                self.satrec.isimp = 1
                tc = argpm = nodem = mm = 0
                inclm = self.satrec.inclo

                # Call dscom function to compute deep-space common variables
                self.ds.dscom(tc)

                # Call dpper function to adjust for perturbations
                if not self.satrec.init:
                    self.ds.dpper(self.satrec.t)
                    self.satrec.ecco = self.ds.ep
                    self.satrec.inclo = self.ds.inclp
                    self.satrec.nodeo = self.ds.nodep
                    self.satrec.argpo = self.ds.argpp
                    self.satrec.mo = self.ds.mp

                # Call dsinit function for further deep-space initialization
                self.ds.dsinit(
                    self.grav_const.xke,
                    self.satrec.argpo,
                    self.satrec.t,
                    tc,
                    self.sgp4init_out.gsto,
                    self.satrec.mo,
                    self.satrec.mdot,
                    self.satrec.no,
                    self.satrec.nodeo,
                    self.satrec.nodedot,
                    xpidot,
                    self.satrec.ecco,
                    self.sgp4init_out.eccsq,
                    inclm,
                    nodem,
                    argpm,
                    mm,
                )

            # Non-deep space initialization
            else:
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
                    3 * self.satrec.d3
                    + self.satrec.cc1 * (12 * self.satrec.d2 + 10 * cc1sq)
                )
                self.satrec.t5cof = 0.2 * (
                    3 * self.satrec.d4
                    + 12 * self.satrec.cc1 * self.satrec.d3
                    + 6 * self.satrec.d2**2
                    + 15 * cc1sq * (2 * self.satrec.d2 + cc1sq)
                )

        # Propagate to zero epoch
        self.propagate(0)

    def propagate(self, t: float):
        """
        Perform the propagation for the satellite at time t.
        """
        pass
