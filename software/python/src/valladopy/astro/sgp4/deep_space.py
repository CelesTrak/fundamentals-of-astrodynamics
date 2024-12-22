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
from .utils import SatRec


@dataclass
class DscomOutput:
    # TODO: Do we need to return all these variables?
    sinim: float = 0.0
    cosim: float = 0.0
    sinomm: float = 0.0
    cosomm: float = 0.0
    snodm: float = 0.0
    cnodm: float = 0.0
    day: float = 0.0
    e3: float = 0.0
    ee2: float = 0.0
    em: float = 0.0
    emsq: float = 0.0
    gam: float = 0.0
    peo: float = 0.0
    pgho: float = 0.0
    pho: float = 0.0
    pinco: float = 0.0
    plo: float = 0.0
    rtemsq: float = 0.0
    se2: float = 0.0
    se3: float = 0.0
    sgh2: float = 0.0
    sgh3: float = 0.0
    sgh4: float = 0.0
    sh2: float = 0.0
    sh3: float = 0.0
    si2: float = 0.0
    si3: float = 0.0
    sl2: float = 0.0
    sl3: float = 0.0
    sl4: float = 0.0
    s1: float = 0.0
    s2: float = 0.0
    s3: float = 0.0
    s4: float = 0.0
    s5: float = 0.0
    s6: float = 0.0
    s7: float = 0.0
    ss1: float = 0.0
    ss2: float = 0.0
    ss3: float = 0.0
    ss4: float = 0.0
    ss5: float = 0.0
    ss6: float = 0.0
    ss7: float = 0.0
    sz1: float = 0.0
    sz2: float = 0.0
    sz3: float = 0.0
    sz11: float = 0.0
    sz12: float = 0.0
    sz13: float = 0.0
    sz21: float = 0.0
    sz22: float = 0.0
    sz23: float = 0.0
    sz31: float = 0.0
    sz32: float = 0.0
    sz33: float = 0.0
    xgh2: float = 0.0
    xgh3: float = 0.0
    xgh4: float = 0.0
    xh2: float = 0.0
    xh3: float = 0.0
    xi2: float = 0.0
    xi3: float = 0.0
    xl2: float = 0.0
    xl3: float = 0.0
    xl4: float = 0.0
    nm: float = 0.0
    z1: float = 0.0
    z2: float = 0.0
    z3: float = 0.0
    z11: float = 0.0
    z12: float = 0.0
    z13: float = 0.0
    z21: float = 0.0
    z22: float = 0.0
    z23: float = 0.0
    z31: float = 0.0
    z32: float = 0.0
    z33: float = 0.0
    zmol: float = 0.0
    zmos: float = 0.0


@dataclass
class DsInitOutput:
    em: float = 0.0
    argpm: float = 0.0
    inclm: float = 0.0
    mm: float = 0.0
    nm: float = 0.0
    nodem: float = 0.0
    irez: int = 0
    atime: float = 0.0
    d2201: float = 0.0
    d2211: float = 0.0
    d3210: float = 0.0
    d3222: float = 0.0
    d4410: float = 0.0
    d4422: float = 0.0
    d5220: float = 0.0
    d5232: float = 0.0
    d5421: float = 0.0
    d5433: float = 0.0
    dedt: float = 0.0
    didt: float = 0.0
    dmdt: float = 0.0
    dndt: float = 0.0
    dnodt: float = 0.0
    domdt: float = 0.0
    del1: float = 0.0
    del2: float = 0.0
    del3: float = 0.0
    xfact: float = 0.0
    xlamo: float = 0.0
    xli: float = 0.0
    xni: float = 0.0


class DeepSpace:
    def __init__(
        self,
        epoch: float,
        ep: float,
        inclp: float,
        nodep: float,
        argpp: float,
        np_: float,
        mp: float = 0,
        use_afspc_mode: bool = True,
    ):
        """Initializes the DeepSpace class.

        Args:
            epoch (float): Epoch time (units = ?)  # TODO: check units
            ep (float): Eccentricity
            inclp (float): Inclination in radians
            nodep (float): RAAN (right ascension of ascending node) in radians
            argpp (float): Argument of perigee in radians
            np_ (float): Mean motion in rad/s
            mp (float): Mean anomaly in radians (default = 0)
            use_afspc_mode (bool): Flag to use AFSPC mode (default = True)
        """
        self.epoch = epoch
        self.ep = ep
        self.inclp = inclp
        self.nodep = nodep
        self.argpp = argpp
        self.np_ = np_
        self.mp = mp
        self.use_afspc_mode = use_afspc_mode

        # Initialize output dataclasses
        self.dscom_out = None
        self.dsinit_out = None

    def dscom(self, tc: float):
        """Computes deep space common terms for SGP4 (used by both the secular and
        periodics subroutines).

        References:
            - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
            - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
            - Hoots, Schumacher, and Glover, 2004
            - Vallado, Crawford, Hujsak, Kelso, 2006

        Args:
            tc (float): Time correction in minutes

        Returns:
            None (sets dscom_out attribute)
        """
        # Initialize output dataclass
        out = DscomOutput()

        # Constants
        zes = 0.01675
        zel = 0.05490
        c1ss = 2.9864797e-6
        c1l = 4.7968065e-7
        zsinis, zcosis = 0.39785416, 0.91744867
        zcosgs, zsings = 0.1945905, -0.98088458

        # Initialize variables
        out.nm, out.em = self.np_, self.ep
        out.snodm, out.cnodm = np.sin(self.nodep), np.cos(self.nodep)
        out.sinomm, out.cosomm = np.sin(self.argpp), np.cos(self.argpp)
        out.sinim, out.cosim = np.sin(self.inclp), np.cos(self.inclp)
        out.emsq = out.em**2
        betasq = 1.0 - out.emsq
        out.rtemsq = np.sqrt(betasq)

        # Initialize lunar solar terms
        out.day = self.epoch + 18261.5 + tc / const.DAY2MIN
        xnodce = np.remainder(4.5236020 - 9.2422029e-4 * out.day, const.TWOPI)
        stem, ctem = np.sin(xnodce), np.cos(xnodce)
        zcosil = 0.91375164 - 0.03568096 * ctem
        zsinil = np.sqrt(1.0 - zcosil**2)
        zsinhl = 0.089683511 * stem / zsinil
        zcoshl = np.sqrt(1.0 - zsinhl**2)
        out.gam = 5.8351514 + 0.0019443680 * out.day
        zy = zcoshl * ctem + 0.91744867 * zsinhl * stem
        zx = np.arctan2(0.39785416 * stem / zsinil, zy)
        zx = out.gam + zx - xnodce
        zcosgl, zsingl = np.cos(zx), np.sin(zx)

        # Solar and lunar terms
        zcosg, zsing = zcosgs, zsings
        zcosi, zsini = zcosis, zsinis
        zcosh, zsinh = out.cnodm, out.snodm
        cc = c1ss

        # Loop for solar and lunar terms
        for lsflg in range(2):
            # Compute intermediate values
            a1 = zcosg * zcosh + zsing * zcosi * zsinh
            a3 = -zsing * zcosh + zcosg * zcosi * zsinh
            a7 = -zcosg * zsinh + zsing * zcosi * zcosh
            a8 = zsing * zsini
            a9 = zsing * zsinh + zcosg * zcosi * zcosh
            a10 = zcosg * zsini
            a2 = out.cosim * a7 + out.sinim * a8
            a4 = out.cosim * a9 + out.sinim * a10
            a5 = -out.sinim * a7 + out.cosim * a8
            a6 = -out.sinim * a9 + out.cosim * a10

            # Further calculations
            x1 = a1 * out.cosomm + a2 * out.sinomm
            x2 = a3 * out.cosomm + a4 * out.sinomm
            x3 = -a1 * out.sinomm + a2 * out.cosomm
            x4 = -a3 * out.sinomm + a4 * out.cosomm
            x5 = a5 * out.sinomm
            x6 = a6 * out.sinomm
            x7 = a5 * out.cosomm
            x8 = a6 * out.cosomm

            # Z matrix calculations
            out.z31 = 12 * x1**2 - 3 * x3**2
            out.z32 = 24 * x1 * x2 - 6 * x3 * x4
            out.z33 = 12 * x2**2 - 3 * x4**2
            out.z1 = 3 * (a1**2 + a2**2) + out.z31 * out.emsq
            out.z2 = 6 * (a1 * a3 + a2 * a4) + out.z32 * out.emsq
            out.z3 = 3 * (a3**2 + a4**2) + out.z33 * out.emsq
            out.z11 = -6 * a1 * a5 + out.emsq * (-24 * x1 * x7 - 6 * x3 * x5)
            out.z12 = -6 * (a1 * a6 + a3 * a5) + out.emsq * (
                -24 * (x2 * x7 + x1 * x8) - 6 * (x3 * x6 + x4 * x5)
            )
            out.z13 = -6 * a3 * a6 + out.emsq * (-24 * x2 * x8 - 6 * x4 * x6)
            out.z21 = 6 * a2 * a5 + out.emsq * (24 * x1 * x5 - 6 * x3 * x7)
            out.z22 = 6 * (a4 * a5 + a2 * a6) + out.emsq * (
                24 * (x2 * x5 + x1 * x6) - 6 * (x4 * x7 + x3 * x8)
            )
            out.z23 = 6 * a4 * a6 + out.emsq * (24 * x2 * x6 - 6 * x4 * x8)

            # Update z values with beta square terms
            out.z1 += out.z1 + betasq * out.z31
            out.z2 += out.z2 + betasq * out.z32
            out.z3 += out.z3 + betasq * out.z33

            # S coefficients
            out.s3 = cc * 1 / out.nm
            out.s2 = -0.5 * out.s3 / out.rtemsq
            out.s4 = out.s3 * out.rtemsq
            out.s1 = -15 * out.em * out.s4
            out.s5 = x1 * x3 + x2 * x4
            out.s6 = x2 * x3 + x1 * x4
            out.s7 = x2 * x4 - x1 * x3

            # Update lunar terms on first iteration
            if lsflg == 0:
                out.ss1 = out.s1
                out.ss2 = out.s2
                out.ss3 = out.s3
                out.ss4 = out.s4
                out.ss5 = out.s5
                out.ss6 = out.s6
                out.ss7 = out.s7
                out.sz1 = out.z1
                out.sz2 = out.z2
                out.sz3 = out.z3
                out.sz11 = out.z11
                out.sz12 = out.z12
                out.sz13 = out.z13
                out.sz21 = out.z21
                out.sz22 = out.z22
                out.sz23 = out.z23
                out.sz31 = out.z31
                out.sz32 = out.z32
                out.sz33 = out.z33
                zcosg, zsing = zcosgl, zsingl
                zcosi, zsini = zcosil, zsinil
                zcosh = zcoshl * out.cnodm + zsinhl * out.snodm
                zsinh = out.snodm * zcoshl - out.cnodm * zsinhl
                cc = c1l

        # Compute additional terms
        out.zmol = np.remainder(4.7199672 + 0.22997150 * out.day - out.gam, const.TWOPI)
        out.zmos = np.remainder(6.2565837 + 0.017201977 * out.day, const.TWOPI)

        # Compute solar terms
        out.se2 = 2 * out.ss1 * out.ss6
        out.se3 = 2 * out.ss1 * out.ss7
        out.si2 = 2 * out.ss2 * out.sz12
        out.si3 = 2 * out.ss2 * (out.sz13 - out.sz11)
        out.sl2 = -2 * out.ss3 * out.sz2
        out.sl3 = -2 * out.ss3 * (out.sz3 - out.sz1)
        out.sl4 = -2 * out.ss3 * (-21 - 9 * out.emsq) * zes
        out.sgh2 = 2 * out.ss4 * out.sz32
        out.sgh3 = 2 * out.ss4 * (out.sz33 - out.sz31)
        out.sgh4 = -18 * out.ss4 * zes
        out.sh2 = -2 * out.ss2 * out.sz22
        out.sh3 = -2 * out.ss2 * (out.sz23 - out.sz21)

        # Compute lunar terms
        out.ee2 = 2 * out.s1 * out.s6
        out.e3 = 2 * out.s1 * out.s7
        out.xi2 = 2 * out.s2 * out.z12
        out.xi3 = 2 * out.s2 * (out.z13 - out.z11)
        out.xl2 = -2 * out.s3 * out.z2
        out.xl3 = -2 * out.s3 * (out.z3 - out.z1)
        out.xl4 = -2 * out.s3 * (-21 - 9 * out.emsq) * zel
        out.xgh2 = 2 * out.s4 * out.z32
        out.xgh3 = 2 * out.s4 * (out.z33 - out.z31)
        out.xgh4 = -18 * out.s4 * zel
        out.xh2 = -2 * out.s2 * out.z22
        out.xh3 = -2 * out.s2 * (out.z23 - out.z21)

        self.dscom_out = out

    def dpper(self, t: float, incl_tol: float = 0.2):
        """Deep space long period periodic contributions to mean elements.

        References:
            - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
            - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
            - Hoots, Schumacher, and Glover, 2004
            - Vallado, Crawford, Hujsak, Kelso, 2006

        Args:
            t (float): Elapsed time since epoch in minutes
            incl_tol (float): Inclination tolerance for periodics (default = 0.2)

        Returns:
            None (updates class attributes)

        Notes:
            - These periodics are zero at epoch by design.
        """
        # Check if dscom_out is set
        if not self.dscom_out:
            raise ValueError("dscom_out not set. Run dscom() first.")

        # Constants
        zns = 1.19459e-5
        zes = 0.01675
        znl = 1.5835218e-4
        zel = 0.05490

        # Calculate time-varying periodics
        zm = self.dscom_out.zmos + zns * t
        zf = zm + 2.0 * zes * np.sin(zm)
        sinzf = np.sin(zf)
        f2 = 0.5 * sinzf**2 - 0.25
        f3 = -0.5 * sinzf * np.cos(zf)
        ses = self.dscom_out.se2 * f2 + self.dscom_out.se3 * f3
        sis = self.dscom_out.si2 * f2 + self.dscom_out.si3 * f3
        sls = (
            self.dscom_out.sl2 * f2
            + self.dscom_out.sl3 * f3
            + self.dscom_out.sl4 * sinzf
        )
        sghs = (
            self.dscom_out.sgh2 * f2
            + self.dscom_out.sgh3 * f3
            + self.dscom_out.sgh4 * sinzf
        )
        shs = self.dscom_out.sh2 * f2 + self.dscom_out.sh3 * f3
        zm = self.dscom_out.zmol + znl * t
        zf = zm + 2.0 * zel * np.sin(zm)
        sinzf = np.sin(zf)
        f2 = 0.5 * sinzf**2 - 0.25
        f3 = -0.5 * sinzf * np.cos(zf)
        sel = self.dscom_out.ee2 * f2 + self.dscom_out.e3 * f3
        sil = self.dscom_out.xi2 * f2 + self.dscom_out.xi3 * f3
        sll = (
            self.dscom_out.xl2 * f2
            + self.dscom_out.xl3 * f3
            + self.dscom_out.xl4 * sinzf
        )
        sghl = (
            self.dscom_out.xgh2 * f2
            + self.dscom_out.xgh3 * f3
            + self.dscom_out.xgh4 * sinzf
        )
        shll = self.dscom_out.xh2 * f2 + self.dscom_out.xh3 * f3

        # Initialize periodics
        pe = ses + sel
        pinc = sis + sil
        pl = sls + sll
        pgh = sghs + sghl
        ph = shs + shll

        # Subtract baseline offsets
        pe -= self.dscom_out.peo
        pinc -= self.dscom_out.pinco
        pl -= self.dscom_out.plo
        pgh -= self.dscom_out.pgho
        ph -= self.dscom_out.pho

        # Update inclination and eccentricity
        self.inclp += pinc
        self.ep += pe
        sinip = np.sin(self.inclp)
        cosip = np.cos(self.inclp)

        # Apply periodics based on inclination
        if self.inclp >= incl_tol:
            # Apply periodics directly
            ph /= sinip
            pgh -= cosip * ph
            self.argpp += pgh
            self.nodep += ph
            self.mp += pl
        else:
            # Apply periodics with Lyddane modification
            sinop, cosop = np.sin(self.nodep), np.cos(self.nodep)
            alfdp = sinip * sinop
            betdp = sinip * cosop
            dalf = ph * cosop + pinc * cosip * sinop
            dbet = -ph * sinop + pinc * cosip * cosop
            alfdp += dalf
            betdp += dbet

            # SGP4 fix for AFSPC-written intrinsic functions
            self.nodep = np.remainder(self.nodep, const.TWOPI)
            if self.nodep < 0.0 and self.use_afspc_mode:
                self.nodep += const.TWOPI

            xls = self.mp + self.argpp + cosip * self.nodep
            dls = pl + pgh - pinc * self.nodep * sinip
            xls += dls
            xnoh = self.nodep
            self.nodep = np.arctan2(alfdp, betdp)

            # SGP4 fix for AFSPC-written intrinsic functions
            if self.nodep < 0.0 and self.use_afspc_mode:
                self.nodep += const.TWOPI

            if np.abs(xnoh - self.nodep) > np.pi:
                self.nodep += const.TWOPI if self.nodep < xnoh else -const.TWOPI

            self.mp += pl
            self.argpp = xls - self.mp - cosip * self.nodep

    def dsinit(
        self,
        satrec: SatRec,
        xke: float,
        tc: float,
        gsto: float,
        xpidot: float,
        eccsq: float,
        inclm: float,
        nodem: float,
        argpm: float,
        mm: float,
        incl_tol: float = np.radians(3),
    ):
        """Deep space initialization for SGP4.

        This function provides deep space contributions to mean motion dot due to
        geopotential resonance with half day and one day orbits.

        References:
            - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
            - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
            - Hoots, Schumacher, and Glover, 2004
            - Vallado, Crawford, Hujsak, Kelso, 2006

        Args:
            satrec (SatRec): Satellite record dataclass
            xke (float): SGP4 constant
            tc (float): Time correction in minutes
            gsto (float): GST at epoch in radians
            xpidot (float): RAAN dot in rad/s
            eccsq (float): Eccentricity squared
            inclm (float): Inclination in radians
            nodem (float): RAAN in radians
            argpm (float): Argument of perigee in radians
            mm (float): Mean anomaly in radians
            incl_tol (float): Inclination tolerance for applying periodics
                              (default = 3 deg in radians)

        Returns:
            None (sets dsinit_out attribute)
        """
        # Initialize output dataclass
        out = DsInitOutput()

        # Constants
        rptim = 4.37526908801129966e-3
        q22, q31, q33 = 1.7891679e-6, 2.1460748e-6, 2.2123015e-7
        root22, root44, root54 = 1.7891679e-6, 7.3636953e-9, 2.1765803e-9
        root32, root52 = 3.7393792e-7, 1.1428639e-7
        x2o3 = 2 / 3
        znl, zns = 1.5835218e-4, 1.19459e-5

        # Initialize outputs
        out.em, out.nm = self.dscom_out.em, self.dscom_out.nm
        if 0.0034906585 < out.nm < 0.0052359877:
            out.irez = 1
        if 8.26e-3 <= out.nm <= 9.24e-3 and out.em >= 0.5:
            out.irez = 2
        emsq = self.dscom_out.emsq

        # Solar terms
        ses = self.dscom_out.ss1 * zns * self.dscom_out.ss5
        sis = self.dscom_out.ss2 * zns * (self.dscom_out.sz11 + self.dscom_out.sz13)
        sls = (
            -zns
            * self.dscom_out.ss3
            * (self.dscom_out.sz1 + self.dscom_out.sz3 - 14 - 6 * self.dscom_out.emsq)
        )
        sghs = (
            self.dscom_out.ss4 * zns * (self.dscom_out.sz31 + self.dscom_out.sz33 - 6)
        )
        shs = -zns * self.dscom_out.ss2 * (self.dscom_out.sz21 + self.dscom_out.sz23)
        if inclm < incl_tol or inclm > np.pi - incl_tol:
            shs = 0
        if self.dscom_out.sinim != 0:
            shs /= self.dscom_out.sinim
        sgs = sghs - self.dscom_out.cosim * shs

        # Lunar terms
        out.dedt = ses + self.dscom_out.s1 * znl * self.dscom_out.s5
        out.didt = sis + self.dscom_out.s2 * znl * (
            self.dscom_out.z11 + self.dscom_out.z13
        )
        out.dmdt = sls - znl * self.dscom_out.s3 * (
            self.dscom_out.z1 + self.dscom_out.z3 - 14 - 6 * self.dscom_out.emsq
        )
        sghl = self.dscom_out.s4 * znl * (self.dscom_out.z31 + self.dscom_out.z33 - 6)
        shll = -znl * self.dscom_out.s2 * (self.dscom_out.z21 + self.dscom_out.z23)
        if inclm < incl_tol or inclm > np.pi - incl_tol:
            shll = 0
        out.domdt, out.dnodt = sgs + sghl, shs
        if self.dscom_out.sinim != 0:
            out.domdt -= self.dscom_out.cosim / self.dscom_out.sinim * shll
            out.dnodt += shll / self.dscom_out.sinim

        # Deep space resonance effects
        theta = np.remainder(gsto + tc * rptim, const.TWOPI)
        out.em += out.dedt * satrec.t
        out.inclm = inclm + out.didt * satrec.t
        out.argpm = argpm + out.domdt * satrec.t
        out.nodem = nodem + out.dnodt * satrec.t
        out.mm = mm + out.dmdt * satrec.t

        # Return if no resonance
        if out.irez == 0:
            return out

        aonv = (out.nm / xke) ** x2o3

        # Geopotential resonance for 12-hour orbits
        if out.irez == 2:
            cosisq = self.dscom_out.cosim**2
            emo = out.em
            em = satrec.ecco
            emsqo = self.dscom_out.emsq
            emsq = eccsq
            eoc = em * emsq
            g201 = -0.306 - (em - 0.64) * 0.440

            if em <= 0.65:
                g211 = 3.616 - 13.2470 * em + 16.2900 * emsq
                g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc
                g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc
                g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc
                g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc
                g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc
            else:
                g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc
                g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc
                g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc
                g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc
                g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc
                if em > 0.715:
                    g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc
                else:
                    g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq

            if em < 0.7:
                g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21 * eoc
                g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc
                g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4 * eoc
            else:
                g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc
                g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc
                g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc

            out.em = em
            sini2 = self.dscom_out.sinim**2
            f220 = 0.75 * (1.0 + 2.0 * self.dscom_out.cosim + cosisq)
            f221 = 1.5 * sini2
            f321 = (
                1.875
                * self.dscom_out.sinim
                * (1 - 2 * self.dscom_out.cosim - 3 * cosisq)
            )
            f322 = (
                -1.875
                * self.dscom_out.sinim
                * (1 + 2 * self.dscom_out.cosim - 3 * cosisq)
            )
            f441 = 35 * sini2 * f220
            f442 = 39.3750 * sini2**2
            f522 = (
                9.84375
                * self.dscom_out.sinim
                * (
                    sini2 * (1 - 2 * self.dscom_out.cosim - 5 * cosisq)
                    + 0.33333333 * (-2 + 4 * self.dscom_out.cosim + 6 * cosisq)
                )
            )
            f523 = self.dscom_out.sinim * (
                4.92187512 * sini2 * (-2 - 4 * self.dscom_out.cosim + 10 * cosisq)
                + 6.56250012 * (1 + 2 * self.dscom_out.cosim - 3 * cosisq)
            )
            f542 = (
                29.53125
                * self.dscom_out.sinim
                * (
                    2
                    - 8 * self.dscom_out.cosim
                    + cosisq * (-12 + 8 * self.dscom_out.cosim + 10 * cosisq)
                )
            )
            f543 = (
                29.53125
                * self.dscom_out.sinim
                * (
                    -2
                    - 8 * self.dscom_out.cosim
                    + cosisq * (12 + 8 * self.dscom_out.cosim - 10 * cosisq)
                )
            )

            xno2 = out.nm**2
            ainv2 = aonv**2
            temp1 = 3 * xno2 * ainv2
            temp = temp1 * root22
            out.d2201 = temp * f220 * g201
            out.d2211 = temp * f221 * g211
            temp1 *= aonv
            temp = temp1 * root32
            out.d3210 = temp * f321 * g310
            out.d3222 = temp * f322 * g322
            temp1 *= aonv
            temp = 2 * temp1 * root44
            out.d4410 = temp * f441 * g410
            out.d4422 = temp * f442 * g422
            temp1 *= aonv
            temp = temp1 * root52
            out.d5220 = temp * f522 * g520
            out.d5232 = temp * f523 * g532
            temp = 2 * temp1 * root54
            out.d5421 = temp * f542 * g521
            out.d5433 = temp * f543 * g533
            out.xlamo = (
                satrec.mo + satrec.nodeo + satrec.nodeo - theta - theta
            ) % const.TWOPI
            out.xfact = (
                satrec.mdot
                + out.dmdt
                + 2 * (satrec.nodedot + out.dnodt - rptim)
                - satrec.no
            )
            out.em, emsq = emo, emsqo

        # Synchronous resonance terms
        if out.irez == 1:
            g200 = 1 + emsq * (-2.5 + 0.8125 * emsq)
            g310 = 1 + 2 * emsq
            g300 = 1 + emsq * (-6 + 6.60937 * emsq)
            f220 = 0.75 * (1 + self.dscom_out.cosim) ** 2
            f311 = 0.9375 * self.dscom_out.sinim**2 * (
                1 + 3 * self.dscom_out.cosim
            ) - 0.75 * (1 + self.dscom_out.cosim)
            f330 = 1 + self.dscom_out.cosim
            f330 = 1.875 * f330**3
            out.del1 = 3 * out.nm**2 * aonv**2
            out.del2 = 2 * out.del1 * f220 * g200 * q22
            out.del3 = 3 * out.del1 * f330 * g300 * q33 * aonv
            out.del1 *= f311 * g310 * q31 * aonv
            out.xlamo = (satrec.mo + satrec.nodeo + satrec.argpo - theta) % const.TWOPI
            out.xfact = (
                satrec.mdot
                + xpidot
                - rptim
                + out.dmdt
                + out.domdt
                + out.dnodt
                - satrec.no
            )

        # Initialize the integrator for SGP4
        out.xli, out.xni = out.xlamo, satrec.no
        out.nm += out.dndt

        self.dsinit_out = out
