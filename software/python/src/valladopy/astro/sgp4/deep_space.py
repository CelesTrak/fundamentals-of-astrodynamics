# -----------------------------------------------------------------------------
# Authors: David Vallado, Jeff Beck
# Date: 28 June 2005
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from ... import constants as const


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


def dscom(
    epoch: float,
    tc: float,
    ep: float,
    inclp: float,
    nodep: float,
    argpp: float,
    np_: float,
) -> DscomOutput:
    """Computes deep space common terms for SGP4 ( used by both the secular and
    periodics subroutines).

    References:
        - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
        - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
        - Hoots, Schumacher, and Glover, 2004
        - Vallado, Crawford, Hujsak, Kelso, 2006

    Args:
        epoch (float): Epoch time (units = ?)  # TODO: check units
        tc (float): Time since epoch (units = ?)  # TODO: check units
        ep (float): Eccentricity
        inclp (float): Inclination in radians
        nodep (float): RAAN (right ascension of ascending node) in radians
        argpp (float): Argument of perigee in radians
        np_ (float): Mean motion in rad/s  # TODO: check units

    Returns:
        DscomOutput: Dataclass containing deep space common terms
    """
    # Initialize output dataclass
    dscom_out = DscomOutput()

    # Constants
    zes = 0.01675
    zel = 0.05490
    c1ss = 2.9864797e-6
    c1l = 4.7968065e-7
    zsinis, zcosis = 0.39785416, 0.91744867
    zcosgs, zsings = 0.1945905, -0.98088458

    # Initialize variables
    dscom_out.nm, dscom_out.em = np_, ep
    dscom_out.snodm, dscom_out.cnodm = np.sin(nodep), np.cos(nodep)
    dscom_out.sinomm, dscom_out.cosomm = np.sin(argpp), np.cos(argpp)
    dscom_out.sinim, dscom_out.cosim = np.sin(inclp), np.cos(inclp)
    dscom_out.emsq = dscom_out.em**2
    betasq = 1.0 - dscom_out.emsq
    dscom_out.rtemsq = np.sqrt(betasq)

    # Initialize lunar solar terms
    dscom_out.day = epoch + 18261.5 + tc / const.DAY2MIN
    xnodce = np.remainder(4.5236020 - 9.2422029e-4 * dscom_out.day, const.TWOPI)
    stem, ctem = np.sin(xnodce), np.cos(xnodce)
    zcosil = 0.91375164 - 0.03568096 * ctem
    zsinil = np.sqrt(1.0 - zcosil**2)
    zsinhl = 0.089683511 * stem / zsinil
    zcoshl = np.sqrt(1.0 - zsinhl**2)
    dscom_out.gam = 5.8351514 + 0.0019443680 * dscom_out.day
    zy = zcoshl * ctem + 0.91744867 * zsinhl * stem
    zx = np.arctan2(0.39785416 * stem / zsinil, zy)
    zx = dscom_out.gam + zx - xnodce
    zcosgl, zsingl = np.cos(zx), np.sin(zx)

    # Solar and lunar terms
    zcosg, zsing = zcosgs, zsings
    zcosi, zsini = zcosis, zsinis
    zcosh, zsinh = dscom_out.cnodm, dscom_out.snodm
    cc = c1ss
    xnoi = 1.0 / dscom_out.nm

    # Loop for solar and lunar terms
    for lsflg in range(2):
        # Compute intermediate values
        a1 = zcosg * zcosh + zsing * zcosi * zsinh
        a3 = -zsing * zcosh + zcosg * zcosi * zsinh
        a7 = -zcosg * zsinh + zsing * zcosi * zcosh
        a8 = zsing * zsini
        a9 = zsing * zsinh + zcosg * zcosi * zcosh
        a10 = zcosg * zsini
        a2 = dscom_out.cosim * a7 + dscom_out.sinim * a8
        a4 = dscom_out.cosim * a9 + dscom_out.sinim * a10
        a5 = -dscom_out.sinim * a7 + dscom_out.cosim * a8
        a6 = -dscom_out.sinim * a9 + dscom_out.cosim * a10

        # Further calculations
        x1 = a1 * dscom_out.cosomm + a2 * dscom_out.sinomm
        x2 = a3 * dscom_out.cosomm + a4 * dscom_out.sinomm
        x3 = -a1 * dscom_out.sinomm + a2 * dscom_out.cosomm
        x4 = -a3 * dscom_out.sinomm + a4 * dscom_out.cosomm
        x5 = a5 * dscom_out.sinomm
        x6 = a6 * dscom_out.sinomm
        x7 = a5 * dscom_out.cosomm
        x8 = a6 * dscom_out.cosomm

        # Z matrix calculations
        dscom_out.z31 = 12 * x1**2 - 3 * x3**2
        dscom_out.z32 = 24 * x1 * x2 - 6 * x3 * x4
        dscom_out.z33 = 12 * x2**2 - 3 * x4**2
        dscom_out.z1 = 3 * (a1**2 + a2**2) + dscom_out.z31 * dscom_out.emsq
        dscom_out.z2 = 6 * (a1 * a3 + a2 * a4) + dscom_out.z32 * dscom_out.emsq
        dscom_out.z3 = 3 * (a3**2 + a4**2) + dscom_out.z33 * dscom_out.emsq
        dscom_out.z11 = -6 * a1 * a5 + dscom_out.emsq * (-24 * x1 * x7 - 6 * x3 * x5)
        dscom_out.z12 = -6 * (a1 * a6 + a3 * a5) + dscom_out.emsq * (
            -24 * (x2 * x7 + x1 * x8) - 6 * (x3 * x6 + x4 * x5)
        )
        dscom_out.z13 = -6 * a3 * a6 + dscom_out.emsq * (-24 * x2 * x8 - 6 * x4 * x6)
        dscom_out.z21 = 6 * a2 * a5 + dscom_out.emsq * (24 * x1 * x5 - 6 * x3 * x7)
        dscom_out.z22 = 6 * (a4 * a5 + a2 * a6) + dscom_out.emsq * (
            24 * (x2 * x5 + x1 * x6) - 6 * (x4 * x7 + x3 * x8)
        )
        dscom_out.z23 = 6 * a4 * a6 + dscom_out.emsq * (24 * x2 * x6 - 6 * x4 * x8)

        # Update z values with beta square terms
        dscom_out.z1 += dscom_out.z1 + betasq * dscom_out.z31
        dscom_out.z2 += dscom_out.z2 + betasq * dscom_out.z32
        dscom_out.z3 += dscom_out.z3 + betasq * dscom_out.z33

        # S coefficients
        dscom_out.s3 = cc * xnoi
        dscom_out.s2 = -0.5 * dscom_out.s3 / dscom_out.rtemsq
        dscom_out.s4 = dscom_out.s3 * dscom_out.rtemsq
        dscom_out.s1 = -15 * dscom_out.em * dscom_out.s4
        dscom_out.s5 = x1 * x3 + x2 * x4
        dscom_out.s6 = x2 * x3 + x1 * x4
        dscom_out.s7 = x2 * x4 - x1 * x3

        # Update lunar terms on first iteration
        if lsflg == 0:
            dscom_out.ss1 = dscom_out.s1
            dscom_out.ss2 = dscom_out.s2
            dscom_out.ss3 = dscom_out.s3
            dscom_out.ss4 = dscom_out.s4
            dscom_out.ss5 = dscom_out.s5
            dscom_out.ss6 = dscom_out.s6
            dscom_out.ss7 = dscom_out.s7
            dscom_out.sz1 = dscom_out.z1
            dscom_out.sz2 = dscom_out.z2
            dscom_out.sz3 = dscom_out.z3
            dscom_out.sz11 = dscom_out.z11
            dscom_out.sz12 = dscom_out.z12
            dscom_out.sz13 = dscom_out.z13
            dscom_out.sz21 = dscom_out.z21
            dscom_out.sz22 = dscom_out.z22
            dscom_out.sz23 = dscom_out.z23
            dscom_out.sz31 = dscom_out.z31
            dscom_out.sz32 = dscom_out.z32
            dscom_out.sz33 = dscom_out.z33
            zcosg, zsing = zcosgl, zsingl
            zcosi, zsini = zcosil, zsinil
            zcosh = zcoshl * dscom_out.cnodm + zsinhl * dscom_out.snodm
            zsinh = dscom_out.snodm * zcoshl - dscom_out.cnodm * zsinhl
            cc = c1l

    # Compute additional terms
    dscom_out.zmol = np.remainder(
        4.7199672 + 0.22997150 * dscom_out.day - dscom_out.gam, const.TWOPI
    )
    dscom_out.zmos = np.remainder(6.2565837 + 0.017201977 * dscom_out.day, const.TWOPI)

    # Compute solar terms
    dscom_out.se2 = 2 * dscom_out.ss1 * dscom_out.ss6
    dscom_out.se3 = 2 * dscom_out.ss1 * dscom_out.ss7
    dscom_out.si2 = 2 * dscom_out.ss2 * dscom_out.sz12
    dscom_out.si3 = 2 * dscom_out.ss2 * (dscom_out.sz13 - dscom_out.sz11)
    dscom_out.sl2 = -2 * dscom_out.ss3 * dscom_out.sz2
    dscom_out.sl3 = -2 * dscom_out.ss3 * (dscom_out.sz3 - dscom_out.sz1)
    dscom_out.sl4 = -2 * dscom_out.ss3 * (-21 - 9 * dscom_out.emsq) * zes
    dscom_out.sgh2 = 2 * dscom_out.ss4 * dscom_out.sz32
    dscom_out.sgh3 = 2 * dscom_out.ss4 * (dscom_out.sz33 - dscom_out.sz31)
    dscom_out.sgh4 = -18 * dscom_out.ss4 * zes
    dscom_out.sh2 = -2 * dscom_out.ss2 * dscom_out.sz22
    dscom_out.sh3 = -2 * dscom_out.ss2 * (dscom_out.sz23 - dscom_out.sz21)

    # Compute lunar terms
    dscom_out.ee2 = 2 * dscom_out.s1 * dscom_out.s6
    dscom_out.e3 = 2 * dscom_out.s1 * dscom_out.s7
    dscom_out.xi2 = 2 * dscom_out.s2 * dscom_out.z12
    dscom_out.xi3 = 2 * dscom_out.s2 * (dscom_out.z13 - dscom_out.z11)
    dscom_out.xl2 = -2 * dscom_out.s3 * dscom_out.z2
    dscom_out.xl3 = -2 * dscom_out.s3 * (dscom_out.z3 - dscom_out.z1)
    dscom_out.xl4 = -2 * dscom_out.s3 * (-21 - 9 * dscom_out.emsq) * zel
    dscom_out.xgh2 = 2 * dscom_out.s4 * dscom_out.z32
    dscom_out.xgh3 = 2 * dscom_out.s4 * (dscom_out.z33 - dscom_out.z31)
    dscom_out.xgh4 = -18 * dscom_out.s4 * zel
    dscom_out.xh2 = -2 * dscom_out.s2 * dscom_out.z22
    dscom_out.xh3 = -2 * dscom_out.s2 * (dscom_out.z23 - dscom_out.z21)

    return dscom_out


def dpper(
    dscom_out: DscomOutput,
    t: float,
    ep: float,
    inclp: float,
    nodep: float,
    argpp: float,
    mp: float,
    use_afspc_mode: bool = True,
) -> Tuple[float, float, float, float, float]:
    """Deep space long period periodic contributions to mean elements.

    References:
        - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
        - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
        - Hoots, Schumacher, and Glover, 2004
        - Vallado, Crawford, Hujsak, Kelso, 2006

    Args:
        dscom_out (DscomOutput): Deep space common terms
        t (float): Time since epoch (units = ?)  # TODO: check units
        ep (float): Eccentricity
        inclp (float): Inclination in radians
        nodep (float): RAAN (right ascension of ascending node) in radians
        argpp (float): Argument of perigee in radians
        mp (float): Mean anomaly in radians
        use_afspc_mode (bool): Flag to use AFSPC mode (default = True)

    Returns:
        tuple: (ep, inclp, nodep, argpp, mp)
            ep (float): Updated eccentricity
            inclp (float): Updated inclination in radians
            nodep (float): Updated RAAN (right ascension of ascending node) in radians
            argpp (float): Updated argument of perigee in radians
            mp (float): Updated mean anomaly in radians

    Notes:
        - These periodics are zero at epoch by design.
    """
    # Constants
    zns = 1.19459e-5
    zes = 0.01675
    znl = 1.5835218e-4
    zel = 0.05490

    # Calculate time-varying periodics
    zm = dscom_out.zmos + zns * t
    zf = zm + 2.0 * zes * np.sin(zm)
    sinzf = np.sin(zf)
    f2 = 0.5 * sinzf**2 - 0.25
    f3 = -0.5 * sinzf * np.cos(zf)
    ses = dscom_out.se2 * f2 + dscom_out.se3 * f3
    sis = dscom_out.si2 * f2 + dscom_out.si3 * f3
    sls = dscom_out.sl2 * f2 + dscom_out.sl3 * f3 + dscom_out.sl4 * sinzf
    sghs = dscom_out.sgh2 * f2 + dscom_out.sgh3 * f3 + dscom_out.sgh4 * sinzf
    shs = dscom_out.sh2 * f2 + dscom_out.sh3 * f3
    zm = dscom_out.zmol + znl * t
    zf = zm + 2.0 * zel * np.sin(zm)
    sinzf = np.sin(zf)
    f2 = 0.5 * sinzf**2 - 0.25
    f3 = -0.5 * sinzf * np.cos(zf)
    sel = dscom_out.ee2 * f2 + dscom_out.e3 * f3
    sil = dscom_out.xi2 * f2 + dscom_out.xi3 * f3
    sll = dscom_out.xl2 * f2 + dscom_out.xl3 * f3 + dscom_out.xl4 * sinzf
    sghl = dscom_out.xgh2 * f2 + dscom_out.xgh3 * f3 + dscom_out.xgh4 * sinzf
    shll = dscom_out.xh2 * f2 + dscom_out.xh3 * f3

    # Initialize periodics
    pe = ses + sel
    pinc = sis + sil
    pl = sls + sll
    pgh = sghs + sghl
    ph = shs + shll

    # Subtract baseline offsets
    pe -= dscom_out.peo
    pinc -= dscom_out.pinco
    pl -= dscom_out.plo
    pgh -= dscom_out.pgho
    ph -= dscom_out.pho

    # Update inclination and eccentricity
    inclp += pinc
    ep += pe
    sinip = np.sin(inclp)
    cosip = np.cos(inclp)

    # Apply periodics based on inclination
    if inclp >= 0.2:
        # Apply periodics directly
        ph /= sinip
        pgh -= cosip * ph
        argpp += pgh
        nodep += ph
        mp += pl
    else:
        # Apply periodics with Lyddane modification
        sinop, cosop = np.sin(nodep), np.cos(nodep)
        alfdp = sinip * sinop
        betdp = sinip * cosop
        dalf = ph * cosop + pinc * cosip * sinop
        dbet = -ph * sinop + pinc * cosip * cosop
        alfdp += dalf
        betdp += dbet

        # SGP4 fix for AFSPC-written intrinsic functions
        nodep = np.remainder(nodep, const.TWOPI)
        if nodep < 0.0 and use_afspc_mode:
            nodep += const.TWOPI

        xls = mp + argpp + cosip * nodep
        dls = pl + pgh - pinc * nodep * sinip
        xls += dls
        xnoh = nodep
        nodep = np.arctan2(alfdp, betdp)

        # SGP4 fix for AFSPC-written intrinsic functions
        if nodep < 0.0 and use_afspc_mode:
            nodep += const.TWOPI

        if np.abs(xnoh - nodep) > np.pi:
            nodep += const.TWOPI if nodep < xnoh else -const.TWOPI

        mp += pl
        argpp = xls - mp - cosip * nodep

    return ep, inclp, nodep, argpp, mp
