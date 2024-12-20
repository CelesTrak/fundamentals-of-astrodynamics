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


@dataclass
class DsinitOutput:
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
    incl_tol: float = 0.2,
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
        incl_tol (float): Inclination tolerance for applying periodics (default = 0.2)

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
    if inclp >= incl_tol:
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


def dsinit(
    dscom_out: DscomOutput,
    xke: float,
    argpo: float,
    t: float,
    tc: float,
    gsto: float,
    mo: float,
    mdot: float,
    no: float,
    nodeo: float,
    nodedot: float,
    xpidot: float,
    ecco: float,
    eccsq: float,
    inclm: float,
    nodem: float,
    argpm: float,
    mm: float,
    incl_tol: float = np.radians(3),
) -> DsinitOutput:
    """Deep space initialization for SGP4.

    This function provides deep space contributions to mean motion dot due to
    geopotential resonance with half day and one day orbits.

    References:
        - Hoots, Roehrich, NORAD SpaceTrack Report #3, 1980
        - Hoots, Roehrich, NORAD SpaceTrack Report #6, 1986
        - Hoots, Schumacher, and Glover, 2004
        - Vallado, Crawford, Hujsak, Kelso, 2006

    Args:
        dscom_out (DscomOutput): Deep space common terms
        xke (float): SGP4 constant
        argpo (float): Argument of perigee in radians
        t (float): Time since epoch (units = ?)  # TODO: check units
        tc (float): Time since epoch (units = ?)  # TODO: check units
        gsto (float): GST at epoch in radians
        mo (float): Mean anomaly in radians
        mdot (float): Mean anomaly dot in rad/s
        no (float): Mean motion in rad/s  # TODO: check units
        nodeo (float): RAAN at epoch in radians
        nodedot (float): RAAN dot in rad/s
        xpidot (float): RAAN dot in rad/s
        ecco (float): Eccentricity
        eccsq (float): Eccentricity squared
        inclm (float): Inclination in radians
        nodem (float): RAAN in radians
        argpm (float): Argument of perigee in radians
        mm (float): Mean anomaly in radians
        incl_tol (float): Inclination tolerance for applying periodics
                          (default = 3 deg in radians)

    Returns:
        DsinitOutput: Dataclass containing deep space initialization terms
    """
    # Initialize output dataclass
    dsinit_out = DsinitOutput()

    # Constants
    rptim = 4.37526908801129966e-3
    q22, q31, q33 = 1.7891679e-6, 2.1460748e-6, 2.2123015e-7
    root22, root44, root54 = 1.7891679e-6, 7.3636953e-9, 2.1765803e-9
    root32, root52 = 3.7393792e-7, 1.1428639e-7
    x2o3 = 2 / 3
    znl, zns = 1.5835218e-4, 1.19459e-5

    # Initialize outputs
    dsinit_out.em, dsinit_out.nm = dscom_out.em, dscom_out.nm
    if 0.0034906585 < dsinit_out.nm < 0.0052359877:
        dsinit_out.irez = 1
    if 8.26e-3 <= dsinit_out.nm <= 9.24e-3 and dsinit_out.em >= 0.5:
        dsinit_out.irez = 2
    emsq = dscom_out.emsq

    # Solar terms
    ses = dscom_out.ss1 * zns * dscom_out.ss5
    sis = dscom_out.ss2 * zns * (dscom_out.sz11 + dscom_out.sz13)
    sls = (
        -zns * dscom_out.ss3 * (dscom_out.sz1 + dscom_out.sz3 - 14 - 6 * dscom_out.emsq)
    )
    sghs = dscom_out.ss4 * zns * (dscom_out.sz31 + dscom_out.sz33 - 6)
    shs = -zns * dscom_out.ss2 * (dscom_out.sz21 + dscom_out.sz23)
    if inclm < incl_tol or inclm > np.pi - incl_tol:
        shs = 0
    if dscom_out.sinim != 0:
        shs /= dscom_out.sinim
    sgs = sghs - dscom_out.cosim * shs

    # Lunar terms
    dsinit_out.dedt = ses + dscom_out.s1 * znl * dscom_out.s5
    dsinit_out.didt = sis + dscom_out.s2 * znl * (dscom_out.z11 + dscom_out.z13)
    dsinit_out.dmdt = sls - znl * dscom_out.s3 * (
        dscom_out.z1 + dscom_out.z3 - 14 - 6 * dscom_out.emsq
    )
    sghl = dscom_out.s4 * znl * (dscom_out.z31 + dscom_out.z33 - 6)
    shll = -znl * dscom_out.s2 * (dscom_out.z21 + dscom_out.z23)
    if inclm < incl_tol or inclm > np.pi - incl_tol:
        shll = 0
    dsinit_out.domdt, dsinit_out.dnodt = sgs + sghl, shs
    if dscom_out.sinim != 0:
        dsinit_out.domdt -= dscom_out.cosim / dscom_out.sinim * shll
        dsinit_out.dnodt += shll / dscom_out.sinim

    # Deep space resonance effects
    theta = np.remainder(gsto + tc * rptim, const.TWOPI)
    dsinit_out.em += dsinit_out.dedt * t
    dsinit_out.inclm = inclm + dsinit_out.didt * t
    dsinit_out.argpm = argpm + dsinit_out.domdt * t
    dsinit_out.nodem = nodem + dsinit_out.dnodt * t
    dsinit_out.mm = mm + dsinit_out.dmdt * t

    # Return if no resonance
    if dsinit_out.irez == 0:
        return dsinit_out

    aonv = (dsinit_out.nm / xke) ** x2o3

    # Geopotential resonance for 12-hour orbits
    if dsinit_out.irez == 2:
        cosisq = dscom_out.cosim**2
        emo = dsinit_out.em
        em = ecco
        emsqo = dscom_out.emsq
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

        dsinit_out.em = em
        sini2 = dscom_out.sinim**2
        f220 = 0.75 * (1.0 + 2.0 * dscom_out.cosim + cosisq)
        f221 = 1.5 * sini2
        f321 = 1.875 * dscom_out.sinim * (1 - 2 * dscom_out.cosim - 3 * cosisq)
        f322 = -1.875 * dscom_out.sinim * (1 + 2 * dscom_out.cosim - 3 * cosisq)
        f441 = 35 * sini2 * f220
        f442 = 39.3750 * sini2**2
        f522 = (
            9.84375
            * dscom_out.sinim
            * (
                sini2 * (1 - 2 * dscom_out.cosim - 5 * cosisq)
                + 0.33333333 * (-2 + 4 * dscom_out.cosim + 6 * cosisq)
            )
        )
        f523 = dscom_out.sinim * (
            4.92187512 * sini2 * (-2 - 4 * dscom_out.cosim + 10 * cosisq)
            + 6.56250012 * (1 + 2 * dscom_out.cosim - 3 * cosisq)
        )
        f542 = (
            29.53125
            * dscom_out.sinim
            * (
                2
                - 8 * dscom_out.cosim
                + cosisq * (-12 + 8 * dscom_out.cosim + 10 * cosisq)
            )
        )
        f543 = (
            29.53125
            * dscom_out.sinim
            * (
                -2
                - 8 * dscom_out.cosim
                + cosisq * (12 + 8 * dscom_out.cosim - 10 * cosisq)
            )
        )

        xno2 = dsinit_out.nm**2
        ainv2 = aonv**2
        temp1 = 3 * xno2 * ainv2
        temp = temp1 * root22
        dsinit_out.d2201 = temp * f220 * g201
        dsinit_out.d2211 = temp * f221 * g211
        temp1 *= aonv
        temp = temp1 * root32
        dsinit_out.d3210 = temp * f321 * g310
        dsinit_out.d3222 = temp * f322 * g322
        temp1 *= aonv
        temp = 2 * temp1 * root44
        dsinit_out.d4410 = temp * f441 * g410
        dsinit_out.d4422 = temp * f442 * g422
        temp1 *= aonv
        temp = temp1 * root52
        dsinit_out.d5220 = temp * f522 * g520
        dsinit_out.d5232 = temp * f523 * g532
        temp = 2 * temp1 * root54
        dsinit_out.d5421 = temp * f542 * g521
        dsinit_out.d5433 = temp * f543 * g533
        dsinit_out.xlamo = (mo + nodeo + nodeo - theta - theta) % const.TWOPI
        dsinit_out.xfact = (
            mdot + dsinit_out.dmdt + 2 * (nodedot + dsinit_out.dnodt - rptim) - no
        )
        dsinit_out.em, emsq = emo, emsqo

    # Synchronous resonance terms
    if dsinit_out.irez == 1:
        g200 = 1 + emsq * (-2.5 + 0.8125 * emsq)
        g310 = 1 + 2 * emsq
        g300 = 1 + emsq * (-6 + 6.60937 * emsq)
        f220 = 0.75 * (1 + dscom_out.cosim) ** 2
        f311 = 0.9375 * dscom_out.sinim**2 * (1 + 3 * dscom_out.cosim) - 0.75 * (
            1 + dscom_out.cosim
        )
        f330 = 1 + dscom_out.cosim
        f330 = 1.875 * f330**3
        dsinit_out.del1 = 3 * dsinit_out.nm**2 * aonv**2
        dsinit_out.del2 = 2 * dsinit_out.del1 * f220 * g200 * q22
        dsinit_out.del3 = 3 * dsinit_out.del1 * f330 * g300 * q33 * aonv
        dsinit_out.del1 *= f311 * g310 * q31 * aonv
        dsinit_out.xlamo = (mo + nodeo + argpo - theta) % const.TWOPI
        dsinit_out.xfact = (
            mdot
            + xpidot
            - rptim
            + dsinit_out.dmdt
            + dsinit_out.domdt
            + dsinit_out.dnodt
            - no
        )

    # Initialize the integrator for SGP4
    dsinit_out.xli, dsinit_out.xni = dsinit_out.xlamo, no
    dsinit_out.nm += dsinit_out.dndt

    return dsinit_out
