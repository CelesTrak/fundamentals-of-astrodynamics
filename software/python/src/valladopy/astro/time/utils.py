# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 25 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import math
import numpy as np
from typing import Tuple

from .data import iau80in
from ...constants import ARCSEC2RAD, DEG2ARCSEC, TWOPI


def fundarg(
    ttt: float, opt: str
) -> Tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
]:
    """Calculates the Delaunay variables and planetary values for several theories.

    References:
        Vallado: 2022, p. 210-212, 226

    Args:
        ttt (float): Julian centuries of TT
        opt (str): Method option ('06', '02', '96', or '80')

    Returns:
        tuple: (l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat,
                lonurn, lonnep, precrate)
            l (float): Delaunay element in radians
            l1 (float): Delaunay element in radians
            f (float): Delaunay element in radians
            d (float): Delaunay element in radians
            omega (float): Delaunay element in radians
            lonmer (float): Longitude of Mercury in radians
            lonven (float): Longitude of Venus in radians
            lonear (float): Longitude of Earth in radians
            lonmar (float): Longitude of Mars in radians
            lonjup (float): Longitude of Jupiter in radians
            lonsat (float): Longitude of Saturn in radians
            lonurn (float): Longitude of Uranus in radians
            lonnep (float): Longitude of Neptune in radians
            precrate (float): Precession rate in radians per Julian century

    TODO:
        - Implement commented out methods (from m-file)?
        - Use enums instead of strings for option/method
    """

    def calc_delunay_elem(ttt, coeffs):
        """Delaunay fundamental arguments formed in arcsec, converted to deg"""
        return (
            (((coeffs[0] * ttt + coeffs[1]) * ttt + coeffs[2]) * ttt + coeffs[3]) * ttt
            + coeffs[4]
        ) / DEG2ARCSEC

    def calc_delunay_elem_80(ttt, coeffs, extra):
        return (
            ((coeffs[0] * ttt + coeffs[1]) * ttt + coeffs[2]) * ttt
        ) / DEG2ARCSEC + extra

    # Determine coefficients from IAU 2006 nutation theory
    if opt == "06":
        # Delaunay fundamental arguments in deg
        l = calc_delunay_elem(  # noqa
            ttt, [-0.0002447, 0.051635, 31.8792, 1717915923.2178, 485868.249036]
        )
        l1 = calc_delunay_elem(
            ttt, [-0.00001149, 0.000136, -0.5532, 129596581.0481, 1287104.793048]
        )
        f = calc_delunay_elem(
            ttt, [0.00000417, -0.001037, -12.7512, 1739527262.8478, 335779.526232]
        )
        d = calc_delunay_elem(
            ttt, [-0.00003169, 0.006593, -6.3706, 1602961601.209, 1072260.703692]
        )
        omega = calc_delunay_elem(
            ttt, [-0.00005939, 0.007702, 7.4722, -6962890.5431, 450160.398036]
        )

        # Planetary arguments in deg (from TN-36)
        lonmer = np.mod((4.402608842 + 2608.7903141574 * ttt), TWOPI)
        lonven = np.mod((3.176146697 + 1021.3285546211 * ttt), TWOPI)
        lonear = np.mod((1.753470314 + 628.3075849991 * ttt), TWOPI)
        lonmar = np.mod((6.203480913 + 334.0612426700 * ttt), TWOPI)
        lonjup = np.mod((0.599546497 + 52.9690962641 * ttt), TWOPI)
        lonsat = np.mod((0.874016757 + 21.3299104960 * ttt), TWOPI)
        lonurn = np.mod((5.481293872 + 7.4781598567 * ttt), TWOPI)
        lonnep = np.mod((5.311886287 + 3.8133035638 * ttt), TWOPI)
        precrate = (0.024381750 + 0.00000538691 * ttt) * ttt

    # Determine coefficients from IAU 2000b theory
    elif opt == "02":
        # Delaunay fundamental arguments in deg
        l = 134.96340251 + (1717915923.2178 * ttt) / DEG2ARCSEC  # noqa
        l1 = 357.52910918 + (129596581.0481 * ttt) / DEG2ARCSEC
        f = 93.27209062 + (1739527262.8478 * ttt) / DEG2ARCSEC
        d = 297.85019547 + (1602961601.2090 * ttt) / DEG2ARCSEC
        omega = 125.04455501 + (-6962890.5431 * ttt) / DEG2ARCSEC

        # Planetary arguments in deg
        (lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate) = (
            0,
        ) * 9

    # Determine coefficients from IAU 1996 theory
    elif opt == "96":
        # Delaunay fundamental arguments in deg
        l = (  # noqa
            calc_delunay_elem(ttt, [-0.00024470, 0.051635, 31.8792, 1717915923.2178, 0])
            + 134.96340251
        )
        l1 = (
            calc_delunay_elem(ttt, [-0.00001149, -0.000136, -0.5532, 129596581.0481, 0])
            + 357.52910918
        )
        f = (
            calc_delunay_elem(ttt, [0.00000417, 0.001037, -12.7512, 1739527262.8478, 0])
            + 93.27209062
        )
        d = (
            calc_delunay_elem(ttt, [-0.00003169, 0.006593, -6.3706, 1602961601.2090, 0])
            + 297.85019547
        )
        omega = (
            calc_delunay_elem(ttt, [-0.00005939, 0.007702, 7.4722, -6962890.2665, 0])
            + 125.04455501
        )

        # Planetary arguments in deg
        lonmer, lonurn, lonnep = (0,) * 3
        lonven = 181.979800853 + 58517.8156748 * ttt
        lonear = 100.466448494 + 35999.3728521 * ttt
        lonmar = 355.433274605 + 19140.299314 * ttt
        lonjup = 34.351483900 + 3034.90567464 * ttt
        lonsat = 50.0774713998 + 1222.11379404 * ttt
        precrate = 1.39697137214 * ttt + 0.0003086 * ttt**2

    # Determine coefficients from IAU 1980 theory
    elif opt == "80":
        # Delaunay fundamental arguments in deg
        l = calc_delunay_elem_80(  # noqa
            ttt, [0.064, 31.31, 1717915922.633], 134.96298139
        )
        l1 = calc_delunay_elem_80(ttt, [-0.012, -0.577, 129596581.224], 357.52772333)
        f = calc_delunay_elem_80(ttt, [0.011, -13.257, 1739527263.137], 93.27191028)
        d = calc_delunay_elem_80(ttt, [0.019, -6.891, 1602961601.328], 297.85036306)
        omega = calc_delunay_elem_80(ttt, [0.008, 7.455, -6962890.539], 125.04452222)

        # Planetary arguments in deg
        lonmer = 252.3 + 149472 * ttt
        lonven = 179.9 + 58517.8 * ttt
        lonear = 98.4 + 35999.4 * ttt
        lonmar = 353.3 + 19140.3 * ttt
        lonjup = 32.3 + 3034.9 * ttt
        lonsat = 48 + 1222.1 * ttt
        lonurn, lonnep, precrate = (0,) * 3
    else:
        raise ValueError(
            "Method must be one of the following: '06', '02', '96', or '80'"
        )

    # Convert units to radians
    twopi_deg = np.degrees(TWOPI)
    l = float(np.radians(np.mod(l, twopi_deg)))  # noqa
    l1 = float(np.radians(np.mod(l1, twopi_deg)))
    f = float(np.radians(np.mod(f, twopi_deg)))
    d = float(np.radians(np.mod(d, twopi_deg)))
    omega = float(np.radians(np.mod(omega, twopi_deg)))
    lonmer = float(np.radians(np.mod(lonmer, twopi_deg)))
    lonven = float(np.radians(np.mod(lonven, twopi_deg)))
    lonear = float(np.radians(np.mod(lonear, twopi_deg)))
    lonmar = float(np.radians(np.mod(lonmar, twopi_deg)))
    lonjup = float(np.radians(np.mod(lonjup, twopi_deg)))
    lonsat = float(np.radians(np.mod(lonsat, twopi_deg)))
    lonurn = float(np.radians(np.mod(lonurn, twopi_deg)))
    lonnep = float(np.radians(np.mod(lonnep, twopi_deg)))
    precrate = float(np.radians(np.mod(precrate, twopi_deg)))

    return (
        l,
        l1,
        f,
        d,
        omega,
        lonmer,
        lonven,
        lonear,
        lonmar,
        lonjup,
        lonsat,
        lonurn,
        lonnep,
        precrate,
    )


def precess(ttt: float, opt: str) -> Tuple[np.ndarray, float, float, float, float]:
    """Calculates the transformation matrix that accounts for the effects of
    precession. Both the 1980 and 2006 IAU theories are handled, as well as the
    FK B1950 theory.

    References:
        Vallado: 2022, p. 219, 227-229

    Args:
        ttt (float): Julian centuries of Terrestrail Time (TT)
        opt (str): Method option ('50', '80', or '06')
                   '50' = FK4 B1950
                   '80' = IAU 1980
                   '06' = IAU 2006

    Returns:
        tuple: (prec, psia, wa, ea, xa)
            prec (np.array): Transformation matrix for MOD to J2000
            psia (float): Canonical precession angle in radians
            wa (float): Canonical precession angle in radians
            ea (float): Canonical precession angle in radians
            xa (float): Canonical precession angle in radians

    TODO:
        - Implement commented out methods (from m-file)?
        - Use enums instead of strings for option/method
    """

    def calc_prec_angle(ttt, coeffs):
        return (
            (((coeffs[0] * ttt + coeffs[1]) * ttt + coeffs[2]) * ttt + coeffs[3]) * ttt
            + coeffs[4]
        ) * ttt

    # Initialize some variables
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    prec = np.eye(3)

    # FK4 B1950 precession angles
    if opt == "50":
        # Commenting these out because they seem important but not used
        # TODO: Decide if these need to be used instead of definitions below
        # psia = 50.3708 + 0.0050 * ttt
        # wa = 0.0
        # ea = 84428.26 - 46.845 * ttt - 0.00059 * ttt2 + 0.00181 * ttt3
        xa = 0.1247 - 0.0188 * ttt

        # GTDS pg 3-17 using days from 1950 - avoids long precession constants
        zeta = 2304.9969 * ttt + 0.302 * ttt2 + 0.01808 * ttt3
        theta = 2004.298 * ttt - 0.425936 * ttt2 - 0.0416 * ttt3
        z = 2304.9969 * ttt + 1.092999 * ttt2 + 0.0192 * ttt3

        # ttt is tropical centuries from 1950 (36524.22 days)
        prec[0, 0] = 1 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3
        prec[0, 1] = 2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3
        prec[0, 2] = 9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3
        prec[1, 0] = -prec[0, 1]
        prec[1, 1] = 1 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3
        prec[1, 2] = -1.0858e-4 * ttt2
        prec[2, 0] = -prec[0, 2]
        prec[2, 1] = prec[1, 2]
        prec[2, 2] = 1 - 4.721e-5 * ttt2

        # Pass these back out for testing
        # TODO: decide if these need to be removed
        psia, wa, ea = zeta, theta, z

    # IAU 80 precession angles
    elif opt == "80":
        psia = 5038.7784 * ttt - 1.07259 * ttt2 - 0.001147 * ttt3
        wa = 84381.448 + 0.05127 * ttt2 - 0.007726 * ttt3
        ea = 84381.448 - 46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3
        xa = 10.5526 * ttt - 2.38064 * ttt2 - 0.001125 * ttt3

        zeta = 2306.2181 * ttt + 0.30188 * ttt2 + 0.017998 * ttt3
        theta = 2004.3109 * ttt - 0.42665 * ttt2 - 0.041833 * ttt3
        z = 2306.2181 * ttt + 1.09468 * ttt2 + 0.018203 * ttt3

    # IAU 06 precession angles
    elif opt == "06":
        oblo = 84381.406
        psia = calc_prec_angle(
            ttt, [-0.0000000951, 0.000132851, -0.00114045, -1.0790069, 5038.481507]
        )
        wa = (
            calc_prec_angle(
                ttt, [0.0000003337, -0.000000467, -0.00772503, 0.0512623, -0.025754]
            )
            + oblo
        )
        ea = (
            calc_prec_angle(
                ttt, [-0.0000000434, -0.000000576, 0.0020034, -0.0001831, -46.836769]
            )
            + oblo
        )
        xa = calc_prec_angle(
            ttt, [-0.000000056, 0.000170663, -0.00121197, -2.3814292, 10.556403]
        )
        zeta = (
            calc_prec_angle(
                ttt, [-0.0000003173, -0.000005971, 0.01801828, 0.2988499, 2306.083227]
            )
            + 2.650545
        )
        theta = calc_prec_angle(
            ttt, [-0.0000001274, -0.000007089, -0.04182264, -0.4294934, 2004.191903]
        )
        z = (
            calc_prec_angle(
                ttt, [0.0000002904, -0.000028596, 0.01826837, 1.0927348, 2306.077181]
            )
            - 2.650545
        )
    else:
        raise ValueError("Method must be one of the following: '50', '80', or '06'")

    # Convert units to radians
    zeta *= ARCSEC2RAD
    theta *= ARCSEC2RAD
    z *= ARCSEC2RAD

    # IAU precession angles
    if opt in ["80", "06"]:
        coszeta = np.cos(zeta)
        sinzeta = np.sin(zeta)
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        cosz = np.cos(z)
        sinz = np.sin(z)

        # Form matrix MOD to J2000
        prec[0, 0] = coszeta * costheta * cosz - sinzeta * sinz
        prec[0, 1] = coszeta * costheta * sinz + sinzeta * cosz
        prec[0, 2] = coszeta * sintheta
        prec[1, 0] = -sinzeta * costheta * cosz - coszeta * sinz
        prec[1, 1] = -sinzeta * costheta * sinz + coszeta * cosz
        prec[1, 2] = -sinzeta * sintheta
        prec[2, 0] = -sintheta * cosz
        prec[2, 1] = -sintheta * sinz
        prec[2, 2] = costheta

    return prec, psia * ARCSEC2RAD, wa * ARCSEC2RAD, ea * ARCSEC2RAD, xa * ARCSEC2RAD


def nutation(
    ttt: float, ddpsi: float, ddeps: float
) -> Tuple[float, float, float, float, np.ndarray]:
    """Calculates the transformation matrix that accounts for the effects of nutation.

    References:
        Vallado: 2022, p. 225-227

    Args:
        ttt (float): Julian centuries of TT
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta eps correction to GCRF in radians

    Returns:
        tuple:
            deltapsi (float): Nutation angle in radians
            trueeps (float): True obliquity of the ecliptic in radians
            meaneps (float): Mean obliquity of the ecliptic in radians
            omega (float): Delaunay element in radians
            nut (np.ndarray): Transformation matrix for TOD - MOD
    """
    # Load nutation coefficients
    iar80, rar80 = iau80in()

    # Calculate powers of ttt
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt

    # Mean obliquity of the ecliptic
    meaneps = -46.815 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3 + 84381.448
    meaneps = float(np.radians(np.remainder(meaneps / DEG2ARCSEC, np.degrees(TWOPI))))

    # Fundamental arguments using the IAU80 theory
    (
        l,
        l1,
        f,
        d,
        omega,
        lonmer,
        lonven,
        lonear,
        lonmar,
        lonjup,
        lonsat,
        lonurn,
        lonnep,
        precrate,
    ) = fundarg(ttt, "80")

    # Calculate nutation parameters
    deltapsi, deltaeps = 0, 0
    for i in range(len(iar80)):
        tempval = (
            iar80[i, 0] * l
            + iar80[i, 1] * l1
            + iar80[i, 2] * f
            + iar80[i, 3] * d
            + iar80[i, 4] * omega
        )
        deltapsi += (rar80[i, 0] + rar80[i, 1] * ttt) * np.sin(tempval)
        deltaeps += (rar80[i, 2] + rar80[i, 3] * ttt) * np.cos(tempval)

    # Add corrections
    deltapsi = math.remainder(deltapsi + ddpsi, TWOPI)
    deltaeps = math.remainder(deltaeps + ddeps, TWOPI)
    trueeps = meaneps + deltaeps

    # Sine/cosine values of psi and eps
    cospsi = np.cos(deltapsi)
    sinpsi = np.sin(deltapsi)
    coseps = np.cos(meaneps)
    sineps = np.sin(meaneps)
    costrueeps = np.cos(trueeps)
    sintrueeps = np.sin(trueeps)

    # Construct nutation rotation matrix
    nut = np.zeros((3, 3))
    nut[0, 0] = cospsi
    nut[0, 1] = costrueeps * sinpsi
    nut[0, 2] = sintrueeps * sinpsi
    nut[1, 0] = -coseps * sinpsi
    nut[1, 1] = costrueeps * coseps * cospsi + sintrueeps * sineps
    nut[1, 2] = sintrueeps * coseps * cospsi - sineps * costrueeps
    nut[2, 0] = -sineps * sinpsi
    nut[2, 1] = costrueeps * sineps * cospsi - sintrueeps * coseps
    nut[2, 2] = sintrueeps * sineps * cospsi + costrueeps * coseps

    return deltapsi, trueeps, meaneps, omega, nut


def polarm(xp: float, yp: float, ttt: float, use_iau80: bool = True) -> np.ndarray:
    """Calculate the transformation matrix that accounts for polar motion.

    References:
        Vallado: 2022, p. 213, 224

    Both the 1980 and 2000 theories are handled. Note that the rotation order
    is different between 1980 and 2000.

    Args:
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ttt (float): Julian centuries of TT (only used in IAU 2000 method)
        use_iau80 (bool, optional): Whether to use the IAU 1980 method instead
                                    of IAU 2000 method (defaults to True)

    Returns:
        pm (np.ndarray): Transformation matrix for ECEF to PEF
    """
    cosxp, sinxp = np.cos(xp), np.sin(xp)
    cosyp, sinyp = np.cos(yp), np.sin(yp)

    # Use IAU 1980 theory
    if use_iau80:
        pm = np.array(
            [
                [cosxp, 0, -sinxp],
                [sinxp * sinyp, cosyp, cosxp * sinyp],
                [sinxp * cosyp, -sinyp, cosxp * cosyp],
            ]
        )
    # Use IAU 2000 theory
    else:
        # −47e-6 corresponds to a constant drift in the Terrestrial
        # Intermediate Origin (TIO) locator 𝑠′, which is approximately −47
        # microarcseconds per century (applied to the IAU 2000 theory)
        # See: https://pyerfa.readthedocs.io/en/latest/api/erfa.pom00.html
        # TODO: consider using pyerfa for this
        sp = -47e-6 * ttt * ARCSEC2RAD
        cossp, sinsp = np.cos(sp), np.sin(sp)

        pm = np.array(
            [
                [
                    cosxp * cossp,
                    -cosyp * sinsp + sinyp * sinxp * cossp,
                    -sinyp * sinsp - cosyp * sinxp * cossp,
                ],
                [
                    cosxp * sinsp,
                    cosyp * cossp + sinyp * sinxp * sinsp,
                    sinyp * cossp - cosyp * sinxp * sinsp,
                ],
                [sinxp, -sinyp * cosxp, cosyp * cosxp],
            ]
        )

    return pm
