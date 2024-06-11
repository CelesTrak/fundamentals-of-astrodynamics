# -----------------------------------------------------------------------------
# Author: David Vallado
# Date: 25 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------


import numpy as np

from ...constants import ARCSEC2RAD, DEG2ARCSEC, TWOPI


def fundarg(ttt, opt):
    """Calculates the Delaunay variables and planetary values for several
    theories.

    References:
        Vallado: 2022, p. 210-212, 226

    Args:
        ttt (float): Julian centuries of TT
        opt (str): Method option ('06', '02', '96', or '80')

    Returns:
        tuple: A tuple containing:
            - l (float): Delaunay element in radians
            - l1 (float): Delaunay element in radians
            - f (float): Delaunay element in radians
            - d (float): Delaunay element in radians
            - omega (float): Delaunay element in radians
            - lonmer (float): Planetary longitude in radians
            - lonven (float): Planetary longitude in radians
            - lonear (float): Planetary longitude in radians
            - lonmar (float): Planetary longitude in radians
            - lonjup (float): Planetary longitude in radians
            - lonsat (float): Planetary longitude in radians
            - lonurn (float): Planetary longitude in radians
            - lonnep (float): Planetary longitude in radians
            - precrate (float): Precession rate in radians

    TODO:
        - Implement commented out methods (from m-file)?
        - Use enums instead of strings for option/method
    """
    def calc_delunay_elem(ttt, coeffs):
        """Delaunay fundamental arguments formed in arcsec, converted to deg"""
        return (
            (
                (
                    (coeffs[0] * ttt + coeffs[1]) * ttt + coeffs[2]
                ) * ttt + coeffs[3]
            ) * ttt + coeffs[4]
        ) / DEG2ARCSEC

    def calc_delunay_elem_80(ttt, coeffs, extra):
        return (
            (
                (coeffs[0] * ttt + coeffs[1]) * ttt + coeffs[2]
            ) * ttt
        ) / DEG2ARCSEC + extra

    # Determine coefficients from IAU 2006 nutation theory
    if opt == '06':
        # Delaunay fundamental arguments in deg
        l = calc_delunay_elem(  # noqa
            ttt,
            [-0.00024470, 0.051635, 31.8792, 1717915923.2178, 485868.249036]
        )

        l1 = calc_delunay_elem(
            ttt,
            [-0.00001149, 0.000136, -0.5532, 129596581.0481, 1287104.793048]
        )

        f = calc_delunay_elem(
            ttt,
            [0.00000417, -0.001037, -12.7512, 1739527262.8478, 335779.526232]
        )

        d = calc_delunay_elem(
            ttt,
            [-0.00003169, 0.006593, -6.3706, 1602961601.2090, 1072260.703692]
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
        precrate = ((0.024381750 + 0.00000538691 * ttt) * ttt)

    # Determine coefficients from IAU 2000b theory
    elif opt == '02':
        # Delaunay fundamental arguments in deg
        l = 134.96340251 + (1717915923.2178 * ttt) / DEG2ARCSEC  # noqa
        l1 = 357.52910918 + (129596581.0481 * ttt) / DEG2ARCSEC
        f = 93.27209062 + (1739527262.8478 * ttt) / DEG2ARCSEC
        d = 297.85019547 + (1602961601.2090 * ttt) / DEG2ARCSEC
        omega = 125.04455501 + (-6962890.5431 * ttt) / DEG2ARCSEC

        # Planetary arguments in deg
        (lonmer, lonven, lonear, lonmar, lonjup,
         lonsat, lonurn, lonnep, precrate) = (0.0,) * 9

    # Determine coefficients from IAU 1996 theory
    elif opt == '96':
        # Delaunay fundamental arguments in deg
        l = calc_delunay_elem(  # noqa
            ttt,
            [-0.00024470, 0.051635, 31.8792, 1717915923.2178, 0.0]
        ) + 134.96340251

        l1 = calc_delunay_elem(
            ttt,
            [-0.00001149, -0.000136, -0.5532, 129596581.0481, 0.0]
        ) + 357.52910918

        f = calc_delunay_elem(
            ttt,
            [0.00000417, 0.001037, -12.7512, 1739527262.8478, 0.0]
        ) + 93.27209062

        d = calc_delunay_elem(
            ttt,
            [-0.00003169, 0.006593, -6.3706, 1602961601.2090, 0.0]
        ) + 297.85019547

        omega = calc_delunay_elem(
            ttt, [-0.00005939, 0.007702, 7.4722, -6962890.2665, 0.0]
        ) + 125.04455501

        # Planetary arguments in deg
        lonmer, lonurn, lonnep = (0.0,) * 3
        lonven = 181.979800853 + 58517.8156748 * ttt
        lonear = 100.466448494 + 35999.3728521 * ttt
        lonmar = 355.433274605 + 19140.299314 * ttt
        lonjup = 34.351483900 + 3034.90567464 * ttt
        lonsat = 50.0774713998 + 1222.11379404 * ttt
        precrate = 1.39697137214 * ttt + 0.0003086 * ttt ** 2

    elif opt == '80':
        # Delaunay fundamental arguments in deg
        l = calc_delunay_elem_80(  # noqa
            ttt, [0.064, 31.310, 1717915922.6330], 134.96298139
        )

        l1 = calc_delunay_elem_80(
            ttt, [-0.012, -0.577, 129596581.2240], 357.52772333
        )

        f = calc_delunay_elem_80(
            ttt, [0.011, -13.257, 1739527263.1370], 93.27191028
        )

        d = calc_delunay_elem_80(
            ttt, [0.019, -6.891, 1602961601.3280], 297.85036306
        )

        omega = calc_delunay_elem_80(
            ttt, [0.008, 7.455, -6962890.5390], 125.04452222
        )

        # Planetary arguments in deg
        lonmer = 252.3 + 149472.0 * ttt
        lonven = 179.9 + 58517.8 * ttt
        lonear = 98.4 + 35999.4 * ttt
        lonmar = 353.3 + 19140.3 * ttt
        lonjup = 32.3 + 3034.9 * ttt
        lonsat = 48.0 + 1222.1 * ttt
        lonurn, lonnep, precrate = (0.0,) * 3
    else:
        raise ValueError(
            "Method must be one of the following: '50', '80', or '06'"
        )

    # Convert units to radians
    twopi_deg = np.degrees(TWOPI)
    l = np.radians(np.mod(l, twopi_deg))  # noqa
    l1 = np.radians(np.mod(l1, twopi_deg))
    f = np.radians(np.mod(f, twopi_deg))
    d = np.radians(np.mod(d, twopi_deg))
    omega = np.radians(np.mod(omega, twopi_deg))
    lonmer = np.radians(np.mod(lonmer, twopi_deg))
    lonven = np.radians(np.mod(lonven, twopi_deg))
    lonear = np.radians(np.mod(lonear, twopi_deg))
    lonmar = np.radians(np.mod(lonmar, twopi_deg))
    lonjup = np.radians(np.mod(lonjup, twopi_deg))
    lonsat = np.radians(np.mod(lonsat, twopi_deg))
    lonurn = np.radians(np.mod(lonurn, twopi_deg))
    lonnep = np.radians(np.mod(lonnep, twopi_deg))
    precrate = np.radians(np.mod(precrate, twopi_deg))

    return (
        l, l1, f, d, omega,
        lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep,
        precrate
    )


def precess(ttt, opt):
    """Calculates the transformation matrix that accounts for the effects of
    precession. Both the 1980 and 2006 IAU theories are handled, as well as the

    References:
        Vallado: 2022, p. 219, 227-229

    Args:
        ttt (float): Julian centuries of Terrestrail Time (TT)
        opt (str): Method option ('50', '80', or '06')
                   '50' = FK4 B1950
                   '80' = IAU 1980
                   '06' = IAU 2006

    Returns:
        tuple: A tuple containing:
            - prec (np.array): Transformation matrix for MOD to J2000
            - psia (float): Canonical precession angle in radians
            - wa (float): Canonical precession angle in radians
            - ea (float): Canonical precession angle in radians
            - xa (float): Canonical precession angle in radians

    TODO:
        - Implement commented out methods (from m-file)?
        - Use enums instead of strings for option/method
    """
    def calc_prec_angle(ttt, coeffs):
        return (
            (
                (
                    (coeffs[0] * ttt + coeffs[1]) * ttt + coeffs[2]
                ) * ttt + coeffs[3]
            ) * ttt + coeffs[4]
        ) * ttt

    # Initialize some variables
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    prec = np.eye(3)

    # FK4 B1950 precession angles
    if opt == '50':
        # Commenting these out because they seem important but not used
        # TODO: Decide if these need to be used instead of definitions below
        # psia = 50.3708 + 0.0050 * ttt
        # wa = 0.0
        # ea = 84428.26 - 46.845 * ttt - 0.00059 * ttt2 + 0.00181 * ttt3
        xa = 0.1247 - 0.0188 * ttt

        # GTDS pg 3-17 using days from 1950 - avoids long precession constants
        zeta = 2304.9969 * ttt + 0.302 * ttt2 + 0.01808 * ttt3
        theta = 2004.2980 * ttt - 0.425936 * ttt2 - 0.0416 * ttt3
        z = 2304.9969 * ttt + 1.092999 * ttt2 + 0.0192 * ttt3

        # ttt is tropical centuries from 1950 (36524.22 days)
        prec[0, 0] = 1.0 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3
        prec[0, 1] = 2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3
        prec[0, 2] = 9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3
        prec[1, 0] = -prec[0, 1]
        prec[1, 1] = 1.0 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3
        prec[1, 2] = -1.0858e-4 * ttt2
        prec[2, 0] = -prec[0, 2]
        prec[2, 1] = prec[1, 2]
        prec[2, 2] = 1.0 - 4.721e-5 * ttt2

        # Pass these back out for testing
        # TODO: decide if these need to be removed
        psia = zeta
        wa = theta
        ea = z

    # IAU 80 precession angles
    elif opt == '80':
        psia = 5038.7784 * ttt - 1.07259 * ttt2 - 0.001147 * ttt3
        wa = 84381.448 + 0.05127 * ttt2 - 0.007726 * ttt3
        ea = 84381.448 - 46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3
        xa = 10.5526 * ttt - 2.38064 * ttt2 - 0.001125 * ttt3

        zeta = 2306.2181 * ttt + 0.30188 * ttt2 + 0.017998 * ttt3
        theta = 2004.3109 * ttt - 0.42665 * ttt2 - 0.041833 * ttt3
        z = 2306.2181 * ttt + 1.09468 * ttt2 + 0.018203 * ttt3

    # IAU 06 precession angles
    elif opt == '06':
        oblo = 84381.406
        psia = calc_prec_angle(
            ttt,
            [-0.0000000951, 0.000132851, -0.00114045, -1.0790069, 5038.481507]
        )
        wa = calc_prec_angle(
            ttt,
            [0.0000003337, -0.000000467, -0.00772503, 0.0512623, -0.025754]
        ) + oblo
        ea = calc_prec_angle(
            ttt,
            [-0.0000000434, -0.000000576, 0.00200340, -0.0001831, -46.836769]
        ) + oblo
        xa = calc_prec_angle(
            ttt,
            [-0.0000000560, 0.000170663, -0.00121197, -2.3814292, 10.556403]
        )
        zeta = calc_prec_angle(
            ttt,
            [-0.0000003173, -0.000005971, 0.01801828, 0.2988499, 2306.083227]
        ) + 2.650545
        theta = calc_prec_angle(
            ttt,
            [-0.0000001274, -0.000007089, -0.04182264, -0.4294934, 2004.191903]
        )
        z = calc_prec_angle(
            ttt,
            [0.0000002904, -0.000028596, 0.01826837, 1.0927348, 2306.077181]
        ) - 2.650545
    else:
        raise ValueError(
            "Method must be one of the following: '50', '80', or '06'"
        )

    # Convert units to radians
    zeta *= ARCSEC2RAD
    theta *= ARCSEC2RAD
    z *= ARCSEC2RAD

    # IAU precession angles
    if opt in ['80', '06']:
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

    return (
        prec,
        psia * ARCSEC2RAD,
        wa * ARCSEC2RAD,
        ea * ARCSEC2RAD,
        xa * ARCSEC2RAD
    )
