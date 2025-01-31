# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 6 Aug 2008
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike

from ... import constants as const
from ..time.data import IAU80Array
from ..time.frame_conversions import eci2ecef
from ..twobody.frame_conversions import rv2razel, rv2tradec
from ..twobody.kepler import kepler
from ..twobody.utils import site
from .utils import finite_diff


@dataclass
class ObservationData:
    obstype: np.ndarray = None
    time: np.ndarray = None
    timef: np.ndarray = None
    latgd: np.ndarray = None
    lon: np.ndarray = None
    alt: np.ndarray = None
    jdut1: np.ndarray = None
    ttt: np.ndarray = None
    xp: np.ndarray = None
    yp: np.ndarray = None
    rng: np.ndarray = None
    az: np.ndarray = None
    el: np.ndarray = None
    trtasc: np.ndarray = None
    tdecl: np.ndarray = None
    noise_rng: np.ndarray = None
    noise_az: np.ndarray = None
    noise_el: np.ndarray = None
    noise_trtasc: np.ndarray = None
    noise_tdecl: np.ndarray = None


def findatwaatwb(
    iau80arr: IAU80Array,
    firstobs: int,
    lastobs: int,
    obsrecarr: ObservationData,
    percentchg: float,
    deltaamtchg: float,
    xnom: ArrayLike,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """Find the a and b matrices for the differential correction problem.

    References:
        Vallado: 2022, p. 782-785

    Args:
        iau80arr (IAU80Array): IAU 1980 data for nutation
        firstobs (int): Index of the first observation
        lastobs (int): Index of the last observation
        obsrecarr (ObservationData): Observation data records
        percentchg (float): Amount to modify the vectors in finite differencing
        deltaamtchg (float): Tolerance for small value in finite differencing
        xnom (np.ndarray): State vector in km and km/s

    Returns:
        tuple (atwa, atwb, atw, b, drng2, daz2, del2):
            atwa (np.ndarray): atwa matrix in km and rad
            atwb (np.ndarray): atwb matrix in km and rad
            atw (np.ndarray): atw matrix in km and rad
            b (np.ndarray): Matrix of observation residuals in km and rad
            drng2 (float): Range residual squared in km^2
            daz2 (float): Azimuth residual squared in rad^2
            del2 (float): Elevation residual squared in rad^2

    Notes:
        - It isn't critical for the propagations to use the highest fidelity techniques
        because we're only trying to find the "slope."

    TODO:
        - Verify units of outputs
        - This function is likely stale and will need to be updated
    """

    def adjust_angle(angle):
        if abs(angle) > np.pi:
            angle -= np.sign(angle) * const.TWOPI
        return angle

    # Initialize variables
    statesize = len(xnom)
    atwa = np.zeros((statesize, statesize))
    atwb = np.zeros((statesize, 1))
    atw = np.zeros((statesize, 3))
    a = np.zeros((3, statesize))
    b = np.zeros((3, 1))
    drng2 = daz2 = del2 = 0

    # Get observation index
    indobs = 1 if obsrecarr.obstype[0] == 0 else 3 if obsrecarr.obstype[0] == 2 else 2

    # Loop through all observations
    aznom = elnom = trtascnom = tdeclnom = 0
    for obsktr in range(firstobs - 1, lastobs):
        # Propagate the nominal vector to the epoch time
        dtsec = (
            obsrecarr.time[obsktr]
            - obsrecarr.time[0]
            + obsrecarr.timef[obsktr]
            - obsrecarr.timef[0]
        ) * const.DAY2SEC
        rnom, vnom = xnom[:3], xnom[3:]
        reci1, veci1 = kepler(rnom, vnom, dtsec)

        # Convert ECI to ECEF
        aeci = np.zeros(3)
        lod = ddpsi = ddeps = 0
        recef1, vecef1, _ = eci2ecef(
            reci1,
            veci1,
            aeci,
            obsrecarr.ttt[obsktr],
            obsrecarr.jdut1[obsktr],
            lod,
            obsrecarr.xp[obsktr],
            obsrecarr.yp[obsktr],
            ddpsi,
            ddeps,
            iau80arr,
        )

        # Get site vectors
        rsecef, vsecef = site(
            obsrecarr.latgd[obsktr], obsrecarr.lon[obsktr], obsrecarr.alt[obsktr]
        )

        # Get range, azimuth, and elevation (or right ascension and declination)
        if obsrecarr.obstype[obsktr] != 3:
            rngnom, aznom, elnom, *_ = rv2razel(
                recef1,
                vecef1,
                obsrecarr.latgd[obsktr],
                obsrecarr.lon[obsktr],
                obsrecarr.alt[obsktr],
            )
        else:
            rngnom, trtascnom, tdeclnom, *_ = rv2tradec(recef1, vecef1, rsecef, vsecef)

        # Calculate b matrix
        if obsrecarr.obstype[obsktr] == 0:
            b[0, 0] = obsrecarr.rng[obsktr] - rngnom
        elif obsrecarr.obstype[obsktr] == 1:
            b[0, 0] = obsrecarr.az[obsktr] - aznom
            b[0, 0] = adjust_angle(b[0, 0])
            b[1, 0] = obsrecarr.el[obsktr] - elnom
        elif obsrecarr.obstype[obsktr] == 2:
            b[0, 0] = obsrecarr.rng[obsktr] - rngnom
            b[1, 0] = obsrecarr.az[obsktr] - aznom
            b[1, 0] = adjust_angle(b[1, 0])
            b[2, 0] = obsrecarr.el[obsktr] - elnom
        else:
            b[0, 0] = obsrecarr.trtasc[obsktr] - trtascnom
            b[0, 0] = adjust_angle(b[0, 0])
            b[1, 0] = obsrecarr.tdecl[obsktr] - tdeclnom

        # Perturb each element in the state (elements or vectors)
        for j in range(statesize):
            # Use finite differencing to perturb state
            deltaamt, xnomp = finite_diff(j, percentchg, deltaamtchg, xnom)

            # Propagate the perturbed vector to the epoch time
            rnomp, vnomp = xnomp[:3], xnomp[3:]
            reci3, veci3 = kepler(rnomp, vnomp, dtsec)

            # Convert ECI to ECEF
            recef3, vecef3, _ = eci2ecef(
                reci3,
                veci3,
                aeci,
                obsrecarr.ttt[obsktr],
                obsrecarr.jdut1[obsktr],
                lod,
                obsrecarr.xp[obsktr],
                obsrecarr.yp[obsktr],
                ddpsi,
                ddeps,
                iau80arr,
            )

            # Get range, azimuth, and elevation (or right ascension and declination)
            rngpert = azpert = elpert = trtascpert = tdeclpert = 0
            if obsrecarr.obstype[obsktr] == 3:
                trrpert, trtascpert, tdeclpert, *_ = rv2tradec(
                    recef3, vecef3, rsecef, vsecef
                )
            else:
                rngpert, azpert, elpert, *_ = rv2razel(
                    recef3,
                    vecef3,
                    obsrecarr.latgd[obsktr],
                    obsrecarr.lon[obsktr],
                    obsrecarr.alt[obsktr],
                )

            # Calculate a matrix
            if obsrecarr.obstype[obsktr] == 0:
                a[0, j] = (rngpert - rngnom) / deltaamt
            elif obsrecarr.obstype[obsktr] == 1:
                a[0, j] = (azpert - aznom) / deltaamt
                a[1, j] = (elpert - elnom) / deltaamt
            elif obsrecarr.obstype[obsktr] == 2:
                a[0, j] = (rngpert - rngnom) / deltaamt
                a[1, j] = (azpert - aznom) / deltaamt
                a[2, j] = (elpert - elnom) / deltaamt
            else:
                a[0, j] = (trtascpert - trtascnom) / deltaamt
                a[1, j] = (tdeclpert - tdeclnom) / deltaamt

        # Form matrix combinations
        at = a.T

        # Assign weights
        w2 = w3 = weight = 1
        if obsrecarr.obstype[obsktr] == 0:
            w1 = 1 / (obsrecarr.noise_rng[obsktr] ** 2)
            drng2 += b[0, 0] ** 2 * w1
        elif obsrecarr.obstype[obsktr] == 1:
            w1 = 1 / (obsrecarr.noise_az[obsktr] ** 2)
            w2 = 1 / (obsrecarr.noise_el[obsktr] ** 2)
            daz2 += b[0, 0] ** 2 * w1
            del2 += b[1, 0] ** 2 * w2
        elif obsrecarr.obstype[obsktr] == 2:
            w1 = 1 / (obsrecarr.noise_rng[obsktr] ** 2)
            w2 = 1 / (obsrecarr.noise_az[obsktr] ** 2)
            w3 = 1 / (obsrecarr.noise_el[obsktr] ** 2)
            drng2 += b[0, 0] ** 2 * w1
            daz2 += b[1, 0] ** 2 * w2
            del2 += b[2, 0] ** 2 * w3
        else:
            w1 = 1 / (obsrecarr.noise_trtasc[obsktr] ** 2)
            w2 = 1 / (obsrecarr.noise_tdecl[obsktr] ** 2)
            drng2 += b[0, 0] ** 2 * w1
            daz2 += b[1, 0] ** 2 * w2

        # Form the atw matrix
        atw = np.zeros((statesize, indobs))
        for rowc in range(statesize):
            for colc in range(indobs):
                if colc == 0:
                    weight = w1
                elif colc == 1:
                    weight = w2
                elif colc == 2:
                    weight = w3
                atw[rowc, colc] = at[rowc, colc] * weight

        # Form the atwa and atwb matrices
        atwaacc, atwbacc = np.dot(atw, a), np.dot(atw, b)
        atwa += atwaacc
        atwb[:, 0] += atwbacc[:, 0]

    return atwa, atwb, atw, b, drng2, daz2, del2
