# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import numpy as np
from numpy.typing import ArrayLike
from typing import Tuple

from . import iau_transform as iau
from .data import IAU80Array, IAU06Array, IAU06xysArray
from .sidereal import gstime, sidereal
from .utils import precess, nutation, polarm
from ... import constants as const


def calc_omegaearth(lod: float) -> np.ndarray:
    """Calculates the Earth's rotation vector.

    Args:
        lod (float): Excess length of day in seconds

    Returns:
        np.ndarray: Earth's rotation vector in rad/s
    """
    return np.array([0, 0, const.EARTHROT * (1 - lod / const.DAY2SEC)])


def calc_orbit_effects(
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
    use_iau80: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates the orbit effects from precession, nutation, sidereal time, and polar
    motion.

    Args:
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)
        use_iau80 (bool, optional): Use IAU 1980 model for precession/nutation
                                    (default True)

    Returns:
        tuple: (prec, nut, st, pm, omegaearth)
            prec (np.ndarray): Transformation matrix for MOD to J2000
            nut (np.ndarray): Transformation matrix for TOD - MOD
            st (np.ndarray): Transformation matrix for PEF to TOD
            pm (np.ndarray): Transformation matrix for ECEF to PEF
            omegaearth (np.ndarray): Earth angular rotation vecctor in rad/s
    """
    # Find matrices that account for various orbit effects
    prec, *_ = precess(ttt, opt="80")
    deltapsi, _, meaneps, omega, nut = nutation(ttt, ddpsi, ddeps, iau80arr)
    st, _ = sidereal(jdut1, deltapsi, meaneps, omega, lod, use_iau80, eqeterms)
    pm = polarm(xp, yp, ttt, use_iau80=use_iau80)

    # Calculate the effects of Earth's rotation
    omegaearth = calc_omegaearth(lod)

    return prec, nut, st, pm, omegaearth


########################################################################################
# ECI <-> ECEF Frame Conversions
########################################################################################


def compute_iau06_matrices(
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float,
    ddy: float,
    use_full_series: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Computes the precession/nutation matrix, sidereal time matrix, polar motion
    matrix, and Earth rotation vector.

    Args:
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians
        ddy (float, optional): EOP correction for y in radians
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data

    Returns:
        tuple: (pnb, st, pm, omegaearth)
            pnb (np.ndarray): Precession/nutation matrix
            st (np.ndarray): Sidereal time matrix
            pm (np.ndarray): Polar motion matrix
            omegaearth (np.ndarray): Earth rotation vector
    """
    # Calculate precession-nutation and sidereal matrices
    *_, pnb = iau.iau06xys(ttt, iau06arr, ddx, ddy, iau06xysarr, use_full_series)
    st, _ = sidereal(jdut1, 0, 0, 0, lod, use_iau80=False)

    # Earth rotation
    omegaearth = calc_omegaearth(lod)

    # Polar motion matrix
    pm = polarm(xp, yp, ttt, use_iau80=False)

    return pnb, st, pm, omegaearth


def eci2ecef(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox frame (J2000) to the
    Earth-fixed frame (ECEF).

    References:
        Vallado: 2022, p. 223-230

    Args:
        reci (array_like): ECi position vector in km
        veci (array_like): ECI velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector
            vecef (np.ndarray): ECEF velocity vector
            aecef (np.ndarray): ECEF acceleration vector

    TODO: The acceleration transformation is not correct and needs to be fixed.
    """
    # Find matrices that account for various orbit effects
    prec, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Position transformation
    rpef = st.T @ nut.T @ prec.T @ reci
    recef = pm.T @ rpef

    # Velocity transformation
    vpef = st.T @ nut.T @ prec.T @ veci - np.cross(omegaearth, rpef)
    vecef = pm.T @ vpef

    # Acceleration transformation
    aecef = (
        pm.T @ (st.T @ nut.T @ prec.T @ aeci)
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return recef, vecef, aecef


def eci2ecef06(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float = 0.0,
    ddy: float = 0.0,
    use_full_series: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI frame (J2000) to the Earth-fixed frame (ECEF)
    using the IAU 2006 xys approach.

    References:
        Vallado, 2022, p. 223-230

    Args:
        reci (ArrayLike): ECI position vector in km
        veci (ArrayLike): ECI velocity vector in km/s
        aeci (ArrayLike): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians (default 0)
        ddy (float, optional): EOP correction for y in radians (default 0)
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data
                                          (default True)

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    # Compute the IAU06 matrices
    pnb, st, pm, omegaearth = compute_iau06_matrices(
        ttt, jdut1, lod, xp, yp, iau06arr, iau06xysarr, ddx, ddy, use_full_series
    )

    # Transform position
    rpef = st.T @ pnb.T @ reci
    recef = pm.T @ rpef

    # Transform velocity
    vpef = st.T @ pnb.T @ veci - np.cross(omegaearth, rpef)
    vecef = pm.T @ vpef

    # Transform acceleration
    aecef = pm.T @ (
        st.T @ pnb.T @ aeci
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return recef, vecef, aecef


def ecef2eci(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed frame (ECEF) to the ECI mean equator,
    mean equinox (J2000) frame.

    References:
        Vallado: 2022, p. 223-230

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Find matrices that account for various orbit effects
    prec, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Position transformation
    rpef = np.dot(pm, recef)
    reci = np.dot(prec, np.dot(nut, np.dot(st, rpef)))

    # Velocity transformation
    vpef = np.dot(pm, vecef)
    veci = np.dot(prec, np.dot(nut, np.dot(st, vpef + np.cross(omegaearth, rpef))))

    # Acceleration transformation
    aeci = (
        np.dot(prec, np.dot(nut, np.dot(st, np.dot(pm, aecef))))
        + np.cross(omegaearth, np.cross(omegaearth, rpef))
        + 2 * np.cross(omegaearth, vpef)
    )

    return reci, veci, aeci


def ecef2eci06(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float = 0.0,
    ddy: float = 0.0,
    use_full_series: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed frame (ECEF) to the ECI (J2000) frame.

    References:
        Vallado, 2022, p. 223-230

    Args:
        recef (ArrayLike): ECEF position vector in km
        vecef (ArrayLike): ECEF velocity vector in km/s
        aecef (ArrayLike): ECEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians (default 0)
        ddy (float, optional): EOP correction for y in radians (default 0)
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Compute the IAU06 matrices
    pnb, st, pm, omegaearth = compute_iau06_matrices(
        ttt, jdut1, lod, xp, yp, iau06arr, iau06xysarr, ddx, ddy, use_full_series
    )

    # Transform position
    rpef = pm @ recef
    reci = pnb @ st @ rpef

    # Transform velocity
    vpef = pm @ vecef
    veci = pnb @ st @ (vpef + np.cross(omegaearth, rpef))

    # Transform acceleration
    temp = np.cross(omegaearth, rpef)
    aeci = (
        pnb
        @ st
        @ (pm @ aecef + np.cross(omegaearth, temp) + 2 * np.cross(omegaearth, vpef))
    )

    return reci, veci, aeci


########################################################################################
# ECI <-> PEF Frame Conversions
########################################################################################


def eci2pef(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox frame (J2000), to the
    pseudo Earth-fixed frame (PEF).

    References:
        Vallado: 2022, p. 224

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECi velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        ddpsi (float): Nutation correction for delta psi in radians
        ddeps (float): Nutation correction for delta epsilon in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool): Add terms for ast calculation (default True)

    Returns:
        tuple: (rpef, vpef, apef)
            rpef (np.ndarray): PEF position vector in km
            vpef (np.ndarray): PEF velocity vector in km/s
            apef (np.ndarray): PEF acceleration vector in km/s²
    """
    # Compute the IAU matrices
    prec, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, 0, 0, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Transform vectors
    rpef = st.T @ nut.T @ prec.T @ reci
    vpef = st.T @ nut.T @ prec.T @ veci - np.cross(omegaearth, rpef)
    apef = (
        st.T @ nut.T @ prec.T @ aeci
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return rpef, vpef, apef


def pef2eci(
    rpef: ArrayLike,
    vpef: ArrayLike,
    apef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the pseudo Earth-fixed (PEF) frame to the ECI mean
    equator, mean equinox (J2000) frame using the IAU 1980 model.

    References:
        Vallado: 2022, p. 224

    Args:
        rpef (array_like): PEf position vector in km
        vpef (array_like): PEF velocity vector in km/s
        apef (array_like): PEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Find matrices that account for various orbit effects
    prec, nut, st, _, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, 0, 0, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Transform vectors
    reci = prec @ nut @ st @ np.asarray(rpef)
    veci = prec @ nut @ st @ (np.asarray(vpef) + np.cross(omegaearth, rpef))
    aeci = (
        prec
        @ nut
        @ st
        @ (
            np.asarray(apef)
            + np.cross(omegaearth, np.cross(omegaearth, rpef))
            + 2 * np.cross(omegaearth, vpef)
        )
    )

    return reci, veci, aeci


########################################################################################
# ECI <-> TOD Frame Conversions
########################################################################################


def eci2tod(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox frame (J2000) to the
    true equator, true equinox of date (TOD) frame.

    References:
        Vallado: 2022, p. 225

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation

    Returns:
        tuple: (rtod, vtod, atod)
            rtod (np.ndarray): TOD position vector in km
            vtod (np.ndarray): TOD velocity vector in km/s
            atod (np.ndarray): TOD acceleration vector in km/s²
    """
    # Precession (IAU 1980 model)
    prec, *_ = precess(ttt, opt="80")

    # Nutation
    *_, nut = nutation(ttt, ddpsi, ddeps, iau80arr)

    # Transform vectors
    rtod = nut.T @ prec.T @ np.asarray(reci)
    vtod = nut.T @ prec.T @ np.asarray(veci)
    atod = nut.T @ prec.T @ np.asarray(aeci)

    return rtod, vtod, atod


def tod2eci(
    rtod: ArrayLike,
    vtod: ArrayLike,
    atod: ArrayLike,
    ttt: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the true equator, true equinox of date (TOD) frame to
    the ECI mean equator, mean equinox (J2000) frame.

    References:
        Vallado: 2022, p. 225

    Args:
        rtod (array_like): TOD position vector in km
        vtod (array_like): TOD velocity vector in km/s
        atod (array_like): TOD acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Precession (IAU 1980 model)
    prec, *_ = precess(ttt, opt="80")

    # Nutation
    *_, nut = nutation(ttt, ddpsi, ddeps, iau80arr)

    # Transform vectors
    reci = prec @ nut @ np.asarray(rtod)
    veci = prec @ nut @ np.asarray(vtod)
    aeci = prec @ nut @ np.asarray(atod)

    return reci, veci, aeci


########################################################################################
# ECI <-> MOD Frame Conversions
########################################################################################


def eci2mod(
    reci: ArrayLike, veci: ArrayLike, aeci: ArrayLike, ttt: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox frame (J2000) to the
    mean equator, mean equinox of date (MOD) frame.

    References:
        Vallado: 2022, p. 227

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT

    Returns:
        tuple: (rmod, vmod, amod)
            rmod (np.ndarray): MOD position vector in km
            vmod (np.ndarray): MOD velocity vector in km/s
            amod (np.ndarray): MOD acceleration vector in km/s²
    """
    # Precession (IAU 1980 model)
    prec, *_ = precess(ttt, opt="80")

    # Transform vectors
    rmod = prec.T @ np.asarray(reci)
    vmod = prec.T @ np.asarray(veci)
    amod = prec.T @ np.asarray(aeci)

    return rmod, vmod, amod


def mod2eci(
    rmod: ArrayLike, vmod: ArrayLike, amod: ArrayLike, ttt: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the mean equator, mean equinox of date (MOD) frame
    to the ECI mean equator, mean equinox (J2000) frame.

    References:
        Vallado: 2022, p. 227

    Args:
        rmod (array_like): MOD position vector in km
        vmod (array_like): MOD velocity vector in km/s
        amod (array_like): MOD acceleration vector in km/s²
        ttt (float): Julian centuries of TT

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Precession (IAU 1980 model)
    prec, *_ = precess(ttt, opt="80")

    # Transform vectors
    reci = prec @ np.asarray(rmod)
    veci = prec @ np.asarray(vmod)
    aeci = prec @ np.asarray(amod)

    return reci, veci, aeci


########################################################################################
# ECI <-> TEME Frame Conversions
########################################################################################


def _get_teme_eci_transform(
    ttt: float, ddpsi: float, ddeps: float, iau80arr: IAU80Array
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generates the transformation matrices for TEME ↔ ECI.

    Args:
        ttt (float): Julian centuries of TT
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation

    Returns:
        tuple: (prec, nut, eqe)
            prec (np.ndarray): Precession matrix
            nut (np.ndarray): Nutation matrix
            eqe (np.ndarray): Equation of equinoxes rotation matrix
    """
    # Precession (IAU 1980 model)
    prec, *_ = precess(ttt, opt="80")

    # Nutation
    deltapsi, _, meaneps, _, nut = nutation(ttt, ddpsi, ddeps, iau80arr)

    # Equation of equinoxes (geometric terms only)
    eqeg = deltapsi * np.cos(meaneps)
    eqeg = np.remainder(eqeg, const.TWOPI)

    # Construct the rotation matrix for the equation of equinoxes
    eqe = np.array(
        [[np.cos(eqeg), np.sin(eqeg), 0], [-np.sin(eqeg), np.cos(eqeg), 0], [0, 0, 1]]
    )

    return prec, nut, eqe


def eci2teme(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox (J2000) frame to the
    true equator, mean equinox (TEME) frame.

    References:
        Vallado: 2022, p. 232-234

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation

    Returns:
        tuple: (rteme, vteme, ateme)
            rteme (np.ndarray): TEME position vector in km
            vteme (np.ndarray): TEME velocity vector in km/s
            ateme (np.ndarray): TEME acceleration vector in km/s²
    """
    # Get the individual transformation matrices
    prec, nut, eqe = _get_teme_eci_transform(ttt, ddpsi, ddeps, iau80arr)

    # Combined transformation matrix
    tm = eqe @ nut.T @ prec.T

    # Transform vectors
    rteme = tm @ np.asarray(reci)
    vteme = tm @ np.asarray(veci)
    ateme = tm @ np.asarray(aeci)

    return rteme, vteme, ateme


def teme2eci(
    rteme: ArrayLike,
    vteme: ArrayLike,
    ateme: ArrayLike,
    ttt: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the true equator, mean equinox (TEME) frame to the ECI
    mean equator, mean equinox (J2000) frame.

    References:
        Vallado: 2022, p. 232-234

    Args:
        rteme (array_like): TEME position vector in km
        vteme (array_like): TEME velocity vector in km/s
        ateme (array_like): TEME acceleration vector in km/s²
                            (set to zeros if not available)
        ttt (float): Julian centuries of TT
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Get the individual transformation matrices
    prec, nut, eqe = _get_teme_eci_transform(ttt, ddpsi, ddeps, iau80arr)

    # Combined transformation matrix
    tm = prec @ nut @ eqe.T

    # Transform vectors
    reci = tm @ np.asarray(rteme)
    veci = tm @ np.asarray(vteme)
    aeci = tm @ np.asarray(ateme)

    return reci, veci, aeci


########################################################################################
# ECI <-> CIRS Frame Conversions
########################################################################################


def eci2cirs(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float = 0.0,
    ddy: float = 0.0,
    use_full_series: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox (J2000) frame to the
    Celestial Intermediate Reference System (CIRS) frame using the XYS approach.

    References:
        Vallado: 2022, p. 214

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians (default 0)
        ddy (float, optional): EOP correction for y in radians (default 0)
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data

    Returns:
        tuple: (rcirs, vcirs, acirs)
            rcirs (np.ndarray): CIRS position vector in km
            vcirs (np.ndarray): CIRS velocity vector in km/s
            acirs (np.ndarray): CIRS acceleration vector in km/s²
    """
    # Compute transformation matrix using XYS approach
    *_, pnb = iau.iau06xys(ttt, iau06arr, ddx, ddy, iau06xysarr, use_full_series)

    # Transform vectors
    rcirs = pnb.T @ np.asarray(reci)
    vcirs = pnb.T @ np.asarray(veci)
    acirs = pnb.T @ np.asarray(aeci)

    return rcirs, vcirs, acirs


def cirs2eci(
    rcirs: ArrayLike,
    vcirs: ArrayLike,
    acirs: ArrayLike,
    ttt: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float = 0.0,
    ddy: float = 0.0,
    use_full_series: bool = True,
):
    """Transforms a vector from the Celestial Intermediate Reference System (CIRS) frame
    to the ECI mean equator, mean equinox (J2000) frame using the XYS approach.

    References:
        Vallado: 2022, p. 214

    Args:
        rcirs (array_like): CIRS position vector in km
        vcirs (array_like): CIRS velocity vector in km/s
        acirs (array_like): CIRS acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians (default 0)
        ddy (float, optional): EOP correction for y in radians (default 0)
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Compute transformation matrix using XYS approach
    *_, pnb = iau.iau06xys(ttt, iau06arr, ddx, ddy, iau06xysarr, use_full_series)

    # Transform vectors
    reci = pnb @ np.asarray(rcirs)
    veci = pnb @ np.asarray(vcirs)
    aeci = pnb @ np.asarray(acirs)

    return reci, veci, aeci


########################################################################################
# ECI <-> TIRS Frame Conversions
########################################################################################


def eci2tirs(
    reci: ArrayLike,
    veci: ArrayLike,
    aeci: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float = 0.0,
    ddy: float = 0.0,
    use_full_series: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the ECI mean equator, mean equinox (J2000) frame to the
    Terrestrial Intermediate Reference System (TIRS) frame using the XYS approach.

    References:
        Vallado: 2022, p. 213

    Args:
        reci (array_like): ECI position vector in km
        veci (array_like): ECI velocity vector in km/s
        aeci (array_like): ECI acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians (default 0)
        ddy (float, optional): EOP correction for y in radians (default 0)
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data

    Returns:
        tuple: (rtirs, vtirs, atirs)
            rtirs (np.ndarray): TIRS position vector in km
            vtirs (np.ndarray): TIRS velocity vector in km/s
            atirs (np.ndarray): TIRS acceleration vector in km/s²
    """
    # Compute transformation matrices using XYS approach
    prec = np.eye(3)
    nut, st, _, omegaearth = compute_iau06_matrices(
        ttt, jdut1, lod, 0, 0, iau06arr, iau06xysarr, ddx, ddy, use_full_series
    )

    # Position transformation
    rtirs = st.T @ nut.T @ prec.T @ reci

    # Velocity transformation
    vtirs = st.T @ nut.T @ prec.T @ veci - np.cross(omegaearth, rtirs)

    # Acceleration transformation
    atirs = (
        st.T @ nut.T @ prec.T @ aeci
        - np.cross(omegaearth, np.cross(omegaearth, rtirs))
        - 2 * np.cross(omegaearth, vtirs)
    )

    return rtirs, vtirs, atirs


def tirs2eci(
    rtirs: ArrayLike,
    vtirs: ArrayLike,
    atirs: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    iau06arr: IAU06Array,
    iau06xysarr: IAU06xysArray,
    ddx: float = 0.0,
    ddy: float = 0.0,
    use_full_series: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Terrestrial Intermediate Reference System (TIRS)
    frame to the ECI mean equator, mean equinox (J2000) frame using the XYS approach.

    References:
        Vallado: 2022, p. 213

    Args:
        rtirs (array_like): TIRS position vector in km
        vtirs (array_like): TIRS velocity vector in km/s
        atirs (array_like): TIRS acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        iau06arr (IAU06Array): IAU 2006 data
        iau06xysarr (IAU06xysArray): IAU 2006 XYS data
        ddx (float, optional): EOP correction for x in radians (default 0)
        ddy (float, optional): EOP correction for y in radians (default 0)
        use_full_series (bool, optional): Use full series for IAU 2006 XYS data

    Returns:
        tuple: (reci, veci, aeci)
            reci (np.ndarray): ECI position vector in km
            veci (np.ndarray): ECI velocity vector in km/s
            aeci (np.ndarray): ECI acceleration vector in km/s²
    """
    # Compute transformation matrices using XYS approach
    prec = np.eye(3)
    nut, st, _, omegaearth = compute_iau06_matrices(
        ttt, jdut1, lod, 0, 0, iau06arr, iau06xysarr, ddx, ddy, use_full_series
    )

    # Position transformation
    reci = prec @ nut @ st @ np.asarray(rtirs)

    # Velocity transformation
    veci = (
        prec @ nut @ st @ (np.asarray(vtirs) + np.cross(omegaearth, np.asarray(rtirs)))
    )

    # Acceleration transformation
    aeci = (
        prec
        @ nut
        @ st
        @ (
            np.asarray(atirs)
            + np.cross(omegaearth, np.cross(omegaearth, np.asarray(rtirs)))
            + 2 * np.cross(omegaearth, np.asarray(vtirs))
        )
    )

    return reci, veci, aeci


########################################################################################
# ECEF <-> PEF Frame Conversions
########################################################################################


def ecef2pef(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    xp: float,
    yp: float,
    ttt: float,
    use_iau80: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed (ECEF) frame to the pseudo Earth-fixed
    (PEF) frame.

    References:
        Vallado: 2022, p. 224

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ttt (float): Julian centuries of TT
        use_iau80 (bool, optional): Use IAU 1980 data for polar motion (default True)

    Returns:
        tuple: (rpef, vpef, apef)
            rpef (np.ndarray): PEF position vector in km
            vpef (np.ndarray): PEF velocity vector in km/s
            apef (np.ndarray): PEF acceleration vector in km/s²
    """
    # Compute polar motion matrix
    pm = polarm(xp, yp, ttt, use_iau80)

    # Transform vectors
    rpef = pm @ np.asarray(recef)
    vpef = pm @ np.asarray(vecef)
    apef = pm @ np.asarray(aecef)

    return rpef, vpef, apef


def pef2ecef(
    rpef: ArrayLike,
    vpef: ArrayLike,
    apef: ArrayLike,
    xp: float,
    yp: float,
    ttt: float,
    use_iau80: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the pseudo Earth-fixed (PEF) frame to the Earth-fixed
    (ECEF) frame.

    References:
        Vallado: 2022, p. 224

    Args:
        rpef (array_like): PEF position vector in km
        vpef (array_like): PEF velocity vector in km/s
        apef (array_like): PEF acceleration vector in km/s²
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ttt (float): Julian centuries of TT
        use_iau80 (bool, optional): Use IAU 1980 data for polar motion (default True)

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    # Compute polar motion matrix
    pm = polarm(xp, yp, ttt, use_iau80)

    # Transform vectors
    recef = pm.T @ np.asarray(rpef)
    vecef = pm.T @ np.asarray(vpef)
    aecef = pm.T @ np.asarray(apef)

    return recef, vecef, aecef


########################################################################################
# ECEF <-> TOD Frame Conversions
########################################################################################


def ecef2tod(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed (ECEF) frame to the true-of-date (TOD)
    frame.

    References:
        Vallado: 2022, p. 225

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (rtod, vtod, atod)
            rtod (np.ndarray): TOD position vector in km
            vtod (np.ndarray): TOD velocity vector in km/s
            atod (np.ndarray): TOD acceleration vector in km/s²
    """
    # Find matrices that account for various orbit effects
    _, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Transform position
    rpef = pm @ recef
    rtod = st @ rpef

    # Transform velocity
    vpef = pm @ vecef
    vtod = st @ (vpef + np.cross(omegaearth, rpef))

    # Transform acceleration
    atod = st @ (
        pm @ aecef
        + np.cross(omegaearth, np.cross(omegaearth, rpef))
        + 2 * np.cross(omegaearth, vpef)
    )

    return rtod, vtod, atod


def tod2ecef(
    rtod: ArrayLike,
    vtod: ArrayLike,
    atod: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the true-of-date (TOD) frame to the Earth-fixed (ECEF)
    frame.

    References:
        Vallado: 2022, p. 225

    Args:
        rtod (array_like): TOD position vector in km
        vtod (array_like): TOD velocity vector in km/s
        atod (array_like): TOD acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    # Find matrices that account for various orbit effects
    _, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Transform position
    rpef = st.T @ rtod
    recef = pm.T @ rpef

    # Transform velocity
    vpef = st.T @ vtod - np.cross(omegaearth, rpef)
    vecef = pm.T @ vpef

    # Transform acceleration
    aecef = (
        pm.T @ (st.T @ atod)
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return recef, vecef, aecef


########################################################################################
# ECEF <-> MOD Frame Conversions
########################################################################################


def ecef2mod(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed (ECEF) frame to the mean-of-date (MOD)
    frame.

    References:
        Vallado: 2022, p. 227

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (rmod, vmod, amod)
            rmod (np.ndarray): MOD position vector in km
            vmod (np.ndarray): MOD velocity vector in km/s
            amod (np.ndarray): MOD acceleration vector in km/s²
    """
    # Find matrices that account for various orbit effects
    _, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Transform position
    rpef = pm @ recef
    rmod = nut @ st @ rpef

    # Transform velocity
    vpef = pm @ vecef
    vmod = nut @ st @ (vpef + np.cross(omegaearth, rpef))

    # Transform acceleration
    temp = np.cross(omegaearth, rpef)
    amod = (
        nut
        @ st
        @ (pm @ aecef + np.cross(omegaearth, temp) + 2 * np.cross(omegaearth, vpef))
    )

    return rmod, vmod, amod


def mod2ecef(
    rmod: ArrayLike,
    vmod: ArrayLike,
    amod: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    ddpsi: float,
    ddeps: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the mean-of-date (MOD) frame to the Earth-fixed (ECEF)
    frame.

    References:
        Vallado: 2022, p. 227

    Args:
        rmod (array_like): MOD position vector in km
        vmod (array_like): MOD velocity vector in km/s
        amod (array_like): MOD acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ddpsi (float): Delta psi correction to GCRF in radians
        ddeps (float): Delta epsilon correction to GCRF in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    # Find matrices that account for various orbit effects
    _, nut, st, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, ddpsi, ddeps, iau80arr, eqeterms=eqeterms
    )

    # Transform position
    rpef = st.T @ nut.T @ rmod
    recef = pm.T @ rpef

    # Transform velocity
    vpef = st.T @ nut.T @ vmod - np.cross(omegaearth, rpef)
    vecef = pm.T @ vpef

    # Transform acceleration
    aecef = (
        pm.T @ (st.T @ nut.T @ amod)
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return recef, vecef, aecef


########################################################################################
# ECEF <-> TEME Frame Conversions
########################################################################################


def get_teme_transform_matrices(
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    iau80arr: IAU80Array,
    eqeterms: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the matrices and Earth rotation vector for TEME transformations.

    Args:
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool): Add terms for ast calculation

    Returns:
        tuple: (st, pm, omegaearth)
            st (np.ndarray): Sidereal time matrix (PEF to TEME)
            pm (np.ndarray): Polar motion matrix (ECEF to PEF)
            omegaearth (np.ndarray): Earth's angular velocity vector
    """
    # Compute Greenwich Mean Sidereal Time (GMST)
    gmst = gstime(jdut1)

    # Compute omega from nutation theory
    omega = (
        125.04452222
        + (-6962890.539 * ttt + 7.455 * ttt**2 + 0.008 * ttt**3) / const.DEG2ARCSEC
    )
    omega = np.remainder(np.radians(omega), const.TWOPI)

    # Adjust GMST for geometric terms (kinematic after 1997)
    if jdut1 > 2450449.5 and eqeterms:
        gmstg = (
            gmst
            + 0.00264 * const.ARCSEC2RAD * np.sin(omega)
            + 0.000063 * const.ARCSEC2RAD * np.sin(2 * omega)
        )
    else:
        gmstg = gmst

    gmstg = np.remainder(gmstg, const.TWOPI)

    # Sidereal time matrix
    st = np.array(
        [
            [np.cos(gmstg), -np.sin(gmstg), 0],
            [np.sin(gmstg), np.cos(gmstg), 0],
            [0, 0, 1],
        ]
    )

    # Polar motion matrix and Earth's rotation vector
    *_, pm, omegaearth = calc_orbit_effects(
        ttt, jdut1, lod, xp, yp, 0, 0, iau80arr, eqeterms=eqeterms
    )

    return st, pm, omegaearth


def ecef2teme(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed (ECEF) frame to the true equator, mean
    mean equinox (TEME) frame.

    Results take into account the effects of sidereal time and polar motion.

    References:
        Vallado: 2022, p. 232-234

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (rteme, vteme, ateme)
            rteme (np.ndarray): TEME position vector in km
            vteme (np.ndarray): TEME velocity vector in km/s
            ateme (np.ndarray): TEME acceleration vector in km/s²
    """
    # Get common matrices and Earth rotation vector
    st, pm, omegaearth = get_teme_transform_matrices(
        ttt, jdut1, lod, xp, yp, iau80arr, eqeterms
    )

    # Transform position
    rpef = pm @ recef
    rteme = st @ rpef

    # Transform velocity
    vpef = pm @ vecef
    vteme = st @ (vpef + np.cross(omegaearth, rpef))

    # Transform acceleration
    ateme = st @ (
        pm @ aecef
        + np.cross(omegaearth, np.cross(omegaearth, rpef))
        + 2 * np.cross(omegaearth, vpef)
    )

    return rteme, vteme, ateme


def teme2ecef(
    rteme: ArrayLike,
    vteme: ArrayLike,
    ateme: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
    iau80arr: IAU80Array,
    eqeterms: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the true equator, mean equinox (TEME) frame to the
    Earth-fixed (ECEF) frame.

    Results take into account the effects of sidereal time and polar motion.

    References:
        Vallado: 2022, p. 232-234

    Args:
        rteme (array_like): TEME position vector in km
        vteme (array_like): TEME velocity vector in km/s
        ateme (array_like): TEME acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        iau80arr (IAU80Array): IAU 1980 data for nutation
        eqeterms (bool, optional): Add terms for ast calculation (default True)

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    # Get common matrices and Earth rotation vector
    st, pm, omegaearth = get_teme_transform_matrices(
        ttt, jdut1, lod, xp, yp, iau80arr, eqeterms
    )

    # Transform position
    rpef = st.T @ rteme
    recef = pm.T @ rpef

    # Transform velocity
    vpef = st.T @ vteme - np.cross(omegaearth, rpef)
    vecef = pm.T @ vpef

    # Transform acceleration
    aecef = pm.T @ (
        st.T @ ateme
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return recef, vecef, aecef


########################################################################################
# ECEF <-> CIRS Frame Conversions
########################################################################################


def ecef2cirs(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed (ECEF) frame to the Celestial
    Intermediate Reference System (CIRS) frame.

    References:
        Vallado: 2022, p. 214

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians

    Returns:
        tuple: (rcirs, vcirs, acirs)
            rcirs (np.ndarray): CIRS position vector in km
            vcirs (np.ndarray): CIRS velocity vector in km/s
            acirs (np.ndarray): CIRS acceleration vector in km/s²
    """
    # Compute transformation matrices
    st, _ = sidereal(jdut1, 0, 0, 0, lod, use_iau80=False)
    omegaearth = calc_omegaearth(lod)
    pm = polarm(xp, yp, ttt, use_iau80=False)

    # Transform position
    rpef = pm @ np.asarray(recef)
    rcirs = st @ rpef

    # Transform velocity
    vpef = pm @ np.asarray(vecef)
    vcirs = st @ (vpef + np.cross(omegaearth, rpef))

    # Transform acceleration
    acirs = st @ (
        pm @ np.asarray(aecef)
        + np.cross(omegaearth, np.cross(omegaearth, rpef))
        + 2 * np.cross(omegaearth, vpef)
    )

    return rcirs, vcirs, acirs


def cirs2ecef(
    rcirs: ArrayLike,
    vcirs: ArrayLike,
    acirs: ArrayLike,
    ttt: float,
    jdut1: float,
    lod: float,
    xp: float,
    yp: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Celestial Intermediate Reference System (CIRS) frame
    to the Earth-fixed (ECEF) frame.

    References:
        Vallado: 2022, p. 214

    Args:
        rcirs (array_like): CIRS position vector in km
        vcirs (array_like): CIRS velocity vector in km/s
        acirs (array_like): CIRS acceleration vector in km/s²
        ttt (float): Julian centuries of TT
        jdut1 (float): Julian date of UT1 (days from 4713 BC)
        lod (float): Excess length of day in seconds
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    # Compute transformation matrices
    st, _ = sidereal(jdut1, 0, 0, 0, lod, use_iau80=False)
    omegaearth = calc_omegaearth(lod)
    pm = polarm(xp, yp, ttt, use_iau80=False)

    # Transform position
    rpef = st.T @ rcirs
    recef = pm.T @ rpef

    # Transform velocity
    vpef = st.T @ vcirs - np.cross(omegaearth, rpef)
    vecef = pm.T @ vpef

    # Transform acceleration
    aecef = (
        pm.T @ (st.T @ acirs)
        - np.cross(omegaearth, np.cross(omegaearth, rpef))
        - 2 * np.cross(omegaearth, vpef)
    )

    return recef, vecef, aecef


########################################################################################
# ECEF <-> TIRS Frame Conversions
########################################################################################


def ecef2tirs(
    recef: ArrayLike,
    vecef: ArrayLike,
    aecef: ArrayLike,
    xp: float,
    yp: float,
    ttt: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Earth-fixed (ECEF) frame to the Terrestrial
    Intermediate Reference System (TIRS) frame.

    References:
        Vallado: 2022, p. 213

    Args:
        recef (array_like): ECEF position vector in km
        vecef (array_like): ECEF velocity vector in km/s
        aecef (array_like): ECEF acceleration vector in km/s²
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ttt (float): Julian centuries of TT

    Returns:
        tuple: (rtirs, vtirs, atirs)
            rtirs (np.ndarray): TIRS position vector in km
            vtirs (np.ndarray): TIRS velocity vector in km/s
            atirs (np.ndarray): TIRS acceleration vector in km/s²
    """
    return ecef2pef(recef, vecef, aecef, xp, yp, ttt, use_iau80=False)


def tirs2ecef(
    rtirs: ArrayLike,
    vtirs: ArrayLike,
    atirs: ArrayLike,
    xp: float,
    yp: float,
    ttt: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transforms a vector from the Terrestrial Intermediate Reference System (TIRS)
    frame to the Earth-fixed (ECEF) frame.

    References:
        Vallado: 2022, p. 213

    Args:
        rtirs (array_like): TIRS position vector in km
        vtirs (array_like): TIRS velocity vector in km/s
        atirs (array_like): TIRS acceleration vector in km/s²
        xp (float): Polar motion coefficient in radians
        yp (float): Polar motion coefficient in radians
        ttt (float): Julian centuries of TT

    Returns:
        tuple: (recef, vecef, aecef)
            recef (np.ndarray): ECEF position vector in km
            vecef (np.ndarray): ECEF velocity vector in km/s
            aecef (np.ndarray): ECEF acceleration vector in km/s²
    """
    return pef2ecef(rtirs, vtirs, atirs, xp, yp, ttt, use_iau80=False)
