# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 21 June 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------


import numpy as np


########################################################################################
# B1950 to J2000 conversion matrices
########################################################################################


def fk4() -> np.ndarray:
    """Returns the B1950 to J2000 conversion matrix (book approach).

    References:
        Vallado, 2022, p. 235

    Returns:
        np.ndarray: B1950 to J2000 conversion matrix

    Notes:
        This process is not exact. There are different secular rates for each system,
        and there are differences in the central location. the matrices are multiplied
        directly for speed.
    """
    return np.array(
        [
            [0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
            [0.0111814832391717, 0.9999374848933135, -0.0000271625947142],
            [0.0048590037723143, -0.0000271702937440, 0.9999881946043742],
        ]
    )


def fk4_stk() -> np.ndarray:
    """Returns the B1950 to J2000 conversion matrix (STK approach).

    Returns:
        np.ndarray: B1950 to J2000 conversion matrix

    Notes:
        This way is formed by multiplying the matrices on pages 173 and 174 and adding
        in the correction to equinox given on page 168 of the supplement to the
        astronomical almanac.
    """
    return np.array(
        [
            [0.999925678612394, -0.011181874556714, -0.0048582848126],
            [0.011181874524964, 0.99993748051788, -0.000027169816135],
            [0.004858284884778, -0.000027156932874, 0.999988198095508],
        ]
    )


def fk4_exp_supp_6x6() -> np.ndarray:
    """Returns the B1950 to J2000 conversion matrix (almanac supplement, 6x6).

    Returns:
        np.ndarray: B1950 to J2000 conversion matrix

    Notes:
        Explanatory Supplement to the Astronomical Almanac, page 185
    """
    return np.array(
        [
            [
                0.9999256782,
                -0.0111820611,
                -0.0048579477,
                0.00000242395018,
                -0.00000002710663,
                -0.00000001177656,
            ],
            [
                0.011182061,
                0.9999374784,
                -0.0000271765,
                0.00000002710663,
                0.00000242397878,
                -0.00000000006587,
            ],
            [
                0.0048579479,
                -0.0000271474,
                0.9999881997,
                0.00000001177656,
                -0.00000000006582,
                0.00000242410173,
            ],
            [-0.000551, -0.238565, 0.435739, 0.99994704, -0.01118251, -0.00485767],
            [0.238514, -0.002667, -0.008541, 0.01118251, 0.99995883, -0.00002718],
            [-0.435623, 0.012254, 0.002117, 0.00485767, -0.00002714, 1.00000956],
        ]
    )


########################################################################################
# J2000 to B1950 conversion matrices
# TODO: Consolidate with the above (or remove) - just need to return the transpose
########################################################################################


def fk4i() -> np.ndarray:
    """Returns the J2000 to B1950 conversion matrix (book approach).

    References:
        Vallado, 2022, p. 235

    Returns:
        np.ndarray: J2000 to B1950 conversion matrix
    """
    return np.array(
        [
            [0.9999256794956877, 0.0111814832391717, 0.0048590037723143],
            [-0.0111814832204662, 0.9999374848933135, -0.000027170293744],
            [-0.0048590038153592, -0.0000271625947142, 0.9999881946043742],
        ]
    )


def fk4i_stk() -> np.ndarray:
    """Returns the J2000 to B1950 conversion matrix (STK approach)."""
    return np.array(
        [
            [0.999925678612394, 0.011181874524964, 0.004858284884778],
            [-0.011181874556714, 0.999937480517880, -0.000027156932874],
            [-0.004858284812600, -0.000027169816135, 0.999988198095508],
        ]
    )


def fk4i_exp_supp_3x3() -> np.ndarray:
    """Returns the J2000 to B1950 FK4 conversion matrix
    (almanac supplement, 3x3).

    Returns:
        np.ndarray: J2000 to B1950 conversion matrix
    """
    return np.array(
        [
            [0.9999256795, 0.0111814828, 0.0048590039],
            [-0.0111814828, 0.9999374849, -0.0000271771],
            [-0.0048590040, -0.0000271557, 0.9999881946],
        ]
    )


def fk4i_exp_supp_6x6() -> np.ndarray:
    """Returns the J2000 to B1950 FK4 conversion matrix (almanac supplement, 6x6).

    Returns:
        np.ndarray: J2000 to B1950 conversion matrix
    """
    return np.array(
        [
            [
                0.9999256795,
                0.0111814828,
                0.0048590039,
                -0.00000242389840,
                -0.00000002710544,
                -0.00000001177742,
            ],
            [
                -0.0111814828,
                0.9999374849,
                -0.0000271771,
                0.00000002710544,
                -0.00000242392702,
                0.00000000006585,
            ],
            [
                -0.0048590040,
                -0.0000271557,
                0.9999881946,
                0.00000001177742,
                0.00000000006585,
                -0.00000242404995,
            ],
            [-0.000551, 0.238509, -0.435614, 0.99990432, 0.01118145, 0.00485852],
            [-0.238560, -0.002667, 0.012254, -0.01118145, 0.99991613, -0.00002717],
            [0.435730, -0.008541, 0.002117, -0.00485852, -0.00002716, 0.99996684],
        ]
    )
