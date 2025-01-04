# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 27 May 2002
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np

from ...constants import ARCSEC2RAD


# Default data directory
ROOT_DIR = Path(__file__).resolve().parents[6]
DATA_DIR = ROOT_DIR / "datalib"

# Conversion factors
CONVRTU = 1e-6 * ARCSEC2RAD  # microarcseconds to radians
CONVRTM = 1e-3 * ARCSEC2RAD  # milliarcseconds to radians


@dataclass
class IAU80Array:
    # fmt: off
    """Data class for IAU 1980 nutation data."""
    iar80: np.ndarray = None
    rar80: np.ndarray = None


@dataclass
class IAU06Array:
    ax0: np.ndarray = None
    ax0i: np.ndarray = None
    ay0: np.ndarray = None
    ay0i: np.ndarray = None
    as0: np.ndarray = None
    as0i: np.ndarray = None
    agst: np.ndarray = None
    agsti: np.ndarray = None
    apn0: np.ndarray = None
    apn0i: np.ndarray = None
    apl0: np.ndarray = None
    apl0i: np.ndarray = None
    aapn0: np.ndarray = None
    aapn0i: np.ndarray = None


def iau80in(data_dir: str = DATA_DIR) -> IAU80Array:
    """Initializes the nutation matrices needed for reduction calculations.

    References:
        Vallado, 2022, Section 3.7.1

    Args:
        data_dir (str, optional): Directory containing the nutation data file
                                  "nut80.dat" (default: DATA_DIR)

    Returns:
        IAU80Array: Data object containing the nutation matrices
            iar80 (np.ndarray): Integers for FK5 1980
            rar80 (np.ndarray): Reals for FK5 1980 in radians
    """
    # Load the nutation data and initialize data object
    nut80 = np.loadtxt(os.path.join(data_dir, "nut80.dat"))
    iau80arr = IAU80Array()

    # Split into integer and real parts
    iau80arr.iar80 = nut80[:, :5].astype(int)
    iau80arr.rar80 = nut80[:, 5:9]

    # Convert from 0.0001 arcseconds to radians
    convrt = 1e-4 * ARCSEC2RAD
    iau80arr.rar80 *= convrt

    return iau80arr


def iau06in() -> (
    Tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]
):
    """Initializes the matrices needed for IAU 2006 reduction calculations.

    References:
        Vallado, 2022, Section 3.7.1

    Returns:
        tuple: (axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti)
            axs0 (np.ndarray): Real coefficients for X in radians
            a0xi (np.ndarray): Integer coefficients for X
            ays0 (np.ndarray): Real coefficients for Y in radians
            a0yi (np.ndarray): Integer coefficients for Y
            ass0 (np.ndarray): Real coefficients for S in radians
            a0si (np.ndarray): Integer coefficients for S
            apn (np.ndarray): Real coefficients for nutation in radians
            apni (np.ndarray): Integer coefficients for nutation
            appl (np.ndarray): Real coefficients for planetary nutation in radians
            appli (np.ndarray): Integer coefficients for planetary nutation
            agst (np.ndarray): Real coefficients for GST in radians
            agsti (np.ndarray): Integer coefficients for GST

    Notes:
        Data files are from the IAU 2006 precession-nutation model:
            - iau06xtab5.2.a.dat (file for X coefficients)
            - iau06ytab5.2.b.dat (file for Y coefficients)
            - iau06stab5.2.d.dat (file for S coefficients)
            - iau03n.dat (file for nutation coefficients)
            - iau03pl.dat (file for planetary nutation coefficients)
            - iau06gsttab5.2.e.dat (file for GST coefficients)

    TODO: Update with latest MATLAB updates
    """
    # Conversion factors
    convrtu = 1e-6 * ARCSEC2RAD  # microarcseconds to radians
    convrtm = 1e-3 * ARCSEC2RAD  # milliarcseconds to radians

    def load_data(
        filename, columns_real, columns_int, conv_factor, convert_exclude_last=False
    ):
        """Helper function to load and process data."""
        filepath = os.path.join(DATA_DIR, filename)
        data = np.loadtxt(filepath)
        reals = data[:, columns_real]
        if convert_exclude_last:
            reals[:, :-1] *= conv_factor  # convert all except the last column
        else:
            reals *= conv_factor  # convert all
        integers = data[:, columns_int].astype(int)
        return reals, integers

    # Load data
    axs0, a0xi = load_data(
        "iau06xtab5.2.a.dat",
        columns_real=[1, 2],
        columns_int=range(3, 17),
        conv_factor=convrtu,
    )
    ays0, a0yi = load_data(
        "iau06ytab5.2.b.dat",
        columns_real=[1, 2],
        columns_int=range(3, 17),
        conv_factor=convrtu,
    )
    ass0, a0si = load_data(
        "iau06stab5.2.d.dat",
        columns_real=[1, 2],
        columns_int=range(3, 17),
        conv_factor=convrtu,
    )
    apn, apni = load_data(
        "iau03n.dat",
        columns_real=range(6, 14),
        columns_int=range(0, 5),
        conv_factor=convrtm,
    )
    appl, appli = load_data(
        "iau03pl.dat",
        columns_real=range(16, 21),  # include column 21 (extra)
        columns_int=range(1, 15),
        conv_factor=convrtm,
        convert_exclude_last=True,
    )
    agst, agsti = load_data(
        "iau06gsttab5.2.e.dat",
        columns_real=[1, 2],
        columns_int=range(3, 17),
        conv_factor=convrtu,
    )

    return axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti


def iau06in2(data_dir: str = DATA_DIR) -> IAU06Array:
    """Initializes the matrices needed for IAU 2006 reduction calculations.

    Args:
        data_dir: Path to the directory containing the IAU 2006 data files.

    Returns:
        IAU06Array: Dataclass containing all coefficients for IAU 2006 calculations.
    """

    def load_data(filename, cols_real, cols_int, conv_factor):
        """Helper function to load and process data."""
        filepath = os.path.join(data_dir, filename)
        data = np.loadtxt(filepath, skiprows=2)
        reals = data[:, cols_real[0] : cols_real[1]]
        if conv_factor:
            reals *= conv_factor
        integers = data[:, cols_int[0] : cols_int[1]].astype(int)
        return reals, integers

    # Load all data files
    ax0, ax0i = load_data(
        "iau06xtab5.2.a.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    ay0, ay0i = load_data(
        "iau06ytab5.2.b.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    as0, as0i = load_data(
        "iau06stab5.2.d.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    agst, agsti = load_data(
        "iau06gsttab5.2.e.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    aapn0, aapn0i = load_data(
        "iau06ansofa.dat", cols_real=(5, 11), cols_int=(0, 5), conv_factor=CONVRTM
    )
    apl0, apl0i = load_data(
        "iau06anplsofa.dat", cols_real=(14, 18), cols_int=(0, 14), conv_factor=CONVRTM
    )
    apn0, apn0i = load_data(
        "iau06nlontab5.3.a.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTM
    )

    return IAU06Array(
        ax0,
        ax0i,
        ay0,
        ay0i,
        as0,
        as0i,
        agst,
        agsti,
        apn0,
        apn0i,
        apl0,
        apl0i,
        aapn0,
        aapn0i,
    )
