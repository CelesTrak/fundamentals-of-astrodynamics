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

import numpy as np

from ...constants import ARCSEC2RAD, JD_TO_MJD_OFFSET


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
class IAU06pnOldArray:
    # fmt: off
    """Data class for IAU 2006 precession-nutation data (older version)."""
    apn: np.ndarray = None
    apni: np.ndarray = None
    appl: np.ndarray = None
    appli: np.ndarray = None


@dataclass
class IAU06Array:
    # fmt: off
    """Data class for IAU 2006 precession-nutation data."""
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


@dataclass
class IAU06xysArray:
    # fmt: off
    """Data class for IAU 2006 XYS data."""
    jd: np.ndarray = None  # Julian Date
    jdf: np.ndarray = None  # fractional Julian Date
    x: np.ndarray = None  # x coordinate in radians
    y: np.ndarray = None  # y coordinate in radians
    s: np.ndarray = None  # s coordinate in radians
    mjd_tt: np.ndarray = None  # modified Julian Date (computed)


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


def iau06in_pnold(data_dir: str = DATA_DIR) -> IAU06pnOldArray:
    """Initializes the nutation matrices needed for IAU 2006 reduction calculations.

    This is an older version of the function that loads the data from the IAU 2006,
    and we're only returning the coefficients for nutation.

    References:
        Vallado, 2022, Section 3.7.1

    Returns:
        IAU06pnOldArray: Data object containing the nutation coefficients
            apn (np.ndarray): Nutation coefficients for IAU 2006
            apni (np.ndarray): Integers for nutation coefficients
            appl (np.ndarray): Planetary nutation coefficients for IAU 2006
            appli (np.ndarray): Integers for planetary nutation coefficients

    Notes:
        Data files are from the IAU 2006 precession-nutation model:
            - iau03n.dat (file for nutation coefficients)
            - iau03pl.dat (file for planetary nutation coefficients)

    TODO:
        - Deprecate when updates to iau06pna/b are complete.
    """

    def load_data(
        filename, columns_real, columns_int, conv_factor, convert_exclude_last=False
    ):
        """Helper function to load and process data."""
        filepath = os.path.join(data_dir, filename)
        data = np.loadtxt(filepath)
        reals = data[:, columns_real]
        if convert_exclude_last:
            reals[:, :-1] *= conv_factor  # convert all except the last column
        else:
            reals *= conv_factor  # convert all
        integers = data[:, columns_int].astype(int)
        return reals, integers

    # Initialize data object
    iau06arr = IAU06pnOldArray()

    # Load all data files
    iau06arr.apn, iau06arr.apni = load_data(
        "iau03n.dat",
        columns_real=range(6, 14),
        columns_int=range(0, 5),
        conv_factor=CONVRTM,
    )
    iau06arr.appl, iau06arr.appli = load_data(
        "iau03pl.dat",
        columns_real=range(16, 21),  # include column 21 (extra)
        columns_int=range(1, 15),
        conv_factor=CONVRTM,
        convert_exclude_last=True,
    )

    return iau06arr


def iau06in(data_dir: str = DATA_DIR) -> IAU06Array:
    """Initializes the matrices needed for IAU 2006 reduction calculations.

    Args:
        data_dir (str, optional): Path to the directory containing the IAU 2006 data
                                  files (default: DATA_DIR)

    Returns:
        IAU06Array: Dataclass containing all coefficients for IAU 2006 calculations

    Notes:
        Data files are from the IAU 2006 precession-nutation model:
            - iau06xtab5.2.a.dat (file for X coefficients)
            - iau06ytab5.2.b.dat (file for Y coefficients)
            - iau06stab5.2.d.dat (file for S coefficients)
            - iau06gsttab5.2.e.dat (file for GST coefficients)
            - iau06ansofa.dat (file for SOFA luni-solar nutation coefficients)
            - iau06anplsofa.dat (file for SOFA planetary nutation coefficients)
            - iau06nlontab5.3.a.dat (file for nutation in obliquity coefficients
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

    # Initialize data object
    iau06arr = IAU06Array()

    # Load all data files
    iau06arr.ax0, iau06arr.ax0i = load_data(
        "iau06xtab5.2.a.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    iau06arr.ay0, iau06arr.ay0i = load_data(
        "iau06ytab5.2.b.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    iau06arr.as0, iau06arr.as0i = load_data(
        "iau06stab5.2.d.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    iau06arr.agst, iau06arr.agsti = load_data(
        "iau06gsttab5.2.e.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTU
    )
    iau06arr.aapn0, iau06arr.aapn0i = load_data(
        "iau06ansofa.dat", cols_real=(5, 11), cols_int=(0, 5), conv_factor=CONVRTM
    )
    iau06arr.apl0, iau06arr.apl0i = load_data(
        "iau06anplsofa.dat", cols_real=(14, 18), cols_int=(0, 14), conv_factor=CONVRTM
    )
    iau06arr.apn0, iau06arr.apn0i = load_data(
        "iau06nlontab5.3.a.dat", cols_real=(1, 3), cols_int=(3, 17), conv_factor=CONVRTM
    )

    return iau06arr


def readxys(data_dir: str = DATA_DIR) -> IAU06xysArray:
    """Initializes the XYS IAU2006 data from the input file into a dataclass.

    Args:
        data_dir (str, optional): Directory containing the XYS data file "xysdata.dat"
                                  (default: DATA_DIR)

    Returns:
        IAU06xysArray: Dataclass containing the XYS data
    """
    # Read the file as an np.ndarray
    data = np.loadtxt(os.path.join(data_dir, "xysdata.dat"))

    # Create instance and extract individual columns
    iau06xysarr = IAU06xysArray()
    iau06xysarr.jd = data[:, 0]
    iau06xysarr.jdf = data[:, 1]
    iau06xysarr.x = data[:, 2]
    iau06xysarr.y = data[:, 3]
    iau06xysarr.s = data[:, 4]

    # Compute the derived column 'mjd_tt'
    iau06xysarr.mjd_tt = iau06xysarr.jd + iau06xysarr.jdf - JD_TO_MJD_OFFSET

    return iau06xysarr
