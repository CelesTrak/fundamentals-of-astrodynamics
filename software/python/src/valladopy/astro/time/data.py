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

# Old data directory
DATA_DIR_OLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@dataclass
class IAU80Array:
    # fmt: off
    """Data class for IAU 1980 nutation data."""
    iar80: np.ndarray = None
    rar80: np.ndarray = None


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
        filepath = os.path.join(DATA_DIR_OLD, filename)
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
