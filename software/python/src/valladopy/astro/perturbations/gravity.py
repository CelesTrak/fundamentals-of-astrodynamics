# --------------------------------------------------------------------------------------
# Author: David Vallado
# Date: 10 Oct 2019
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from dataclasses import dataclass

import numpy as np


@dataclass
class GravityFieldData:
    c: np.ndarray = None
    s: np.ndarray = None
    c_unc: np.ndarray = None
    s_unc: np.ndarray = None
    normalized: bool = False


def read_gravity_field(filename: str, normalized: bool) -> GravityFieldData:
    """Reads and stores gravity field coefficients.

    References:
        Vallado: 2022, p. 550-551

    Args:
        filename (str): The filename of the gravity field data
        normalized (bool): True if the gravity field data is normalized

    Returns:
        GravityFieldData: A dataclass containing gravity field data:
            - c (np.ndarray): Cosine coefficients
            - s (np.ndarray): Sine coefficients
            - normalized (bool): True if the gravity field data is normalized
    """
    # Load gravity field data
    file_data = np.loadtxt(filename)

    # Get the maximum degree of the gravity field
    max_degree = int(np.max(file_data[:, 0]))
    size = max_degree + 1

    # Initialize gravity field data
    gravarr = GravityFieldData(
        c=np.zeros((size, size)), s=np.zeros((size, size)), normalized=normalized
    )

    # Check if uncertainties are included in the data (columns 5 and 6)
    has_uncertainty = file_data.shape[1] >= 6
    if has_uncertainty:
        gravarr.c_unc = np.zeros((size, size))
        gravarr.s_unc = np.zeros((size, size))

    # Store gravity field coefficients
    for row in file_data:
        n, m = int(row[0]), int(row[1])
        c_value, s_value = row[2], row[3]
        gravarr.c[n, m] = c_value
        gravarr.s[n, m] = s_value

        if has_uncertainty:
            gravarr.c_unc[n, m] = row[4]
            gravarr.s_unc[n, m] = row[5]

    return gravarr
