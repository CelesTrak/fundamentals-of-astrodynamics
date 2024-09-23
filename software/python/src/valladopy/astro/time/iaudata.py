import numpy as np
import os

from ...constants import ARCSEC2RAD


def iau80in():
    """Initializes the nutation matrices needed for reduction calculations

    Returns:
        tuple: A tuple containing:
            - iar80 (np.ndarray): Integers for FK5 1980
            - rar80 (np.ndarray): Reals for FK5 1980 (in radians)
    """
    # Define the path to the nut80.dat file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "data", "nut80.dat")

    # Load the nutation data
    nut80 = np.loadtxt(file_path)

    # Split into integer and real parts
    iar80 = nut80[:, :5].astype(int)
    rar80 = nut80[:, 5:9]

    # Convert from 0.0001 arcseconds to radians
    convrt = 0.0001 * ARCSEC2RAD
    rar80 *= convrt

    return iar80, rar80
