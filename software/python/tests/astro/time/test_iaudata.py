import os
import pytest

import numpy as np
import scipy

from src.valladopy.astro.time.iaudata import iau80in
from ...conftest import custom_allclose


def load_matlab_data(file_path, keys):
    """Load MATLAB .mat file data

    Args:
        file_path (str): Path to the .mat file
        keys (list): list of variable names to grab

    Returns:
        (dict [str, np.ndarray]): dictionary of the input keys and their
                                  associated matlab data as numpy arrays
    """
    # Load the .m data file
    data = scipy.io.loadmat(file_path)

    # Grab data for each key
    return {key: data[key] for key in keys}


@pytest.fixture()
def iau80_mat_data():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "data", "nutation_data.mat")
    return load_matlab_data(file_path, keys=["iar80", "rar80"])


def test_iau80in(iau80_mat_data):
    matlab_iar80 = iau80_mat_data["iar80"]
    matlab_rar80 = iau80_mat_data["rar80"]

    # Load Python data using iau80in
    iar80, rar80 = iau80in()

    # Check that they are the same
    assert np.array_equal(iar80, matlab_iar80)
    assert custom_allclose(rar80, matlab_rar80, rtol=1e-12)
