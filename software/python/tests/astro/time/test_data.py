import os
import pytest

import numpy as np
import scipy

from src.valladopy.astro.time.data import iau80in, iau06in
from ...conftest import custom_allclose


def load_matlab_data(file_path: str, keys: list) -> dict:
    """Load MATLAB .mat file data

    Args:
        file_path (str): Path to the .mat file
        keys (list): list of variable names to grab

    Returns:
        dict [str, np.ndarray]: Dictionary of the input keys and their associated matlab
                                data as numpy arrays
    """
    # Load the .m data file
    data = scipy.io.loadmat(file_path)

    # Grab data for each key
    return {key: data[key] for key in keys}


@pytest.fixture()
def iau80_mat_data():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "data", "iau80in_data.mat")
    return load_matlab_data(file_path, keys=["iar80", "rar80"])


@pytest.fixture()
def iau06_mat_data():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "data", "iau06in_data.mat")
    keys = [
        "axs0",
        "a0xi",
        "ays0",
        "a0yi",
        "ass0",
        "a0si",
        "apn",
        "apni",
        "appl",
        "appli",
        "agst",
        "agsti",
    ]
    return load_matlab_data(file_path, keys=keys)


def test_iau80in(iau80_mat_data):
    matlab_iar80 = iau80_mat_data["iar80"]
    matlab_rar80 = iau80_mat_data["rar80"]

    # Load Python data using iau80in
    iau80arr = iau80in()

    # Check that they are the same
    assert np.array_equal(iau80arr.iar80, matlab_iar80)
    assert custom_allclose(iau80arr.rar80, matlab_rar80)


def test_iau06in(iau06_mat_data):
    # Load MATLAB data
    matlab_data = iau06_mat_data

    # Load Python data using iau06in
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()

    # Check that they are the same
    assert custom_allclose(axs0, matlab_data["axs0"])
    assert np.array_equal(a0xi, matlab_data["a0xi"])
    assert custom_allclose(ays0, matlab_data["ays0"])
    assert np.array_equal(a0yi, matlab_data["a0yi"])
    assert custom_allclose(ass0, matlab_data["ass0"])
    assert np.array_equal(a0si, matlab_data["a0si"])
    assert custom_allclose(apn, matlab_data["apn"])
    assert np.array_equal(apni, matlab_data["apni"])
    assert custom_allclose(appl, matlab_data["appl"])
    assert np.array_equal(appli, matlab_data["appli"])
    assert custom_allclose(agst, matlab_data["agst"])
    assert np.array_equal(agsti, matlab_data["agsti"])
