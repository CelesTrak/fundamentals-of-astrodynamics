import os
import pytest

import numpy as np
import scipy

from src.valladopy.astro.time.data import iau80in, iau06in, iau06in2
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


def load_matlab_data2(file_path: str, keys: list) -> dict:
    """Load MATLAB .mat file data and handle structures.

    Args:
        file_path (str): Path to the .mat file
        keys (list): List of variable names to grab

    Returns:
        dict [str, np.ndarray or dict]: Dictionary of the input keys and their
                                        associated MATLAB data as numpy arrays or dicts
    """

    def unpack_structure(struct):
        """Recursively unpack MATLAB structure arrays into dictionaries."""
        if isinstance(struct, np.ndarray) and struct.dtype.names:
            return {name: struct[name] for name in struct.dtype.names}
        return struct

    # Load the .mat file
    data = scipy.io.loadmat(file_path, struct_as_record=False, squeeze_me=True)

    # Process keys and unpack structures
    result = {}
    for key in keys:
        if key in data:
            value = data[key]
            if isinstance(value, np.ndarray) and value.dtype.names:
                # If the key is a structure, unpack its fields
                result[key] = unpack_structure(value)
            else:
                result[key] = value
    return result


@pytest.fixture()
def current_dir():
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture()
def iau80_mat_data(current_dir):
    file_path = os.path.join(current_dir, "data", "iau80in_data.mat")
    return load_matlab_data2(file_path, keys=["iau80arr"])


@pytest.fixture()
def iau06_mat_data(current_dir):
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


@pytest.fixture()
def iau06_mat_data2(current_dir):
    file_path = os.path.join(current_dir, "data", "iau06in2_data.mat")
    return load_matlab_data2(file_path, keys=["iau06arr"])


def test_iau80in(iau80_mat_data):
    matlab_data = iau80_mat_data["iau80arr"]

    # Load Python data using iau80in
    iau80arr = iau80in()

    # Check that they are the same
    assert np.array_equal(iau80arr.iar80, matlab_data.iar80)
    assert custom_allclose(iau80arr.rar80, matlab_data.rar80)


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


def test_iau06in2(iau06_mat_data2):
    # Load MATLAB data
    matlab_data = iau06_mat_data2["iau06arr"]

    # Load Python data using iau06in
    iau06arr = iau06in2()

    # Check that they are the same
    assert custom_allclose(iau06arr.ax0, matlab_data.ax0)
    assert np.array_equal(iau06arr.ax0i, matlab_data.a0xi)
    assert custom_allclose(iau06arr.ay0, matlab_data.ay0)
    assert np.array_equal(iau06arr.ay0i, matlab_data.a0yi)
    assert custom_allclose(iau06arr.as0, matlab_data.as0)
    assert np.array_equal(iau06arr.as0i, matlab_data.a0si)
    assert custom_allclose(iau06arr.agst, matlab_data.agst)
    assert np.array_equal(iau06arr.agsti, matlab_data.agsti[:, :14])
    assert custom_allclose(iau06arr.apn0, matlab_data.apn0)
    assert np.array_equal(iau06arr.apn0i, matlab_data.apn0i)
    assert custom_allclose(iau06arr.apl0, matlab_data.apl0)
    assert np.array_equal(iau06arr.apl0i, matlab_data.apl0i)
    assert custom_allclose(iau06arr.aapn0[:, :5], matlab_data.aapn0[:, :5])
    assert np.array_equal(iau06arr.aapn0i, matlab_data.aapn0i)
