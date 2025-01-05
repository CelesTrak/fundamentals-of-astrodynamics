import os
import pytest

import numpy as np
import scipy

from src.valladopy.astro.time.data import iau80in
from ...conftest import custom_allclose


def load_matlab_data(file_path: str, keys: list) -> dict:
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
    return load_matlab_data(file_path, keys=["iau80arr"])


@pytest.fixture()
def iau06_pnold_mat_data(current_dir):
    file_path = os.path.join(current_dir, "data", "iau06in_pnold_data.mat")
    return load_matlab_data(file_path, keys=["apn", "apni", "appl", "appli"])


@pytest.fixture()
def iau06_mat_data(current_dir):
    file_path = os.path.join(current_dir, "data", "iau06in_data.mat")
    return load_matlab_data(file_path, keys=["iau06arr"])


def test_iau80in(iau80_mat_data):
    matlab_data = iau80_mat_data["iau80arr"]

    # Load Python data using iau80in
    iau80arr = iau80in()

    # Check that they are the same
    assert np.array_equal(iau80arr.iar80, matlab_data.iar80)
    assert custom_allclose(iau80arr.rar80, matlab_data.rar80)


def test_iau06in_pnold(iau06data_old, iau06_pnold_mat_data):
    # Check that the data is the same
    assert custom_allclose(iau06data_old.apn, iau06_pnold_mat_data["apn"])
    assert np.array_equal(iau06data_old.apni, iau06_pnold_mat_data["apni"])
    assert custom_allclose(iau06data_old.appl, iau06_pnold_mat_data["appl"])
    assert np.array_equal(iau06data_old.appli, iau06_pnold_mat_data["appli"])


def test_iau06in(iau06arr, iau06_mat_data):
    # Load MATLAB data
    matlab_data = iau06_mat_data["iau06arr"]

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
