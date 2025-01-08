import os
import pytest

import numpy as np
import scipy

from ...conftest import custom_isclose, custom_allclose


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
    struct_name = "iau80arr"
    file_path = os.path.join(current_dir, "data", "iau80in_data.mat")
    return load_matlab_data(file_path, keys=[struct_name])[struct_name]


@pytest.fixture()
def iau06_pnold_mat_data(current_dir):
    file_path = os.path.join(current_dir, "data", "iau06in_pnold_data.mat")
    return load_matlab_data(file_path, keys=["apn", "apni", "appl", "appli"])


@pytest.fixture()
def iau06_mat_data(current_dir):
    struct_name = "iau06arr"
    file_path = os.path.join(current_dir, "data", "iau06in_data.mat")
    return load_matlab_data(file_path, keys=[struct_name])[struct_name]


@pytest.fixture()
def xys_data(current_dir):
    struct_name = "xys06table_struct"
    file_path = os.path.join(current_dir, "data", "xys_data.mat")
    return load_matlab_data(file_path, keys=[struct_name])[struct_name]


def test_iau80in(iau80arr, iau80_mat_data):
    # Check that they are the same
    assert np.array_equal(iau80arr.iar80, iau80_mat_data.iar80)
    assert custom_allclose(iau80arr.rar80, iau80_mat_data.rar80)


def test_iau06in_pnold(iau06data_old, iau06_pnold_mat_data):
    # Check that the data is the same
    assert custom_allclose(iau06data_old.apn, iau06_pnold_mat_data["apn"])
    assert np.array_equal(iau06data_old.apni, iau06_pnold_mat_data["apni"])
    assert custom_allclose(iau06data_old.appl, iau06_pnold_mat_data["appl"])
    assert np.array_equal(iau06data_old.appli, iau06_pnold_mat_data["appli"])


def test_iau06in(iau06arr, iau06_mat_data):
    # Check that they are the same
    assert custom_allclose(iau06arr.ax0, iau06_mat_data.ax0)
    assert np.array_equal(iau06arr.ax0i, iau06_mat_data.a0xi)
    assert custom_allclose(iau06arr.ay0, iau06_mat_data.ay0)
    assert np.array_equal(iau06arr.ay0i, iau06_mat_data.a0yi)
    assert custom_allclose(iau06arr.as0, iau06_mat_data.as0)
    assert np.array_equal(iau06arr.as0i, iau06_mat_data.a0si)
    assert custom_allclose(iau06arr.agst, iau06_mat_data.agst)
    assert np.array_equal(iau06arr.agsti, iau06_mat_data.agsti[:, :14])
    assert custom_allclose(iau06arr.apn0, iau06_mat_data.apn0)
    assert np.array_equal(iau06arr.apn0i, iau06_mat_data.apn0i)
    assert custom_allclose(iau06arr.apl0, iau06_mat_data.apl0)
    assert np.array_equal(iau06arr.apl0i, iau06_mat_data.apl0i)
    assert custom_allclose(iau06arr.aapn0[:, :5], iau06_mat_data.aapn0[:, :5])
    assert np.array_equal(iau06arr.aapn0i, iau06_mat_data.aapn0i)


def test_readxys(iau06xysarr, xys_data):
    # Check that the data is the same
    assert custom_allclose(iau06xysarr.jd, xys_data.jd)
    assert custom_allclose(iau06xysarr.jdf, xys_data.jdf)
    assert custom_allclose(iau06xysarr.x, xys_data.x)
    assert custom_allclose(iau06xysarr.y, xys_data.y)
    assert custom_allclose(iau06xysarr.s, xys_data.s)
    assert custom_allclose(iau06xysarr.mjd_tt, xys_data.mjd_tt)


def test_readeop(eoparr):
    # Check that the first line is correct
    assert int(eoparr.mjd[0]) == 37665
    assert custom_isclose(eoparr.xp[0], -0.0127)
    assert custom_isclose(eoparr.yp[0], 0.213)
    assert custom_isclose(eoparr.dut1[0], 0.0326338)
    assert custom_isclose(eoparr.lod[0], 0.001723)
    assert custom_isclose(eoparr.ddpsi[0], 0.064261)
    assert custom_isclose(eoparr.ddeps[0], 0.006067)
    assert custom_isclose(eoparr.dx[0], 0)
    assert custom_isclose(eoparr.dy[0], 0)
    assert eoparr.dat[0] == 2

    # Check that the last line is correct
    assert int(eoparr.mjd[-1]) == 60126
    assert custom_isclose(eoparr.xp[-1], 0.203662)
    assert custom_isclose(eoparr.yp[-1], 0.492913)
    assert custom_isclose(eoparr.dut1[-1], -0.0114449)
    assert custom_isclose(eoparr.lod[-1], -0.0009071)
    assert custom_isclose(eoparr.ddpsi[-1], -0.113661)
    assert custom_isclose(eoparr.ddeps[-1], -0.009266)
    assert custom_isclose(eoparr.dx[-1], 0.000121)
    assert custom_isclose(eoparr.dy[-1], -0.000211)
    assert eoparr.dat[-1] == 37


def test_readspw(spwarr):
    # Check that the first line is correct
    assert custom_isclose(spwarr.mjd[0], 36112)
    assert custom_allclose(spwarr.kparray[0], [43, 40, 30, 20, 37, 23, 43, 37])
    assert custom_isclose(spwarr.sumkp[0], 273)
    assert custom_allclose(spwarr.aparray[0], [32, 27, 15, 7, 22, 9, 32, 22])
    assert custom_isclose(spwarr.avgap[0], 21)
    assert custom_isclose(spwarr.adjf107[0], 268.0)
    assert custom_isclose(spwarr.adjctrf81[0], 265.2)
    assert custom_isclose(spwarr.obsf107[0], 269.3)
    assert custom_isclose(spwarr.obsctrf81[0], 266.6)
    assert custom_isclose(spwarr.obslstf81[0], 230.9)

    # Check that the last line is correct
    assert custom_isclose(spwarr.mjd[-1], 59467)
    assert custom_allclose(spwarr.kparray[-1], [13, 13, 13, 13, 13, 13, 13, 13])
    assert custom_isclose(spwarr.sumkp[-1], 104)
    assert custom_allclose(spwarr.aparray[-1], [5, 5, 5, 5, 5, 5, 5, 5])
    assert custom_isclose(spwarr.avgap[-1], 5)
    assert custom_isclose(spwarr.adjf107[-1], 80.0)
    assert custom_isclose(spwarr.adjctrf81[-1], 69.0)
    assert custom_isclose(spwarr.obsf107[-1], 78.9)
    assert custom_isclose(spwarr.obsctrf81[-1], 68.3)
    assert custom_isclose(spwarr.obslstf81[-1], 78.9)
