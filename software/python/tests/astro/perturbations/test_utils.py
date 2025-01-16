import os

import numpy as np
import pytest

import src.valladopy.astro.perturbations.utils as utils
import src.valladopy.astro.time.data as data

from ...conftest import DEFAULT_TOL, load_matlab_data, custom_allclose


@pytest.fixture()
def gravarr_norm(test_data_dir):
    struct_name = "gravarr_norm"
    file_path = os.path.join(test_data_dir, "gravarr_norm.mat")
    return load_matlab_data(file_path, keys=[struct_name])[struct_name]


@pytest.mark.parametrize(
    "x, expected_result",
    [
        # Test scalar input
        (0.5, -0.4375),
        # Test vector input
        ([0.5, 0.6, 0.7], [-0.4375, -0.36, -0.1925]),
        # Test matrix input
        (
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]],
            [
                [-0.1475, -0.28, -0.3825],
                [-0.44, -0.4375, -0.36],
                [-0.1925, 0.08, 0.4725],
            ],
        ),
    ],
)
def test_leg_poly(x, expected_result):
    # Call leg_poly method
    result = utils.leg_poly(x, order=3)

    # Check results
    assert np.allclose(result, expected_result, rtol=DEFAULT_TOL)


def test_read_gravity_field(gravarr_norm):
    # Read gravity field data
    filepath = os.path.join(data.DATA_DIR, "EGM-08norm100.txt")
    gravity_field_data = utils.read_gravity_field(filepath, normalized=True)

    # Check results
    assert custom_allclose(gravarr_norm.cNor, gravity_field_data.c)
    assert custom_allclose(gravarr_norm.sNor, gravity_field_data.s)
    assert gravity_field_data.normalized
