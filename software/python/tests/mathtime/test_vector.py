import numpy as np
import pytest

import src.valladopy.mathtime.vector as vec


###############################################################################
# Axes Rotations
###############################################################################


@pytest.mark.parametrize(
    "func, angle, vector, expected",
    [
        (vec.rot1, np.pi / 2, [1, 0, 0], [1, 0, 0]),
        (vec.rot2, np.pi / 2, [1, 0, 0], [0, 0, 1]),
        (vec.rot3, np.pi / 2, [1, 0, 0], [0, -1, 0]),
        (vec.rot1, np.pi, [0, 1, 0], [0, -1, 0]),
        (vec.rot2, np.pi, [0, 0, 1], [0, 0, -1]),
        (vec.rot3, np.pi, [0, 1, 0], [0, -1, 0]),
    ],
)
def test_rotations(func, angle, vector, expected):
    assert np.allclose(func(vector, angle), expected, rtol=1e-12)


###############################################################################
# Rotation Matrices
###############################################################################


@pytest.mark.parametrize(
    "func, angle, expected",
    [
        (vec.rot1mat, np.pi / 2, np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])),
        (vec.rot2mat, np.pi / 2, np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])),
        (vec.rot3mat, np.pi / 2, np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])),
        (
            vec.rot1mat,
            np.pi / 4,
            np.array(
                [
                    [1, 0, 0],
                    [0, 1 / np.sqrt(2), 1 / np.sqrt(2)],
                    [0, -1 / np.sqrt(2), 1 / np.sqrt(2)],
                ]
            ),
        ),
        (
            vec.rot2mat,
            np.pi / 4,
            np.array(
                [
                    [1 / np.sqrt(2), 0, -1 / np.sqrt(2)],
                    [0, 1, 0],
                    [1 / np.sqrt(2), 0, 1 / np.sqrt(2)],
                ]
            ),
        ),
        (
            vec.rot3mat,
            np.pi / 4,
            np.array(
                [
                    [1 / np.sqrt(2), 1 / np.sqrt(2), 0],
                    [-1 / np.sqrt(2), 1 / np.sqrt(2), 0],
                    [0, 0, 1],
                ]
            ),
        ),
    ],
)
def test_rotation_matrices(func, angle, expected):
    assert np.allclose(func(angle), expected, rtol=1e-12)


###############################################################################
# Vector Math
###############################################################################


@pytest.mark.parametrize(
    "v1, v2, expected_angle",
    [
        (np.array([1, 0, 0]), np.array([1, 0, 0]), 0),
        (np.array([1, 0, 0]), np.array([-1, 0, 0]), np.pi),
        (np.array([1, 0, 0]), np.array([0, 1, 0]), np.pi / 2),
        (np.array([1, 0, 0]), np.array([1, 1, 0]), np.arccos(1 / np.sqrt(2))),
        (np.array([0, 0, 0]), np.array([1, 0, 0]), np.nan),
        (np.array([1, 0, 0]), np.array([0, 0, 0]), np.nan),
        (np.array([1e-10, 0, 0]), np.array([1e-10, 0, 0]), np.nan),
    ],
)
def test_angle_between_vectors(v1, v2, expected_angle):
    angle = vec.angle(v1, v2)
    assert (
        np.isnan(angle)
        if np.isnan(expected_angle)
        else abs(angle - expected_angle) < 1e-12
    )


@pytest.mark.parametrize(
    "v, expected_v_unit",
    [
        (np.array([1, 0, 0]), np.array([1, 0, 0])),
        (np.array([1, 1, 0]), np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0])),
        (np.array([0, 0, 0]), np.array([0, 0, 0])),
    ],
)
def test_unit_vector(v, expected_v_unit):
    v_unit = vec.unit(v)
    assert np.allclose(v_unit, expected_v_unit, rtol=1e-12)
