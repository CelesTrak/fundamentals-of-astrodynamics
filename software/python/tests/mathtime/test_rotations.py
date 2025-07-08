import pytest

import src.valladopy.mathtime.rotations as rot

from ..conftest import custom_allclose


@pytest.fixture
def quats():
    qa = [0.3, -0.5, 0.4, 0.7]
    qb = [-0.2, 0.6, 0.1, 0.75]
    return qa, qb


@pytest.mark.parametrize(
    "dir_a, dir_b, expected",
    [
        (1, 1, [-0.205, -0.065, 0.45, 0.845]),
        (-1, 1, [-0.075, 0.905, -0.31, 0.205]),
        (1, -1, [0.655, -0.685, 0.15, 0.205]),
        (-1, -1, [-0.375, -0.155, -0.29, 0.845]),
    ],
)
def test_quat_multiply(quats, dir_a, dir_b, expected):
    qa, qb = quats
    q_out = rot.quat_multiply(qa, qb, dir_a, dir_b)
    assert custom_allclose(q_out, expected)


def test_quat_multiply_invalid_direction(quats):
    qa, qb = quats
    with pytest.raises(ValueError):
        rot.quat_multiply(qa, qb, dir_a=2)
    with pytest.raises(ValueError):
        rot.quat_multiply(qa, qb, dir_b=-3)


def test_quat_transform(quats):
    qi, qf = quats
    assert custom_allclose(
        rot.quat_transform(qi, qf),
        [
            -0.07643616004405838,
            0.9223296645316378,
            -0.31593612818210803,
            0.20892550412042618,
        ],
    )


def test_quat2body(quats):
    q = quats[0]
    x_axis, y_axis, z_axis = rot.quat2body(q)
    assert custom_allclose(x_axis, [0.18, 0.26, 0.94])
    assert custom_allclose(y_axis, [-0.86, 0.5, 0.02])
    assert custom_allclose(z_axis, [-0.46, -0.82, 0.32])


@pytest.mark.parametrize(
    "direction, expected", [(1, [3.52, 0.2, -1.14]), (-1, [-2.92, -1.2, 1.94])]
)
def test_vec_by_quat(quats, direction, expected):
    rotated_vec = rot.vec_by_quat(quats[0], vec=[1, 2, 3], direction=direction)
    assert custom_allclose(rotated_vec, expected)


def test_vec_by_quat_invalid_direction(quats):
    with pytest.raises(ValueError):
        rot.vec_by_quat(quats[0], vec=[1, 2, 3], direction=3)
