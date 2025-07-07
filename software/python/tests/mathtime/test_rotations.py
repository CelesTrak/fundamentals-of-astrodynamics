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
        rot.quat_multiply(qa, qb, dir_a=2)  # Invalid direction sign
    with pytest.raises(ValueError):
        rot.quat_multiply(qa, qb, dir_b=-3)  # Invalid direction sign
