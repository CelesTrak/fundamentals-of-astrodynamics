import numpy as np
import pytest

import src.valladopy.mathtime.rotations as rot

from ..conftest import custom_allclose


@pytest.fixture
def quats():
    # Normalized quaternions for testing
    qa = np.array([0.3, -0.5, 0.4, 0.7])
    qb = np.array([-0.2, 0.6, 0.1, 0.75])
    return qa / np.linalg.norm(qa), qb / np.linalg.norm(qb)


@pytest.mark.parametrize(
    "dir_a, dir_b, expected",
    [
        (
            1,
            1,
            [
                -0.20892550412042618,
                -0.06624467203818406,
                0.4586169602643503,
                0.861180736496391,
            ],
        ),
        (
            -1,
            1,
            [
                -0.07643616004405843,
                0.9223296645316378,
                -0.31593612818210803,
                0.20892550412042607,
            ],
        ),
        (
            1,
            -1,
            [
                0.6675424643847765,
                -0.6981169284023998,
                0.1528723200881168,
                0.20892550412042607,
            ],
        ),
        (
            -1,
            -1,
            [
                -0.38218080022029194,
                -0.1579680640910539,
                -0.29555315217035905,
                0.861180736496391,
            ],
        ),
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
    x_axis, y_axis, z_axis = rot.quat2body(quats[0])
    assert custom_allclose(
        x_axis, [0.1717171717171717, 0.2626262626262626, 0.9494949494949494]
    )
    assert custom_allclose(
        y_axis, [-0.8686868686868687, 0.4949494949494949, 0.020202020202020166]
    )
    assert custom_allclose(
        z_axis, [-0.46464646464646453, -0.8282828282828283, 0.31313131313131315]
    )


@pytest.mark.parametrize(
    "direction, expected",
    [
        (1, [3.545454545454545, 0.18181818181818143, -1.1818181818181817]),
        (-1, [-2.9595959595959593, -1.2323232323232327, 1.929292929292929]),
    ],
)
def test_vec_by_quat(quats, direction, expected):
    rotated_vec = rot.vec_by_quat(quats[0], vec=[1, 2, 3], direction=direction)
    assert custom_allclose(rotated_vec, expected)


def test_vec_by_quat_invalid_direction(quats):
    with pytest.raises(ValueError):
        rot.vec_by_quat(quats[0], vec=[1, 2, 3], direction=3)


class TestQuatRV:
    @pytest.fixture
    def orbit_state(self):
        # Get orbit quantities from test state vectors
        ro = [-6518.1083, -2403.8479, -22.1722]
        vo = [2.604057, -7.105717, -0.263218]
        rmag = np.linalg.norm(ro)
        dot = np.dot(ro, vo) / rmag
        omega = np.linalg.norm(np.cross(ro, vo)) / rmag**2
        return rmag, dot, omega

    @pytest.fixture
    def rv(self):
        # Position and velocity vectors
        r = [3228.0295171104126, 5754.3134870229105, -2175.4111963135388]
        v = [1.307894892614813, 2.0022297314596504, 7.184851961675554]
        return r, v

    def test_quat2rv(self, quats, orbit_state, rv):
        # Compute position and velocity vectors
        r, v = rot.quat2rv(q=[*quats[0], *orbit_state])

        # Compare with expected values
        assert custom_allclose(r, rv[0])
        assert custom_allclose(v, rv[1])

    def test_rv2quat(self, orbit_state, rv):
        # Convert position and velocity vectors to 7-element quaternion
        q = rot.rv2quat(*rv)

        # Expected quaternion components
        q_expected = [
            0.3015113445777637,
            -0.502518907629606,
            0.40201512610368495,
            0.7035264706814484,
            *orbit_state,
        ]

        assert custom_allclose(q, q_expected)


class TestQuatDCM:
    @pytest.fixture
    def dcm(self):
        return [
            [0.17171717171717157, 0.2626262626262626, 0.9494949494949494],
            [-0.8686868686868687, 0.4949494949494947, 0.020202020202020166],
            [-0.46464646464646453, -0.8282828282828283, 0.3131313131313129],
        ]

    def test_quat2dcm(self, quats, dcm):
        dcm_out = rot.quat2dcm(quats[0])
        assert custom_allclose(dcm_out, dcm)

    def test_dcm2quat(self, dcm, quats):
        q_out = rot.dcm2quat(dcm)
        assert custom_allclose(q_out, quats[0])
