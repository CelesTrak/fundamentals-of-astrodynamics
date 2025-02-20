import numpy as np
import pytest

import src.valladopy.astro.iod.lambert as lambert
import src.valladopy.constants as const
from ...conftest import DEFAULT_TOL, custom_allclose


@pytest.fixture
def lambert_inputs():
    r1 = np.array([2.5 * const.RE, 0, 0])
    r2 = np.array([1.9151111 * const.RE, 1.6069690 * const.RE, 0])
    v1 = np.array([0, 4.999792554221911, 0])
    nrev = 1
    dtsec = 92854.234
    return r1, r2, v1, nrev, dtsec


def test_calculate_mag_and_angle(lambert_inputs):
    # Unpack inputs
    r1, r2, *_ = lambert_inputs

    # Compute magnitude and angle
    magr1, magr2, angle = lambert.calculate_mag_and_angle(r1, r2)

    # Check results
    assert np.isclose(magr1, 15945.34075, rtol=DEFAULT_TOL)
    assert np.isclose(magr2, 15945.340612620314, rtol=DEFAULT_TOL)
    assert np.isclose(angle, 0.7660444465999809, rtol=DEFAULT_TOL)


class TestLambert:
    @pytest.mark.parametrize(
        "dm, v_exp, aminenergy_exp, tminenergy_exp, tminabs_exp",
        [
            (
                lambert.DirectionOfMotion.LONG,
                [-2.0474089759890735, -2.924003076447717, 0],
                10699.484172968232,
                15554.50821587732,
                1534.8915813389815,
            ),
            (
                lambert.DirectionOfMotion.SHORT,
                [2.0474089759890735, 2.924003076447717, 0],
                10699.484172968232,
                17488.265508772805,
                1534.8915813389815,
            ),
        ],
    )
    def test_min_energy(
        self, lambert_inputs, dm, v_exp, aminenergy_exp, tminenergy_exp, tminabs_exp
    ):
        # Unpack inputs
        r1, r2, _, nrev, _ = lambert_inputs

        # Compute Lambert minimum energy
        out = lambert.min_energy(r1, r2, dm, nrev)
        v, aminenergy, tminenergy, tminabs = out

        # Check results
        assert np.allclose(v, v_exp, rtol=DEFAULT_TOL)
        assert np.isclose(aminenergy, aminenergy_exp, rtol=DEFAULT_TOL)
        assert np.isclose(tminenergy, tminenergy_exp, rtol=DEFAULT_TOL)
        assert np.isclose(tminabs, tminabs_exp, rtol=DEFAULT_TOL)

    @pytest.mark.parametrize(
        "dm, de, tmin_exp, tminp_exp, tminenergy_exp",
        [
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.LOW,
                12468.267989702517,
                3139.7046602121873,
                17488.26550877280,
            ),
            (
                lambert.DirectionOfMotion.SHORT,
                lambert.DirectionOfEnergy.LOW,
                15048.526832075213,
                1534.8915813389815,
                15554.50821587732,
            ),
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.HIGH,
                21648.968497701877,
                3139.7046602121873,
                17488.265508772805,
            ),
            (
                lambert.DirectionOfMotion.SHORT,
                lambert.DirectionOfEnergy.HIGH,
                16972.34235386572,
                1534.8915813389815,
                15554.50821587732,
            ),
        ],
    )
    def test_min_time(
        self, lambert_inputs, dm, de, tmin_exp, tminp_exp, tminenergy_exp
    ):
        # Unpack inputs
        r1, r2, _, nrev, _ = lambert_inputs

        # Compute Lambert minimum time
        tmin, tminp, tminenergy = lambert.min_time(r1, r2, dm, de, nrev)

        # Check results
        assert np.allclose(tmin, tmin_exp, rtol=DEFAULT_TOL)
        assert np.allclose(tminp, tminp_exp, rtol=DEFAULT_TOL)
        assert np.allclose(tminenergy, tminenergy_exp, rtol=DEFAULT_TOL)

    @pytest.mark.parametrize(
        "dm, tmaxrp_exp, v1t_exp",
        [
            (lambert.DirectionOfMotion.LONG, 37850.19485964372, np.zeros(3)),
            (
                lambert.DirectionOfMotion.SHORT,
                22264.820236293686,
                np.array([-1.1835179322811704e-07, 4.999792624745, 0]),
            ),
        ],
    )
    def test_tmax_rp(self, lambert_inputs, dm, tmaxrp_exp, v1t_exp):
        # Unpack inputs
        r1, r2, _, nrev, _ = lambert_inputs

        # Compute Lambert maximum radius of perigee
        tmaxrp, v1t = lambert.tmax_rp(r1, r2, dm, nrev)

        # Check results
        assert np.isclose(tmaxrp, tmaxrp_exp, rtol=DEFAULT_TOL)
        assert custom_allclose(v1t, v1t_exp)


class TestLambertBattin:
    @pytest.fixture
    def v(self):
        return 0.123

    def test_seebatt(self, v):
        # Check nominal case
        assert np.isclose(lambert.seebatt(v), 5.153421950753984, rtol=DEFAULT_TOL)

        # Check when `v` is less than -1
        v = -1.5
        with pytest.raises(ValueError):
            lambert.seebatt(v)

    def test_kbatt(self, v):
        assert np.isclose(lambert.kbatt(v), 0.327568960337347, rtol=DEFAULT_TOL)

    def test_hodograph(self):
        r1 = [6888, 0, 0]
        v1 = [0, 7.6072, 0]
        r2 = [8000 * np.sin(np.pi / 4), 8000 * np.cos(np.pi / 4), 0]
        p = 7000
        ecc = 1e-3
        dnu = np.radians(2.1)
        dtsec = 1234

        # Check nominal case
        v1t, v2t = lambert.hodograph(r1, r2, v1, p, ecc, dnu, dtsec)
        assert custom_allclose(v1t, [0, 7.668753340719345, 0])
        assert custom_allclose(v2t, [-4.665702986379401, 4.672061552209318, 0])

        # Check when r1 and v1 are parallel and dnu is zero
        new_v1 = [7.6072, 0, 0]
        new_dnu = 0
        with pytest.raises(ValueError):
            lambert.hodograph(r1, r2, new_v1, p, ecc, new_dnu, dtsec)

        # Check when r1 and r2 are parallel
        new_r2 = [8000, 0, 0]
        with pytest.raises(ValueError):
            lambert.hodograph(r1, new_r2, v1, p, ecc, dnu, dtsec)

    @pytest.mark.parametrize(
        "dm, de, v1t_exp, v2t_exp",
        [
            (
                lambert.DirectionOfMotion.SHORT,
                lambert.DirectionOfEnergy.HIGH,
                [-0.8696153795282852, 6.3351545812502374, 0],
                [-3.405994961791248, 5.41198791828363, 0],
            ),
            (
                lambert.DirectionOfMotion.SHORT,
                lambert.DirectionOfEnergy.LOW,
                [5.832522716212579, 1.4319944881331306, 0],
                [-5.388439978490882, -2.652101898141935, 0],
            ),
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.HIGH,
                [-6.241103309400493, -1.351339299630816, 0],
                [5.6495867154901545, 2.976517897853268, 0],
            ),
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.LOW,
                [0.6411191586146303, -5.957501823796459, 0],
                [3.33828270226307, -4.975814585231199, 0],
            ),
        ],
    )
    def test_battin(self, lambert_inputs, dm, de, v1t_exp, v2t_exp):
        # Unpack inputs and set dtsec
        r1, r2, v1, nrev, dtsec = lambert_inputs

        # Compute Lambert minimum time
        v1t, v2t = lambert.battin(r1, r2, v1, dm, de, nrev, dtsec)

        # Check results
        assert np.allclose(v1t, v1t_exp, rtol=DEFAULT_TOL)
        assert np.allclose(v2t, v2t_exp, rtol=DEFAULT_TOL)

    def test_battin_bad_inputs(self, lambert_inputs):
        # Unpack inputs
        r1, r2, v1, nrev, _ = lambert_inputs
        dm = lambert.DirectionOfMotion.LONG
        df = lambert.DirectionOfFlight.DIRECT

        # Check when dtsec is 0
        with pytest.raises(ValueError):
            lambert.battin(r1, r2, v1, dm, df, nrev, dtsec=0)

        # Check when dtsec is positive but within a range without a solution
        # (hits the error in a different part of the algorithm)
        with pytest.raises(ValueError):
            lambert.battin(r1, r2, v1, dm, df, nrev, dtsec=4560)


class TestLambertUniversal:
    @pytest.mark.parametrize(
        "psi, c2dot_exp, c3dot_exp, c2ddot_exp, c3ddot_exp",
        [
            (5, 0.08926, -0.00134, -0.036254, 0.009596),
            (0, -0.041666666667, -0.008333333333, 0, 0),
        ],
    )
    def test_calculate_c2dot_c3dot(
        self, psi, c2dot_exp, c3dot_exp, c2ddot_exp, c3ddot_exp
    ):
        # Check nominal case
        c2, c3 = 0.0232, 0.0122
        out = lambert._calculate_c2dot_c3dot(psi, c2, c3)
        c2dot, c3dot, c2ddot, c3ddot = out
        assert np.isclose(c2dot, c2dot_exp, rtol=DEFAULT_TOL)
        assert np.isclose(c3dot, c3dot_exp, rtol=DEFAULT_TOL)
        assert np.isclose(c2ddot, c2ddot_exp, rtol=DEFAULT_TOL)
        assert np.isclose(c3ddot, c3ddot_exp, rtol=DEFAULT_TOL)

    def test_calculate_dtdpsi(self):
        # Check nominal case
        x, y, vara = 1346.4822, 42071.6835, -21190.1786
        c2, c3, c2dot, c3dot = 0.0232, 0.0122, 1.756e-04, -8.67e-05
        dtdpsi = lambert._calculate_dtdpsi(x, c2, c3, c2dot, c3dot, vara, y)
        assert np.isclose(dtdpsi, -2162.3620971040746, rtol=DEFAULT_TOL)

    def test_get_kbiu(self, lambert_inputs):
        # Unpack inputs
        r1, r2, *_ = lambert_inputs

        # Expected values
        kbi_exp = np.array(
            [
                [57.81592755194815, 15048.626881962446],
                [195.80115373119972, 26269.95060241361],
                [412.98000230103014, 37370.33586904963],
            ]
        )
        kbil_exp = np.array(
            [
                [113.60281509526207, 16971.386078436695],
                [291.6899864823016, 28199.561813444412],
                [548.316527089472, 39301.8313984678],
            ]
        )

        # Compute universal variable Lambert problem
        kbi, kbil = lambert.get_kbiu(r1, r2, order=3)

        # Check results
        assert np.allclose(kbi, kbi_exp, rtol=DEFAULT_TOL)
        assert np.allclose(kbil, kbil_exp, rtol=DEFAULT_TOL)

    @pytest.mark.parametrize(
        "dm, kbi_exp, tof_exp",
        [
            (lambert.DirectionOfMotion.LONG, 113.602815095262, 16971.386078436695),
            (lambert.DirectionOfMotion.SHORT, 57.8159275519482, 15048.626881962446),
        ],
    )
    def test_universal_min(self, lambert_inputs, dm, kbi_exp, tof_exp):
        # Unpack inputs
        r1, r2, _, nrev, _ = lambert_inputs

        # Compute universal variable Lambert problem
        kbi, tof = lambert.universal_min(r1, r2, dm, nrev)

        # Check results
        assert np.isclose(kbi, kbi_exp, rtol=DEFAULT_TOL)
        assert np.isclose(tof, tof_exp, rtol=DEFAULT_TOL)

    @pytest.mark.parametrize(
        "dm, de, nrev, kbi, v1dv_exp, v2dv_exp",
        [
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.HIGH,
                1,
                113,
                [-6.241103311949874, -1.3513393001017955, 0],
                [5.649586717745847, 2.976517899131199, 0],
            ),
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.LOW,
                1,
                113,
                [0.6411191590247107, -5.957501826310465, 0],
                [3.338282703564938, -4.975814587420603, 0],
            ),
            (
                lambert.DirectionOfMotion.SHORT,
                lambert.DirectionOfEnergy.HIGH,
                1,
                113,
                [-0.8696153799562063, 6.335154583804734, 0],
                [-3.40599496310541, 5.411987920515576, 0],
            ),
            (
                lambert.DirectionOfMotion.SHORT,
                lambert.DirectionOfEnergy.LOW,
                1,
                57,  # does not converge with 113
                [5.8325227187105835, 1.4319944886088742, 0],
                [-5.388439980710266, -2.6521018993831795, 0],
            ),
            (
                lambert.DirectionOfMotion.LONG,
                lambert.DirectionOfEnergy.HIGH,
                0,  # test when nrev is 0
                113,
                [0.8790611303752992, -6.351161592594172, 0],
                [3.409048205683391, -5.43032161186189, 0],
            ),
        ],
    )
    def test_universal(
        self, lambert_inputs, dm, de, nrev, kbi, v1dv_exp, v2dv_exp, caplog
    ):
        # Unpack inputs
        r1, r2, v1, _, dtsec = lambert_inputs

        # Compute universal variable Lambert problem
        v1dv, v2dv = lambert.universal(r1, v1, r2, dtsec, dm, de, nrev, kbi)

        # Check results
        assert np.allclose(v1dv, v1dv_exp, rtol=DEFAULT_TOL)
        assert np.allclose(v2dv, v2dv_exp, rtol=DEFAULT_TOL)

        # Check that nothing was logged (results are OK)
        assert not caplog.records

    @pytest.mark.parametrize(
        "dtsec, kbi",
        [(4560, 113), (92854.234, 20)],  # dtsec is too small  # kbi is too small
    )
    def test_universal_bad_inputs(self, lambert_inputs, dtsec, kbi, caplog):
        # Unpack inputs
        r1, r2, v1, nrev, _ = lambert_inputs

        # Define other inputs
        dm = lambert.DirectionOfMotion.LONG
        de = lambert.DirectionOfEnergy.HIGH

        # Compute universal variable Lambert problem
        v1dv, v2dv = lambert.universal(r1, v1, r2, dtsec, dm, de, nrev, kbi)

        # Check results
        assert np.allclose(v1dv, np.zeros(3), rtol=DEFAULT_TOL)
        assert np.allclose(v2dv, np.zeros(3), rtol=DEFAULT_TOL)

        # Check logged error
        assert caplog.records[0].levelname == "ERROR"
        assert caplog.records[0].message == (
            "Lambert did not converge! Try different values of dtsec or psi"
        )

    def test_universal_bad_orbit(self, lambert_inputs, caplog):
        # Unpack inputs
        *_, v1, nrev, dtsec = lambert_inputs

        # Define bad position vectors
        r1 = np.array([0.1, 0, 0])
        r2 = np.array([0.15, 0, 0])

        # Define other inputs
        dm = lambert.DirectionOfMotion.LONG
        de = lambert.DirectionOfEnergy.HIGH
        kbi = 113

        # Compute universal variable Lambert problem
        v1dv, v2dv = lambert.universal(r1, v1, r2, dtsec, dm, de, nrev, kbi)

        # Check results
        v1dv_exp = np.array([0, 2823.473599000267, 0])
        v2dv_exp = np.array([0, 2305.3561303158463, 0])
        assert np.allclose(v1dv, v1dv_exp, rtol=DEFAULT_TOL)
        assert np.allclose(v2dv, v2dv_exp, rtol=DEFAULT_TOL)

        # Check logged error
        assert caplog.records[0].levelname == "ERROR"
        assert caplog.records[0].message == "Orbit is not possible"
