import os

import numpy as np
import pytest

import src.valladopy.astro.perturbations.gravity as gravity

from ...conftest import custom_allclose


@pytest.mark.parametrize(
    "filename, shape_exp, has_uncertainties",
    # NOTE: These files are assumed to exist in the datalib directory!
    [
        ("EGM-08norm100.txt", (101, 101), False),  # 100x100 with no uncertainties
        ("egm2008 norm 486.txt", (487, 487), True),  # 486x486 with uncertainties
    ],
)
def test_read_gravity_field(data_dir, filename, shape_exp, has_uncertainties):
    # Read gravity field data
    filepath = os.path.join(data_dir, filename)
    gravity_field_data = gravity.read_gravity_field(filepath, normalized=True)

    # Check the second row of the gravity field data
    c_exp = [-0.484165143790815e-03, -0.206615509074176e-09, 0.243938357328313e-05]
    s_exp = [0, 0.138441389137979e-08, -0.140027370385934e-05]
    assert custom_allclose(gravity_field_data.c[2, :3], c_exp)
    assert custom_allclose(gravity_field_data.s[2, :3], s_exp)

    if has_uncertainties:
        c_unc_exp = [0.748123949e-11, 0.7063781502e-11, 0.7230231722e-11]
        s_unc_exp = [0, 0.7348347201e-11, 0.7425816951e-11]
        assert custom_allclose(gravity_field_data.c_unc[2, :3], c_unc_exp)
        assert custom_allclose(gravity_field_data.s_unc[2, :3], s_unc_exp)

    # Check the last 3 elements of the 100th row
    c_exp = [0.599897072379349e-09, 0.580871480377766e-10, 0.995655505739113e-09]
    s_exp = [0.495325263424430e-09, 0.138141678432454e-08, -0.801941613138099e-09]
    assert custom_allclose(gravity_field_data.c[100, 98:101], c_exp)
    assert custom_allclose(gravity_field_data.s[100, 98:101], s_exp)

    if has_uncertainties:
        c_unc_exp = [0.1189226918e-09, 0.1194318875e-09, 0.119352849e-09]
        s_unc_exp = [0.1189221187e-09, 0.119431496e-09, 0.119350785e-09]
        assert custom_allclose(gravity_field_data.c_unc[100, 98:101], c_unc_exp)
        assert custom_allclose(gravity_field_data.s_unc[100, 98:101], s_unc_exp)

    # Check if uncertainties are included in the data
    if not has_uncertainties:
        assert not gravity_field_data.c_unc
        assert not gravity_field_data.s_unc

    # Structural checks
    assert gravity_field_data.c.shape == shape_exp
    assert gravity_field_data.s.shape == shape_exp
    assert gravity_field_data.normalized

    # Ensure all values are finite (not NaN or inf)
    assert np.isfinite(gravity_field_data.c).all()
    assert np.isfinite(gravity_field_data.s).all()


@pytest.fixture
def gravarr():
    gravarr = gravity.GravityFieldData()
    # fmt: off
    gravarr.c = np.array([
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [-0.000484165143790815, -2.06615509074176e-10, 2.43938357328313e-06,
         0, 0, 0],
        [9.57161207093473e-07, 2.03046201047864e-06, 9.04787894809528e-07,
         7.21321757121568e-07, 0, 0],
        [5.39965866638991e-07, -5.36157389388867e-07, 3.50501623962649e-07,
         9.90856766672321e-07, -1.88519633023033e-07, 0],
        [6.86702913736681e-08, -6.29211923042529e-08, 6.52078043176164e-07,
         -4.51847152328843e-07, -2.95328761175629e-07, 1.74811795496002e-07]
    ])
    gravarr.s = np.array([
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 1.38441389137979e-09, -1.40027370385934e-06, 0, 0, 0],
        [0, 2.48200415856872e-07, -6.19005475177618e-07, 1.41434926192941e-06, 0, 0],
        [0, -4.73567346518086e-07, 6.62480026275829e-07,
         -2.00956723567452e-07, 3.08803882149194e-07, 0],
        [0, -9.43698073395769e-08, -3.23353192540522e-07,
         -2.14955408306046e-07, 4.98070550102351e-08, -6.69379935180165e-07]
    ])
    # fmt: on
    gravarr.normalized = True
    return gravarr


def test_get_norm():
    norm_arr = gravity.get_norm(degree=5)
    # fmt: off
    assert custom_allclose(
        norm_arr,
        np.array(
            [
                [1, 0, 0, 0, 0, 0],
                [1.7320508075688772, 1.7320508075688772, 0, 0, 0, 0],
                [2.23606797749979, 1.2909944487358056, 0.6454972243679028, 0, 0, 0],
                [2.6457513110645907, 1.0801234497346435, 0.3415650255319866,
                 0.13944333775567927, 0, 0],
                [3, 0.9486832980505138, 0.22360679774997896,
                 0.05976143046671968, 0.021128856368212913, 0],
                [3.3166247903554, 0.8563488385776752, 0.1618347187425374,
                 0.033034373632170495, 0.007786276535852612, 0.0024622368345219954]
            ]
        )
    )


def test_get_norm_gott():
    norms = gravity.get_norm_gott(degree=5)
    norm1, norm2, norm11, normn10, norm1m, norm2m, normn1 = norms

    # Expected normalization arrays
    # fmt: off
    assert custom_allclose(
        norm1,
        np.array([
            0, 1.2909944487358056, 1.1832159566199232,
            1.1338934190276817, 1.1055415967851334, 1.087114613009218
        ])
    )
    assert custom_allclose(
        norm2,
        np.array([
            0, 2.23606797749979, 1.5275252316519468,
            1.3416407864998738, 1.2535663410560174, 1.2018504251546631
        ])
    )
    assert custom_allclose(
        norm11,
        np.array([
            0, 0.37267799624996495, 0.2160246899469287,
            0.1515228816828316, 0.11653431646335018, 0.09462118179391511
        ])
    )
    assert custom_allclose(
        normn10,
        np.array([
            0, 1.7320508075688772, 2.449489742783178,
            3.1622776601683795, 3.872983346207417, 4.58257569495584
        ])
    )
    assert custom_allclose(
        norm1m,
        np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [0.7453559924999299, 0, 0, 0, 0, 0],
                [0.8366600265340756, 0.5291502622129182, 0, 0, 0, 0],
                [0.8783100656536799, 0.6546536707079771, 0.4285714285714286, 0, 0, 0],
                [0.90267093384844, 0.7237468644557459, 0.5527707983925667,
                 0.3685138655950444, 0, 0],
                [0.9187795405622853, 0.7687061147858074, 0.6276459144608478,
                 0.4861724348043977, 0.3277773886785445, 0]
            ]
        )
    )
    assert custom_allclose(
        norm2m,
        np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0.6236095644623235, 0, 0, 0, 0, 0],
                [0.7348469228349535, 0.34641016151377546, 0, 0, 0, 0],
                [0.7928249671720918, 0.47380354147934284, 0.23690177073967142, 0, 0, 0],
                [0.8293555858801988, 0.5563486402641868, 0.3469443332443555,
                 0.17916128329552333, 0, 0]
            ]
        )
    )
    assert custom_allclose(
        normn1,
        np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [2, 0, 0, 0, 0, 0],
                [3.1622776601683795, 2.449489742783178, 0, 0, 0, 0],
                [4.242640687119285, 3.7416573867739413, 2.8284271247461903, 0, 0, 0],
                [5.291502622129181, 4.898979485566356, 4.242640687119285,
                 3.1622776601683795, 0, 0],
                [6.324555320336759, 6, 5.477225575051661,
                 4.69041575982343, 3.4641016151377544, 0]
            ]
        )
    )


def test_accel_gott(gravarr):
    # Test acceleration calculation
    recef = np.array([-2110.289523, -5511.916033, 3491.913394])
    degree = 5
    order = 5
    leg_gott_n, accel = gravity.accel_gott(recef, gravarr, degree, order)

    # Expected results
    # fmt: off
    leg_gott_n_exp = np.array(
        [
            [1, 0, 0, 0, 0, 0],
            [0.8819538072013969, 1.7320508075688772, 0, 0, 0, 0],
            [-0.24837961568968236, 1.972108665917067, 1.9364916731037085, 0, 0, 0],
            [-1.147547834722164, 0.4802296852291552, 2.608854544205961,
             2.0916500663351894, 0, 0],
            [-0.9095609259590189, -1.4311242208979975, 1.3667393225859021,
             3.1951812232228924, 2.218529918662356, 0],
            [0.22917490142604802, -1.955966243255054, -0.9611189873469996,
             2.312742992243308, 3.746682091224139, 2.3268138086232857]
        ]
    )
    accel_exp = np.array(
        [0.0026070624549907433, 0.006809494130544516, -0.004326141565980874]
    )
    # fmt: on

    # Check the acceleration values
    assert custom_allclose(leg_gott_n, leg_gott_n_exp)
    assert custom_allclose(accel, accel_exp)

    # Check that we get an error if the gravity field data is not normalized
    gravarr.normalized = False
    with pytest.raises(ValueError):
        gravity.accel_gott(recef, gravarr, degree, order)
