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
