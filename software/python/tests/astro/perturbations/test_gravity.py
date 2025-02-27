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
