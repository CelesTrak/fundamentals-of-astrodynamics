import os

import numpy as np
import pytest

import src.valladopy.astro.perturbations.utils as utils

from ...conftest import custom_allclose


def test_legpolyn():
    # Define input values
    latgc = np.radians(30.6103084177511)
    order = 5

    # Call leg_polyn method
    legarr_mu, legarr_gu, legarr_mn, legarr_gn = utils.legpolyn(latgc, order)

    # Expected results
    # fmt: off
    legarr_mu_exp = np.array([
        [1, 0, 0, 0, 0, 0],
        [0.5091962686273478, 0.8606504284644177, 0, 0, 0, 0],
        [-0.11107874002397877, 1.314719960299829, 2.2221574800479575, 0, 0, 0],
        [-0.43373231232496057, 0.38265060248979266, 5.6575714857138495,
         9.562503936593428, 0, 0],
        [-0.3031869762652333, -1.2983233427895509, 4.52745631494301,
         34.084339262733884, 57.609811771552714, 0],
        [0.06909883124077135, -1.9657914067577456, -4.39905138560079,
         44.631518892380164, 264.0123107135865, 446.2371826644717]
    ])
    legarr_mn_exp = np.array([
        [1, 0, 0, 0, 0, 0],
        [0.8819538082870567, 1.490690269656295, 0, 0, 0, 0],
        [-0.24837961354864316, 1.697296170389238, 1.4343964854793299, 0., 0, 0],
        [-1.147547833984841, 0.4133098888043146, 1.9324285489668902,
         1.33342746622041, 0, 0],
        [-0.9095609287956998, -1.2316976707735587, 1.0123700085373268,
         2.036928870853955, 1.2172294383209188, 0],
        [0.22917489667772642, -1.6834031880629696, -0.7119192437226732,
         1.4743742708621626, 2.0556728600854277, 1.0987416280897822]
    ])
    # fmt: on

    # Check results
    assert custom_allclose(legarr_mu, legarr_mu_exp)
    assert custom_allclose(legarr_gu, legarr_mu_exp)
    assert custom_allclose(legarr_mn, legarr_mn_exp)
    assert custom_allclose(legarr_gn, legarr_mn_exp)


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
    gravity_field_data = utils.read_gravity_field(filepath, normalized=True)

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
        s_unc_exp = [0.1189221187e-09, 0.1194314960e-09, 0.119350785e-09]
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
