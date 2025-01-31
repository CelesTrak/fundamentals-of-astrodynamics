import os

import numpy as np
import pytest

import src.valladopy.astro.perturbations.utils as utils
import src.valladopy.astro.time.data as data

from ...conftest import load_matlab_data, custom_allclose


@pytest.fixture()
def gravarr_norm(test_data_dir):
    struct_name = "gravarr_norm"
    file_path = os.path.join(test_data_dir, "gravarr_norm.mat")
    return load_matlab_data(file_path, keys=[struct_name])[struct_name]


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


def test_read_gravity_field(gravarr_norm):
    # Read gravity field data
    filepath = os.path.join(data.DATA_DIR, "EGM-08norm100.txt")
    gravity_field_data = utils.read_gravity_field(filepath, normalized=True)

    # Check results
    assert custom_allclose(gravarr_norm.cNor, gravity_field_data.c)
    assert custom_allclose(gravarr_norm.sNor, gravity_field_data.s)
    assert gravity_field_data.normalized
