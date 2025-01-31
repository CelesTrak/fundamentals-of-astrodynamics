import pytest

import src.valladopy.astro.sgp4.utils as utils

from ...conftest import custom_isclose


@pytest.mark.parametrize(
    "wgs_model, expected",
    [
        (
            utils.WGSModel.WGS_72_LOW_PRECISION,
            {
                "tumin": 13.446839702957643,
                "mu": 398600.79964,
                "radiusearthkm": 6378.135,
                "xke": 0.0743669161,
                "j2": 0.001082616,
                "j3": -0.00000253881,
                "j4": -0.00000165597,
                "j3oj2": -0.002345069720011528,
            },
        ),
        (
            utils.WGSModel.WGS_72,
            {
                "tumin": 13.446839696959309,
                "mu": 398600.8,
                "radiusearthkm": 6378.135,
                "xke": 60 / ((6378.135**3) / 398600.8) ** 0.5,
                "j2": 0.001082616,
                "j3": -0.00000253881,
                "j4": -0.00000165597,
                "j3oj2": -0.002345069720011528,
            },
        ),
        (
            utils.WGSModel.WGS_84,
            {
                "tumin": 13.446851082044981,
                "mu": 398600.5,
                "radiusearthkm": 6378.137,
                "xke": 60 / ((6378.137**3) / 398600.5) ** 0.5,
                "j2": 0.00108262998905,
                "j3": -0.00000253215306,
                "j4": -0.00000161098761,
                "j3oj2": -0.0023388905587420003,
            },
        ),
    ],
)
def test_getgravc(wgs_model, expected):
    grav_const = utils.getgravc(wgs_model)
    assert custom_isclose(grav_const.mu, expected["mu"])
    assert custom_isclose(grav_const.radiusearthkm, expected["radiusearthkm"])
    assert custom_isclose(grav_const.xke, expected["xke"])
    assert custom_isclose(grav_const.j2, expected["j2"])
    assert custom_isclose(grav_const.j3, expected["j3"])
    assert custom_isclose(grav_const.j4, expected["j4"])
    assert custom_isclose(grav_const.j3oj2, expected["j3oj2"])
