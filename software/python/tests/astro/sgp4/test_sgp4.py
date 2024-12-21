import numpy as np
import pytest

import src.valladopy.astro.sgp4.sgp4 as sgp4
from src.valladopy.astro.sgp4.utils import WGSModel

from ...conftest import DEFAULT_TOL, custom_isclose


def set_jd_startstop(sgp4_obj):
    start_ymdhms = (2024, 12, 12, 2, 3, 4)
    stop_ymdhms = (2024, 12, 14, 5, 6, 7)
    sgp4_obj.set_jd_from_from_ymdhms(start_ymdhms, stop_ymdhms)


@pytest.mark.parametrize(
    "typerun, startmfe_exp, stopmfe_exp, deltamin_exp",
    [
        (sgp4.TypeRun.Catalog, -1440, 1440, 20),
        (sgp4.TypeRun.Verification, 0, 4320, 360),
        (sgp4.TypeRun.FromJD, 12863952.7377735, 12867015.7877738, 1),
    ],
)
def test_twoline2rv(typerun, startmfe_exp, stopmfe_exp, deltamin_exp):
    # TODO: Add more test cases
    # Initialize SGP4 class
    sgp4_obj = sgp4.SGP4()

    # Inputs
    tle_line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    tle_line2 = (
        "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667     0.00"
        "      4320.0        360.00"
    )

    # Expected sgp4_obj.satrec attributes
    expected = {
        "error": 0,
        "satnum": 5,
        "classification": sgp4.Classification.Unclassified,
        "intldesg": "58002B_",
        "epochyr": 0.0,
        "epochdays": 179.78495062,
        "ndot": 6.96919666594958e-13,
        "nddot": 0.0,
        "bstar": 2.8098e-05,
        "elnum": 475,
        "inclo": 0.5980929187319208,
        "nodeo": 6.08638547138321,
        "ecco": 0.1859667,
        "argpo": 5.790416027488515,
        "mo": 0.3373093125574321,
        "no_kozai": 0.04722944544077857,
        "revnum": 41366,
        "jdsatepoch": 2451722.5,
        "jdsatepochf": 0.7849506199999999,
    }

    # Set JD start and stop if typerun is `FromJD`
    if typerun == sgp4.TypeRun.FromJD:
        set_jd_startstop(sgp4_obj)

    # Call method
    startmfe, stopmfe, deltamin = sgp4_obj.twoline2rv(tle_line1, tle_line2, typerun)

    # Check results
    assert np.isclose(startmfe, startmfe_exp, rtol=DEFAULT_TOL)
    assert np.isclose(stopmfe, stopmfe_exp, rtol=DEFAULT_TOL)
    assert np.isclose(deltamin, deltamin_exp, rtol=DEFAULT_TOL)
    for key in expected:
        # Non-float comparisons
        if key in ["error", "satnum", "elnum", "revnum", "classification", "intldesg"]:
            assert getattr(sgp4_obj.satrec, key) == expected[key]
        # Float comparisons
        else:
            assert custom_isclose(getattr(sgp4_obj.satrec, key), expected[key])


def test_initl(epoch, oe_params):
    # Inputs
    sgp4_obj = sgp4.SGP4(wgs_model=WGSModel.WGS_72)
    ecco, inclo, *_ = oe_params
    no_kozai = 0.00874808688806747

    # Expected outputs
    expected = {
        "ainv": 0.24008817584663986,
        "ao": 4.165136398215487,
        "con41": -0.43002188663163776,
        "con42": 0.05003647771939623,
        "cosio": 0.43588152571096744,
        "cosio2": 0.18999270445612076,
        "eccsq": 0.47295137105315993,
        "omeosq": 0.5270486289468401,
        "posq": 4.819032241803303,
        "rp": 1.3007112861712828,
        "rteosq": 0.7259811491676903,
        "sinio": 0.9000040530708066,
        "gsto": 0.574180126902192,
        "no_unkozai": 0.008748547019630244,
    }

    # Call method
    sgp4init_out = sgp4_obj.initl(epoch, ecco, inclo, no_kozai)

    # Check results
    for key in expected:
        assert custom_isclose(getattr(sgp4init_out, key), expected[key])
