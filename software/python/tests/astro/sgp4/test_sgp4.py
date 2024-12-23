import numpy as np
import pytest

import src.valladopy.astro.sgp4.sgp4 as sgp4
from src.valladopy.astro.sgp4.utils import WGSModel

from ...conftest import DEFAULT_TOL, custom_isclose


def set_jd_startstop(sgp4_obj):
    start_ymdhms = (2024, 12, 12, 2, 3, 4)
    stop_ymdhms = (2024, 12, 14, 5, 6, 7)
    sgp4_obj.set_jd_from_from_ymdhms(start_ymdhms, stop_ymdhms)


def test_initl(epoch, oe_params):
    # Inputs
    sgp4_obj = sgp4.SGP4(wgs_model=WGSModel.WGS_72)
    sgp4_obj.satrec.ecco, sgp4_obj.satrec.inclo, *_ = oe_params
    sgp4_obj.satrec.no_kozai = 0.00874808688806747

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
    sgp4_obj.initl(epoch)

    # Check results
    for key in expected:
        assert custom_isclose(getattr(sgp4_obj.sgp4init_out, key), expected[key])


@pytest.mark.parametrize(
    "typerun, startmfe_exp, stopmfe_exp, deltamin_exp",
    [
        (sgp4.TypeRun.Catalog, -1440, 1440, 20),
        (sgp4.TypeRun.Verification, 0, 4320, 360),
        (sgp4.TypeRun.FromJD, 12863952.7377735, 12867015.7877738, 1),
    ],
)
def test_twoline2rv(typerun, startmfe_exp, stopmfe_exp, deltamin_exp, monkeypatch):
    # TODO: Add more test cases
    # Patch SGP4 initialization method
    monkeypatch.setattr(sgp4.SGP4, "sgp4init", lambda *args: None)

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


def test_sgp4(oe_params, monkeypatch):
    # Patch SGP4 propagation method
    monkeypatch.setattr(sgp4.SGP4, "propagate", lambda *args: None)

    # Initialize SGP4 class
    sgp4_obj = sgp4.SGP4(wgs_model=WGSModel.WGS_72)

    # Inputs
    ecc, incl, node, argp, _, m = oe_params
    epoch = 20630.332154440228

    # Update `satrec` attributes
    sgp4_obj.satrec.satnum = 8195
    sgp4_obj.satrec.jdsatepoch = 2453911
    sgp4_obj.satrec.jdsatepochf = 0.8321544402
    sgp4_obj.satrec.no_kozai = 0.00874808688806747
    sgp4_obj.satrec.ecco = ecc
    sgp4_obj.satrec.inclo = incl
    sgp4_obj.satrec.nodeo = node
    sgp4_obj.satrec.argpo = argp
    sgp4_obj.satrec.mo = m
    sgp4_obj.satrec.bstar = 0.00011873
    sgp4_obj.satrec.ndot = 2.99978465186525e-12
    sgp4_obj.satrec.elnum = 813
    sgp4_obj.satrec.revnum = 22565
    sgp4_obj.satrec.classification = sgp4.Classification.Unclassified
    sgp4_obj.satrec.intldesg = "75081A"

    # Call method
    sgp4_obj.sgp4init(epoch)

    # Check results (OEs do not change)
    satrec_expected = {
        "ecco": ecc,
        "inclo": incl,
        "nodeo": node,
        "argpo": argp,
        "mo": m,
        "a": 4.165136398215488,
        "alta": 6.0295615102596924,
        "altp": 0.30071128617128307,
        "aycof": 0.0010552861263719983,
        "cc1": 1.9425158347208908e-13,
        "cc4": 1.2056623400361753e-09,
        "cc5": 2.0662701271427242e-07,
        "delmo": 6.3571512831144235,
        "eta": 0.9085028530232736,
        "argpdot": -7.437009265331242e-08,
        "omgcof": -3.510477365586775e-21,
        "sinmao": 0.3457191249346113,
        "t2cof": 2.9137737520813363e-13,
        "x1mth2": 0.8100072955438793,
        "x7thm1": 0.3299489311928454,
        "xlcof": 0.0019032757631354295,
        "xmcof": -2.4105286259805916e-15,
        "mdot": 0.008748086886633134,
        "nodedot": -1.2845672158012387e-06,
        "nodecf": -4.604617513547723e-19,
    }

    # Check results for default case
    assert sgp4_obj.use_deep_space
    for key in satrec_expected:
        assert custom_isclose(getattr(sgp4_obj.satrec, key), satrec_expected[key])

    # Check results for forced non-deep space case with mismatched inputs
    sgp4_obj.satrec.no_kozai = 1
    with pytest.raises(ValueError):
        sgp4_obj.sgp4init(epoch)
