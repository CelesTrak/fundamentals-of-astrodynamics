import numpy as np
import pytest

import src.valladopy.astro.sgp4.sgp4 as sgp4
from src.valladopy.astro.sgp4.utils import WGSModel

from ...conftest import DEFAULT_TOL, custom_isclose, custom_allclose


@pytest.fixture
def satrec_coeffs_nonds():
    """Coefficients for the satellite record when not using deep space."""
    return {
        "d2": -2.9337430867525674e-24,
        "d3": 8.72425270068883e-36,
        "d4": -3.0009425550003753e-47,
        "t3cof": 4.083102888579937e-24,
        "t4cof": 6.486675392955838e-36,
        "t5cof": 6.711692282239235e-48,
    }


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


def test_sgp4init(oe_params, monkeypatch):
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

    # Expected values (OEs do not change)
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


def test_initialize_non_deep_space(satrec_coeffs_nonds):
    # Initialize SGP4 class and sgp4init_out
    sgp4_obj = sgp4.SGP4()
    sgp4_obj.sgp4init_out = sgp4.SGP4InitOutput()

    # Update some fields
    sgp4_obj.satrec.cc1 = 1.8730784787793202e-12
    sgp4_obj.sgp4init_out.ao = 0.1734465142323644

    # Call method
    sgp4_obj._initialize_non_deep_space(
        tsi=-1.2052706021039807, sfour=1.003135712869044
    )

    # Check results
    for key in satrec_coeffs_nonds:
        assert custom_isclose(getattr(sgp4_obj.satrec, key), satrec_coeffs_nonds[key])


def test_adjust_perigee():
    # Initialize SGP4 class and sgp4init_out
    sgp4_obj = sgp4.SGP4(wgs_model=WGSModel.WGS_72)
    sgp4_obj.sgp4init_out = sgp4.SGP4InitOutput()

    # Update periapsis radius
    sgp4_obj.sgp4init_out.rp = 0.054164814075659616

    # Call method
    sfour, qzms24 = sgp4_obj._adjust_perigee(
        ss=1.0122292801892716, qzms2t=1.880279159015271e-09
    )

    # Check results
    assert np.isclose(sfour, 1.003135712869044, rtol=DEFAULT_TOL)
    assert np.isclose(qzms24, 6.042618427427583e-08, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "use_deep_space, isimp, r_expected, v_expected",
    [
        (
            True,
            True,
            [15223.917136637867, -17852.958817081857, 25280.395582370667],
            [1.0790417322823393, 0.8751873723939548, 2.485682812729583],
        ),
        (
            False,
            False,
            [15258.747732593562, -17842.408055528413, 25305.446672098853],
            [1.0779326175214992, 0.8749235951487931, 2.482444210397105],
        ),
    ],
)
def test_propagate(
    ds,
    dscom_data,
    dsinit_data,
    satrec_coeffs_nonds,
    use_deep_space,
    isimp,
    r_expected,
    v_expected,
):
    # Initialize SGP4 class and sgp4init_out
    sgp4_obj = sgp4.SGP4(wgs_model=WGSModel.WGS_72)
    sgp4_obj.sgp4init_out = sgp4.SGP4InitOutput()

    # Set SGP4 object to use deep space
    sgp4_obj.ds = ds
    sgp4_obj.ds.dscom_out = dscom_data
    sgp4_obj.ds.dsinit_out = dsinit_data
    sgp4_obj.use_deep_space = use_deep_space

    # Update satellite record fields
    sgp4_obj.satrec.ecco = 0.6877146
    sgp4_obj.satrec.inclo = 1.11977881347003
    sgp4_obj.satrec.nodeo = 4.87072001413786
    sgp4_obj.satrec.argpo = 4.62102273937204
    sgp4_obj.satrec.no = 0.00874854701963024
    sgp4_obj.satrec.mo = 0.353005058520617
    sgp4_obj.satrec.nodedot = -1.28456721580123e-06
    sgp4_obj.satrec.argpdot = -7.43700926533354e-08
    sgp4_obj.satrec.mdot = 0.00874808688663313
    sgp4_obj.satrec.nodecf = -4.60461751354763e-19
    sgp4_obj.satrec.bstar = 0.00011873
    sgp4_obj.satrec.delmo = 6.35715128311442
    sgp4_obj.satrec.eta = 0.908502853023273
    sgp4_obj.satrec.sinmao = 0.345719124934611
    sgp4_obj.satrec.cc1 = 1.94251583472087e-13
    sgp4_obj.satrec.cc4 = 1.20566234003616e-09
    sgp4_obj.satrec.cc5 = 2.0662701271427e-07
    sgp4_obj.satrec.t2cof = 2.91377375208131e-13
    sgp4_obj.satrec.aycof = 0.001055286126372
    sgp4_obj.satrec.omgcof = -3.51047736558681e-21
    sgp4_obj.satrec.xlcof = 0.00190327576313543
    sgp4_obj.satrec.xmcof = -2.41052862598059e-15
    sgp4_obj.satrec.x1mth2 = 0.810007295543882
    sgp4_obj.satrec.x7thm1 = 0.329948931192823
    sgp4_obj.satrec.isimp = isimp

    # Update other fields
    ds.dsinit_out.argpm = 4.62101381496092
    ds.dsinit_out.nodem = 4.87056586607196
    ds.dsinit_out.mm = 1.40277548491659
    sgp4_obj.grav_const.xke = 0.0743669161331734
    sgp4_obj.sgp4init_out.gsto = 0.574180126924752
    sgp4_obj.sgp4init_out.con41 = -0.430021886631647

    # Add non-deep space fields to satellite record
    if not use_deep_space:
        for key in satrec_coeffs_nonds:
            setattr(sgp4_obj.satrec, key, satrec_coeffs_nonds[key])

    # Call method
    r, v = sgp4_obj.propagate(t=120)

    # Check results
    assert custom_allclose(r, r_expected)
    assert custom_allclose(v, v_expected)
