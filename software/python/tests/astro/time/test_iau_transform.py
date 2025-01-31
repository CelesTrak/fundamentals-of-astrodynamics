import os

import numpy as np
import pytest

import src.valladopy.astro.time.iau_transform as iau_transform
from src.valladopy.astro.time.utils import FundArgs
from src.valladopy.constants import ARCSEC2RAD
from ...conftest import custom_isclose, custom_allclose, DEFAULT_TOL


ROTATION_MATRIX_TOL = 1e-9


@pytest.fixture()
def ttt():
    """Terrestrial Time (TT) in Julian centuries of TT"""
    return 0.1


@pytest.fixture()
def fundargs():
    return FundArgs(
        l=5.844239313494585,
        l1=6.23840254543787,
        f=3.0276889929096353,
        d=3.2212027489393993,
        omega=5.089920270731961,
        lonmer=1.387857356197415,
        lonven=4.778037244236629,
        lonear=1.7523757421141397,
        lonmar=1.9104933369224852,
        lonjup=5.896456123410001,
        lonsat=3.0070078066000003,
        lonurn=6.22910985767,
        lonnep=5.69321664338,
        precrate=0.0024382288691000005,
    )


@pytest.fixture()
def jdtt_jdttf():
    return 2453101.5, 0.328154745474537


@pytest.fixture()
def xys_test_data(test_data_dir):
    # Load the data from the file
    file_path = os.path.join(test_data_dir, "xys_data.dat")
    data = np.genfromtxt(
        file_path, skip_header=1, names=["ttt", "jdtt", "jdftt", "x", "y", "s"]
    )
    return data


def mock_iau06xys_series_closure(xys_test_data):
    """Create a closure for the mock function with xys_test_data."""

    def mock_iau06xys_series(ttt, *args):
        # Find the closest row based on ttt
        closest_index = np.abs(xys_test_data["ttt"] - ttt).argmin()
        row = xys_test_data[closest_index]

        # Return x, y, s values
        return row["x"], row["y"], row["s"]

    return mock_iau06xys_series


def test_iau06era():
    # Expected values
    era_exp = np.array(
        [
            [-0.8884015255896265, -0.4590672383540609, 0],
            [0.4590672383540609, -0.8884015255896265, 0],
            [0, 0, 1],
        ]
    )

    # Call function
    era = iau_transform.iau06era(2448855.009722)

    # Check that they are the same
    assert np.allclose(era, era_exp, rtol=DEFAULT_TOL, atol=DEFAULT_TOL)


def test_iau06gst(ttt, fundargs, iau06arr):
    # Definitions
    judt1 = 2448855.009722  # Julian date of UT1
    deltapsi = -5.978331920752922e-05  # change in longitude

    # Call function
    gst, st = iau_transform.iau06gst(judt1, ttt, deltapsi, fundargs, iau06arr)

    # Check against expected values
    st_exp = np.array(
        [
            [-0.8894007799222967, -0.4571282672003672, 0],
            [0.4571282672003672, -0.8894007799222967, 0],
            [0, 0, 1],
        ]
    )
    assert custom_allclose(gst, 2.666828984331911)
    assert custom_allclose(st, st_exp)


def test_iau06pna(ttt, fundargs, iau06data_old):
    # Call function
    deltapsi, pnb, prec, nut, fundargs_out = iau_transform.iau06pna(ttt, iau06data_old)

    # Check against expected values
    pnb_exp = np.array(
        [
            [0.9999968301876362, 0.0023093262377770324, 0.0010033080322202625],
            [-0.002309339928413813, 0.9999973333930268, 1.248721984662418e-05],
            [-0.0010032765197276387, -1.4804159563682448e-05, 0.9999994966084041],
        ]
    )
    prec_exp = np.array(
        [
            [0.9999970278755955, 0.0022361037218361396, 0.0009716378547293649],
            [-0.0022361037097214137, 0.9999974999163708, -1.0988114751975253e-06],
            [-0.0009716378826099, -1.073874801980047e-06, 0.9999995279592244],
        ]
    )
    nut_exp = np.array(
        [
            [0.9999999968203361, 7.316525989077928e-05, 3.1720221013211795e-05],
            [-7.316569417457402e-05, 0.9999999972296809, 1.3690126666354718e-05],
            [-3.171921928366118e-05, -1.3692447454805023e-05, 0.999999999403204],
        ]
    )
    assert custom_allclose(pnb, pnb_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_allclose(prec, prec_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_allclose(nut, nut_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_isclose(deltapsi, 7.974539290449716e-05)
    assert custom_isclose(fundargs_out.l, fundargs.l)
    assert custom_isclose(fundargs_out.l1, fundargs.l1)
    assert custom_isclose(fundargs_out.f, fundargs.f)
    assert custom_isclose(fundargs_out.d, fundargs.d)
    assert custom_isclose(fundargs_out.omega, fundargs.omega)
    assert custom_isclose(fundargs_out.lonmer, fundargs.lonmer)
    assert custom_isclose(fundargs_out.lonven, fundargs.lonven)
    assert custom_isclose(fundargs_out.lonear, fundargs.lonear)
    assert custom_isclose(fundargs_out.lonmar, fundargs.lonmar)
    assert custom_isclose(fundargs_out.lonjup, fundargs.lonjup)
    assert custom_isclose(fundargs_out.lonsat, fundargs.lonsat)
    assert custom_isclose(fundargs_out.lonurn, fundargs.lonurn)
    assert custom_isclose(fundargs_out.lonnep, fundargs.lonnep)
    assert custom_isclose(fundargs_out.precrate, fundargs.precrate)


def test_iau06pnb(ttt, iau06data_old):
    # Call function
    deltapsi, pnb, prec, nut, fundargs = iau_transform.iau06pnb(ttt, iau06data_old)

    # Check against expected values
    pnb_exp = np.array(
        [
            [0.9999968301909938, 0.0023093250156609305, 0.0010033074986805486],
            [-0.0023093387049374255, 0.9999973333958689, 1.2485872490362154e-05],
            [-0.0010032759893189437, -1.4802809752234849e-05, 0.9999994966089566],
        ]
    )
    prec_exp = np.array(
        [
            [0.9999970277955893, 0.002236174504540456, 0.0009715572925536376],
            [-0.0022361744601857273, 0.9999974997581255, -1.1319406437726443e-06],
            [-0.0009715573946422197, -1.0406343248994258e-06, 0.9999995280374617],
        ]
    )
    nut_exp = np.array(
        [
            [0.9999999968204425, 7.316403646053429e-05, 3.171969048747087e-05],
            [-7.316447069429023e-05, 0.9999999972297888, 1.3688778120468648e-05],
            [-3.1718688873339126e-05, -1.3691098831258852e-05, 0.9999999994032394],
        ]
    )
    assert custom_allclose(pnb, pnb_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_allclose(prec, prec_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_allclose(nut, nut_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_isclose(deltapsi, 7.974405939816953e-05)
    assert custom_isclose(fundargs.l, 5.844237767697117)
    assert custom_isclose(fundargs.l1, 6.2384025722571055)
    assert custom_isclose(fundargs.f, 3.02768961111037)
    assert custom_isclose(fundargs.d, 3.2212030577628026)
    assert custom_isclose(fundargs.omega, 5.08991990843217)
    assert custom_isclose(fundargs.lonmer, 0)
    assert custom_isclose(fundargs.lonven, 0)
    assert custom_isclose(fundargs.lonear, 0)
    assert custom_isclose(fundargs.lonmar, 0)
    assert custom_isclose(fundargs.lonjup, 0)
    assert custom_isclose(fundargs.lonsat, 0)
    assert custom_isclose(fundargs.lonurn, 0)
    assert custom_isclose(fundargs.lonnep, 0)
    assert custom_isclose(fundargs.precrate, 0)


@pytest.mark.parametrize(
    "interp, x_exp, y_exp, s_exp",
    [
        (None, 0.000414038342, 4.7434075e-05, -8.602e-09),
        (
            iau_transform.InterpolationMode.LINEAR,
            0.00041412564757004343,
            4.74339906642304e-05,
            -8.603968928473071e-09,
        ),
        (
            iau_transform.InterpolationMode.SPLINE,
            0.00041412564757004343,
            4.74339906642304e-05,
            -8.603907496739735e-09,
        ),
    ],
)
def test_findxysparam(iau06xysarr, jdtt_jdttf, interp, x_exp, y_exp, s_exp):
    # Call function
    x, y, s = iau_transform.findxysparam(*jdtt_jdttf, iau06xysarr, interp=interp)

    # Check against expected values
    assert custom_isclose(x, x_exp)
    assert custom_isclose(y, y_exp)
    assert custom_isclose(s, s_exp)


def test_create_xys(iau06arr, xys_test_data, monkeypatch):
    # Mock the iau06xys_series function
    mock_function = mock_iau06xys_series_closure(xys_test_data)
    monkeypatch.setattr(iau_transform, "iau06xys_series", mock_function)

    # Call function
    xys_data = iau_transform.create_xys(iau06arr)

    # Check against expected values
    assert custom_allclose(
        xys_data[0], [2435839.5, 0, -0.004148973464, -4.713574e-05, -6.7608e-08]
    )
    assert custom_allclose(
        xys_data[-1], [2436204.5, 0, -0.004060043156, -5.7176757e-05, -9.0369e-08]
    )


def test_iau06xys_series(ttt, iau06arr):
    # Call function
    x, y, s = iau_transform.iau06xys_series(ttt, iau06arr)

    # Check against expected values
    assert custom_isclose(x, 0.001003308023544622)
    assert custom_isclose(y, 1.2487218189628587e-05)
    assert custom_isclose(s, 7.766920280134978e-09)


@pytest.mark.parametrize(
    "use_full_series, x_exp, y_exp, s_exp, pn_exp",
    [
        (
            True,
            0.0010033070296765758,
            1.2486558843022279e-05,
            7.766920280134978e-09,
            [
                [0.9999994966873754, 1.5029886625084027e-09, 0.0010033070296765758],
                [-1.4030847987967036e-08, 0.9999999999220429, 1.2486558843022279e-05],
                [-0.0010033070295795937, -1.2486566635627995e-05, 0.9999994966094183],
            ],
        ),
        (
            False,
            0.0009716708951319538,
            4.653719865339369e-05,
            -2.0661e-08,
            [
                [0.9999995279277234, -4.327041633163345e-08, 0.0009716708951319538],
                [-1.948426107491633e-09, 0.9999999989171445, 4.653719865339369e-05],
                [-0.0009716708960934587, -4.653717857770131e-05, 0.9999995268448684],
            ],
        ),
    ],
)
def test_iau06xys(
    ttt, iau06arr, iau06xysarr, use_full_series, x_exp, y_exp, s_exp, pn_exp
):
    # Define EOP corrections
    ddx, ddy = -0.000205 * ARCSEC2RAD, -0.000136 * ARCSEC2RAD

    # Call function
    x, y, s, pn = iau_transform.iau06xys(
        ttt, iau06arr, ddx, ddy, iau06xysarr, use_full_series
    )

    assert custom_isclose(x, x_exp)
    assert custom_isclose(y, y_exp)
    assert custom_isclose(s, s_exp)
    assert custom_allclose(pn, np.array(pn_exp), rtol=ROTATION_MATRIX_TOL)


def test_iau06xys_bad(ttt, iau06arr):
    """Check for raised error when full series is not used without the XYS data."""
    with pytest.raises(ValueError):
        iau_transform.iau06xys(ttt, iau06arr, use_full_series=False)


@pytest.mark.parametrize(
    "interp, eop_params_exp",
    # eop_params_exp in tuple (dut1, dat, lod, xp, yp, ddpsi, ddeps, dx, dy)
    [
        (
            None,
            (
                -0.4399623,
                32,
                0.0015308,
                -6.829327918949478e-07,
                1.6164366867345482e-06,
                -2.532763632852438e-07,
                -1.9154988540637765e-08,
                -5.914726909536339e-10,
                -5.187506387872035e-10,
            ),
        ),
        (
            iau_transform.InterpolationMode.LINEAR,
            (
                -0.4404350725418587,
                32,
                0.0014723556398243643,
                -6.817554969598743e-07,
                1.6210488191896547e-06,
                -2.541195610090128e-07,
                -1.9110442245797138e-08,
                -5.994273864608887e-10,
                -4.630677702364201e-10,
            ),
        ),
        (
            iau_transform.InterpolationMode.SPLINE,
            (
                -0.4404508444133234,
                32,
                0.0014873148833432495,
                -6.817301763697862e-07,
                1.6207593674283049e-06,
                -2.5409384906011997e-07,
                -1.9005342733934667e-08,
                -6.29147015448284e-10,
                -4.747972687881268e-10,
            ),
        ),
    ],
)
def test_findeopparam(eoparr, jdtt_jdttf, interp, eop_params_exp):
    eop_params = iau_transform.findeopparam(*jdtt_jdttf, eoparr, interp=interp)
    for i, param in enumerate(eop_params):
        if i == 1:  # dat
            assert int(param) == eop_params_exp[i]
        else:
            assert custom_isclose(param, eop_params_exp[i])
