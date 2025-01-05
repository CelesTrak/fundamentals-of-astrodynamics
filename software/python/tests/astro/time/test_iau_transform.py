import numpy as np
import pytest

import src.valladopy.astro.time.iau_transform as iau_transform
from src.valladopy.constants import ARCSEC2RAD
from ...conftest import custom_isclose, custom_allclose, DEFAULT_TOL


ROTATION_MATRIX_TOL = 1e-9


@pytest.fixture()
def ttt():
    """Terrestrial Time (TT) in Julian centuries of TT"""
    return 0.1


@pytest.fixture()
def delunay_elems():
    """Delunay elements in radians"""
    return [
        5.844239313494585,  # l
        6.23840254543787,  # l1
        3.0276889929096353,  # f
        3.2212027489393993,  # d
        5.089920270731961,  # omega
    ]


@pytest.fixture()
def planet_lon():
    """Planet longitudes in radians"""
    return [
        0.02422268041366862,  # mercury
        0.08339248169484563,  # venus
        0.030584726431970796,  # earth
        0.03334439906671072,  # mars
        0.1029125735528856,  # jupiter
        0.05248218685834288,  # saturn
        0.10871847648477687,  # uranus
        0.09936537545632083,  # neptune
    ]


@pytest.fixture()
def precrate():
    """Precession rate in radians per Julian century"""
    return 4.2555121682972836e-05


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


def test_iau06gst(ttt, delunay_elems, planet_lon, precrate, iau06arr):
    # Definitions
    judt1 = 2448855.009722  # Julian date of UT1
    deltapsi = -5.978331920752922e-05  # change in longitude

    # Call function
    gst, st = iau_transform.iau06gst(
        judt1, ttt, deltapsi, *delunay_elems, *planet_lon, precrate, iau06arr
    )

    # Check against expected values
    st_exp = np.array(
        [
            [-0.8894007799234658, -0.4571282671980926, 0],
            [0.4571282671980926, -0.8894007799234658, 0],
            [0, 0, 1],
        ]
    )
    assert custom_allclose(gst, 2.6668289843344684)
    assert custom_allclose(st, st_exp)


def test_iau06pna(ttt, iau06data_old):
    # Call function
    (
        deltapsi,
        pnb,
        prec,
        nut,
        l,
        l1,
        f,
        d,
        omega,
        lonmer,
        lonven,
        lonear,
        lonmar,
        lonjup,
        lonsat,
        lonurn,
        lonnep,
        precrate,
    ) = iau_transform.iau06pna(ttt, iau06data_old)

    # Check against expected values
    pnb_exp = np.array(
        [
            [0.9999968301767307, 0.0023093302141092642, 0.0010033097493433538],
            [-0.0023093439022892153, 0.9999973333838805, 1.2484743714435523e-05],
            [-0.0010032782425055283, -1.4801691391685055e-05, 0.9999994966067123],
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
            [0.9999999968199907, 7.31692338336742e-05, 3.172194366990678e-05],
            [-7.316966806263472e-05, 0.9999999972294239, 1.368765438014563e-05],
            [-3.172094206683474e-05, -1.3689975420727585e-05, 0.9999999994031832],
        ]
    )
    assert custom_allclose(pnb, pnb_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_allclose(prec, prec_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_allclose(nut, nut_exp, rtol=ROTATION_MATRIX_TOL)
    assert custom_isclose(deltapsi, 7.97497241593155e-05)
    assert custom_isclose(l, 5.844239313494585)
    assert custom_isclose(l1, 6.23840254543787)
    assert custom_isclose(f, 3.0276889929096353)
    assert custom_isclose(d, 3.2212027489393993)
    assert custom_isclose(omega, 5.089920270731961)
    assert custom_isclose(lonmer, 0.02422268041366862)
    assert custom_isclose(lonven, 0.08339248169484563)
    assert custom_isclose(lonear, 0.030584726431970796)
    assert custom_isclose(lonmar, 0.03334439906671072)
    assert custom_isclose(lonjup, 0.1029125735528856)
    assert custom_isclose(lonsat, 0.05248218685834288)
    assert custom_isclose(lonurn, 0.10871847648477687)
    assert custom_isclose(lonnep, 0.09936537545632083)
    assert custom_isclose(precrate, 4.2555121682972836e-05)


def test_iau06pnb(ttt, iau06data_old):
    # Call function
    (
        deltapsi,
        pnb,
        prec,
        nut,
        l,
        l1,
        f,
        d,
        omega,
        lonmer,
        lonven,
        lonear,
        lonmar,
        lonjup,
        lonsat,
        lonurn,
        lonnep,
        precrate,
    ) = iau_transform.iau06pnb(ttt, iau06data_old)

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
    assert custom_isclose(l, 5.844237767697117)
    assert custom_isclose(l1, 6.2384025722571055)
    assert custom_isclose(f, 3.02768961111037)
    assert custom_isclose(d, 3.2212030577628026)
    assert custom_isclose(omega, 5.08991990843217)
    assert custom_isclose(lonmer, 0.0)
    assert custom_isclose(lonven, 0.0)
    assert custom_isclose(lonear, 0.0)
    assert custom_isclose(lonmar, 0.0)
    assert custom_isclose(lonjup, 0.0)
    assert custom_isclose(lonsat, 0.0)
    assert custom_isclose(lonurn, 0.0)
    assert custom_isclose(lonnep, 0.0)
    assert custom_isclose(precrate, 0.0)


def test_iau06xys_series(ttt, delunay_elems, planet_lon, precrate, iau06arr):
    # Call function
    x, y, s = iau_transform.iau06xys_series(
        ttt, *delunay_elems, *planet_lon, precrate, iau06arr
    )

    # Check against expected values
    assert custom_isclose(x, 0.0010033097412569085)
    assert custom_isclose(y, 1.248474017229823e-05)
    assert custom_isclose(s, 7.765595105574456e-09)


def test_iau06xys(ttt, iau06arr):
    # Define EOP corrections
    ddx = -0.000205 * ARCSEC2RAD
    ddy = -0.000136 * ARCSEC2RAD

    # Call function
    x, y, s, nut = iau_transform.iau06xys(ttt, iau06arr, ddx, ddy)

    # Check against expected values
    nut_exp = np.array(
        [
            [0.999999496685652, 1.502895872980454e-09, 0.0010033087473888622],
            [-1.402829042902788e-08, 0.9999999999220739, 1.2484080825691922e-05],
            [-0.001003308747291916, -1.248408861698142e-05, 0.9999994966077259],
        ]
    )
    assert custom_isclose(x, 0.0010033087473888622)
    assert custom_isclose(y, 1.2484080825691922e-05)
    assert custom_isclose(s, 7.765595105574456e-09)
    assert custom_allclose(nut, nut_exp, rtol=ROTATION_MATRIX_TOL)
