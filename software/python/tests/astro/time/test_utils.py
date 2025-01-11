import numpy as np
import pytest

import src.valladopy.astro.time.utils as utils
from src.valladopy.constants import ARCSEC2RAD
from ...conftest import custom_isclose, custom_allclose

DEFAULT_TOL = 1e-12


@pytest.mark.parametrize(
    "ttt, opt, fundargs_expected",
    [
        (
            0.1,
            "06",
            {
                "l": 5.844239313494585,
                "l1": 6.23840254543787,
                "f": 3.0276889929096353,
                "d": 3.2212027489393993,
                "omega": 5.089920270731961,
                "lonmer": 1.387857356197415,
                "lonven": 4.778037244236629,
                "lonear": 1.7523757421141397,
                "lonmar": 1.9104933369224852,
                "lonjup": 5.896456123410001,
                "lonsat": 3.0070078066000003,
                "lonurn": 6.22910985767,
                "lonnep": 5.69321664338,
                "precrate": 0.0024382288691000005,
            },
        ),
        (
            1.3,
            "96",
            {
                "l": 3.7263866519539266,
                "l1": 6.218507382213416,
                "f": 0.9754435804880215,
                "d": 4.626729788935919,
                "omega": 2.2806396540189735,
                "lonmer": 0,
                "lonven": 5.151167862645337,
                "lonear": 1.739240779528124,
                "lonmar": 0.6601251108969298,
                "lonjup": 0.34433358093148264,
                "lonsat": 3.4701579303306653,
                "lonurn": 0,
                "lonnep": 0,
                "precrate": 0.03170537748539087,
            },
        ),
        (
            -0.15,
            "02",
            {
                "l": 3.4057180143687678,
                "l1": 6.242546458897546,
                "f": 5.8114135943580365,
                "d": 1.8811765778266971,
                "omega": 0.9628107824807162,
                "lonmer": 0,
                "lonven": 0,
                "lonear": 0,
                "lonmar": 0,
                "lonjup": 0,
                "lonsat": 0,
                "lonurn": 0,
                "lonnep": 0,
                "precrate": 0,
            },
        ),
        (
            -1.1,
            "80",
            {
                "l": 1.679344348277534,
                "l1": 6.258264795030512,
                "f": 5.0797486481605745,
                "d": 1.8155859432858112,
                "omega": 1.6161209896111113,
                "lonmer": 6.162757588791571,
                "lonven": 4.368908183592073,
                "lonear": 1.7289231570255421,
                "lonmar": 3.1236157622942082,
                "lonjup": 5.129697204536531,
                "lonsat": 2.5078636021906533,
                "lonurn": 0,
                "lonnep": 0,
                "precrate": 0,
            },
        ),
    ],
)
def test_fundarg(ttt, opt, fundargs_expected):
    """
    Some values are a different representation of the MATLAB outputs due to
    `np.mod` always returning a positive angle (unlike the `rem` function) but
    they are still equivalent
    """
    fundargs = utils.fundarg(ttt, opt).__dict__
    for farg in fundargs_expected:
        assert custom_isclose(fundargs[farg], fundargs_expected[farg])


@pytest.mark.parametrize(
    "ttt, opt, prec, psia, wa, ea, xa",
    [
        (
            0.1,
            "80",
            np.array(
                [
                    [0.99999702751987, 0.00223623766396228, 9.716956978347777e-04],
                    [-0.0022362376639804, 0.99999749961684, -1.08645419889999e-06],
                    [-9.7169569779302e-04, -1.0864915479894e-06, 0.999999527903034],
                ]
            ),
            0.0024428166982077833,
            0.409092806670512,
            0.40907010765003343,
            5.0006227129429165e-06,
        ),
        (
            3.4,
            "50",
            np.array(
                [
                    [0.99656203288, 0.07597927776, 0.03297579896],
                    [-0.07597927776, 0.9971069944000001, -0.0012551847999999997],
                    [-0.03297579896, -0.0012551847999999997, 0.9994542524],
                ]
            ),
            0.03801516768171585,
            0.03300637882588775,
            0.03805971221205611,
            2.9466975537837596e-07,
        ),
        (
            -0.2,
            "06",
            np.array(
                [
                    [0.999988113087504, -0.0044717893242489, -0.00194339495019456],
                    [0.00447178937247984, 0.999990001490586, -4.320434382393e-06],
                    [0.0019433948392142, -4.370069859043412e-06, 0.999998111596918],
                ]
            ),
            -0.004885658734881702,
            0.4090926358130436,
            0.409138014700149,
            -1.0697548647819031e-05,
        ),
    ],
)
def test_precess(ttt, opt, prec, psia, wa, ea, xa):
    prec_out, psia_out, wa_out, ea_out, xa_out = utils.precess(ttt, opt)
    assert np.allclose(prec, prec_out)
    assert custom_isclose(psia, psia_out)
    assert custom_isclose(wa, wa_out)
    assert custom_isclose(ea, ea_out)
    assert custom_isclose(xa, xa_out)


def test_precess_bad():
    # Test invalid option
    with pytest.raises(ValueError):
        _ = utils.precess(0.5, "25")


def test_nutation(iau80arr):
    # Inputs
    ttt = 0.042623631888994
    ddpsi = -0.052195 * ARCSEC2RAD
    ddeps = -0.003875 * ARCSEC2RAD

    # Expected nutation transformation matrix
    nut_exp = np.array(
        [
            [0.999999998212977, -5.484951020917109e-05, -2.37818519260497e-05],
            [5.4850353003009e-05, 0.9999999978677465, 3.5439324690356155e-05],
            [2.377990804573938e-05, -3.544062907001599e-05, 0.9999999990892389],
        ]
    )

    # Call the function with test inputs
    deltapsi, trueeps, meaneps, omega, nut = utils.nutation(ttt, ddpsi, ddeps, iau80arr)

    # Check if the outputs are close to the expected values
    assert custom_isclose(deltapsi, -5.978331920752922e-05)
    assert custom_isclose(trueeps, 0.4091185700997511)
    assert custom_isclose(meaneps, 0.40908313012283176)
    assert custom_isclose(omega, 0.7435907904484494)
    assert custom_allclose(nut, nut_exp)


@pytest.mark.parametrize(
    "use_eutelsat_approx, nut_expected",
    [
        (
            False,
            [
                [0.9999999982280666, -5.4617220076710384e-05, -2.3681763640372747e-05],
                [5.4618059767382264e-05, 0.9999999978797962, 3.545807249405453e-05],
                [2.367982696881368e-05, -3.545936588322807e-05, 0.9999999990909495],
            ],
        ),
        (
            True,
            [
                [0.9999999997196329, 8.396563360840841e-10, -2.3679826965700458e-05],
                [0, 0.9999999993713397, 3.545871922005356e-05],
                [2.3679826980587025e-05, -3.5458719210112105e-05, 0.9999999990909726],
            ],
        ),
    ],
)
def test_nutation_qmod(iau80arr, use_eutelsat_approx, nut_expected):
    # Inputs
    ttt = 0.042623631888994

    # Call the function with test inputs
    deltapsi, trueeps, meaneps, omega, nut = utils.nutation_qmod(
        ttt, iau80arr, use_eutelsat_approx
    )

    # Check if the outputs are close to the expected values
    assert custom_isclose(deltapsi, -5.953038436137153e-05)
    assert custom_isclose(trueeps, 0.40912826294155646)
    assert custom_isclose(meaneps, 0.40909280422232897)
    assert custom_isclose(omega, 0.7435914192041119)
    assert custom_allclose(nut, np.array(nut_expected))


@pytest.mark.parametrize(
    "xp, yp, ttt, use_iau80, pm_expected",
    [
        (
            0.1,
            0.2,
            0.042623631889,
            True,
            np.array(
                [
                    [0.9950041652780257, 0, -0.09983341664682815],
                    [0.019833838076209875, 0.9800665778412416, 0.19767681165408385],
                    [0.09784339500725571, -0.19866933079506122, 0.9751703272018158],
                ]
            ),
        ),
        (
            0.1,
            0.2,
            0.042623631889,
            False,
            np.array(
                [
                    [0.9950041652780257, 0.0198338380857286, -0.09784339500532617],
                    [-9.663803175648039e-12, 0.980066577841049, 0.1986693307960115],
                    [0.09983341664682815, -0.19767681165408385, 0.9751703272018158],
                ]
            ),
        ),
        (0, 0, 0.042623631889, True, np.identity(3)),
        (
            0,
            0,
            0.042623631889,
            False,
            np.array(
                [[1, 9.71232434283103e-12, 0], [-9.71232434283103e-12, 1, 0], [0, 0, 1]]
            ),
        ),
    ],
)
def test_polarm(xp, yp, ttt, use_iau80, pm_expected):
    pm = utils.polarm(xp, yp, ttt, use_iau80)
    assert custom_allclose(pm, pm_expected)


@pytest.mark.parametrize(
    "kp, ap", [(0, 0), (3, 15), (9, 400), (2.5, 10.437984357072379)]
)
def test_kp_ap_conversions(kp, ap):
    assert custom_isclose(utils.kp2ap(kp), ap)
    assert custom_isclose(utils.ap2kp(ap), kp)


def test_kp_ap_conversions_bad():
    assert utils.kp2ap(-0.5) is None
    assert utils.ap2kp(900) is None
