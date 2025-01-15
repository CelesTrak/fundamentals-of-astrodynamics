import numpy as np
import pytest

from src.valladopy.astro.time.sidereal import gstime, gstime0, lstime, sidereal
from ...conftest import custom_isclose, custom_allclose


@pytest.mark.parametrize(
    "jdut1, lon, gstime_exp, lstime_exp",
    [
        # J2000 epoch
        (2451545, np.radians(5), np.radians(280.460618375), np.radians(285.460618375)),
        # Vallado, Example 3-5 (2013)
        (
            2448855.009722,
            np.radians(-104),
            np.radians(152.57870762835046),
            np.radians(48.57870762835046),
        ),
    ],
)
def test_gstime_lstime(jdut1, lon, gstime_exp, lstime_exp):
    assert custom_isclose(gstime(jdut1), float(gstime_exp))
    assert custom_allclose(lstime(lon, jdut1), [float(lstime_exp), float(gstime_exp)])


def test_gstime0():
    assert custom_isclose(gstime0(1989), 1.7561909422996962)


@pytest.mark.parametrize(
    "use_iau80, st_exp, stdot_exp",
    [
        (
            True,
            np.array(
                [
                    [0.17979642915209632, 0.9837038396103551, 0],
                    [-0.9837038396103551, 0.17979642915209632, 0],
                    [0, 0, 1],
                ]
            ),
            np.array(
                [
                    [7.173281395194778e-05, -1.3110962143545667e-05, 0],
                    [1.3110962143545667e-05, 7.173281395194778e-05, 0],
                    [0, 0, 0],
                ]
            ),
        ),
        (
            False,
            np.array(
                [
                    [0.678868449558204, 0.7342599186898602, 0],
                    [-0.7342599186898602, 0.678868449558204, 0],
                    [0, 0, 1],
                ]
            ),
            np.array(
                [
                    [5.354307670549993e-05, -4.950386714897323e-05, 0],
                    [4.950386714897323e-05, 5.354307670549993e-05, 0],
                    [0, 0, 0],
                ]
            ),
        ),
    ],
)
def test_sidereal(use_iau80, st_exp, stdot_exp):
    # Definitions
    jdut1 = 2453101.8274067831  # Julian date of UT1
    deltapsi = 6.230930776908879  # nutation angle, rad
    meaneps = 0.40908313012283176  # mean obliquity of the ecliptic, rad
    omega = 0.7435907904482468  # moon longitude of AN, rad
    lod = 0.001556  # length of day, sec

    # Compute transformation matrix and matrix rate
    st, stdot = sidereal(jdut1, deltapsi, meaneps, omega, lod, use_iau80)

    # Compare results
    assert custom_allclose(st, st_exp)
    assert custom_allclose(stdot, stdot_exp)
