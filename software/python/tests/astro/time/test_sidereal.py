import numpy as np

from src.valladopy.astro.time.sidereal import gstime, sidereal
from ...conftest import custom_isclose, custom_allclose


DEFAULT_TOL = 1e-12


def test_gstime():
    # Test Julian Date of J2000 epoch
    assert custom_isclose(gstime(2451545.0), np.radians(280.460618375))

    # Vallado, Example 3-5 (2013)
    assert custom_isclose(
        gstime(2448855.009722), np.radians(152.57870762835), rtol=DEFAULT_TOL
    )


def test_sidereal():
    # Definitions
    jdut1 = 2453101.82740678310  # Julian date of UT1
    deltapsi = 6.230930776908879  # Nutation angle, rad
    meaneps = 0.40908313012283176  # Mean obliquity of the ecliptic, rad
    omega = 0.7435907904482468  # Moon longitude of AN, rad
    lod = 0.001556  # Length of day, sec
    eqeterms = True  # Add terms for ast calculation

    # Compute transformation matrix and matrix rate
    st, stdot = sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    # Compare matrices
    st_exp = np.array(
        [
            [0.17979642915209632, 0.9837038396103551, 0.0],
            [-0.9837038396103551, 0.17979642915209632, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    stdot_exp = np.array(
        [
            [7.173281395194778e-05, -1.3110962143545667e-05, 0.0],
            [1.3110962143545667e-05, 7.173281395194778e-05, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    assert custom_allclose(st, st_exp, rtol=DEFAULT_TOL)
    assert custom_allclose(stdot, stdot_exp, rtol=DEFAULT_TOL)
