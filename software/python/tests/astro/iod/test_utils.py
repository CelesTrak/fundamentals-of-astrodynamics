import numpy as np

import src.valladopy.constants as const
from src.valladopy.astro.iod.utils import gibbs, hgibbs

from ...conftest import DEFAULT_TOL, custom_isclose


def test_gibbs():
    # Vallado 2022, Example 7-3
    r1 = [0, 0, const.RE]
    r2 = [0, -4464.696, -5102.509]
    r3 = [0, 5740.323, 3189.068]

    # Call Gibbs method
    v2, theta12, theta23, copa = gibbs(r1, r2, r3)

    # Check results
    v2_expected = [0, 5.5311472050176125, -5.191806413494606]
    assert np.allclose(v2, v2_expected, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(theta12), 138.81407085944375, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(theta23), 160.24053069723146, rtol=DEFAULT_TOL)
    assert custom_isclose(float(np.degrees(copa)), 0)


def test_hgibbs():
    # Vallado 2022, Example 7-4
    r1 = [3419.85564, 6019.82602, 2784.60022]
    r2 = [2935.91195, 6326.18324, 2660.59584]
    r3 = [2434.95202, 6597.38674, 2521.52311]
    jd1 = 0
    jd2 = (60 + 16.48) / const.DAY2SEC
    jd3 = (120 + 33.04) / const.DAY2SEC

    # Call Herrick-Gibbs method
    v2, theta12, theta23, copa = hgibbs(r1, r2, r3, jd1, jd2, jd3)

    # Check results
    v2_expected = [-6.441557227511062, 3.777559606719521, -1.7205675602414345]
    assert np.allclose(v2, v2_expected, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(theta12), 4.499996147374992, rtol=DEFAULT_TOL)
    assert np.isclose(np.degrees(theta23), 4.499998402168982, rtol=DEFAULT_TOL)
    assert custom_isclose(float(np.degrees(copa)), -2.6712072741804056e-06)
