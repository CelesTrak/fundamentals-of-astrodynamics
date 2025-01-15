import numpy as np
import pytest

import src.valladopy.astro.twobody.kepler as kepler

from ...conftest import DEFAULT_TOL, custom_isclose


@pytest.mark.parametrize(
    # TODO: find test case where np.abs(alpha) < SMALL
    "ro, vo, dtsec, r_expected, v_expected",
    [
        # alpha >= SMALL
        (
            [-6518.1083, -2403.8479, -22.1722],
            [2.604057, -7.105717, -0.263218],
            120,
            [-6150.8180772417, -3233.5647020192746, -53.47890911642744],
            [3.5085906394995146, -6.703227854017405, -0.2578170063354219],
        ),
        # alpha < SMALL
        (
            [15912.1061, 1352024060, 0],
            [-0.000783895533, 31889.3576, 0],
            100000,
            [15833.716546696181, 4540959819.99955, 0],
            [-0.000783895533049524, 31889.357599993506, 0],
        ),
        # alpha < SMALL, not converged
        (
            [15912.1061, 1352024060, 0],
            [-0.000783895533, 31889.3576, 0],
            5485549291.841488,
            [0, 0, 0],
            [0, 0, 0],
        ),
    ],
)
def test_kepler(ro, vo, dtsec, r_expected, v_expected):
    r, v = kepler.kepler(ro, vo, dtsec)
    assert np.allclose(r, r_expected, rtol=DEFAULT_TOL)
    assert np.allclose(v, v_expected, rtol=DEFAULT_TOL)


def test_pkepler():
    # TODO: test other orbit types
    # Input values
    ro = [-6518.1083, -2403.8479, -22.1722]
    vo = [2.604057, -7.105717, -0.263218]
    dtsec = 12345
    ndot, nddot = 5e-10, 1e-15

    # Expected values
    r_expected = [-1918.8942755619782, -6636.821817502705, -210.30326936057827]
    v_expected = [7.295164512561593, -2.115661226327176, -0.13016935223789428]

    # Compute the position and velocity vectors
    r, v = kepler.pkepler(ro, vo, dtsec, ndot, nddot)

    assert np.allclose(r, r_expected, rtol=DEFAULT_TOL)
    assert np.allclose(v, v_expected, rtol=DEFAULT_TOL)


def test_pkepler_bad_inputs(caplog):
    # Input values
    ro = [15912.1061, 1352024060, 0]
    vo = [-0.000783895533, 31889.3576, 0]

    # Compute the position and velocity vectors
    r, v = kepler.pkepler(ro, vo, 100000, 0, 0)

    # Check that the returned values are all zeros and that an error was logged
    assert np.allclose(r, [0, 0, 0], rtol=DEFAULT_TOL)
    assert np.allclose(v, [0, 0, 0], rtol=DEFAULT_TOL)
    assert caplog.records[0].levelname == "ERROR"
    assert caplog.records[0].message == "Negative semi-major axis encountered"


@pytest.mark.parametrize(
    "method, f_expected, g_expected, fdot_expected, gdot_expected",
    [
        (
            kepler.FGMethod.PQW,
            0.9977990947097023,
            59.95591110155602,
            -7.344342475157822e-05,
            0.9977958834883551,
        ),
        (
            kepler.FGMethod.SERIES,
            0.9977958909411967,
            59.955911263407685,
            -7.344342494983264e-05,
            0.9977958834764249,
        ),
        (
            kepler.FGMethod.C2C3,
            0.9999993072243452,
            59.99999974349895,
            -1.210288978296214e-06,
            0.9999993072219989,
        ),
    ],
)
def test_findfandg(method, f_expected, g_expected, fdot_expected, gdot_expected):
    # Input values
    r1 = [4938.49830042171, -1922.24810472241, 4384.68293292613]
    v1 = [0.738204644165659, 7.20989453238397, 2.32877392066299]
    r2 = [4971.8730437207, -1485.73546345938, 4514.6423760934]
    v2 = [0.373877325819783, 7.33517956912611, 2.00161489967163]
    x, z, c2, c3 = 0.1, 0.57483, 0.47650299902524496, 0.16194145738041812
    dtsec = 60

    # Compute f and g values
    f, g, fdot, gdot = kepler.findfandg(r1, v1, r2, v2, dtsec, x, z, c2, c3, method)

    # Expected values
    assert custom_isclose(f, f_expected)
    assert custom_isclose(g, g_expected)
    assert custom_isclose(fdot, fdot_expected)
    assert custom_isclose(gdot, gdot_expected)
