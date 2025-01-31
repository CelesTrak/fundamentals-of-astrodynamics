import os
import pytest

import numpy as np

import src.valladopy.astro.time.data as data
from src.valladopy.astro.od.diffcorrection import findatwaatwb

from ...conftest import DEFAULT_TOL, load_matlab_data


@pytest.fixture()
def iau80arr():
    """Load the IAU 1980 data"""
    return data.iau80in()


@pytest.fixture()
def obsrecarr(test_data_dir):
    struct_name = "obsrecarr_array"
    file_path = os.path.join(test_data_dir, "obsrecarr_struct.mat")
    return load_matlab_data(file_path, keys=[struct_name])[struct_name]


def test_findatwaatwb(iau80arr, obsrecarr):
    # Inputs (Vallado 2022, Example 10-4)
    # TODO: Test cases where obstype != 2
    firstobs, lastobs = 1, 10
    percentchg = deltaamtchg = 0.01
    xnom = [5975.2904, 2568.64, 3120.5845, 3.983846, -2.071159, -5.917095]  # km, km/s

    # Update keys
    obsrecarr.noise_rng = obsrecarr.noiserng
    obsrecarr.noise_az = obsrecarr.noiseaz
    obsrecarr.noise_el = obsrecarr.noiseel

    # Expected results
    # TODO: This doesn't quite match MATLAB results - need to investigate differences
    # fmt: off
    atwa_exp = np.array(
        [
            [136.80298278376176, 240.4270437000345, 73.73103119321861,
             9321.489531070378, 16864.492360822725, 3947.3446276205896],
            [240.42704370003452, 1006.2675602784025, 300.20096730776476,
             15677.86347349234, 56210.7953424631, 11706.150843524285],
            [73.73103119321861, 300.2009673077647, 145.40111224999106,
             3930.117946336082, 12460.548674625494, 5638.488343409961],
            [9321.489531070378, 15677.86347349234, 3930.1179463360822,
             771118.3881618785, 1329485.5098113604, 277196.09059560933],
            [16864.492360822725, 56210.7953424631, 12460.548674625494,
             1329485.5098113604, 4361530.689572721, 753241.0605554599],
            [3947.3446276205896, 11706.150843524285, 5638.488343409961,
             277196.09059560933, 753241.0605554599, 362539.3495064649]
        ]
    )
    atwb_exp = np.array(
        [
            [38720.582925101255],
            [215578.1737859935],
            [79313.40690559476],
            [2282667.0524854274],
            [10508125.9772133],
            [2972546.3660315783]
        ]
    )
    atw_exp = np.array(
        [
            [46.43622813245533, 2046.8140647745468, 10004.162211833363],
            [108.0787983583951, -372.6380776727631, -4301.827857133288],
            [11.49080512614076, -4936.397149585619, 4502.476819242878],
            [4821.961994557116, 225428.55787438713, 1113812.5280228828],
            [11883.310125538548, -41879.880488808216, -488234.25549887645],
            [1055.1987643446184, -543083.8082404276, 502369.8912974823],
        ]
    )
    bt_exp = np.array(
        [
            [129.43510400356035],
            [-0.29027620242991525],
            [-0.06843658979127826]
        ]
    )
    # fmt: on

    # Call function
    atwa, atwb, atw, b, drng2, daz2, del2 = findatwaatwb(
        iau80arr, firstobs, lastobs, obsrecarr, percentchg, deltaamtchg, xnom
    )
    assert np.allclose(atwa, atwa_exp, rtol=DEFAULT_TOL)
    assert np.allclose(atwb, atwb_exp, rtol=DEFAULT_TOL)
    assert np.allclose(atw, atw_exp, rtol=DEFAULT_TOL)
    assert np.allclose(b, bt_exp, rtol=DEFAULT_TOL)
    assert np.isclose(drng2, 48157540.09644246, atol=DEFAULT_TOL)
    assert np.isclose(daz2, 3209189.821954687, atol=DEFAULT_TOL)
    assert np.isclose(del2, 1162318.2830344548, atol=DEFAULT_TOL)
