import numpy as np
import pytest

import src.valladopy.astro.iod.angles as angles

from ...conftest import DEFAULT_TOL


@pytest.fixture
def obs_data():
    """Observation data for three sightings (Vallado 2022, Example 7-2)"""
    # Right Ascension (radians)
    rtasc1, rtasc2, rtasc3 = np.radians([0.939913, 45.025748, 67.886655])

    # Declination (radians)
    decl1, decl2, decl3 = np.radians([18.667717, 35.664741, 36.996583])

    # Julian dates (integer and fractional parts)
    jd1 = jd2 = jd3 = 2456159.5
    jdf1 = 0.4864351851851852
    jdf2 = 0.49199074074074073
    jdf3 = 0.4947685185185185

    return decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3


@pytest.fixture
def site_data():
    """Site position vectors for three sightings (Vallado 2022, Example 7-2)"""
    # Site position vectors (ECI, in km)
    rseci1 = [4054.881, 2748.195, 4074.237]
    rseci2 = [3956.224, 2888.232, 4074.364]
    rseci3 = [3905.073, 2956.935, 4074.430]

    return rseci1, rseci2, rseci3


@pytest.fixture
def range_data():
    range1, range2 = 12742.21211926773, 12997.05636165308
    return range1, range2


@pytest.fixture
def time_intervals():
    tau12, tau13, tau32 = -480, -720, 240
    return tau12, tau13, tau32


@pytest.fixture
def los_vectors():
    los1 = [0.9472633016837996, 0.015540847420215897, 0.3200792391651927]
    los2 = [0.5742253604222706, 0.5747416917454863, 0.5830413563525745]
    los3 = [0.3006519052747138, 0.739921912374617, 0.6017673931367287]
    return los1, los2, los3


def test_calculate_time_intervals(obs_data, time_intervals):
    # Input data
    *_, jd1, jdf1, jd2, jdf2, jd3, jdf3 = obs_data

    # Calculate time intervals
    tau12, tau13, tau32 = angles.calculate_time_intervals(
        jd1, jdf1, jd2, jdf2, jd3, jdf3
    )

    # Check results
    tau12_expected, tau13_expected, tau32_expected = time_intervals
    assert np.isclose(tau12, tau12_expected, rtol=DEFAULT_TOL)
    assert np.isclose(tau13, tau13_expected, rtol=DEFAULT_TOL)
    assert np.isclose(tau32, tau32_expected, rtol=DEFAULT_TOL)


def test_calculation_los_vectors(obs_data, los_vectors):
    # Input data
    decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, *_ = obs_data

    # Calculate line-of-sight unit vectors
    los1, los2, los3 = angles.calculate_los_vectors(
        [decl1, decl2, decl3], [rtasc1, rtasc2, rtasc3]
    )

    # Check results
    los1_expected, los2_expected, los3_expected = los_vectors
    assert np.allclose(los1, los1_expected, rtol=DEFAULT_TOL)
    assert np.allclose(los2, los2_expected, rtol=DEFAULT_TOL)
    assert np.allclose(los3, los3_expected, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "diffsites, reci_expected, veci_expected",
    [
        (
            True,
            [6375.656382985365, 5309.839885514231, 6530.941565891183],
            [22.32079637714886, 26.24161577925484, 28.464611900869997],
        ),
        (
            False,
            [6375.826207838072, 5310.009863069861, 6531.113998039411],
            [5.836959439606365, 14.208716188023974, 11.488650926772648],
        ),
    ],
)
def test_angles_laplace(obs_data, site_data, diffsites, reci_expected, veci_expected):
    # Calculate position and velocity vectors
    reci, veci = angles.laplace(*obs_data, diffsites, *site_data)

    # Expected results
    assert np.allclose(reci, reci_expected, rtol=DEFAULT_TOL)
    assert np.allclose(veci, veci_expected, rtol=DEFAULT_TOL)


def test_angles_gauss(obs_data, site_data):
    # Expected results
    reci_expected = [6313.378130210396, 5247.50563344895, 6467.707164431651]
    veci_expected = [-4.185488280436629, 4.7884929168898145, 1.721714659663034]

    # Calculate position and velocity vectors
    reci, veci = angles.gauss(*obs_data, *site_data)

    # Expected results
    assert np.allclose(reci, reci_expected, rtol=DEFAULT_TOL)
    assert np.allclose(veci, veci_expected, rtol=DEFAULT_TOL)


def test_doubler_iter(site_data, range_data, time_intervals, los_vectors):
    # Input data
    tau12, _, tau32 = time_intervals
    n12, n13, n23 = 0.12, 1.23, 2.34

    # Calculate outputs
    r2, r3, f1, f2, q1, magr1, magr2, a, deltae32 = angles.doubler_iter(
        *range_data, *los_vectors, *site_data, tau12, tau32, n12, n13, n23
    )

    # Expected results
    r2_expected = [7779.700168197198, 6715.14615726086, 7956.5414249888745]
    r3_expected = [6027.317367172247, 8179.90243554486, 8322.191075272185]
    assert np.allclose(r2, r2_expected, rtol=DEFAULT_TOL)
    assert np.allclose(r3, r3_expected, rtol=DEFAULT_TOL)
    assert np.isclose(f1, 2167.210768046717, rtol=DEFAULT_TOL)
    assert np.isclose(f2, -36502.69736491624, rtol=DEFAULT_TOL)
    assert np.isclose(q1, 36566.97577634498, rtol=DEFAULT_TOL)
    assert np.isclose(magr1, 12742.21211926773, rtol=DEFAULT_TOL)
    assert np.isclose(magr2, 12997.05636165308, rtol=DEFAULT_TOL)
    assert np.isclose(a, 13404.16181471764, rtol=DEFAULT_TOL)
    assert np.isclose(deltae32, 14.875483923982403, rtol=DEFAULT_TOL)


def test_angles_doubler(obs_data, site_data, range_data):
    # Expected results
    reci_expected = [6356.3974225567745, 5290.563607834945, 6511.386925528436]
    veci_expected = [-4.172963733537418, 4.77645053907524, 1.7201932997557194]

    # Calculate position and velocity vectors
    reci, veci = angles.doubler(*obs_data, *site_data, *range_data)

    # Expected results
    assert np.allclose(reci, reci_expected, rtol=DEFAULT_TOL)
    assert np.allclose(veci, veci_expected, rtol=DEFAULT_TOL)
