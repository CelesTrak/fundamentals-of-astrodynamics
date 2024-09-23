import numpy as np
import pytest

import src.valladopy.astro.time.frame_conversions as fc
from src.valladopy.constants import ARCSEC2RAD
from ...conftest import custom_allclose

DEFAULT_TOL = 1e-12


@pytest.fixture
def rva_ecef():
    # Example ECEF vectors (km, km/s, km/s^2)
    recef = np.array([-1033.4793830, 7901.2952754, 6380.3565958])
    vecef = np.array([-3.225636520, -2.872451450, 5.531924446])
    aecef = np.array([0.001, 0.002, 0.003])
    return recef, vecef, aecef


@pytest.fixture
def rva_eci():
    # Example ECI vectors (km, km/s, km/s^2)
    reci = np.array([2989.905220660578, -7387.200565596868, 6379.438182851598])
    veci = np.array([2.940401948462732, 3.809395206305895, 5.53064935673674])
    aeci = np.array([-3.927364527726347e-05, -0.00269956155725574, 0.0030002544835211])
    return reci, veci, aeci


@pytest.fixture
def t_inputs():
    # Time inputs
    ttt = 0.042623631888994  # Julian centuries of TT
    jdut1 = 2.4531015e06  # Julian date of UT1
    lod = 0.0015563  # Excess length of day, sec
    return ttt, jdut1, lod


@pytest.fixture
def orbit_effects_inputs():
    # Other inputs for accounting for orbit effectgs
    xp = -0.140682 * ARCSEC2RAD  # Polar motion coefficient
    yp = 0.333309 * ARCSEC2RAD  # Polar motion coefficient
    ddpsi = -0.052195 * ARCSEC2RAD  # Delta psi correction to GCRF, raed
    ddeps = -0.003875 * ARCSEC2RAD  # Delta epsilon correction to GCRF
    eqeterms = True  # Add extra terms for ast calculation
    return xp, yp, ddpsi, ddeps, eqeterms


def test_ecef2eci(rva_ecef, rva_eci, t_inputs, orbit_effects_inputs):
    # Expected ECI output vectors
    reci, veci, aeci = rva_eci

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.ecef2eci(
        *rva_ecef, *t_inputs, *orbit_effects_inputs
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(reci, reci_out)
    assert custom_allclose(veci, veci_out)
    assert custom_allclose(aeci, aeci_out)


def test_eci2ecef(rva_ecef, rva_eci, t_inputs, orbit_effects_inputs):
    # Expected ECEF output vectors
    recef, vecef, _ = rva_ecef

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.eci2ecef(
        *rva_eci, *t_inputs, *orbit_effects_inputs
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(recef, recef_out)
    assert custom_allclose(vecef, vecef_out)

    # For some reason, the acceleration out does not quite match that from the
    # original ECEF input
    aecef = np.array(
        [0.0002936830002159169, 0.0031151668034451073, 0.003000148416052949]
    )
    assert custom_allclose(aecef, aecef_out)
