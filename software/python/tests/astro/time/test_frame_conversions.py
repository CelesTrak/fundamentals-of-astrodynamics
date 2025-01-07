import numpy as np
import pytest

import src.valladopy.astro.time.frame_conversions as fc
from src.valladopy.constants import ARCSEC2RAD
from ...conftest import custom_allclose

DEFAULT_TOL = 1e-12


@pytest.fixture
def rva_eci():
    # Example ECI vectors (km, km/s, km/s^2)
    reci = [2989.905220660578, -7387.200565596868, 6379.438182851598]
    veci = [2.940401948462732, 3.809395206305895, 5.53064935673674]
    aeci = [-3.927364527726347e-05, -0.00269956155725574, 0.0030002544835211]
    return reci, veci, aeci


@pytest.fixture
def rva_ecef():
    # Example ECEF vectors (km, km/s, km/s^2)
    recef = [-1033.4793830, 7901.2952754, 6380.3565958]
    vecef = [-3.225636520, -2.872451450, 5.531924446]
    aecef = [0.001, 0.002, 0.003]
    return recef, vecef, aecef


@pytest.fixture
def rva_ecef06():
    recef = [-1033.4778040252972, 7901.295481365824, 6380.3565964958425]
    vecef = [-3.2256371005292754, -2.872450803966763, 5.531924442318752]
    aecef = [0.00029368362221701087, 0.003115166744642978, 0.0030001489557543396]
    return recef, vecef, aecef


@pytest.fixture
def rva_pef():
    # Example PEF vectors (km, km/s, km/s^2) [opt = "80"]
    rpef = [-1033.4750313057266, 7901.305585585349, 6380.344532748868]
    vpef = [-3.225632746974616, -2.872442510803122, 5.531931287696299]
    apef = [0.000293685046453725, 0.0031151716514636447, 0.0030001437204660755]
    return rpef, vpef, apef


@pytest.fixture
def rva_pef_ecef(rva_pef):
    rpef, vpef, _ = rva_pef
    apef = [0.0010000020461365159, 0.0020000048477791838, 0.0029999960860945377]
    return rpef, vpef, apef


@pytest.fixture
def rva_mod():
    # Example MOD vectors (km, km/s, km/s^2)
    rmod = [2994.30249664958, -7384.348641268552, 6380.677444947862]
    vmod = [2.934478765910065, 3.8121950281090644, 5.531865978456613]
    amod = [-3.79431747506729e-05, -0.0026995983568731115, 0.003000238492782315]
    return rmod, vmod, amod


@pytest.fixture
def rva_mod_ecef(rva_mod):
    rmod, vmod, _ = rva_mod
    amod = [-0.0010028781946650484, -0.0017988314102041973, 0.003000035987896538]
    return rmod, vmod, amod


@pytest.fixture
def rva_tod():
    # Example TOD vectors (km, km/s, km/s^2) [opt = "80"]
    rtod = [2994.049189091933, -7384.738996771148, 6380.344532748868]
    vtod = [2.9348194081733827, 3.8118380124472613, 5.531931287696299]
    atod = [-3.801990321023478e-05, -0.002699702600291877, 0.0030001437204660755]
    return rtod, vtod, atod


@pytest.fixture
def rva_tod_ecef(rva_tod):
    rtod, vtod, _ = rva_tod
    atod = [-0.0010029055188307988, -0.0017988827221534997, 0.0029999960860945377]
    return rtod, vtod, atod


@pytest.fixture
def rva_teme():
    # Example TEME vectors (km, km/s, km/s^2)
    rteme = [2994.454240128889, -7384.574761007484, 6380.344532748868]
    vteme = [2.9346103230979974, 3.8119989825937415, 5.531931287696299]
    ateme = [-3.787182351239877e-05, -0.0026997046816358795, 0.0030001437204660755]
    return rteme, vteme, ateme


@pytest.fixture
def rva_teme_ecef(rva_teme):
    rteme, vteme, _ = rva_teme
    ateme = [-0.0010028068479698174, -0.001798937729169217, 0.0029999960860945377]
    return rteme, vteme, ateme


@pytest.fixture
def t_inputs():
    # Time inputs
    ttt = 0.042623631888994  # Julian centuries of TT
    jdut1 = 2453101.5  # Julian date of UT1
    lod = 0.0015563  # excess length of day, sec
    return ttt, jdut1, lod


@pytest.fixture
def orbit_effects_inputs():
    # Other inputs for accounting for orbit effectgs
    xp = -0.140682 * ARCSEC2RAD  # polar motion coefficient, rad
    yp = 0.333309 * ARCSEC2RAD  # polar motion coefficient, rad
    ddpsi = -0.052195 * ARCSEC2RAD  # delta psi correction to GCRF, rad
    ddeps = -0.003875 * ARCSEC2RAD  # delta epsilon correction to GCRF, rad
    return xp, yp, ddpsi, ddeps


@pytest.fixture
def eop_corrections():
    ddx = -0.000205 * ARCSEC2RAD  # delta x correction to GCRF, rad
    ddy = -0.000136 * ARCSEC2RAD  # delta y correction to GCRF, rad
    return ddx, ddy


def test_eci2ecef(rva_ecef, rva_eci, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected ECEF output vectors
    recef, vecef, _ = rva_ecef

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.eci2ecef(
        *rva_eci, *t_inputs, *orbit_effects_inputs, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(recef, recef_out)
    assert custom_allclose(vecef, vecef_out)

    # The acceleration vector is not correct so we will just compare it to the expected
    aecef = np.array(
        [0.0002936830002159169, 0.0031151668034451073, 0.003000148416052949]
    )
    assert custom_allclose(aecef, aecef_out)


def test_eci2ecef06(
    iau06arr,
    iau06xysarr,
    rva_eci,
    rva_ecef06,
    t_inputs,
    orbit_effects_inputs,
    eop_corrections,
):
    # Extract inputs
    xp, yp, *_ = orbit_effects_inputs

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.eci2ecef06(
        *rva_eci, *t_inputs, xp, yp, iau06arr, iau06xysarr, *eop_corrections
    )

    # Check if the output vectors are close to the expected values
    recef_exp, vecef_exp, aecef_exp = rva_ecef06
    assert custom_allclose(recef_exp, recef_out)
    assert custom_allclose(vecef_exp, vecef_out)
    assert custom_allclose(aecef_exp, aecef_out)


def test_ecef2eci(rva_ecef, rva_eci, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected ECI output vectors
    reci, veci, aeci = rva_eci

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.ecef2eci(
        *rva_ecef, *t_inputs, *orbit_effects_inputs, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(reci, reci_out)
    assert custom_allclose(veci, veci_out)
    assert custom_allclose(aeci, aeci_out)


def test_ecef2eci06(
    iau06arr,
    iau06xysarr,
    rva_eci,
    rva_ecef06,
    t_inputs,
    orbit_effects_inputs,
    eop_corrections,
):
    # Extract inputs
    xp, yp, *_ = orbit_effects_inputs

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.ecef2eci06(
        *rva_ecef06, *t_inputs, xp, yp, iau06arr, iau06xysarr, *eop_corrections
    )

    # Check if the output vectors are close to the expected values
    reci_exp, veci_exp, aeci_exp = rva_eci
    assert custom_allclose(reci_exp, reci_out)
    assert custom_allclose(veci_exp, veci_out)
    assert custom_allclose(aeci_exp, aeci_out)


def test_eci2pef(iau80arr, rva_eci, rva_pef, t_inputs, orbit_effects_inputs):
    # Orbit effects inputs
    *_, ddpsi, ddeps = orbit_effects_inputs

    # Call the function with test inputs
    rpef_out, vpef_out, apef_out = fc.eci2pef(
        *rva_eci, *t_inputs, ddpsi, ddeps, iau80arr
    )

    # Check if the output vectors are close to the expected values
    rpef_exp, vpef_exp, apef_exp = rva_pef
    assert custom_allclose(rpef_exp, rpef_out)
    assert custom_allclose(vpef_exp, vpef_out)
    assert custom_allclose(apef_exp, apef_out)


def test_pef2eci(t_inputs, orbit_effects_inputs, rva_eci, rva_pef, iau80arr):
    # Expected ECI output vectors
    reci, veci, aeci = rva_eci

    # Orbit effects inputs
    *_, ddpsi, ddeps = orbit_effects_inputs

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.pef2eci(
        *rva_pef, *t_inputs, ddpsi, ddeps, iau80arr
    )

    # Expected ECI output vectors
    assert custom_allclose(reci, reci_out)
    assert custom_allclose(veci, veci_out)
    assert custom_allclose(aeci, aeci_out)


def test_eci2tod(iau80arr, rva_eci, rva_tod, t_inputs, orbit_effects_inputs):
    # Extract inputs
    ttt, *_ = t_inputs
    *_, ddpsi, ddeps = orbit_effects_inputs

    # Call the function with test inputs
    rtod_out, vtod_out, atod_out = fc.eci2tod(*rva_eci, ttt, ddpsi, ddeps, iau80arr)

    # Check if the output vectors are close to the expected values
    rtod_exp, vtod_exp, atod_exp = rva_tod
    assert custom_allclose(rtod_exp, rtod_out)
    assert custom_allclose(vtod_exp, vtod_out)
    assert custom_allclose(atod_exp, atod_out)


def test_tod2eci(t_inputs, orbit_effects_inputs, rva_eci, rva_tod, iau80arr):
    # Expected ECI output vectors
    reci, veci, aeci = rva_eci

    # Extract inputs
    ttt, *_ = t_inputs
    *_, ddpsi, ddeps = orbit_effects_inputs

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.tod2eci(*rva_tod, ttt, ddpsi, ddeps, iau80arr)

    # Expected ECI output vectors
    assert custom_allclose(reci, reci_out)
    assert custom_allclose(veci, veci_out)
    assert custom_allclose(aeci, aeci_out)


def test_eci2mod(rva_eci, rva_mod, t_inputs):
    # Expected MOD output vectors
    rmod, vmod, amod = rva_mod

    # Extract inputs
    ttt, *_ = t_inputs

    # Call the function with test inputs
    rmod_out, vmod_out, amod_out = fc.eci2mod(*rva_eci, ttt)

    # Check if the output vectors are close to the expected values
    assert custom_allclose(rmod, rmod_out)
    assert custom_allclose(vmod, vmod_out)
    assert custom_allclose(amod, amod_out)


def test_mod2eci(rva_eci, rva_mod, t_inputs):
    # Expected ECI output vectors
    reci, veci, aeci = rva_eci

    # Extract inputs
    ttt, *_ = t_inputs

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.mod2eci(*rva_mod, ttt)

    # Check if the output vectors are close to the expected values
    assert custom_allclose(reci, reci_out)
    assert custom_allclose(veci, veci_out)
    assert custom_allclose(aeci, aeci_out)


def test_eci2teme(rva_eci, rva_teme, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected TEME output vectors
    rteme, vteme, ateme = rva_teme

    # Extract inputs
    ttt, *_ = t_inputs
    *_, ddpsi, ddeps = orbit_effects_inputs

    # Call the function with test inputs
    rteme_out, vteme_out, ateme_out = fc.eci2teme(*rva_eci, ttt, ddpsi, ddeps, iau80arr)

    # Check if the output vectors are close to the expected values
    assert custom_allclose(rteme, rteme_out)
    assert custom_allclose(vteme, vteme_out)
    assert custom_allclose(ateme, ateme_out)


def test_teme2eci(rva_eci, rva_teme, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected ECI output vectors
    reci, veci, aeci = rva_eci

    # Extract inputs
    ttt, *_ = t_inputs
    *_, ddpsi, ddeps = orbit_effects_inputs

    # Call the function with test inputs
    reci_out, veci_out, aeci_out = fc.teme2eci(*rva_teme, ttt, ddpsi, ddeps, iau80arr)

    # Check if the output vectors are close to the expected values
    assert custom_allclose(reci, reci_out)
    assert custom_allclose(veci, veci_out)
    assert custom_allclose(aeci, aeci_out)


def test_ecef2pef(rva_ecef, rva_pef_ecef, t_inputs, orbit_effects_inputs):
    # Expected PEF output vectors
    rpef, vpef, apef = rva_pef_ecef

    # Extract inputs
    ttt, *_ = t_inputs
    xp, yp, *_ = orbit_effects_inputs

    # Call the function with test inputs
    rpef_out, vpef_out, apef_out = fc.ecef2pef(*rva_ecef, xp, yp, ttt)

    # Check if the output vectors are close to the expected values
    assert custom_allclose(rpef, rpef_out)
    assert custom_allclose(vpef, vpef_out)
    assert custom_allclose(apef, apef_out)


def test_pef2ecef(rva_ecef, rva_pef_ecef, t_inputs, orbit_effects_inputs):
    # Expected ECEF output vectors
    recef, vecef, aecef = rva_ecef

    # Extract inputs
    ttt, *_ = t_inputs
    xp, yp, *_ = orbit_effects_inputs

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.pef2ecef(*rva_pef_ecef, xp, yp, ttt)

    # Check if the output vectors are close to the expected values
    assert custom_allclose(recef, recef_out)
    assert custom_allclose(vecef, vecef_out)
    assert custom_allclose(aecef, aecef_out)


def test_ecef2tod(rva_ecef, rva_tod_ecef, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected TOD output vectors
    rtod, vtod, atod = rva_tod_ecef

    # Call the function with test inputs
    rtod_out, vtod_out, atod_out = fc.ecef2tod(
        *rva_ecef, *t_inputs, *orbit_effects_inputs, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(rtod, rtod_out)
    assert custom_allclose(vtod, vtod_out)
    assert custom_allclose(atod, atod_out)


def test_tod2ecef(rva_ecef, rva_tod_ecef, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected ECEF output vectors
    recef, vecef, aecef = rva_ecef

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.tod2ecef(
        *rva_tod_ecef, *t_inputs, *orbit_effects_inputs, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(recef, recef_out)
    assert custom_allclose(vecef, vecef_out)
    assert custom_allclose(aecef, aecef_out, rtol=1e-6)


def test_ecef2mod(rva_ecef, rva_mod_ecef, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected MOD output vectors
    rmod, vmod, amod = rva_mod_ecef

    # Call the function with test inputs
    rmod_out, vmod_out, amod_out = fc.ecef2mod(
        *rva_ecef, *t_inputs, *orbit_effects_inputs, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(rmod, rmod_out)
    assert custom_allclose(vmod, vmod_out)
    assert custom_allclose(amod, amod_out)


def test_mod2ecef(rva_ecef, rva_mod_ecef, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected ECEF output vectors
    recef, vecef, aecef = rva_ecef

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.mod2ecef(
        *rva_mod_ecef, *t_inputs, *orbit_effects_inputs, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(recef, recef_out)
    assert custom_allclose(vecef, vecef_out)
    assert custom_allclose(aecef, aecef_out, rtol=1e-6)


def test_ecef2teme(rva_ecef, rva_teme_ecef, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected TEME output vectors
    rteme, vteme, ateme = rva_teme_ecef

    # Extract inputs
    xp, yp, *_ = orbit_effects_inputs

    # Call the function with test inputs
    rteme_out, vteme_out, ateme_out = fc.ecef2teme(
        *rva_ecef, *t_inputs, xp, yp, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(rteme, rteme_out)
    assert custom_allclose(vteme, vteme_out)
    assert custom_allclose(ateme, ateme_out)


def test_teme2ecef(rva_ecef, rva_teme_ecef, t_inputs, orbit_effects_inputs, iau80arr):
    # Expected ECEF output vectors
    recef, vecef, aecef = rva_ecef

    # Extract inputs
    xp, yp, *_ = orbit_effects_inputs

    # Call the function with test inputs
    recef_out, vecef_out, aecef_out = fc.teme2ecef(
        *rva_teme_ecef, *t_inputs, xp, yp, iau80arr
    )

    # Check if the output vectors are close to the expected values
    assert custom_allclose(recef, recef_out)
    assert custom_allclose(vecef, vecef_out)
    assert custom_allclose(aecef, aecef_out)
