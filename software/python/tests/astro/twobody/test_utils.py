import numpy as np
import pytest

import src.valladopy.astro.twobody.utils as utils
from ...conftest import DEFAULT_TOL, custom_isclose, custom_allclose


@pytest.mark.parametrize(
    "ecc, incl, expected_orbit_type",
    [
        (0, 0, utils.OrbitType.CIR_EQUATORIAL),  # circular equatorial
        (0, np.pi, utils.OrbitType.CIR_EQUATORIAL),  # circular equatorial
        (0, 0.1, utils.OrbitType.CIR_INCLINED),  # circular inclined
        (0.1, 0, utils.OrbitType.EPH_EQUATORIAL),  # elliptical equatorial
        (0.1, np.pi, utils.OrbitType.EPH_EQUATORIAL),  # elliptical equatorial
        (0.1, 0.1, utils.OrbitType.EPH_INCLINED),  # elliptical inclined
        (1.1, 0, utils.OrbitType.EPH_EQUATORIAL),  # hyperbolic equatorial
        (1.1, np.pi, utils.OrbitType.EPH_EQUATORIAL),  # hyperbolic equatorial
        (1.1, 0.1, utils.OrbitType.EPH_INCLINED),  # hyperbolic inclined
    ],
)
def test_determine_orbit_type(ecc, incl, expected_orbit_type):
    assert utils.determine_orbit_type(ecc, incl) == expected_orbit_type


@pytest.mark.parametrize(
    "inc, expected",
    [(0, True), (np.pi, True), (np.pi / 2, False), (np.pi / 4, False)]  # fmt: skip
)
def test_is_equatorial(inc, expected):
    assert utils.is_equatorial(inc) == expected


def test_elliptic12():
    """Test the elliptic integrals of the first and second kind and their inverses."""
    # Inputs
    u = de = 1.45388206137796e-06
    m = 6.66133814775094e-16**2  # target ecc

    # Expected
    f_exp = 1.4538820613779603e-06
    e_exp = inv_e_exp = 1.45388206137796e-06
    z_exp = 0

    # Test scalar inputs
    f, e, z = utils.elliptic12(u, m)
    inv_e = utils.inverse_elliptic2(de, m)
    assert custom_isclose(f, f_exp)
    assert custom_isclose(e, e_exp)
    assert custom_isclose(z, z_exp)
    assert custom_isclose(inv_e, inv_e_exp)

    # Test vector inputs
    f, e, z = utils.elliptic12([u, 0.1], [m, 0.01])
    inv_e = utils.inverse_elliptic2([de, 0.1], [m, 0.01])
    assert custom_allclose(f, [f_exp, 0.1000016634111543])
    assert custom_allclose(e, [e_exp, 0.09999833663861171])
    assert custom_allclose(z, [z_exp, 0.0004973097001265137])
    assert custom_allclose(inv_e, [inv_e_exp, 0.10000166344428708])

    # Check bad inputs
    for m in [-0.01, 1.01]:
        with pytest.raises(ValueError):
            utils.elliptic12(u, m)
        with pytest.raises(ValueError):
            utils.inverse_elliptic2(de, m)


@pytest.mark.parametrize(
    "a, b, theta0, theta1, s_exp",
    [
        (5, 5, 0, 2 * np.pi, 5 * 2 * np.pi),  # full circle
        (5, 10, 0, 2 * np.pi, 48.44224110273838),  # full ellipse
        (5, 10, np.pi / 10, np.pi / 2, 9.007389063529619),  # partial ellipse
        (1, 0.5, np.pi * 0.002, np.pi * 0.001, -0.0015708505862375155),  # negative
        (1, 2, np.pi * 0.001, np.pi * 0.002, 0.006283131046402917),  # partial, small
    ],
)
def test_arclength_ellipse(a, b, theta0, theta1, s_exp):
    assert custom_isclose(utils.arclength_ellipse(a, b, theta0, theta1), s_exp)


@pytest.mark.parametrize(
    "latgd, lon, alt, rsecef_exp",
    [
        (0, 0, 0, np.array([6378.1363, 0, 0])),
        (
            np.pi / 4,
            np.pi / 4,
            0,
            np.array([3194.4187944733417, 3194.418794473341, 4487.347916379808]),
        ),
        (
            np.radians(39.007),
            np.radians(-104.883),
            2.187,  # Ex. 7-1
            np.array([-1275.1217701148173, -4797.988500532293, 3994.2970129784385]),
        ),
    ],
)
def test_site(latgd, lon, alt, rsecef_exp):
    rsecef, vsecef = utils.site(latgd, lon, alt)
    assert np.allclose(rsecef, rsecef_exp, rtol=DEFAULT_TOL)
    assert np.allclose(vsecef, [0, 0, 0], rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "znew, c2new_exp, c3new_exp",
    [
        (-39.47842, 6.756775284482481, 1.054067771837564),
        (0, 0.5, 0.16666666666666666),
        (0.57483, 0.47650299902524496, 0.16194145738041812),
        (39.47842, 4.612053269571722e-16, 0.025330293604932358),
    ],
)
def test_findc2c3(znew, c2new_exp, c3new_exp):
    c2new, c3new = utils.findc2c3(znew)
    assert np.isclose(c2new, c2new_exp, rtol=DEFAULT_TOL)
    assert np.isclose(c3new, c3new_exp, rtol=DEFAULT_TOL)


def test_lon_nu_conv():
    # Inputs
    jdut1 = 2456864.5
    incl = np.radians(0.070273056)
    raan = np.radians(19.90450011)
    argp = np.radians(352.5056022)
    lon = np.radians(7.020438698)
    nu = np.radians(298.2221136844455)

    # Convert lon to nu
    nu_out = utils.lon2nu(jdut1, lon, incl, raan, argp)

    # Convert nu to lon
    lon_out = utils.nu2lon(jdut1, nu, incl, raan, argp)

    # Check that the original and output values are the same
    assert np.isclose(nu_out, nu, rtol=DEFAULT_TOL)
    assert np.isclose(lon_out, lon, rtol=DEFAULT_TOL)


def test_gc2gd():
    latgd = utils.gc2gd(np.radians(34.25))
    assert np.isclose(latgd, 0.6009038529757992, rtol=DEFAULT_TOL)


def test_gd2gc():
    latgc = utils.gd2gc(np.radians(125.79))
    assert np.isclose(latgc, -0.9429532487023382, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "altpad, r1, v1, r2, v2, nrev, expected_hitearth, expected_hitearthstr",
    # TODO: Add test cases for parabolic, elliptical, and hyperbolic impacts
    [
        # Case 1: Impact at initial/final position
        (
            100,
            [6378, 0, 0],
            [0, 7.8, 0],
            [7000, 0, 0],
            [0, 7, 0],
            0,
            True,
            "Impact at initial/final radii",
        ),
        # Case 2: Impact during multi-revolution transfer
        (
            100,
            [7000, 0, 0],
            [0, 7, 0],
            [8000, 0, 0],
            [0, 6.5, 0],
            1,
            True,
            "Impact during nrev",
        ),
        # Case 3: No impact
        (
            100,
            [8000, 0, 0],
            [0, 7, 0],
            [9000, 0, 0],
            [0, 6.5, 0],
            0,
            False,
            "No impact",
        ),
    ],
)
def test_checkhitearth(
    altpad, r1, v1, r2, v2, nrev, expected_hitearth, expected_hitearthstr
):
    hitearth, hitearthstr = utils.checkhitearth(altpad, r1, v1, r2, v2, nrev)
    assert hitearth == expected_hitearth
    assert hitearthstr == expected_hitearthstr


@pytest.mark.parametrize(
    "r, p, tof_expected",
    [
        ([5000, 5000, 0], 8000, 679.4405031710606),  # elliptical case
        ([20000, 20000, 0], 8000, 2919.2991312782533),  # hyperbolic case
        ([14000, 0, 0], 0, 799.557259683841),  # parabolic case
    ],
)
def test_findtof(r, p, tof_expected):
    ro = [7000, 0, 0]
    tof = utils.findtof(ro, r, p)
    assert np.isclose(tof, tof_expected, rtol=DEFAULT_TOL)


def test_findtof_bad_inputs():
    # Nonsensical orbit (a = 0, p = nonzero)
    ro = [7000, 0, 0]
    r1 = [14000, 0, 0]
    p = 14000
    assert np.isnan(utils.findtof(ro, r1, p))
