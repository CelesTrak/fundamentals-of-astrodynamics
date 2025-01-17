import numpy as np
import pytest

import src.valladopy.astro.twobody.frame_conversions as fc
from src.valladopy.astro.time.data import iau80in
from src.valladopy.astro.twobody.utils import OrbitType
from src.valladopy.constants import ARCSEC2RAD
from ...conftest import custom_isclose, custom_allclose


DEFAULT_TOL = 1e-12


@pytest.fixture()
def iau80arr():
    """Load the IAU 1980 data"""
    return iau80in()


# ECI position/velocity vectors for a number of tests
@pytest.fixture
def rv():
    reci = [1525.9870698051157, -5867.209915411114, 3499.601587508083]
    veci = [1.4830443958075603, -7.093267951700349, 0.9565730381487033]
    return np.array(reci), np.array(veci)


# Input parameters related to ECEF conversions
@pytest.fixture
def ecef_inputs():
    ttt = 0.042623631888994
    jdut1 = 2.4531015e06
    lod = 0.0015563
    xp = -0.140682 * ARCSEC2RAD
    yp = 0.333309 * ARCSEC2RAD
    ddpsi = -0.052195 * ARCSEC2RAD
    ddeps = -0.003875 * ARCSEC2RAD
    return ttt, jdut1, lod, xp, yp, ddpsi, ddeps


class TestSpherical:
    @pytest.fixture
    def rv(self):
        # Position and velocity in km and km/s
        r = np.array([4.286607049870562e03, 4.286607049870561e03, 3.5e03])
        v = np.array([4.059474712855235, 4.860051329924127, 4.018776695238445])
        return r, v

    @pytest.fixture
    def rvmag(self, rv):
        rmag = np.linalg.norm(rv[0])
        vmag = np.linalg.norm(rv[1])
        return rmag, vmag

    @pytest.fixture
    def adbarv(self, rvmag):
        rtasc = np.radians(45)
        decl = np.radians(30)
        fpav = np.radians(5)
        az = np.radians(60)
        return rvmag[0], rvmag[1], rtasc, decl, fpav, az

    def test_adbar2rv(self, rv, adbarv):
        # Unpack variables
        rmag, vmag, rtasc, decl, fpav, az = adbarv
        expected_r, expected_v = rv  # expected values

        # Call the function with test inputs
        r, v = fc.adbar2rv(rmag, vmag, rtasc, decl, fpav, az)

        # Check if the output is close to the expected values
        assert np.allclose(r, expected_r, rtol=DEFAULT_TOL)
        assert np.allclose(v, expected_v, rtol=DEFAULT_TOL)

    def test_rv2adbar(self, rv, adbarv):
        expected_elems = adbarv

        # Call the function with test inputs
        out = fc.rv2adbar(rv[0], rv[1])

        # Check if the output is close to the expected values
        assert np.allclose(out, expected_elems, rtol=DEFAULT_TOL)


class TestClassical:
    @pytest.fixture
    def coe(self):
        # Vallado, 2007, Ex. 2-6
        p = 11067.790  # semi-latus rectum, km
        ecc = 0.83285  # eccentricity
        incl = np.radians(87.87)  # inclination, rad
        raan = np.radians(227.89)  # RAAN, rad
        argp = np.radians(53.38)  # arg. of periapsis, rad
        nu = np.radians(92.335)  # true anomaly, rad
        return p, ecc, incl, raan, argp, nu

    @pytest.fixture
    def rv(self):
        # Vallado, 2007, Ex. 2-5
        # Position and velocity in km and km/s
        r = np.array([6524.834, 6862.875, 6448.296])
        v = np.array([4.901327, 5.533756, -1.976341])
        return r, v

    def test_coe2rv(self, coe):
        # Vallado, 2007, Ex. 2-6
        p, ecc, incl, raan, _, nu = coe
        r_exp = np.array([-4.049198890323112e03, -4.479765179366826e03, 0])
        v_exp = np.array([0.303279847002191, -0.274130533804499, 10.9917080783198])

        # Call the function with test inputs
        r_out, v_out = fc.coe2rv(p, ecc, incl, raan, nu)

        # Check if the output is close to the expected values
        assert np.allclose(r_out, r_exp, rtol=DEFAULT_TOL)
        assert np.allclose(v_out, v_exp, rtol=DEFAULT_TOL)

    def test_rv2coe(self, rv):
        # Vallado, 2007, Ex. 2-5
        # TODO: add tests for other orbit type cases
        # Call the function with test inputs
        (p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper, orbit_type) = (
            fc.rv2coe(*rv)
        )

        # Check if the output is close to the expected values
        # TODO: lonper is not `nan` in the book example (but is in matlab)
        assert abs(p - 11067.798350991814) < DEFAULT_TOL
        assert abs(a - 36127.337763974785) < DEFAULT_TOL
        assert abs(ecc - 0.8328533990836885) < DEFAULT_TOL
        assert abs(incl - 1.5336055626394494) < DEFAULT_TOL
        assert abs(raan - 3.9775750028016947) < DEFAULT_TOL
        assert abs(argp - 0.9317428111437854) < DEFAULT_TOL
        assert abs(nu - 1.6115524999414736) < DEFAULT_TOL
        assert abs(m - 0.1327277817291186) < DEFAULT_TOL
        assert abs(arglat - 2.5432953110852594) < DEFAULT_TOL
        assert np.isnan(truelon)
        assert np.isnan(lonper)
        assert orbit_type is OrbitType.EPH_INCLINED


class TestEquinoctial:
    # TODO: validate with more test cases
    @pytest.fixture
    def eq(self):
        a = 7000  # semimajor axis, km
        af = 0.001  # eccentricity component
        ag = 0.001  # eccentricity component
        chi = 0.001  # EQW node vector component
        psi = 0.001  # EQW node vector component
        meanlon = np.pi / 4  # mean longitude, rad
        fr = 1  # retrograde factor (+1 = prograde)
        return a, af, ag, chi, psi, meanlon, fr

    @pytest.fixture
    def rv(self):
        r = np.array([4942.747468305833, 4942.747468305833, 0])
        v = np.array([-5.34339547370427, 5.343395473704271, 0.021373624642066363])
        return r, v

    def test_eq2rv(self, eq, rv):
        # Expected outputs
        r_exp, v_exp = rv

        # Call the function with test inputs
        r_out, v_out = fc.eq2rv(*eq)

        # Check if the outputs are close to the expected values
        assert np.allclose(r_out, r_exp, rtol=DEFAULT_TOL)
        assert np.allclose(v_out, v_exp, rtol=DEFAULT_TOL)

    def test_rv2eq(self, rv, eq):
        # expected outputs
        a_exp, af_exp, ag_exp, chi_exp, psi_exp, meanlon_exp, fr_exp = eq

        # Call the function with test inputs
        a, n, af, ag, chi, psi, meanlon, truelon, fr = fc.rv2eq(*rv)

        # Check if the outputs are close to the expected values
        assert abs(a - a_exp) < DEFAULT_TOL
        assert abs(n - 0.001078007612466833) < DEFAULT_TOL
        assert abs(af - af_exp) < DEFAULT_TOL
        assert abs(ag - ag_exp) < DEFAULT_TOL
        assert abs(chi - chi_exp) < DEFAULT_TOL
        assert abs(psi - psi_exp) < DEFAULT_TOL
        assert abs(meanlon - meanlon_exp) < DEFAULT_TOL
        assert abs(truelon - np.pi / 4) < DEFAULT_TOL
        assert int(fr) == int(fr_exp)


class TestTopocentric:
    @pytest.fixture
    def rvseci(self):
        # ECI site position and velocity vector in km and km/s
        rseci = [-2.968655122428691e03, 3.980613919662232e03, 3.992860345291290e03]
        vseci = [-0.290278922351514, -0.216325537609299, -0.000157672327972]
        return rseci, vseci

    @pytest.fixture
    def tradec(self):
        # Topocentric coordinates
        rho = 4.437731184421759e09  # range, km
        trtasc = 5.148532095674960  # right ascension, rad
        tdecl = -0.363438990548242  # declination, rad
        drho = -25.599038196399519  # range rate, km/s
        tdrtasc = -2.051513501139983e-09  # right ascension rate, rad/s
        tddecl = -3.189648164446254e-10  # declination rate, rad/s
        return rho, trtasc, tdecl, drho, tdrtasc, tddecl

    @pytest.fixture
    def rveci(self):
        r_eci = [1752246215.6652846, -3759563434.243893, -1577568101.96675]
        v_eci = [-18.323497062513614, 18.332049491766764, 7.777041227057346]
        return r_eci, v_eci

    def test_tradec2rv(self, rvseci, tradec, rveci):
        # Expected outputs
        r_eci_exp, v_eci_exp = rveci

        # Call the function with test inputs
        r_eci, v_eci = fc.tradec2rv(*tradec, *rvseci)

        # Check if the output is close to the expected values
        assert np.allclose(r_eci, r_eci_exp, rtol=DEFAULT_TOL)
        assert np.allclose(v_eci, v_eci_exp, rtol=DEFAULT_TOL)

    def test_rv2tradec(self, rvseci, tradec, rveci):
        # Call the function with test inputs
        tradec_out = fc.rv2tradec(*rveci, *rvseci)

        # Check if the output is close to the expected values
        assert np.allclose(tradec_out, np.array(tradec), rtol=DEFAULT_TOL)


class TestFlight:
    @pytest.fixture
    def rvmag(self):
        rmag = 7000  # position magnitude, km
        vmag = 7.546  # velocity magnitude, km
        return rmag, vmag

    @pytest.fixture
    def flight(self):
        latgc = np.pi / 6  # 30 degrees
        lon = np.pi / 2  # 90 degrees
        fpa = -np.pi / 6  # -30 degrees
        az = np.pi / 4  # 45 degrees
        return latgc, lon, fpa, az

    def test_flt2rv(self, iau80arr, rv, rvmag, flight, ecef_inputs):
        # Expected outputs
        reci_exp, veci_exp = rv

        # Call the function with test inputs
        reci, veci = fc.flt2rv(*rvmag, *flight, *ecef_inputs, iau80arr)

        # Check if the output is close to the expected values
        assert np.allclose(reci, reci_exp, rtol=DEFAULT_TOL)
        assert np.allclose(veci, veci_exp, rtol=DEFAULT_TOL)

    def test_rv2flt(self, iau80arr, rv, rvmag, flight, ecef_inputs):
        # Expected outputs
        rmag_exp, _ = rvmag  # vmag will not match original
        latgc_exp, lon_exp, _, _ = flight  # fpa and az will not match original

        # Call the function with test inputs
        lon, latgc, rtasc, decl, fpa, az, rmag, vmag = fc.rv2flt(
            *rv, *ecef_inputs, iau80arr
        )

        # Check if the output is close to the expected values
        # Some values differ from the orignal ones inputted in the test above
        assert np.isclose(lon, lon_exp, rtol=DEFAULT_TOL)
        assert np.isclose(latgc, latgc_exp, rtol=DEFAULT_TOL)
        assert np.isclose(rtasc, -1.3163464523837871, rtol=DEFAULT_TOL)
        assert np.isclose(decl, 0.5235330558280773, rtol=DEFAULT_TOL)
        assert np.isclose(fpa, 1.1758912957606626, rtol=DEFAULT_TOL)
        assert np.isclose(az, -3.016745055466374, rtol=DEFAULT_TOL)
        assert np.isclose(rmag, rmag_exp, rtol=DEFAULT_TOL)
        assert np.isclose(vmag, 7.309507705165137, rtol=DEFAULT_TOL)


class TestEcliptic:
    @pytest.fixture
    def ell(self):
        rr = 7000
        ecllon = -1.205591131641763
        ecllat = 0.9142332891895931
        drr = 6.746917593386583
        decllon = -0.0001879622667434771
        decllat = -0.0003849993355596977
        return rr, ecllon, ecllat, drr, decllon, decllat

    def test_ell2rv(self, rv, ell):
        # Expected outputs
        reci_exp, veci_exp = rv

        # Call function with test inputs
        reci, veci = fc.ell2rv(*ell)

        # Check if output values are close
        assert np.allclose(reci, reci_exp, rtol=DEFAULT_TOL)
        assert np.allclose(veci, veci_exp, rtol=DEFAULT_TOL)

    def test_rv2ell(self, rv, ell):
        # Expected outputs
        rr_exp, ecllon_exp, ecllat_exp, drr_exp, decllon_exp, decllat_exp = ell

        # Call function with test inputs
        rr, ecllon, ecllat, drr, decllon, decllat = fc.rv2ell(*rv)

        # Check if output values are close
        assert np.isclose(rr, rr_exp, rtol=DEFAULT_TOL)
        assert np.isclose(ecllon, ecllon_exp, rtol=DEFAULT_TOL)
        assert np.isclose(ecllat, ecllat_exp, rtol=DEFAULT_TOL)
        assert np.isclose(drr, drr_exp, rtol=DEFAULT_TOL)
        assert custom_isclose(decllon, decllon_exp)
        assert custom_isclose(decllat, decllat_exp)


class TestCelestial:
    @pytest.fixture
    def radec(self):
        rr = 7000
        rtasc = 4.9668388547958
        decl = 0.5235330558280773
        drr = 6.746917593386583
        drtasc = -5.776166833354767e-05
        ddecl = -0.0003986042868175495
        return rr, rtasc, decl, drr, drtasc, ddecl

    def test_radec2rv(self, rv, radec):
        # Expected outputs
        reci_exp, veci_exp = rv

        # Call function with test inputs
        reci, veci = fc.radec2rv(*radec)

        # Check if output values are close
        assert np.allclose(reci, reci_exp, rtol=DEFAULT_TOL)
        assert np.allclose(veci, veci_exp, rtol=DEFAULT_TOL)

    def test_rv2radec(self, rv, radec):
        # Expected outputs
        rr_exp, rtasc_exp, decl_exp, drr_exp, drtasc_exp, ddecl_exp = radec

        # Call function with test inputs
        rr, rtasc, decl, drr, drtasc, ddecl = fc.rv2radec(*rv)

        # Check if output values are close
        assert np.isclose(rr, rr_exp, rtol=DEFAULT_TOL)
        assert np.isclose(rtasc, rtasc_exp, rtol=DEFAULT_TOL)
        assert np.isclose(decl, decl_exp, rtol=DEFAULT_TOL)
        assert np.isclose(drr, drr_exp, rtol=DEFAULT_TOL)
        assert custom_isclose(drtasc, drtasc_exp)
        assert custom_isclose(ddecl, ddecl_exp)


class TestAzEl:
    @pytest.fixture
    def rvecef(self):
        # ECEF position and velocity vector in km and km/s
        recef = [-5225.545532658024, -3070.9865969614993, 3501.8159870845006]
        vecef = [-5.373672578283977, -3.3021038634286506, 4.020523711114683]
        return recef, vecef

    @pytest.fixture
    def rvsez(self):
        rhosez = [-29.206599480225734, -4261.469869512799, -818.4151925779145]
        drhosez = [-0.24683193316601043, -4.345266040330323, 6.0829758219664605]
        return np.array(rhosez), np.array(drhosez)

    @pytest.fixture
    def lla(self):
        # Site geodetic coordinates
        latgd = np.radians(39.007)
        lon = np.radians(-104.883)
        alt = 2.102
        return latgd, lon, alt

    @pytest.fixture
    def azel(self):
        rho = 4339.444884044615
        az = -1.5639427896150881
        el = -0.18973540227923288
        drho = 3.1216042513457634
        daz = 5.093097833129788e-05
        del_el = 0.0015655515448868324
        return rho, az, el, drho, daz, del_el

    def test_raz2rvs(self, azel, rvsez):
        # Expected outputs
        rhosez_exp, drhosez_exp = rvsez

        # Call function with test inputs
        rhosez, drhosez = fc.raz2rvs(*azel)

        # Check if output values are close
        assert np.allclose(rhosez, rhosez_exp, rtol=DEFAULT_TOL)
        assert np.allclose(drhosez, drhosez_exp, rtol=DEFAULT_TOL)

    def test_rvs2raz(self, azel, rvsez):
        # Expected outputs
        rho_exp, az_exp, el_exp, drho_exp, daz_exp, del_el_exp = azel

        # Call function with test inputs
        rho, az, el, drho, daz, del_el = fc.rvs2raz(*rvsez)

        # Check if output values are close
        assert np.isclose(rho, rho_exp, rtol=DEFAULT_TOL)
        assert np.isclose(az, az_exp, rtol=DEFAULT_TOL)
        assert np.isclose(el, el_exp, rtol=DEFAULT_TOL)
        assert np.isclose(drho, drho_exp, rtol=DEFAULT_TOL)
        assert custom_isclose(daz, daz_exp)
        assert custom_isclose(del_el, del_el_exp)

    def test_razel2rv(self, iau80arr, rvecef, ecef_inputs, lla, azel):
        # Expected outputs
        recef_exp, vecef_exp = rvecef

        # Call function with test inputs
        recef, vecef = fc.razel2rv(*azel, *lla)

        # Check if output values are close
        assert np.allclose(recef, recef_exp, rtol=DEFAULT_TOL)
        assert np.allclose(vecef, vecef_exp, rtol=DEFAULT_TOL)

    def test_rv2razel(self, iau80arr, rvecef, lla, azel):
        # Expected outputs
        rho_exp, az_exp, el_exp, drho_exp, daz_exp, del_el_exp = azel

        # Call function with test inputs
        rho, az, el, drho, daz, del_el = fc.rv2razel(*rvecef, *lla)

        # Check if output values are close
        assert np.isclose(rho, rho_exp, rtol=DEFAULT_TOL)
        assert np.isclose(az, az_exp, rtol=DEFAULT_TOL)
        assert np.isclose(el, el_exp, rtol=DEFAULT_TOL)
        assert np.isclose(drho, drho_exp, rtol=DEFAULT_TOL)
        assert custom_isclose(daz, daz_exp)
        assert custom_isclose(del_el, del_el_exp)


class TestSatCoord:
    def test_rv2rsw(self, rv):
        # Expected outputs
        rrsw_exp = np.array(
            [6999.999999999998, -4.547473508864641e-13, 2.2737367544323206e-13]
        )
        vrsw_exp = np.array(
            [6.746917593386584, 2.8121176860009087, 2.220446049250313e-16]
        )
        transmat_exp = np.array(
            [
                [0.21799815282930227, -0.8381728450587307, 0.4999430839297262],
                [0.0043486171359898, -0.5114241285968169, -0.8593174327441467],
                [0.9759394934584914, 0.18950367409403857, -0.10784462254949785],
            ]
        )

        # Call function with test inputs
        rrsw, vrsw, transmat = fc.rv2rsw(*rv)

        # Check if output values are close
        assert np.allclose(rrsw, rrsw_exp, rtol=DEFAULT_TOL)
        assert np.allclose(vrsw, vrsw_exp, rtol=DEFAULT_TOL)
        assert custom_allclose(transmat, transmat_exp)

    def test_rv2ntw(self, rv):
        # Expected outputs
        rntw_exp = np.array(
            [2693.043717307586, 6461.231735255294, 2.2737367544323206e-13]
        )
        vntw_exp = np.array(
            [2.220446049250313e-16, 7.309507705165138, 2.220446049250313e-16]
        )
        transmat_exp = np.array(
            [
                [0.07985444754543984, 0.1495990993516855, 0.9855168068989881],
                [0.2028925141921107, -0.9704166460742649, 0.13086695804051995],
                [0.9759394934584914, 0.18950367409403857, -0.10784462254949785],
            ]
        )

        # Call function with test inputs
        rntw, vntw, transmat = fc.rv2ntw(*rv)

        # Check if output values are close
        assert np.allclose(rntw, rntw_exp, rtol=DEFAULT_TOL)
        assert np.allclose(vntw, vntw_exp, rtol=DEFAULT_TOL)
        assert custom_allclose(transmat, transmat_exp)

    def test_rv2pqw(self, rv):
        # Expected outputs
        rpqw_exp = np.array([-6528.341050907748, 2526.017245197032, 0])
        vpqw_exp = np.array([-7.307090980337026, -0.18794759095509114, 0])

        # Call function with test inputs
        rpqw, vpqw = fc.rv2pqw(*rv)

        # Check if output values are close
        assert np.allclose(rpqw, rpqw_exp, rtol=DEFAULT_TOL)
        assert np.allclose(vpqw, vpqw_exp, rtol=DEFAULT_TOL)


class TestGeodetic:
    @pytest.fixture
    def recef(self):
        # ECEF position vector in km (Example 3-3, Vallado, 2001)
        return np.array([6524.834, 6862.875, 6448.296])

    def test_ecef2ll(self, recef):
        # Astronomical Almanac method
        latgc, latgd, lon, hellp = fc.ecef2ll(recef)
        assert np.isclose(latgc, 0.597826066235814, rtol=DEFAULT_TOL)
        assert np.isclose(latgd, 0.5995641464668334, rtol=DEFAULT_TOL)
        assert np.isclose(lon, 0.8106428999047803, rtol=DEFAULT_TOL)
        assert np.isclose(hellp, 5085.219430346959, rtol=DEFAULT_TOL)

    def test_ecef2llb(self, recef):
        # Borkowski method
        latgc, latgd, lon, hellp = fc.ecef2llb(recef)
        assert np.isclose(latgc, 0.597826066235814, rtol=DEFAULT_TOL)
        assert np.isclose(latgd, 0.5995641464669065, rtol=DEFAULT_TOL)
        assert np.isclose(lon, 0.8106428999047803, rtol=DEFAULT_TOL)
        assert np.isclose(hellp, 5085.2194303451715, rtol=DEFAULT_TOL)


def test_perifocal_transform():
    # From Vallado Example 2-6, 2022
    i = np.radians(87.87)
    raan = np.radians(227.89)
    w = np.radians(53.38)

    # Call the function with test inputs
    p_, q_, w_ = fc.perifocal_transform(i, raan, w)

    # Check if the output is close to the expected values
    assert custom_allclose(
        p_, [-0.37786007180211706, -0.46252560096516104, 0.8020547578498088]
    )
    assert custom_allclose(
        q_, [0.5546417852962037, 0.580556376955298, 0.596092931664164]
    )
    assert custom_allclose(
        w_, [-0.7413462457860962, 0.6700928007584878, 0.03716695078301055]
    )
