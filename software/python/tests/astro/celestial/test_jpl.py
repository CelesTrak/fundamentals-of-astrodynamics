import numpy as np
import pytest
from pathlib import Path

import src.valladopy.astro.celestial.jpl as jpl

from ...conftest import DEFAULT_TOL


@pytest.fixture
def sunmooneph_filepath(test_data_dir):
    filepath = Path(test_data_dir) / "sunmooneph_430t.txt"
    assert filepath.exists()
    return filepath


@pytest.fixture
def sunmooneph_filepath_12hr(test_data_dir):
    filepath = Path(test_data_dir) / "sunmooneph_430t12.txt"
    assert filepath.exists()
    return filepath


@pytest.mark.parametrize(
    "include_hr, jdstart, jdstartf, yr_start_stop, rsun0_exp, mjd0_exp",
    [
        (
            False,
            2435839.5,
            0,
            (1957, 2098),
            [27869796.6511, -132514359.7840, -57466883.5187],
            35839,
        ),
        (
            True,
            2435838.5,
            0.5,
            (1956, 2098),
            [26583662.5476, -132737262.3623, -57563602.5894],
            35838.5,
        ),
    ],
)
def test_read_jplde(
    sunmooneph_filepath,
    sunmooneph_filepath_12hr,
    include_hr,
    jdstart,
    jdstartf,
    yr_start_stop,
    rsun0_exp,
    mjd0_exp,
):
    # Get ephem filepath
    filepath = sunmooneph_filepath_12hr if include_hr else sunmooneph_filepath

    # Call the function
    jpldearr, jdjpldestart, jdjpldestart_frac = jpl.read_jplde(filepath, include_hr)

    # Check start Julian date and fractional part
    assert np.isclose(jdjpldestart, jdstart, rtol=DEFAULT_TOL)
    assert np.isclose(jdjpldestart_frac, jdstartf, rtol=DEFAULT_TOL)

    # Validate the structure of jpldearr
    expected_keys = {
        "year",
        "month",
        "day",
        "hour",
        "rsun1",
        "rsun2",
        "rsun3",
        "rsmag",
        "rmoon1",
        "rmoon2",
        "rmoon3",
        "mjd",
    }
    assert set(jpldearr.keys()) == expected_keys

    # Ensure data is present and that all arrays have the same length
    n_pts = len(jpldearr["year"])
    assert n_pts > 0  # ensure data is present
    for keys in expected_keys:
        assert len(jpldearr[keys]) == n_pts

    # Spot-check some values
    rsun = np.array([jpldearr["rsun1"][0], jpldearr["rsun2"][0], jpldearr["rsun3"][0]])
    assert jpldearr["year"][0], jpldearr["year"][-1] == yr_start_stop
    assert np.allclose(rsun, rsun0_exp, rtol=DEFAULT_TOL)
    assert np.isclose(jpldearr["mjd"][0], mjd0_exp, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "include_hr, interp, rsun_exp, rmoon_exp",
    [
        (
            False,
            None,
            [96576094.2145, 106598001.2476, 46210616.7776],
            [-252296.5509, -302841.7334, -93212.772],
        ),
        (
            False,
            jpl.JPLInterp.LINEAR,
            [96262266.54879406, 106839759.58852758, 46315422.3039953],
            [-241011.64921426764, -309676.4021824469, -96114.65287251282],
        ),
        (
            False,
            jpl.JPLInterp.SPLINE,
            [96264116.08394983, 106841825.88313565, 46316319.46430736],
            [-241686.85556172315, -310637.0152811595, -96414.10549215478],
        ),
        (
            True,
            None,
            [96576094.2145, 106598001.2476, 46210616.7776],
            [-252296.5509, -302841.7334, -93212.772],
        ),
        (
            True,
            jpl.JPLInterp.LINEAR,
            [96419730.22037412, 106719496.94497885, 46263287.23218732],
            [-246849.20349792443, -306551.71943246725, -94755.44468661583],
        ),
        (
            True,
            jpl.JPLInterp.SPLINE,
            [96264116.11254415, 106841825.90246898, 46316319.48198082],
            [-241690.58697019098, -310643.5099432287, -96416.15534852046],
        ),
    ],
)
def test_find_jplde_param(
    sunmooneph_filepath,
    sunmooneph_filepath_12hr,
    include_hr,
    interp,
    rsun_exp,
    rmoon_exp,
):
    jd, jd_frac = 2457884.5, 0.1609116400462963  # (2017, 5, 11, 3, 51, 42.7657)

    # Get the data
    filepath = sunmooneph_filepath_12hr if include_hr else sunmooneph_filepath
    jpldearr, *_ = jpl.read_jplde(filepath, include_hr=include_hr)

    # Call the function
    rsun_out, rmoon_out = jpl.find_jplde_param(jd, jd_frac, jpldearr, interp)

    # Check the outputs
    assert np.allclose(rsun_out, rsun_exp, rtol=DEFAULT_TOL)
    assert np.allclose(rmoon_out, rmoon_exp, rtol=DEFAULT_TOL)


def test_sunmoon(monkeypatch):
    def mock_find_jplde_param(*args):
        return rsun_exp, rmoon_exp

    # Patch the main dependent function
    rsun_exp = [95604355.9737, 107353047.2919, 46537942.1006]
    rmoon_exp = [-218443.5158, -325897.7785, -102799.8515]
    monkeypatch.setattr(jpl, "find_jplde_param", mock_find_jplde_param)

    # Call the function with dummy input values (all used in patched function)
    rsun, rtascs, decls, rmoon, rtascm, declm = jpl.sunmoon(0, 0, {})

    # Check the outputs
    assert np.allclose(rsun, rsun_exp, rtol=DEFAULT_TOL)
    assert np.isclose(rtascs, 0.8432211071356744, rtol=DEFAULT_TOL)
    assert np.isclose(decls, 0.3130881243121989, rtol=DEFAULT_TOL)
    assert np.allclose(rmoon, rmoon_exp, rtol=DEFAULT_TOL)
    assert np.isclose(rtascm, -2.1612978926663886, rtol=DEFAULT_TOL)
    assert np.isclose(declm, -0.25625963750585706, rtol=DEFAULT_TOL)
