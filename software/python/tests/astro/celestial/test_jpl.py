import numpy as np
import pytest
from pathlib import Path

import src.valladopy.astro.celestial.jpl as jpl

from ...conftest import DEFAULT_TOL


@pytest.fixture
def sunmooneph_filepath(data_dir):
    filepath = Path(data_dir) / "sunmooneph_430t.txt"
    assert filepath.exists()
    return filepath


@pytest.fixture
def sunmooneph_filepath_12hr(data_dir):
    filepath = Path(data_dir) / "sunmooneph_430t12.txt"
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
            35838,
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
    "interp, rsun_exp, rmoon_exp",
    [
        (
            jpl.JPLInterp.NONE,
            [98498799.3378, 105065115.0082, 45546069.1442],
            [-311173.3387, -246845.7203, -71047.7913],
        ),
        (
            jpl.JPLInterp.LINEAR,
            [98189413.7028805, 105311774.2471488, 45653002.59384974],
            [-301699.37820820103, -255856.13060991265, -74614.39469838185],
        ),
        (
            jpl.JPLInterp.SPLINE,
            [98191299.56964143, 105313814.22568269, 45653888.30512573],
            [-302576.59863313846, -256669.07629020023, -74853.18748977917],
        ),
    ],
)
def test_find_jplde_param(sunmooneph_filepath, interp, rsun_exp, rmoon_exp):
    jd, jd_frac = 2457884.5, 0.1609116400462963

    # Get the data
    jpldearr, jdjpldestart, _ = jpl.read_jplde(sunmooneph_filepath, include_hr=False)

    # Call the function
    rsun_out, rmoon_out = jpl.find_jplde_param(
        jd, jd_frac, jpldearr, jdjpldestart, interp
    )

    # Check the outputs
    assert np.allclose(rsun_out, rsun_exp, rtol=DEFAULT_TOL)
    assert np.allclose(rmoon_out, rmoon_exp, rtol=DEFAULT_TOL)
