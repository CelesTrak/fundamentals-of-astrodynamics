import numpy as np
import pytest

import src.valladopy.mathtime.julian_date as julian_date

from ..conftest import DEFAULT_TOL


def test_jday():
    year, month, day, hour, minute, second = (1992, 8, 20, 12, 14, 33)
    jd, jd_frac = julian_date.jday(year, month, day, hour, minute, second)
    assert np.isclose(jd, 2448854.5, rtol=DEFAULT_TOL)
    assert np.isclose(jd_frac, 0.5101041666666667, rtol=DEFAULT_TOL)


@pytest.mark.parametrize(
    "ymdhms, calendar_type, expected_jd, expected_jd_frac",
    [
        (
            (1989, 1, 24, 20, 35, 23),
            julian_date.CalendarType.JULIAN,
            2447563.5,
            0.8579050925925926,
        ),
        ((2024, 12, 2, 12, 0, 0), julian_date.CalendarType.GREGORIAN, 2460646.5, 0.5),
    ],
)
def test_jday_all(ymdhms, calendar_type, expected_jd, expected_jd_frac):
    jd, jd_frac = julian_date.jday_all(*ymdhms, calendar_type)
    assert np.isclose(jd, expected_jd, rtol=DEFAULT_TOL)
    assert np.isclose(jd_frac, expected_jd_frac, rtol=DEFAULT_TOL)


def test_invjday():
    jd, jd_frac = 2447550.5, 0.8579050925925926
    ymdhms = julian_date.invjday(jd, jd_frac)
    assert ymdhms == (1989, 1, 24, 20, 35, 23)


def test_day_of_week():
    jd = 2448854.125
    dow = julian_date.day_of_week(jd)
    assert dow == 4


@pytest.mark.parametrize(
    "year, exp_startday, exp_stopday, exp_jdstartdst, exp_jdstopdst",
    [
        (2000, 2, 29, 2451636.7916666665, 2451846.7916666665),  # before DST rule change
        (2024, 10, 3, 2460379.7916666665, 2460617.7916666665),  # after DST rule change
    ],
)
def test_daylight_savings(
    year, exp_startday, exp_stopday, exp_jdstartdst, exp_jdstopdst
):
    lon = np.radians(-75)
    startday, stopday, jdstartdst, jdstopdst = julian_date.daylight_savings(year, lon)
    assert startday == exp_startday
    assert stopday == exp_stopday
    assert np.isclose(jdstartdst, exp_jdstartdst, rtol=DEFAULT_TOL)
    assert np.isclose(jdstopdst, exp_jdstopdst, rtol=DEFAULT_TOL)


def test_convtime():
    (
        ut1,
        tut1,
        jdut1,
        jdut1frac,
        utc,
        tai,
        tt,
        ttt,
        jdtt,
        jdttfrac,
        tdb,
        ttdb,
        jdtdb,
        jdtdbfrac,
    ) = julian_date.convtime(
        year=2024,
        month=5,
        day=10,
        hour=12,
        minute=42,
        second=23.4,
        timezone=3,
        dut1=0.1,
        dat=35.7,
    )

    assert np.isclose(ut1, 56543.5, rtol=DEFAULT_TOL)
    assert np.isclose(tut1, 0.24356343432326408, rtol=DEFAULT_TOL)
    assert np.isclose(jdut1, 2460440.5, rtol=DEFAULT_TOL)
    assert np.isclose(jdut1frac, 0.6544386574074074, rtol=DEFAULT_TOL)
    assert np.isclose(utc, 56543.4, rtol=DEFAULT_TOL)
    assert np.isclose(tai, 56579.1, rtol=DEFAULT_TOL)
    assert np.isclose(tt, 56611.284, rtol=DEFAULT_TOL)
    assert np.isclose(ttt, 0.24356345580272443, rtol=DEFAULT_TOL)
    assert np.isclose(jdtt, 2460440.5, rtol=DEFAULT_TOL)
    assert np.isclose(jdttfrac, 0.6552231944444444, rtol=DEFAULT_TOL)
    assert np.isclose(tdb, 56611.28533755204, rtol=DEFAULT_TOL)
    assert np.isclose(ttdb, 0.24356345580314515, rtol=DEFAULT_TOL)
    assert np.isclose(jdtdb, 2460440.5, rtol=DEFAULT_TOL)
    assert np.isclose(jdtdbfrac, 0.6552232099253709, rtol=DEFAULT_TOL)
