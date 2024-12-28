import pytest

import src.valladopy.mathtime.calendar as cal


class TestInitializeTime:
    def test_non_leap_year(self):
        """Test initialize_time for a non-leap year."""
        time_data = cal.initialize_time(2023)  # 2023 is a non-leap year
        assert time_data["lmonth"] == [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        assert time_data["monthtitle"] == [
            "jan",
            "feb",
            "mar",
            "apr",
            "may",
            "jun",
            "jul",
            "aug",
            "sep",
            "oct",
            "nov",
            "dec",
        ]
        assert time_data["daytitle"] == [
            "sun",
            "mon",
            "tue",
            "wed",
            "thu",
            "fri",
            "sat",
        ]

    def test_leap_year(self):
        """Test initialize_time for a leap year."""
        time_data = cal.initialize_time(2024)  # 2024 is a leap year
        assert time_data["lmonth"] == [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


class TestGetIntMonth:
    @pytest.mark.parametrize(
        "month_str, month_int_exp",
        [
            ("jan", 1),
            ("feb", 2),
            ("mar", 3),
            ("apr", 4),
            ("may", 5),
            ("jun", 6),
            ("jul", 7),
            ("aug", 8),
            ("sep", 9),
            ("oct", 10),
            ("nov", 11),
            ("dec", 12),
        ],
    )
    def test_valid_months(self, month_str, month_int_exp):
        assert cal.get_int_month(month_str) == month_int_exp

    @pytest.mark.parametrize(
        "month_str, month_int_exp", [("JAN", 1), ("Feb", 2), ("DeC", 12)]
    )
    def test_case_insensitivity(self, month_str, month_int_exp):
        """Test that the function is case-insensitive."""
        assert cal.get_int_month(month_str) == month_int_exp

    @pytest.mark.parametrize("month_str", ["xyz", ""])
    def test_invalid_months(self, month_str):
        """Test invalid month abbreviations."""
        with pytest.raises(ValueError, match="Invalid month abbreviation:"):
            cal.get_int_month(month_str)


class TestGetIntDay:
    @pytest.mark.parametrize(
        "day_str, day_int_exp",
        [
            ("sun", 1),
            ("mon", 2),
            ("tue", 3),
            ("wed", 4),
            ("thu", 5),
            ("fri", 6),
            ("sat", 7),
        ],
    )
    def test_valid_days(self, day_str, day_int_exp):
        """Test valid day abbreviations."""
        assert cal.get_int_day(day_str) == day_int_exp

    @pytest.mark.parametrize(
        "day_str, day_int_exp", [("SUN", 1), ("Mon", 2), ("Fri", 6)]
    )
    def test_case_insensitivity(self, day_str, day_int_exp):
        """Test case-insensitivity for day abbreviations."""
        assert cal.get_int_day(day_str) == day_int_exp

    @pytest.mark.parametrize("day_str", ["xyz", ""])
    def test_invalid_day(self, day_str):
        """Test invalid day abbreviations."""
        with pytest.raises(ValueError, match="Invalid day abbreviation:"):
            cal.get_int_day(day_str)


@pytest.mark.parametrize(
    "year, days, mdh_exp",
    [
        (2023, 1.0, (1, 1, 0, 0, 0.0)),  # Jan 1, 2023, 00:00:00
        (2023, 365.0, (12, 31, 0, 0, 0.0)),  # Dec 31, 2023, 00:00:00
        (2023, 34.567, (2, 3, 13, 36, 28.8)),  # Feb 3, 2023, 13:36:28.800
        (2024, 60.0, (2, 29, 0, 0, 0.0)),  # Feb 29, 2024, 00:00:00 (leap year)
    ],
)
def test_days2mdh(year, days, mdh_exp):
    assert cal.days2mdh(year, days) == mdh_exp


@pytest.mark.parametrize(
    "ymdhms, days_frac_exp",
    [
        ((2023, 1, 1, 0, 0, 0.0), 1.0),  # Jan 1, 2023, 00:00:00
        ((2023, 2, 3, 13, 36, 28.8), 34.567),  # non-leap year, fraction
        ((2024, 2, 29, 6, 30, 30.5), 60.27118634259259),  # leap year, fraction
    ],
)
def test_find_days(ymdhms, days_frac_exp):
    assert cal.find_days(*ymdhms) == days_frac_exp
