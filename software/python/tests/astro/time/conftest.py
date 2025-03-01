import pytest

import src.valladopy.astro.time.data as data


@pytest.fixture()
def iau80arr():
    """Load the IAU 1980 data"""
    return data.iau80in()


@pytest.fixture()
def iau06arr():
    """Load the IAU 2006 data"""
    return data.iau06in()


@pytest.fixture()
def iau06data_old():
    """Load the IAU 2006 data"""
    return data.iau06in_pnold()


@pytest.fixture()
def iau06xysarr():
    """Load the IAU 2006 XYS data"""
    return data.readxys()


@pytest.fixture()
def eoparr(test_data_dir):
    """Load the Earth Orientation Parameters"""
    return data.readeop(test_data_dir / "EOP-All-v1.1_2023-01-01.txt")


@pytest.fixture()
def spwarr(test_data_dir):
    """Load the Space Weather data"""
    return data.readspw(test_data_dir / "SpaceWeather-All-v1.2_2025-01-10.txt")
