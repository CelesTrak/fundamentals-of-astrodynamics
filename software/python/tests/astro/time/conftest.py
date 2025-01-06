import pytest

from src.valladopy.astro.time.data import iau80in, iau06in_pnold, iau06in, readxys


@pytest.fixture()
def iau80arr():
    """Load the IAU 1980 data"""
    return iau80in()


@pytest.fixture()
def iau06arr():
    """Load the IAU 2006 data"""
    return iau06in()


@pytest.fixture()
def iau06data_old():
    """Load the IAU 2006 data"""
    return iau06in_pnold()


@pytest.fixture()
def iau06xysarr():
    """Load the IAU 2006 XYS data"""
    return readxys()
