import pytest

from src.valladopy.astro.time.data import iau06in, iau06in2


@pytest.fixture()
def iau06arr():
    """Load the IAU 2006 data"""
    return iau06in2()


@pytest.fixture()
def iau06data_old():
    """Load the IAU 2006 data"""
    return iau06in()
