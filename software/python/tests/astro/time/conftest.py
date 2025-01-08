import os
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
def eoparr():
    """Load the Earth Orientation Parameters"""
    return data.readeop(os.path.join(data.DATA_DIR, "EOP-All-v1.1_2023-01-01.txt"))
