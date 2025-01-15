import numpy as np
import pytest

import src.valladopy.astro.twobody.hills as hills
import src.valladopy.constants as const
from ...conftest import DEFAULT_TOL, custom_allclose


@pytest.fixture
def alt_tgt():
    return 590  # target satellite altitude in km


@pytest.fixture
def dts():
    return 300  # desired time in seconds


@pytest.mark.parametrize(
    "r, v, rint_exp, vint_exp",
    [
        (
            [500 + const.RE, 0, 0],
            [0, np.sqrt(const.MU / (500 + const.RE)), 0],
            [8699.5464562552, 1886.98432911071, 0],
            [12.035252716871558, 3.658653928010125, 0],
        ),
        (
            # Vallado 2007, Example 6-14
            [0, 0, 0],
            [-0.1, -0.04, -0.02],
            [-33.345724432187964, -1.473560307487475, -5.894530311607677],
            [-0.120337093866688, 0.032387584447866, -0.018949031717993],
        ),
    ],
)
def test_hillsr(alt_tgt, dts, r, v, rint_exp, vint_exp):
    rint, vint = hills.hillsr(r, v, alt_tgt, dts)
    assert np.allclose(rint, rint_exp, rtol=DEFAULT_TOL)
    assert np.allclose(vint, vint_exp, rtol=DEFAULT_TOL)


def test_hillsv(alt_tgt, dts):
    r = [-70.933, 20.357, -11.17]
    vint_exp = [0.274190838333, 0.013466718841345, 0.035907981315]
    vint = hills.hillsv(r, alt_tgt, dts)
    assert custom_allclose(vint, vint_exp)
