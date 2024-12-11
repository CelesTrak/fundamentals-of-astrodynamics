import numpy as np
import pytest

import src.valladopy.astro.perturbations.utils as utils

from ...conftest import DEFAULT_TOL


@pytest.mark.parametrize(
    "x, expected_result",
    [
        # Test scalar input
        (0.5, -0.4375),
        # Test vector input
        ([0.5, 0.6, 0.7], [-0.4375, -0.36, -0.1925]),
        # Test matrix input
        (
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]],
            [
                [-0.1475, -0.28, -0.3825],
                [-0.44, -0.4375, -0.36],
                [-0.1925, 0.08, 0.4725],
            ],
        ),
    ],
)
def test_leg_poly(x, expected_result):
    # Call leg_poly method
    result = utils.leg_poly(x, order=3)

    # Check results
    assert np.allclose(result, expected_result, rtol=DEFAULT_TOL)
