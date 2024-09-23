import numpy as np
import pytest

from src.valladopy.mathtime.polynomial import quadric


@pytest.mark.parametrize(
    "a, b, c, option, expected",
    [
        (1, -3, 2, "include_imaginary", (2, 0, 1, 0)),
        (1, -3, 2, "real_only", (2, 0, 1, 0)),
        (1, 1, 1, "include_imaginary", (-0.5, np.sqrt(3) / 2, -0.5, -np.sqrt(3) / 2)),
        (1, 1, 1, "real_only", (np.nan, 0, np.nan, 0)),
        (1, 2, 1, "unique_real", (-1, 0, np.nan, 0)),  # Double-root case
    ],
)
def test_quadric(a, b, c, option, expected):
    result = quadric(a, b, c, option)
    for x, y in zip(expected, result):
        if np.isnan(x):
            assert np.isnan(y)
        else:
            assert abs(x - y) < 1e-12
