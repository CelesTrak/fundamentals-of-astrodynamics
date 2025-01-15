import pytest
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike


DEFAULT_TOL = 1e-12  # default tolerance


def custom_isclose(
    result: float, expected: float, rtol: float = DEFAULT_TOL, atol: float = DEFAULT_TOL
) -> bool | np.ndarray:
    """Compare if result value is close to expected

    Custom function to compare result and expected values with a tolerance
    that scales with the magnitude of the expected value.

    Argss:
        result (float): result value
        expected (float): expected value
        rtol (float, optional): Relative tolerance (defaults to 1e-12)
        atol (float, optional): Absolute tolerance (defaults to 1e-12)

    Returns:
        bool: True if the two values are close within the scaled tolerances
    """
    scale_factor = 10 ** np.floor(np.log10(abs(expected))) if expected != 0 else 1
    scaled_atol = atol * scale_factor
    return np.isclose(result, expected, rtol=rtol, atol=scaled_atol)


def custom_allclose(
    a: ArrayLike, b: ArrayLike, rtol: float = DEFAULT_TOL, atol: float = DEFAULT_TOL
) -> bool:
    """Compare if result array is close to expected

    Custom function to compare arrays `a` and `b` with a tolerance that
    scales with the magnitude of the expected values in `b`

    Args:
        a (array_like): Array of results
        b (array_like): Array of expected values
        rtol (float, optional): Relative tolerance (defaults to 1e-12)
        atol (float, optional): Absolute tolerance (defaults to 1e-12)

    Returns:
        bool: True if all elements are close within the scaled tolerances
    """
    a = np.asarray(a)
    b = np.asarray(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale_factors = np.where(b != 0, 10 ** np.floor(np.log10(np.abs(b))), 1)
    scaled_atol = atol * scale_factors
    return np.all(np.isclose(a, b, rtol=rtol, atol=scaled_atol))


@pytest.fixture
def data_dir() -> Path:
    """Fixture providing the path to the data directory at the root of the repository.

    Returns:
        Path: The path to the `datalib` directory
    """
    return Path(__file__).resolve().parents[3] / "datalib"
