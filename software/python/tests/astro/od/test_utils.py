import numpy as np

from src.valladopy.astro.od.utils import finite_diff

from ...conftest import DEFAULT_TOL


def test_finite_diff():
    # Inputs
    xnom = np.array([5975.2904, 2568.64, 3120.5845, 3.983846, -2.071159, -5.917095])
    percentchg = deltaamtchg = 0.01
    pertelem = 0  # element index to perturb

    # Call finite_diff method
    deltaamt, xnomp = finite_diff(pertelem, percentchg, deltaamtchg, xnom)

    # Check results
    assert np.isclose(deltaamt, 59.752904, rtol=DEFAULT_TOL)
    assert np.allclose(xnomp, np.array([6035.043304, *xnom[1:]]), rtol=DEFAULT_TOL)
