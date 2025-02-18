# Python Code Overview

This directory contains the **Python implementation** of the routines from the project. The Python scripts are structured with additional granularity, using modules and submodules to better organize functionality. The last tab of the [spreadsheet](../../Astro%20Software.xlsx) provides a detailed mapping between the Python modules and their equivalent MATLAB routines.

## Key Differences
- **Function Modeling & Validation**: Python functions were designed based on MATLAB implementations and validated against them.
- **Use of NumPy & SciPy**: Lower-level math and linear algebra functions (e.g., vector/matrix operations, interpolation) were excluded in favor of `numpy` and `scipy`. Some routines were also intentionally omitted based on their utility. Functions that won't be implemented are greyed out on the spreadsheet.
- **Differences in Arguments & Returns**: Some function signatures differ (e.g. using enums and booleans instead of strings) for clarity. Refer to the docstrings for details on function usage.

## Installation
You can install this package via pip:
```bash
pip install valladopy
```

## Development & Contribution Guidelines
- **Code Style**: Follow **PEP8** rules, use **Google-style docstrings**, and format code with `flake8` and `black` (see [Contribution Guidelines](../../CONTRIBUTING.md) for full instructions).
- **Testing**: Use `pytest` and ensure all tests pass before submitting changes:
  ```bash
  python -m pytest -v --disable-warnings tests/
  ```
- **Questions?** If you have any issues or need clarification, feel free to reach out to [Samira](https://github.com/samotiwala), the original author of the Python implementation.
