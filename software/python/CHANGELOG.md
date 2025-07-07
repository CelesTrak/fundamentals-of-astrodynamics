# Changelog

## [0.3.0] - 2025-07-02
- Bugfix ([Issue #127](https://github.com/CelesTrak/fundamentals-of-astrodynamics/issues/127)) + updates for Lambert minimum time function
- Minor updates to time data function signatures
- Added GPS time output to `convtime` routine
- Added gravity acceleration functions (Gottlieb, Lear, GTDS, Montenbruck, and Pines)

## [0.2.0] - 2025-03-10
- Updated functions to match current MATLAB routines:
  - IOD angles (`laplace` and `gauss`)
  - Minor name refactoring in IAU functions
  - Add `rv2sez` and `sez2rv` frame conversions, and refactor other function names
  - Add Lambert function `tmax_rp` (`lamberttmaxrp` in MATLAB) and update `battin` function
  - Add `checkhitearthc` (canonical version of `checkhitearth`)
  - Update gravity functions
  - Bugfix for `rv2coe` ([Issue #125](https://github.com/CelesTrak/fundamentals-of-astrodynamics/issues/125))

## [0.1.1] - 2025-02-09
- Added separate description for PyPI.

## [0.1.0] - 2025-02-01
- Initial release
