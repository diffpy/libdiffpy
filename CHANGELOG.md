# Release notes

## Unreleased - Version 1.4.0

Notable differences from version 1.3.4.

### Added

- Test coverage analysis at https://codecov.io/gh/diffpy/libdiffpy.
- Make scons scripts compatible with Python 3 and Python 2.
- Atomic numbers and possible ions for elements 99-103.
- Define libdiffpy version requirements for client Anaconda packages.
- Support free-standing attribute getter and setter functions.
- `uisowidth` and `bisowidth` attributes to `ConstantPeakWidth` class.
  This simplifies PDF modeling of structures with uniform isotropic
  atom displacements.


### Changed

- Build Anaconda package with Anaconda C++ compiler.
- Adopt language standard c++11.
- Replaced boost unordered set and map types with their STL version.
- Switch to direct serialization of `PeakWidthModel`, `PeakProfile`,
  `PDFEnvelope`, `PDFBaseline`, and `ScatteringFactorTable` classes and
  derived types.

### Deprecated

- `StructureAdapter::customPQConfig` function because it permits unexpected
  changes in calculators setup and complicates parallel evaluation.

### Removed

- Virtual class inheritance of `BaseDebyeSum` and `DebyePDFCalculator`.

### Fixed

- Incorrect results for crystals with off-origin inversion center.
- Incomplete `scons install` when shared library fails to compile.
- Use of invalid iterators when removing `map` items.
