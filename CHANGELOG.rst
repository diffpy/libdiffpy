=============
Release notes
=============

.. current developments

1.4.1
=====

Released 2025-06-25. Notable differences from version 1.4.0.

**Added:**

* ``qbroad_seperable`` attribute to ``JeongPeakWidth`` for a Q-dependent
  peak broadening term, see https://doi.org/10.1107/S0108767391011327.
* Build scripts for Windows.
* GitHub release with the source tarball SHA256 created automatically
  after a tag push.

**Changed:**

* Require libobjcryst 2017.2.1 or later.
* Add the conda environment include and library paths to the build when
  ``CONDA_PREFIX`` is set.
* Do not create ``OverlapCalculator`` self-neighborhoods at sites excluded
  by the pair mask.

**Removed:**

* Runtime check for the ObjCryst inversion center bug, which is fixed in
  libobjcryst 2017.2.1.

**Fixed:**

* ``OverlapCalculator`` now respects its ``rmin`` value.
* Build with ``REAL=double`` to match the floating-point type of libobjcryst.
* Read the installation prefix from the ``PREFIX`` environment variable.

1.4.0
=====

Released 2019-03-09. Notable differences from version 1.3.4.

**Added:**

* Test coverage analysis at https://codecov.io/gh/diffpy/libdiffpy.
* Make scons scripts compatible with Python 3 and Python 2.
* Atomic numbers and possible ions for elements 99-103.
* Define libdiffpy version requirements for client Anaconda packages.
* Support free-standing attribute getter and setter functions.
* ``uisowidth`` and ``bisowidth`` attributes to ``ConstantPeakWidth`` class.
  This simplifies PDF modeling of structures with uniform isotropic
  atom displacements.

**Changed:**

* Build Anaconda package with Anaconda C++ compiler.
* Adopt language standard c++11.
* Replaced boost unordered set and map types with their STL version.
* Switch to direct serialization of ``PeakWidthModel``, ``PeakProfile``,
  ``PDFEnvelope``, ``PDFBaseline``, and ``ScatteringFactorTable`` classes and
  derived types.

**Deprecated:**

* ``StructureAdapter::customPQConfig`` function because it permits unexpected
  changes in calculators setup and complicates parallel evaluation.

**Fixed:**

* Incorrect results for crystals with off-origin inversion center.
* Incomplete ``scons install`` when shared library fails to compile.
* Use of invalid iterators when removing ``map`` items.

**Removed:**

* Virtual class inheritance of ``BaseDebyeSum`` and ``DebyePDFCalculator``.
