**Added:**

* ``diffpy/hash.hpp`` with ``hash_combine``, ``hash_range``, ``hash_values``
  and ``PairHash``.
* SIMD calculation kernels for ``PDFCalculator`` and ``DebyePDFCalculator``:
  AVX-512, AVX2 and SSE2 on x86-64, NEON on arm64. They are built when the
  header-only xsimd library is found, unless the ``LIBDIFFPY_SIMD`` CMake
  option is ``OFF``. The best kernel the CPU supports is used at runtime;
  ``availableCalculationKernels``, ``getCalculationKernel`` and
  ``setCalculationKernel`` select it by name.
* ``PeakProfile::addRDFPeak`` to add a peak to an RDF array, and protected
  ``GaussianProfile::kernelShape``, which derived profiles that keep
  ``operator()`` can override to use the calculation kernels.

**Changed:**

* ``R3::Vector`` and ``R3::Matrix`` are fixed-size ``std::array`` types,
  zero-initialized by default. Existing archives and pickles still load.
* ``R3::row`` and ``R3::column`` return copies, so code such as
  ``R3::row(A, 0) = v`` no longer compiles; use the new ``R3::setrow`` and
  ``R3::setcolumn``.
* ``R3::inverse``, ``R3::floor``, ``R3::mxvecproduct`` and the ``Lattice``
  methods ``cartesian``, ``fractional``, ``ucvCartesian``, ``ucvFractional``,
  ``cartesianMatrix`` and ``fractionalMatrix`` return by value instead of a
  reference to a shared static result, so later calls no longer overwrite
  earlier results. They are thread-safe, as are ``Lattice::distance`` and
  ``Lattice::ucMaxDiagonal``.
* Faster Gaussian peak evaluation and bond generation for periodic and
  crystal structures. ``CrystalStructureAdapter`` expands symmetry faster
  and only after its atoms, lattice or symmetry operations change.

**Removed:**

* Dependencies on Boost.Functional/Hash, Boost.StringAlgo and Boost.uBLAS.
  ``R3::Vector`` and ``R3::Matrix`` no longer accept uBLAS expressions.
* Protected ``updater1`` method and ``msphere``, ``mrcsphere``, ``msymidx``,
  ``mcstructure`` and ``mpuc1`` members of ``PeriodicStructureBondGenerator``
  and ``CrystalStructureBondGenerator``, which no longer override
  ``iterateSymmetry``.

**Fixed:**

* ``R3::inverse`` throws ``invalid_argument`` for a singular matrix instead
  of aborting the program.
* ``Lattice::setLatPar`` and ``Lattice::setLatBase`` leave the lattice
  unchanged when they throw for degenerate parameters or base vectors.
* ``CroppedGaussianProfile`` peaks of zero width no longer produce NaN in
  ``PDFCalculator`` results.
* Better accuracy for narrow Gaussian peaks at large pair distances.
* Calculators give the correct result on the first evaluation after
  a ``CrystalStructureAdapter`` atom moves in place to or from a special
  position.
* After an evaluation fails with an exception, the ``OPTIMIZED`` evaluator
  recalculates fully instead of updating the incomplete result.
