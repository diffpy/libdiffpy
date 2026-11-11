**Added:**

* ``diffpy/hash.hpp`` with ``hash_combine``, ``hash_range``, ``hash_values``
  and ``PairHash``.

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

**Removed:**

* Dependency on Boost.Functional/Hash, Boost.StringAlgo and Boost.uBLAS.
  ``R3::Vector`` and ``R3::Matrix`` no longer accept uBLAS expressions.

**Fixed:**

* ``R3::inverse`` throws ``invalid_argument`` for a singular matrix instead
  of aborting the program.
* ``Lattice::setLatPar`` and ``Lattice::setLatBase`` leave the lattice
  unchanged when they throw for degenerate parameters or base vectors.
