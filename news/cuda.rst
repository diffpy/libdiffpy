**Added:**

* Build with ``LIBDIFFPY_CUDA=ON`` or either CUDA preset (``debug-cuda``,
  ``release-cuda``) to enable the ``cuda`` and ``cuda-mixed`` kernels for full
  PDF and Debye sums. Query build support with ``libdiffpy_HAS_CUDA`` (CMake)
  or ``DIFFPY_HAS_CUDA`` (C++).
  CUDA-enabled builds also run on CPU-only machines.
* Derived calculators can implement ``addAllPairContributions`` or reuse
  ``addAllPairsOfPDF`` and ``addAllPairsOfDebyeSum``. The constant, Debye-Waller
  and Jeong peak width models now provide ``calculateFromMSD``.

**Changed:**

* ``PDFCalculator`` and ``DebyePDFCalculator`` accelerate full sums with
  batched pairs and multiple CPU threads, up to the hardware concurrency
  reported by the system. Changing the summation order can affect rounding.
* ``OPTIMIZED`` can choose a full calculation when it costs less than updating
  changed sites. ``getEvaluatorTypeUsed`` then returns ``BASIC``.
* ``PDFCalculator::getQmax`` and ``getQstep`` now return ``double`` by value.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Separate calculators can run concurrently and share a read-only structure.
* ``DebyePDFCalculator`` uses each site's own occupancy when several sites
  share an element.
* The ``OPTIMIZED`` evaluator recalculates when a change of adapter class
  changes the structure, for example when a cluster becomes periodic.
* Copies of a calculator now have independent evaluators.

**Security:**

* <news item>
