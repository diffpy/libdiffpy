**Added:**

* Ninja presets for release, debug and coverage builds.
* An installed CMake package for consumers using ``find_package(libdiffpy)``
  and the ``libdiffpy::diffpy`` target.

**Changed:**

* Building from source requires CMake 3.22 or later. Use ``cmake --install``
  to install the shared library, headers and runtime data.
* Building the unit tests requires Catch2 3.5 or later. CTest runs each
  test case separately.

**Removed:**

* The SCons build scripts and their ``sconsvars.py`` and ``sconscript.local``
  configuration files.
