:tocdepth: -1

.. index:: development

=================
Development guide
=================

This page describes how to build and test libdiffpy from sources when
modifying the library. See the
`README <https://github.com/diffpy/libdiffpy#support-and-contribute>`_
for how to fork the repository and set up the pre-commit hooks.

Set up the build environment
----------------------------

From the root of your libdiffpy clone, create a conda environment with the
build and test dependencies listed in the ``requirements`` folder ::

        conda create -n libdiffpy_env --file requirements/build.txt --file requirements/conda.txt --file requirements/tests.txt
        conda activate libdiffpy_env

Build the library
-----------------

Build the shared library with ::

        cmake --preset release
        cmake --build --preset release --parallel 4

The CMake presets require CMake 3.22 or later and use Ninja. On Windows,
run them in a command prompt with the conda build environment activated.
The build presets are

* ``release`` - optimized build with C assertion checks disabled
* ``debug`` - unoptimized build with debug information and C assertion checks
* ``coverage`` - unoptimized build instrumented for coverage analysis

The library is built in the ``build/<preset>`` directory, for example
``build/release/``. Use ``cmake --install build/release`` to install it.
Pass ``-DCMAKE_INSTALL_PREFIX=<prefix>`` when configuring to choose the
installation directory.

To use the same options for every build, set them in a ``CMakeUserPresets.json``
file in the repository root. To use another generator, configure a separate
directory with ``cmake -S . -B <build-dir> -G "<generator>"``. To disable the
ObjCryst adapter, pass ``-DLIBDIFFPY_ENABLE_OBJCRYST=OFF``. For gprof
instrumentation, use ``-DLIBDIFFPY_PROFILE=ON``.

Run the unit tests
------------------

The unit tests use the CxxTest framework. Build and run all of them with ::

        cmake --preset debug
        cmake --build --preset debug --parallel 4
        ctest --preset debug

To run only some of the tests, pass fixed-string patterns of the test file
names with the ``LIBDIFFPY_TESTS`` option. Patterns are case-insensitive and can be
separated by commas ::

        cmake --preset debug -DLIBDIFFPY_TESTS=PDFCalculator,Lattice
        cmake --build --preset debug --target check

Set ``LIBDIFFPY_TESTS`` to an empty string to build all suites again.

To run the tests against the installed library instead of the one in the
``build`` directory, configure a separate directory with
``-DLIBDIFFPY_TEST_INSTALLED=ON``. If CMake cannot locate the installation,
add its prefix to ``CMAKE_PREFIX_PATH``.

Measure test coverage
---------------------

Build with the ``coverage`` preset, run the tests and summarize the results
with ``gcovr`` ::

        cmake --preset coverage
        cmake --build --preset coverage --parallel 4
        ctest --preset coverage
        gcovr --root . --filter src/diffpy/

Use ``cmake --build --preset coverage --target zerocounters`` to reset the
coverage counts between runs.

Test changes with diffpy.srreal
-------------------------------

`diffpy.srreal <https://github.com/diffpy/diffpy.srreal>`_ compiles the
libdiffpy sources from the commit pinned in its ``extern/libdiffpy`` submodule.
See the diffpy.srreal documentation for its build instructions.

Create a source distribution
----------------------------

From a Git checkout, run ``cmake -P cmake/MakeSdist.cmake`` to package
``HEAD`` as ``libdiffpy-<version>.tar.gz`` in the repository root. To write
the archive into a configured build directory, build the ``sdist`` target.
The archive includes the version header and can be built without Git.
