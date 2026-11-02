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

        scons -j4

The ``build`` option selects the compiler settings, for example
``scons -j4 build=develop``. On Linux and macOS the build variants are

* ``fast`` - optimized build with C assertion checks disabled (default)
* ``debug`` - build with debug information and C assertion checks
* ``develop`` - unoptimized build with C assertion checks
* ``coverage`` - unoptimized build instrumented for coverage analysis

The library is built in the ``build/<variant>-<machine>`` directory,
for example ``build/fast-x86_64/``. Run ``scons --help`` for the full list
of build targets and options.

To use the same options for every build, set them in a ``sconsvars.py``
file in the repository root, for example ::

        build = 'develop'

The SCons construction environment can be further customized in a
``sconscript.local`` file in the repository root. Both files are ignored
by git. For example, the following ``sconscript.local`` builds with the
conda-forge compiler and records the conda library directory as the runtime
path, so that the unit tests find the conda libraries ::

        Import('env')
        import os
        prefix = os.environ['CONDA_PREFIX']
        env.Replace(CXX=os.path.join(prefix, 'bin', 'c++'))
        env.Prepend(CCFLAGS=['-isystem' + os.path.join(prefix, 'include')])
        env.Append(LIBPATH=os.path.join(prefix, 'lib'))
        env.Append(LINKFLAGS='-Wl,-rpath,' + os.path.join(prefix, 'lib'))

Run the unit tests
------------------

The unit tests use the CxxTest framework. Build and run all of them with ::

        scons -j4 test

To run only some of the tests, pass fixed-string patterns of the test file
names with the ``tests`` option. Patterns are case-insensitive and can be
separated by commas ::

        scons -j4 test tests=PDFCalculator,Lattice

To run the tests against the installed library instead of the one in the
``build`` directory, use ``scons test test_installed=True``.

Measure test coverage
---------------------

Build with the ``coverage`` variant, run the tests and summarize the results
with ``gcovr`` ::

        scons -j4 build=coverage test
        gcovr --root . --filter 'build/coverage-.*/diffpy/'

Use ``scons zerocounters`` to reset the coverage counts between runs.

Test changes with diffpy.srreal
-------------------------------

`diffpy.srreal <https://github.com/diffpy/diffpy.srreal>`_ compiles the
libdiffpy sources from the commit pinned in its
``cmake/LibdiffpyVersion.cmake`` file. To build diffpy.srreal against a local
libdiffpy clone, point the CMake ``FETCHCONTENT_SOURCE_DIR_LIBDIFFPY``
variable to it. diffpy.srreal checks that the clone is at the pinned commit,
so update the pin in ``cmake/LibdiffpyVersion.cmake`` when testing new
libdiffpy commits. See the diffpy.srreal documentation for its build
instructions.
