|Icon| |title|_
===============

.. |title| replace:: libdiffpy
.. _title: https://diffpy.github.io/libdiffpy

.. |Icon| image:: https://avatars.githubusercontent.com/diffpy
        :target: https://diffpy.github.io/libdiffpy
        :height: 100px

|Forge| |PR|

|CI| |Tracking|

.. |CI| image:: https://github.com/diffpy/libdiffpy/actions/workflows/matrix-and-codecov.yml/badge.svg
        :target: https://github.com/diffpy/libdiffpy/actions/workflows/matrix-and-codecov.yml

.. |Forge| image:: https://img.shields.io/conda/vn/conda-forge/libdiffpy
        :target: https://anaconda.org/conda-forge/libdiffpy

.. |PR| image:: https://img.shields.io/badge/PR-Welcome-29ab47ff
        :target: https://github.com/diffpy/libdiffpy/pulls

.. |Tracking| image:: https://img.shields.io/badge/issue_tracking-github-blue
        :target: https://github.com/diffpy/libdiffpy/issues

C++ calculators of PDF, bond valence sum and other pair quantities.

libdiffpy is a C++ library for calculating atomic pair distribution function
(PDF), bond valence sums, atom overlaps for a hard-sphere model, bond
distances and directions up to a specified maximum distance. The atomic
structure models are represented by classes for non-periodic, periodic or
structures with space group symmetries. libdiffpy supports Crystal and
Molecule classes from the `ObjCryst crystallographic library
<https://sourceforge.net/projects/objcryst>`_.
Calculators support two evaluation models - BASIC, which performs a full
pair-summation every time, and OPTIMIZED, which updates only pair
contributions that have changed since the last evaluation. libdiffpy supports
object serialization and parallel computations using parallel map function.
PDF calculations can be performed in two modes - either as a real-space
summation of peak profiles (PDFCalculator) or as a reciprocal-space Debye
summation and Fourier transform of the total scattering structure function
(DebyePDFCalculator).

The calculator objects in libdiffpy share a common procedure for iteration
over atom pairs and only specialize the processing of pair contributions.
New calculator classes can thus be readily defined for any quantity that is
obtained by iteration over atom pairs.

A Python interface to libdiffpy is provided by the `diffpy.srreal
<https://github.com/diffpy/diffpy.srreal>`_ package.

For more information about the libdiffpy library, please consult our `online documentation <https://diffpy.github.io/libdiffpy>`_.

Citation
--------

If you use libdiffpy in a scientific publication, we would like you to cite this package as

        libdiffpy Package, https://github.com/diffpy/libdiffpy

Installation
------------

The preferred method is to use `Miniconda Python
<https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html>`_
and install from the "conda-forge" channel of Conda packages.

To add "conda-forge" to the conda channels, run the following in a terminal. ::

        conda config --add channels conda-forge

We want to install our packages in a suitable conda environment.
The following creates and activates a new environment named ``libdiffpy_env`` ::

        conda create -n libdiffpy_env libdiffpy
        conda activate libdiffpy_env

When compiling with the conda version of libdiffpy it is essential to
specify the header path, library path and runtime path of the active conda
environment. To verify the installation, compile and run the included
test code `examples/testlib.cpp <examples/testlib.cpp>`_ ::

        cd examples
        c++ testlib.cpp -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -ldiffpy
        ./a.out

Note the conda package for libdiffpy is built with C++ compilers
provided by conda-forge. This may cause incompatibility with system C++.
In such case please use the conda-forge C++ compiler to link with libdiffpy.

Requirements
^^^^^^^^^^^^

To build libdiffpy from sources, a C++ compiler and the following software are required:

* ``scons`` - software constructions tool (1.0 or later)
* ``libboost-dev`` - Boost C++ libraries development files (1.56 or later)
* ``GSL`` - GNU Scientific Library for C

Recommended software:

* ``libobjcryst`` - C++ library for free objects for crystallography,
  https://github.com/diffpy/libobjcryst
* ``cxxtest`` - CxxTest Unit Testing Framework, http://cxxtest.com

The required software is commonly available in the system package manager.
For example, on Ubuntu Linux the required software can be installed using ::

        sudo apt-get install build-essential scons libboost-dev libgsl0-dev

Alternatively, all of the above can be installed into a conda environment from
the ``requirements`` folder ::

        conda create -n libdiffpy_env --file requirements/build.txt --file requirements/conda.txt --file requirements/tests.txt
        conda activate libdiffpy_env

Installation from sources
^^^^^^^^^^^^^^^^^^^^^^^^^

Use sources from the git repository or extract the latest source
bundle from https://github.com/diffpy/libdiffpy/releases/latest. ::

        tar xzf libdiffpy-VERSION.tar.gz
        cd libdiffpy-VERSION

To build and install the libdiffpy library use ::

        sudo scons -j4 install

This installs libdiffpy for all users under the ``/usr/local`` directory.
If administrator (root) access is not available, see the output from
``scons --help`` for options to install to a user-writable location.

To verify the libdiffpy installation, compile and run the included
test code `examples/testlib.cpp <examples/testlib.cpp>`_ ::

        cd examples
        c++ testlib.cpp -ldiffpy
        ./a.out

If compilation fails because of missing header files or missing libdiffpy
library, adjust the ``CPATH`` and ``LIBRARY_PATH`` environment variables or
use the ``-I`` or ``-L`` compiler options. If the shared library libdiffpy
is unavailable at runtime, add a ``-Wl,-rpath,SomePath`` option to the
c++ command or adjust the ``LD_LIBRARY_PATH`` environment variable.

Getting Started
---------------

You may consult our `online documentation <https://diffpy.github.io/libdiffpy>`_ for tutorials and API references.

Support and Contribute
----------------------

If you see a bug or want to request a feature, please `report it as an issue <https://github.com/diffpy/libdiffpy/issues>`_ and/or `submit a fix as a PR <https://github.com/diffpy/libdiffpy/pulls>`_.

Feel free to fork the project and contribute. When developing it is
preferable to compile with the ``build=develop`` option, which compiles the
library with debug information and C-assertions checks ::

        scons -j4 build=develop

The build script checks for a presence of ``sconsvars.py`` file, which
can be used to permanently set the ``build`` variable. The SCons
construction environment can be further customized in a ``sconscript.local``
script. The library integrity can be verified by executing unit tests with
``scons test`` (requires the CxxTest framework).

See also the `instructions for compiling and modifying source code in development mode <docs/development_instructions.md>`_.

To ensure code quality and to prevent accidental commits into the default branch, please set up the use of our pre-commit
hooks.

1. Install pre-commit in your working environment by running ``conda install pre-commit``.

2. Initialize pre-commit (one time only) ``pre-commit install``.

Thereafter your code will be linted by black and isort and checked against flake8 before you can commit.
If it fails by black or isort, just rerun and it should pass (black and isort will modify the files so should
pass after they are modified). If the flake8 test fails please see the error messages and fix them manually before
trying to commit again.

Improvements and fixes are always appreciated.

Before contributing, please read our `Code of Conduct <https://github.com/diffpy/libdiffpy/blob/main/CODE-OF-CONDUCT.rst>`_.

Contact
-------

For more information on libdiffpy please visit the project `web-page <https://diffpy.github.io/>`_ or email the maintainers ``Simon Billinge(sbillinge@ucsb.edu)``.

Acknowledgements
----------------

libdiffpy is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.

``libdiffpy`` is built and maintained with `scikit-package <https://scikit-package.github.io/scikit-package/>`_.
