[![Build Status](https://travis-ci.org/diffpy/libdiffpy.svg?branch=master)](https://travis-ci.org/diffpy/libdiffpy)
[![codecov](https://codecov.io/gh/diffpy/libdiffpy/branch/master/graph/badge.svg)](https://codecov.io/gh/diffpy/libdiffpy)

# libdiffpy

C++ calculators of PDF, bond valence sum and other pair quantities

libdiffpy is a C++ library for calculating atomic pair distribution function
(PDF), bond valence sums, atom overlaps for a hard-sphere model, bond
distances and directions up to a specified maximum distance.   The atomic
structure models are represented by classes for non-periodic, periodic or
structures with space group symmetries.  libdiffpy supports Crystal and
Molecule classes from the [ObjCryst crystallographic library](
https://sourceforge.net/projects/objcryst).
Calculators support two evaluation models - BASIC, which performs a full
pair-summation every time, and OPTIMIZED, which updates only pair
contributions that have changed since the last evaluation.  libdiffpy supports
object serialization and parallel computations using parallel map function.
PDF calculations can be performed in two modes - either as a real-space
summation of peak profiles (PDFCalculator) or as a reciprocal-space Debye
summation and Fourier transform of the total scattering structure function
(DebyePDFCalculator).

The calculator objects in libdiffpy share a common procedure for iteration
over atom pairs and only specialize the processing of pair contributions.
New calculator classes can thus be readily defined for any quantity that is
obtained by iteration over atom pairs.

For more information see user manual at
http://www.diffpy.org/doc/libdiffpy.


## REQUIREMENTS

libdiffpy library requires C++ compiler and the following software:

* `scons` - software constructions tool (1.0 or later)
* `libboost-dev` - Boost C++ libraries development files (1.43 or later)
* `GSL` - GNU Scientific Library for C

Recommended software:

* `libobjcryst` - C++ library of free objects for crystallography,
  https://github.com/diffpy/libobjcryst
* `cxxtest` - CxxTest Unit Testing Framework, http://cxxtest.com

The required software is commonly available in the system package manager.
For example, on Ubuntu Linux the required software can be installed using

```sh
sudo apt-get install \
    build-essential scons libboost-dev libgsl0-dev
```

libdiffpy is also available as a pre-compiled package for
[Anaconda Python](https://www.continuum.io/downloads).  A Python
interface to libdiffpy is provided by Anaconda package diffpy.srreal.


## INSTALLATION

### Installation from sources

Use sources from the git repository or extract the latest source
bundle from https://github.com/diffpy/libdiffpy/releases/latest.

```sh
tar xzf libdiffpy-VERSION.tar.gz
cd libdiffpy-VERSION
```

To build and install the libdiffpy library use

```sh
sudo scons -j4 install
```

This installs libdiffpy for all users under the `/usr/local` directory.
If administrator (root) access is not available, see the output from
`scons --help` for options to install to a user-writable location.

To verify libdiffpy installation, compile and run the included
test code [examples/testlib.cpp](examples/testlib.cpp)

```sh
cd examples
c++ testlib.cpp -ldiffpy
./a.out
```

If compilation fails because of missing header files or missing libdiffpy
library, adjust the `CPATH` and `LIBRARY_PATH` environment variables or
use the `-I` or `-L` compiler options.  If the shared library libdiffpy
is unavailable at runtime, add a `-Wl,-rpath,SomePath` option to the
c++ command or adjust the `LD_LIBRARY_PATH` environment variable.

### Installation for Anaconda Python

The libdiffpy library can be installed from the "diffpy" channel
of Anaconda packages

```sh
conda config --add channels diffpy
conda install libdiffpy
```

libdiffpy is also included in the "diffpy-cmi" collection of packages
for structure analysis

```sh
conda install diffpy-cmi
```

When compiling with the Anaconda version of libdiffpy it is essential to
specify header path, library path and runtime path of the active Anaconda
environment

```sh
# resolve prefix directory P of the active Anaconda environment
P="$(conda info --json | grep default_prefix | cut -d\" -f4)"
cd examples
c++ testlib.cpp -I$P/include -L$P/lib -Wl,-rpath,$P/lib -ldiffpy
./a.out
```

On Mac OS X the libdiffpy package is built for OS X version
10.7 which may be incompatible with codes emitted on newer OS.
To fix this add `-mmacosx-version-min=10.7` option to the
c++ compiler or set it with an environment variable as
`export MACOSX_DEPLOYMENT_TARGET=10.7`.


## DEVELOPMENT

libdiffpy is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.  The libdiffpy sources are hosted at
https://github.com/diffpy/libdiffpy.

Feel free to fork the project and contribute.  When developing it is
preferable to compile with `build=develop` option, which compiles the
library with debug information and C-assertions checks:

```sh
scons -j4 build=develop
```

The build script checks for a presence of `sconsvars.py` file, which
can be used to permanently set the `build` variable.  The SCons
construction environment can be further customized in a `sconscript.local`
script.  The library integrity can be verified by executing unit tests with
`scons -j4 test` (requires the CxxTest framework).


## CONTACTS

For more information about libdiffpy please visit the project web-page

http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu.
