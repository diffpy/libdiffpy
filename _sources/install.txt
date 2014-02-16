Installation
========================================================================

.. index:: Requirements

libdiffpy library requires C++ compiler and the following software:

* `scons` - software constructions tool (1.0 or later)
* `libboost-dev` - Boost C++ libraries development files (1.43 or later)
* `GSL` - GNU Scientific Library for C

Recommended software:

* `libObjCryst` - C++ library of free objects for crystallography, FIXME-URL
* `cxxtest` - CxxTest Unit Testing Framework, http://cxxtest.com/

Required software is usually available in the system package manager,
for example, on Ubuntu Linux the dependencies can be installed as

    sh
    sudo apt-get install \
    build-essential scons libboost-dev libgsl0-dev

For Mac OS X machine with the MacPorts package manager one could do

    sh
    sudo port install scons boost gsl


For other packages see their project pages for installation instructions.


.. index:: Installation

Installation
------------------------------------------------------------------------


To build and install the libdiffpy library from sources, run

    sh
    sudo scons -j4 install

This installs libdiffpy for all users under the `/usr/local` directory.
If administrator (root) access is not available, see the usage info from
`scons --help` for options to install to a user-writable location.


.. index:: Contribution

Contribution
------------------------------------------------------------------------


libdiffpy is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.  The libdiffpy sources are hosted at
https://github.com/diffpy/libdiffpy.

Feel free to fork the project and contribute.  When developing it is
preferable to compile with `build=develop` option, which compiles the
library with debug information and C-assertions checks:

    sh
    scons -j4 build=develop

The build script checks for a presence of `sconsvars.py` file, which
can be used to permanently set the `build` variable.  The SCons
construction environment can be further customized in a `sconscript.local`
script.  The library integrity can be verified by executing unit tests with
`scons -j4 test` (requires the CxxTest framework).