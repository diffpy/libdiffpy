libdiffpy - C++ calculators of PDF, bond valence sum and other pair quantities

libdiffpy is a C++ library that provides calculators for atomic pair
distribution function (PDF), bond valence sums, atom overlaps for a
hard-sphere model, bond distances and directions up to specified maximum
distance.   The atomic structure models are represented by classes for
non-periodic, periodic or structures with space group symmetries.  libdiffpy
supports Crystal and Molecule classes from the ObjCryst crystallographic
library.  Calculators support two evaluation models - BASIC, which performs
a full pair-summation every time, and OPTIMIZED, which updates only
pair contributions that have changed since the last evaluation.  libdiffpy
provides support for object serialization and an entry points for parallel
computations using a parallel map function.  PDF calculations can be done in
two modes - either as a real-space summation of peak profiles (PDFCalculator)
or as a reciprocal-space Debye summation and Fourier transform of the total
scattering structure function (DebyePDFCalculator).

The calculator classes in libdiffpy share common procedure for iteration
over atom pairs and only specialize the processing of pair contributions.
A new calculator class can thus be readily defined for any quantity that is
obtained by iteration over atom pairs.

For more information about the libdiffpy library, see users manual at

    http://www.diffpy.org/doc/libdiffpy/


REQUIREMENTS

libdiffpy library requires C++ compiler and the following software:

    scons        -- software constructions tool (1.0 or later)
    libboost-dev -- Boost C++ libraries development files (1.43 or later)
    GSL          -- GNU Scientific Library for C

Recommended software:

    libObjCryst  -- C++ library of free objects for crystallography, FIXME-URL
    cxxtest      -- CxxTest Unit Testing Framework, http://cxxtest.com/

Some of the required software may be available in the system package manager,
for example, on Ubuntu Linux the dependencies can be installed as:

    sudo apt-get install \
        build-essential scons libboost-dev libgsl0-dev

For Mac OS X machine with the MacPorts package manager one could do

    sudo port install \
        scons boost gsl

For other required packages see their respective web pages for installation
instructions.


INSTALLATION

To build and install the libdiffpy library from sources, run

    sudo scons -j4 install

This installs libdiffpy for all users under the /usr/local directory.
If administrator (root) access is not available, see the usage info from
"scons --help" for options to install to a user-writable location.
The library integrity can be verified by executing unit tests as
(requires the CxxTest framework):

    scons -j4 test


DEVELOPMENT

libdiffpy is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.  The libdiffpy sources are hosted at

    https://github.com/diffpy/libdiffpy

Feel free to fork the project and contribute.  When developing it is
preferable to compile libdiffpy with the "develop" build mode, which
includes debugging information and preserves all C-assertion checks:

    scons -j4 build=develop

The build process checks for a presence of sconsvars.py file, which
can be used to permanently set SCons variables.  The SCons construction
environment can be customized in a sconscript.local script.


CONTACTS

For more information about libdiffpy please visit the project web-page

    http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu.
