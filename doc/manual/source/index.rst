.. libdiffpy documentation master file, created by
   sphinx-quickstart on Sun Feb  9 11:31:20 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


####################################################
libdiffpy's documentation
####################################################

Software version |release|.

Last updated |today|.

libdiffpy - C++ calculators of PDF, bond valence sum and other pair quantities

libdiffpy is a C++ library for calculating atomic pair distribution function
(PDF), bond valence sums, atom overlaps for a hard-sphere model, bond
distances and directions up to specified maximum distance.   The atomic
structure models are represented by classes for non-periodic, periodic or
structures with space group symmetries.  libdiffpy supports Crystal and
Molecule classes from the ObjCryst crystallographic library.  Calculators
support two evaluation models - BASIC, which performs a full pair-summation
every time, and OPTIMIZED, which updates only pair contributions that have
changed since the last evaluation.  libdiffpy supports object serialization
and parallel computations using parallel map function.  PDF calculations can
be performed in two modes - either as a real-space summation of peak profiles
(PDFCalculator) or as a reciprocal-space Debye summation and Fourier transform
of the total scattering structure function (DebyePDFCalculator).

The calculator objects in libdiffpy share common procedure for iteration
over atom pairs and only specialize the processing of pair contributions.
New calculator class can thus be readily defined for any quantity that is
obtained by iteration over atom pairs.


================
Authors
================

This code was written by members of the Billinge Group at Columbia
University and Brookhaven National Laboratory including Pavol Juhás,
Christopher Farrow and Simon Billinge.
For a full list of contributors see
https://github.com/diffpy/libdiffpy/graphs/contributors.

======================================
Installation
======================================

See the `README.rst <https://github.com/diffpy/libdiffpy#requirements>`_
file included with the distribution.

======================================
API and Index
======================================

.. toctree::
   :maxdepth: 3
   
* :ref:`genindex`
* :ref:`search`
