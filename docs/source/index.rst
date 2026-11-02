#######
|title|
#######

.. |title| replace:: libdiffpy documentation

``libdiffpy`` - DiffPy C++ library for calculation of PDF and other real-space quantities

| Software version |release|
| Last updated |today|.

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

=======
Authors
=======

``libdiffpy`` is developed by Pavol Juhas, Chris Farrow, Simon Billinge, Billinge Group members. This project is maintained by Simon Billinge. For a detailed list of contributors see
https://github.com/diffpy/libdiffpy/graphs/contributors.

============
Installation
============

See the `README <https://github.com/diffpy/libdiffpy#installation>`_
file included with the distribution.

================
Acknowledgements
================

``libdiffpy`` is built and maintained with `scikit-package <https://scikit-package.github.io/scikit-package/>`_.

=================
Table of contents
=================
.. toctree::
   :maxdepth: 2

   development
   release
   license

=======
Indices
=======

* :ref:`genindex`
* :ref:`search`
