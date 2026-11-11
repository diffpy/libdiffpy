/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 ******************************************************************************
 *
 * Internal interface of the calculation kernels, not installed.
 *
 * No standard library includes, see CalculationKernelSIMD.ipp.
 *
 *****************************************************************************/

#ifndef CALCULATIONKERNELIMPL_HPP_INCLUDED
#define CALCULATIONKERNELIMPL_HPP_INCLUDED

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

/// Add one Gaussian peak to the PDFCalculator RDF array.
///
/// For each i in [ilo, ihi) with x = (rlosteps + i) * rstep - dist:
///     out[i] += scale * G(x, fwhm) * (x / dist + 1)
/// where G(x, fwhm) = peakheight * exp(-4 * ln(2) * (x / fwhm)^2), or zero
/// for |x / fwhm| >= halfboundrel when cropped.  peakheight is the profile
/// value at x = 0.  Nothing is added for fwhm <= 0.
struct GaussianWindowArgs {
  double* out;
  int ilo;
  int ihi;
  int rlosteps;
  double rstep;
  double dist;
  double fwhm;
  double scale;
  double peakheight;
  bool cropped;
  double halfboundrel;
};

/// Add one pair to the BaseDebyeSum F(Q) array.
///
/// For kq = kqlo, kqlo + 1, ... while kq < kqhi:
///     q = kq * qstep
///     t = scale * exp(-0.5 * (dwsigma * q)^2) * sf0[kq] * sf1[kq] / dist
///     stop at the first kq where |t| <= prec
///     out[kq] += t * sin(q * dist)
struct DebyeSumArgs {
  double* out;
  int kqlo;
  int kqhi;
  double qstep;
  double dist;
  double dwsigma;
  double scale;
  double prec;
  const double* sf0;
  const double* sf1;
};

/// Select candidate bonds from n points by their squared distance.
///
/// For each i in [0, n) with
///     d2 = (x[i] + ox)^2 + (y[i] + oy)^2 + (z[i] + oz)^2
/// append i to idx when d2lo <= d2 <= d2hi or d2 is NaN.
/// Return the count.  idx must hold n values.
/// Rounding of d2 differs between kernels.
struct BondCandidatesArgs {
  const double* x;
  const double* y;
  const double* z;
  int n;
  double ox;
  double oy;
  double oz;
  double d2lo;
  double d2hi;
  int* idx;
};

struct CalculationKernel {
  /// Name for setCalculationKernel, a string literal.
  const char* name;
  void (*gaussianWindow)(const GaussianWindowArgs&);
  void (*debyeSum)(const DebyeSumArgs&);
  int (*bondCandidates)(const BondCandidatesArgs&);
};

/// The active kernel, exported for the unit tests.
DLL_EXPORT const CalculationKernel& activeCalculationKernel();

// Kernel tables per instruction set.  The SIMD ones exist only when CMake
// compiles them, see cmake/CalculationKernels.cmake.
const CalculationKernel& scalarCalculationKernel();
const CalculationKernel& sse2CalculationKernel();
const CalculationKernel& avx2CalculationKernel();
const CalculationKernel& avx512CalculationKernel();
const CalculationKernel& neonCalculationKernel();

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy

#endif  // CALCULATIONKERNELIMPL_HPP_INCLUDED
