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
class PairQuantity;
namespace kernels {

/// Add one Gaussian peak to the PDFCalculator RDF array.
///
/// For each i in [ilo, ihi) with x = xlo + (i - ilo) * rstep, where
/// xlo = (rlosteps + ilo) * rstep - dist:
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

/// Add count equal pairs to the BaseDebyeSum F(Q) array.
///
/// For kq = kqlo, kqlo + 1, ... while kq < kqhi:
///     q = kq * qstep
///     t = scale * exp(-0.5 * (dwsigma * q)^2) * sf0[kq] * sf1[kq] / dist
///     stop at the first kq where |t| <= prec
///     out[kq] += count * t * sin(q * dist)
/// The sum stops where it does for one pair, so that count pairs add
/// count times what one pair adds.
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
  double count = 1.0;
  double occ0 = 1.0, occ1 = 1.0;
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

/// Positions and pairs of a structure for a sum over all its pairs, as
/// PairTable gives them.  The arrays stay with the caller.
struct PairGeometryArgs {
  int cntsites;
  int cntpositions;
  int cnttrans;
  /// number of entries in the list of pairs, unused without translations
  long long cntpairs;
  /// The pairs are those of the list below, of positions with lattice
  /// translations.  Otherwise they are all pairs of two sites.
  int translations;
  /// Positions of site i are [sitebegin[i], sitebegin[i + 1]).  A pair
  /// starts at the first position of its site0.
  const int* sitebegin;
  /// site of each position
  const int* possite;
  /// Cartesian coordinates x, y, z of each position
  const double* pos;
  /// Displacement tensors of the positions and of the sites, each as its
  /// elements xx, yy, zz, xy, xz, yz.
  const double* posuij;
  const double* siteuij;
  const unsigned char* anisotropy;
  const int* multiplicity;
  /// Cartesian coordinates of the lattice translations
  const double* trans;
  /// Site of the first position, second position and translation of each
  /// pair in the list.
  const int* pairsite0;
  const int* pairpos;
  const int* pairtrans;
  /// changes when the list of pairs does
  int pairversion;
  /// The pairs are those of rmin <= distance <= rmax and of a distance
  /// above zeroeps.
  double rmin;
  double rmax;
  double zeroeps;
};

/// Built-in peak width model and its parameters.
struct PairWidthArgs {
  enum Kind { Jeong, DebyeWaller, Constant } kind;
  double delta1;
  double delta2;
  double qbroad;
  double qbroadsep;
  double width;
};

/// Counts used to choose between a CPU batch and an incremental update.
struct PairSumStats {
  long long pairs = 0;
  long long peaks = 0;
};

/// Owned by the calculator and reused by kernels with the same backend.
struct PairSumWorkspace {
  virtual ~PairSumWorkspace() = default;
};
using MakeWorkspace = PairSumWorkspace* (*)();

struct PairSumArgs {
  PairGeometryArgs g;
  PairWidthArgs w;
  double* out;
  PairSumWorkspace* workspace;
  /// A null mask includes every pair.
  const PairQuantity* mask;
};

/// Add the peaks of PDFCalculator for all pairs of a structure, what
/// gaussianWindow adds for each of them.
struct PDFPairSumArgs : PairSumArgs {
  /// scattering factor of each site with its occupancy
  const double* sfsite;
  /// The peaks are added to the npts values at out, whose first one is at
  /// r = rlosteps * rstep.
  int npts;
  int rlosteps;
  double rstep;
  /// height of the peak profile for the width of one
  double height;
  int cropped;
  double halfboundrel;
};

/// Add the terms of BaseDebyeSum for all pairs of a structure, what
/// debyeSum adds for each of them.
struct DebyePairSumArgs : PairSumArgs {
  /// Each site type maps to one scattering curve and an occupancy.
  /// sf contains ncurves consecutive rows of nq unweighted factors.
  const int* typeofsite;
  const int* curve;
  const double* occupancy;
  int ntypes, ncurves;
  int nq;
  const double* sf;
  /// The terms are added to the values of kqlo <= kq < kqhi at out.
  int kqlo;
  int kqhi;
  double qstep;
  double prec;
};

using GaussianWindow = void (*)(const GaussianWindowArgs&);
using DebyeSum = void (*)(const DebyeSumArgs&);

struct CalculationKernel {
  /// Name for setCalculationKernel, a string literal.
  const char* name;
  GaussianWindow gaussianWindow;
  DebyeSum debyeSum;
  int (*bondCandidates)(const BondCandidatesArgs&);
  // Sums over all pairs of a structure.
  // The output is updated only after the sum succeeds.
  PairSumStats (*pdfPairSum)(const PDFPairSumArgs&);
  PairSumStats (*debyePairSum)(const DebyePairSumArgs&);
  MakeWorkspace makeWorkspace;
};

/// The active kernel, exported for the unit tests.
DLL_EXPORT const CalculationKernel& activeCalculationKernel();

PairSumStats cpuPDFPairSum(const PDFPairSumArgs&, GaussianWindow);
PairSumStats cpuDebyePairSum(const DebyePairSumArgs&, DebyeSum);
PairSumWorkspace* makeCPUWorkspace();

// Kernel tables per instruction set.  The SIMD ones exist only when CMake
// compiles them, see cmake/CalculationKernels.cmake.
const CalculationKernel& scalarCalculationKernel();
const CalculationKernel& sse2CalculationKernel();
const CalculationKernel& avx2CalculationKernel();
const CalculationKernel& avx512CalculationKernel();
const CalculationKernel& neonCalculationKernel();

/// Bind both the per-pair and batch entry points to the same instruction set.
template <GaussianWindow Gaussian, DebyeSum Debye>
constexpr CalculationKernel makeCPUKernel(
  const char* name, int (*candidates)(const BondCandidatesArgs&)) {
  return {name,
          Gaussian,
          Debye,
          candidates,
          [](const PDFPairSumArgs& a) { return cpuPDFPairSum(a, Gaussian); },
          [](const DebyePairSumArgs& a) { return cpuDebyePairSum(a, Debye); },
          makeCPUWorkspace};
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy

#endif  // CALCULATIONKERNELIMPL_HPP_INCLUDED
