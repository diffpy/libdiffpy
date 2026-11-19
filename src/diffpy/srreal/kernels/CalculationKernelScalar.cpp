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
 * Scalar calculation kernel, the reference and fallback implementation
 *
 *****************************************************************************/

#include <cmath>
#include <numbers>

#include <diffpy/srreal/kernels/PairGeometry.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

namespace {

void gaussianWindowScalar(const GaussianWindowArgs& a) {
  if (a.fwhm <= 0.0) return;
  int first = a.ilo, last = a.ihi;
  if (a.cropped)
    cropGaussianWindow(first, last, a.rlosteps, a.rstep, a.dist, a.fwhm,
                       a.halfboundrel);
  const double xlo = (double(a.rlosteps) + a.ilo) * a.rstep - a.dist;
  for (int i = first; i < last; ++i) {
    const double x = xlo + (i - a.ilo) * a.rstep;
    const double xrel = x / a.fwhm;
    const double y = a.peakheight * exp(-4 * std::numbers::ln2 * xrel * xrel);
    // Normalized by dist instead of r, see PeakProfile::addRDFPeak.
    const double yrdf = y * (x / a.dist + 1);
    a.out[i] += a.scale * yrdf;
  }
}

void debyeSumScalar(const DebyeSumArgs& a) {
  for (int kq = a.kqlo; kq < a.kqhi; ++kq) {
    const double q = kq * a.qstep;
    const double dwq = a.dwsigma * q;
    const double dwscale = exp(-0.5 * dwq * dwq);
    const double sinescale = (a.scale / a.dist) * dwscale *
                             (a.sf0[kq] * a.occ0) * (a.sf1[kq] * a.occ1);
    if (fabs(sinescale) <= a.prec) break;
    a.out[kq] += a.count * sinescale * sin(q * a.dist);
  }
}

int bondCandidatesScalar(const BondCandidatesArgs& a) {
  int m = 0;
  for (int i = 0; i < a.n; ++i) {
    const double dx = a.x[i] + a.ox;
    const double dy = a.y[i] + a.oy;
    const double dz = a.z[i] + a.oz;
    const double d2 = dx * dx + dy * dy + dz * dz;
    a.idx[m] = i;
    // NaN fails both comparisons and stays a candidate.
    m += int(!(d2 < a.d2lo || d2 > a.d2hi));
  }
  return m;
}

}  // namespace

const CalculationKernel& scalarCalculationKernel() {
  static constexpr auto kernel =
    makeCPUKernel<gaussianWindowScalar, debyeSumScalar>("scalar",
                                                        bondCandidatesScalar);
  return kernel;
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
