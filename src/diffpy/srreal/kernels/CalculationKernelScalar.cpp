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

#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

namespace {

void gaussianWindowScalar(const GaussianWindowArgs& a) {
  if (a.fwhm <= 0.0) return;
  for (int i = a.ilo; i < a.ihi; ++i) {
    const double x = (a.rlosteps + i) * a.rstep - a.dist;
    const double xrel = x / a.fwhm;
    if (a.cropped && fabs(xrel) >= a.halfboundrel) continue;
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
    const double sinescale = a.scale * dwscale * a.sf0[kq] * a.sf1[kq] / a.dist;
    if (fabs(sinescale) <= a.prec) break;
    a.out[kq] += sinescale * sin(q * a.dist);
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
  static constexpr CalculationKernel kernel = {
    "scalar", gaussianWindowScalar, debyeSumScalar, bondCandidatesScalar};
  return kernel;
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
