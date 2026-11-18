/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 *****************************************************************************/

#ifndef DIFFPY_KERNELS_PAIRGEOMETRY_HPP
#define DIFFPY_KERNELS_PAIRGEOMETRY_HPP

#include <float.h>
#include <math.h>

#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

constexpr double PAIR_SIGMA_TO_FWHM = 2.3548200450309493;

/// Use the per-bond rounding order for the Debye precision cutoff.
inline double fwhmToSigma(double fwhm) {
  constexpr double scale = 1.0 / PAIR_SIGMA_TO_FWHM;
  return fwhm * scale;
}

inline long long pairCount(const PairGeometryArgs& g) {
  return g.translations ? g.cntpairs
                        : (long long)g.cntsites * (g.cntsites - 1) / 2;
}

/// Clip before converting to an integer, including for very wide peaks.
inline int gridIndex(double r, double step, int lo, int offset, int npts) {
  const double i = floor(r / step) - lo + offset;
  return i > 0 ? (i < npts ? int(i) : npts) : 0;
}

struct Pair {
  bool ok;
  int site0;
  int site1;
  double dist;
  double msd;
  unsigned weight;
};

/// As meanSquareDisplacement for the tensor xx, yy, zz, xy, xz, yz.
inline double siteMSD(const double* U, double rx, double ry, double rz,
                      int anisotropy) {
  if (!anisotropy) return U[0];
  const double sUs = U[0] * rx * rx + U[1] * ry * ry + U[2] * rz * rz +
                     2 * (U[3] * rx * ry + U[4] * rx * rz + U[5] * ry * rz);
  return sUs / (rx * rx + ry * ry + rz * rz);
}

/// Coordinates of a pair in the strict lower triangle of a cluster.
struct ClusterPair {
  int site0, site1;

  static ClusterPair at(long long k) {
    // k = i0 * (i0 - 1) / 2 + i1 with i1 < i0.
    long long row = (long long)((1.0 + sqrt(1.0 + 8.0 * (double)k)) / 2);
    while (row * (row - 1) / 2 > k) --row;
    while ((row + 1) * row / 2 <= k) ++row;
    return {int(row), int(k - row * (row - 1) / 2)};
  }

  void next() {
    if (++site1 == site0) {
      ++site0;
      site1 = 0;
    }
  }
};

inline Pair pairAt(const PairGeometryArgs& g, long long k,
                   ClusterPair cluster) {
  Pair rv{};
  int p1 = cluster.site1;
  rv.site0 = cluster.site0;
  double tx = 0.0, ty = 0.0, tz = 0.0;
  if (g.translations) {
    rv.site0 = g.pairsite0[k];
    p1 = g.pairpos[k];
    const double* t = g.trans + 3LL * g.pairtrans[k];
    tx = t[0];
    ty = t[1];
    tz = t[2];
  }
  rv.site1 = g.possite[p1];
  const double* r0 = g.pos + 3LL * g.sitebegin[rv.site0];
  const double* r1 = g.pos + 3LL * p1;
  const double rx = r1[0] + tx - r0[0];
  const double ry = r1[1] + ty - r0[1];
  const double rz = r1[2] + tz - r0[2];
  rv.dist = sqrt(rx * rx + ry * ry + rz * rz);
  rv.ok = !(rv.dist < g.rmin || rv.dist > g.rmax || rv.dist <= g.zeroeps);
  if (!rv.ok) return rv;
  rv.msd =
    siteMSD(g.siteuij + 6LL * rv.site0, rx, ry, rz, g.anisotropy[rv.site0]) +
    siteMSD(g.posuij + 6LL * p1, rx, ry, rz, g.anisotropy[rv.site1]);
  const unsigned m = g.multiplicity[rv.site0];
  rv.weight = (g.translations && rv.site0 == rv.site1) ? m : 2 * m;
  return rv;
}

/// Look up a pair directly from its index.
inline Pair pairAt(const PairGeometryArgs& g, long long k) {
  return pairAt(g, k, g.translations ? ClusterPair{} : ClusterPair::at(k));
}

inline double debyeWallerWidth(double msd) {
  return msd < 0.0 ? 0.0 : PAIR_SIGMA_TO_FWHM * sqrt(msd);
}

inline double sharpeningRatio(double r, double delta1, double delta2,
                              double qbroad) {
  if (r < DBL_EPSILON) return 0.0;
  const double qr = qbroad * r;
  return 1.0 - delta1 / r - delta2 / (r * r) + qr * qr;
}

inline double sharpenedWidth(double r, double corr, double msd,
                             double qbroadsep) {
  const double qr = qbroadsep * r;
  return sqrt(corr) * debyeWallerWidth(msd) + qr * qr;
}

/// FWHM of the built-in peak width models for one pair.
inline double widthOf(const PairWidthArgs& w, double dist, double msd) {
  if (w.kind == PairWidthArgs::Constant) return w.width;
  if (w.kind == PairWidthArgs::DebyeWaller) return debyeWallerWidth(msd);
  const double corr = sharpeningRatio(dist, w.delta1, w.delta2, w.qbroad);
  return corr <= 0.0 ? 0.0 : sharpenedWidth(dist, corr, msd, w.qbroadsep);
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy

#endif
