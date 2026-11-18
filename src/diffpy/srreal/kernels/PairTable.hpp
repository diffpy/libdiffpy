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
 * class PairTable -- all pairs of a structure for a sum over them
 * class PairGroups -- pairs that add the same peak
 *
 * Internal header, not installed.
 *
 * PairTable remembers the pairs of a periodic structure and only updates
 * their distances while the structure stays close to that of the search,
 * as in a refinement.  PairGroups merges pairs of the same distance and
 * width, so that each kind of peak is added once.
 *
 *****************************************************************************/

#ifndef KERNELS_PAIRTABLE_HPP_INCLUDED
#define KERNELS_PAIRTABLE_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <diffpy/Export.hpp>
#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/forwardtypes.hpp>
#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

namespace diffpy {
namespace srreal {

class PeriodicStructureAdapter;
class PeakWidthModel;
class PeakProfile;

namespace kernels {

// Exported for the unit tests.
class DLL_EXPORT PairTable {
 public:
  /// True when the table can give the pairs of the structure.  These are
  /// the structures of the classes AtomicStructureAdapter,
  /// PeriodicStructureAdapter and CrystalStructureAdapter, not of classes
  /// derived from them, also within the nometa and nosymmetry wrappers.
  static bool supports(const StructureAdapter& stru);

  /// Number of searches for the pairs of periodic structures so far.
  int countSearches() const { return msearches; }

  /// The geometry for a kernel that adds all pairs of the structure in
  /// the range.  Its arrays are those of this table, they are valid until
  /// the next call.
  PairGeometryArgs geometryArgs(const StructureAdapter& stru, double rmin,
                                double rmax);

 private:
  /// Sources of the sites and lattice translations.
  struct View {
    /// The structure within any nometa and nosymmetry wrappers.
    const StructureAdapter* source;
    /// The outermost nosymmetry wrapper or null.  Its pairs are those of
    /// the sites alone, without symmetry positions and lattice
    /// translations.
    const StructureAdapter* cluster;
    /// The source when it has a lattice, also for a cluster of its sites.
    const PeriodicStructureAdapter* periodic;
    /// True for the pairs of symmetry positions and lattice translations.
    bool hasTranslations() const { return periodic && !cluster; }
    /// The structure whose sites are those of the pairs.
    const StructureAdapter& sites() const {
      return cluster ? *cluster : *source;
    }
  };
  static View view(const StructureAdapter& stru);

  /// Contiguous storage for the arrays read by both CPU and CUDA kernels.
  struct Geometry {
    std::vector<double> pos, posuij, siteuij, trans;
    std::vector<unsigned char> anisotropy;
    std::vector<int> multiplicity;
  };

  /// Positions and displacement tensors of a periodic structure as its
  /// bond generator has them.
  struct Cell {
    /// Positions of site i are [sitebegin[i], sitebegin[i + 1]), the anchor
    /// of the site is the first one.
    std::vector<int> sitebegin;
    std::vector<R3::Vector> pos;
    std::vector<R3::Matrix> uij;
  };

  // methods
  void preparePeriodic(const PeriodicStructureAdapter& stru, double rmin,
                       double rmax);
  static void fillCell(const PeriodicStructureAdapter& stru, Cell& cell);
  bool hasPairsFor(const PeriodicStructureAdapter& stru, const Cell& cell,
                   double rmin, double rmax) const;
  void searchPairs(const PeriodicStructureAdapter& stru, const Cell& cell,
                   double rmin, double rmax);

  // data
  Cell mcell;
  // Pairs of a periodic structure and what they were searched for.
  /// Positions and lattice basis at the search.
  std::vector<int> msitebegin0;
  std::vector<R3::Vector> mpos0;
  R3::Matrix mbase0;
  /// Range of the distances at the search.
  double mrlo0 = 0.0;
  double mrhi0 = -1.0;
  /// Lattice translations of the pairs in lattice coordinates.
  std::vector<std::array<int, 3>> mmno;
  /// Anchor site, position of the other site and lattice translation of
  /// each pair, ordered by the anchor site.
  std::vector<int> mpairsite0;
  std::vector<int> mpairpos;
  std::vector<int> mpairtrans;
  int msearches = 0;
  // work arrays
  std::vector<int> mpossite;
  Geometry mgeometry;
};

/// Sums of the scales of the pairs that have the same distance and width.
///
/// Even adjacent floating-point values can straddle a peak or Debye cutoff.
class PairGroups {
 public:
  struct Group {
    double dist;
    double width;
    /// Tells apart pairs that the distance and the width do not, for
    /// example by their scattering factors.
    std::int64_t tag;
    /// Sum of the scales of the pairs.
    double scale;
    /// Multiplicity scale of a single Debye pair, before its precision cutoff.
    unsigned weight;
  };

  void clear() {
    mgroups.clear();
    mslots.assign(mslots.size(), -1);
  }

  void add(double dist, double width, std::int64_t tag, double scale,
           unsigned weight = 1) {
    const std::uint64_t kd = key(dist);
    const std::uint64_t kw = key(width);
    if (2 * (mgroups.size() + 1) > mslots.size()) this->grow();
    const std::size_t mask = mslots.size() - 1;
    std::size_t i = hash(kd, kw, tag, weight) & mask;
    for (;; i = (i + 1) & mask) {
      const int g = mslots[i];
      if (g < 0) break;
      Group& group = mgroups[g];
      if (key(group.dist) == kd && key(group.width) == kw && group.tag == tag &&
          group.weight == weight) {
        group.scale += scale;
        return;
      }
    }
    mslots[i] = int(mgroups.size());
    mgroups.push_back({dist, width, tag, scale, weight});
  }

  const std::vector<Group>& groups() const { return mgroups; }

 private:
  static std::uint64_t key(double x) { return std::bit_cast<std::uint64_t>(x); }

  static std::size_t hash(std::uint64_t kd, std::uint64_t kw, std::int64_t tag,
                          unsigned weight) {
    std::uint64_t h = kd * 0x9e3779b97f4a7c15ull;
    h ^= (h >> 29) ^ (kw * 0xbf58476d1ce4e5b9ull);
    h ^= (h >> 31) ^ (std::uint64_t(tag) * 0x94d049bb133111ebull);
    return std::size_t((h ^ (h >> 32)) + weight);
  }

  void grow() {
    const std::size_t n = std::max<std::size_t>(1024, 2 * mslots.size());
    mslots.assign(n, -1);
    const std::size_t mask = n - 1;
    for (std::size_t g = 0; g < mgroups.size(); ++g) {
      const Group& group = mgroups[g];
      std::size_t i =
        hash(key(group.dist), key(group.width), group.tag, group.weight) & mask;
      while (mslots[i] >= 0) i = (i + 1) & mask;
      mslots[i] = int(g);
    }
  }

  std::vector<Group> mgroups;
  /// Open addressing table of the group indices, -1 for a free slot.
  std::vector<int> mslots;
};

/// Peak width of the built-in models from the values of a bond.
bool assignPairWidth(PairWidthArgs&, const PeakWidthModel&);

/// Add a peak with the bounds and values supplied by the profile.
void addProfilePeak(const PeakProfile&, double* out, int npts, int rlosteps,
                    double rstep, double dist, double fwhm, double scale);

/// What a calculator keeps for its sums over all pairs.
struct PairSumState {
  PairTable table;
  PairSumWorkspace* workspaceFor(const CalculationKernel& kernel);

 private:
  MakeWorkspace mfactory = nullptr;
  std::unique_ptr<PairSumWorkspace> mworkspace;
};

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy

#endif  // KERNELS_PAIRTABLE_HPP_INCLUDED
