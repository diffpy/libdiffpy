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
 *
 *****************************************************************************/

#include <cmath>
#include <map>
#include <typeinfo>
#include <utility>

#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/NoMetaStructureAdapter.hpp>
#include <diffpy/srreal/NoSymmetryStructureAdapter.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/kernels/PairTable.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

namespace {

/// The pairs of a periodic structure are searched for distances up to this
/// much beyond the r-range.  They stay valid while the positions and the
/// lattice change the distances by less.
const double SKIN = 0.25;

bool isPeriodicClass(const StructureAdapter& s) {
  return typeid(s) == typeid(PeriodicStructureAdapter) ||
         typeid(s) == typeid(CrystalStructureAdapter);
}

}  // namespace

// Public Methods ------------------------------------------------------------

PairTable::View PairTable::view(const StructureAdapter& stru) {
  View rv = {&stru, nullptr, nullptr};
  while (true) {
    const StructureAdapter& s = *rv.source;
    if (typeid(s) == typeid(NoMetaStructureAdapter)) {
      const auto& nm = static_cast<const NoMetaStructureAdapter&>(s);
      rv.source = nm.getSourceStructure().get();
    } else if (typeid(s) == typeid(NoSymmetryStructureAdapter)) {
      const auto& ns = static_cast<const NoSymmetryStructureAdapter&>(s);
      if (!rv.cluster) rv.cluster = &s;
      rv.source = ns.getSourceStructure().get();
    } else {
      break;
    }
  }
  if (isPeriodicClass(*rv.source)) {
    rv.periodic = static_cast<const PeriodicStructureAdapter*>(rv.source);
  }
  return rv;
}

bool PairTable::supports(const StructureAdapter& stru) {
  const View v = view(stru);
  // The empty structure has no pairs whatever its class.
  if (stru.countSites() == 0) return true;
  if (typeid(*v.source) == typeid(AtomicStructureAdapter)) return true;
  return v.periodic;
}

// Positions, translations and pairs of the structure for the range.

void PairTable::preparePeriodic(const PeriodicStructureAdapter& stru,
                                double rmin, double rmax) {
  Cell& cell = mcell;
  fillCell(stru, cell);
  if (!this->hasPairsFor(stru, cell, rmin, rmax)) {
    this->searchPairs(stru, cell, rmin, rmax);
  }
  // site of each position
  const int cntsites = stru.countSites();
  mpossite.resize(cell.pos.size());
  for (int i = 0; i < cntsites; ++i) {
    for (int p = cell.sitebegin[i]; p < cell.sitebegin[i + 1]; ++p) {
      mpossite[p] = i;
    }
  }
}

PairGeometryArgs PairTable::geometryArgs(const StructureAdapter& stru,
                                         double rmin, double rmax) {
  PairGeometryArgs rv = {};
  Geometry& g = mgeometry;
  const View v = view(stru);
  const StructureAdapter& sites = v.sites();
  const int cntsites = sites.countSites();
  auto tensor = [](std::vector<double>& out, const R3::Matrix& U) {
    for (double u : {U(0, 0), U(1, 1), U(2, 2), U(0, 1), U(0, 2), U(1, 2)}) {
      out.push_back(u);
    }
  };
  rv.translations = v.hasTranslations();
  g.pos.clear();
  g.posuij.clear();
  g.trans.clear();
  if (rv.translations) {
    this->preparePeriodic(*v.periodic, rmin, rmax);
    for (std::size_t p = 0; p < mcell.pos.size(); ++p) {
      g.pos.insert(g.pos.end(), mcell.pos[p].begin(), mcell.pos[p].end());
      tensor(g.posuij, mcell.uij[p]);
    }
    for (const auto& mno : mmno) {
      const auto t =
        v.periodic->getLattice().cartesian(R3::Vector(mno[0], mno[1], mno[2]));
      g.trans.insert(g.trans.end(), t.begin(), t.end());
    }
    rv.cntpairs = mpairsite0.size();
    rv.pairsite0 = mpairsite0.data();
    rv.pairpos = mpairpos.data();
    rv.pairtrans = mpairtrans.data();
    rv.pairversion = msearches;
  } else {
    mcell.sitebegin.resize(cntsites + 1);
    mpossite.resize(cntsites);
    for (int i = 0; i < cntsites; ++i) {
      mcell.sitebegin[i] = mpossite[i] = i;
      const R3::Vector r = sites.siteCartesianPosition(i);
      g.pos.insert(g.pos.end(), r.begin(), r.end());
    }
    mcell.sitebegin[cntsites] = cntsites;
  }
  g.siteuij.clear();
  g.anisotropy.resize(cntsites);
  g.multiplicity.resize(cntsites);
  for (int i = 0; i < cntsites; ++i) {
    tensor(g.siteuij, sites.siteCartesianUij(i));
    g.anisotropy[i] = sites.siteAnisotropy(i);
    g.multiplicity[i] = sites.siteMultiplicity(i);
  }
  rv.cntsites = cntsites;
  rv.cntpositions = int(mpossite.size());
  rv.cnttrans = int(g.trans.size() / 3);
  rv.sitebegin = mcell.sitebegin.data();
  rv.possite = mpossite.data();
  rv.pos = g.pos.data();
  rv.posuij = rv.translations ? g.posuij.data() : g.siteuij.data();
  rv.siteuij = g.siteuij.data();
  rv.anisotropy = g.anisotropy.data();
  rv.multiplicity = g.multiplicity.data();
  rv.trans = g.trans.data();
  rv.rmin = rmin;
  rv.rmax = rmax;
  // as eps_eq(distance, 0.0)
  rv.zeroeps = mathutils::SQRT_DOUBLE_EPS;
  return rv;
}

// The positions of PeriodicStructureBondGenerator and of
// CrystalStructureBondGenerator.

void PairTable::fillCell(const PeriodicStructureAdapter& stru, Cell& cell) {
  const int cntsites = stru.countSites();
  cell.sitebegin.assign(1, 0);
  cell.pos.clear();
  cell.uij.clear();
  if (typeid(stru) == typeid(CrystalStructureAdapter)) {
    const auto& cstru = static_cast<const CrystalStructureAdapter&>(stru);
    for (int i = 0; i < cntsites; ++i) {
      for (const Atom& a : cstru.getEquivalentAtoms(i)) {
        cell.pos.push_back(a.xyz_cartn);
        cell.uij.push_back(a.uij_cartn);
      }
      cell.sitebegin.push_back(int(cell.pos.size()));
    }
    return;
  }
  const Lattice& L = stru.getLattice();
  for (int i = 0; i < cntsites; ++i) {
    const Atom& a = stru[i];
    cell.pos.push_back(L.ucvCartesian(a.xyz_cartn));
    cell.uij.push_back(a.uij_cartn);
    cell.sitebegin.push_back(int(cell.pos.size()));
  }
}

// Whether the remembered pairs include all pairs of the cell in the range.
// Keep the same cell images as a fresh bond generator at cutoff boundaries.

bool PairTable::hasPairsFor(const PeriodicStructureAdapter& stru,
                            const Cell& cell, double rmin, double rmax) const {
  if (cell.sitebegin != msitebegin0) return false;
  const Lattice& L = stru.getLattice();
  // largest change of a position and of a lattice translation
  double dmax = 0.0;
  for (std::size_t i = 0; i < mpos0.size(); ++i) {
    const double d = R3::distance(cell.pos[i], mpos0[i]);
    // NaN makes the pairs invalid
    if (!(d <= SKIN)) return false;
    if (d > dmax) dmax = d;
  }
  double tmax = 0.0;
  if (L.base() != mbase0) {
    // An in-range bond has translation t = r01 - (p1 - p0), even if the
    // previous search missed it. Twice the maximum distance from the
    // first position bounds |p1 - p0|.
    double radius = 0.0;
    for (const auto& p : cell.pos) {
      radius = std::max(radius, R3::distance(p, cell.pos.front()));
    }
    // t * strain is the change back to its old translation. The Frobenius
    // norm bounds this change for every direction, including in skew cells.
    const R3::Matrix strain = R3::prod(L.recbase(), mbase0 - L.base());
    double norm2 = 0.0;
    for (double x : strain.data()) norm2 += x * x;
    tmax = (rmax + 2 * radius) * std::sqrt(norm2);
  }
  // No pair that was out of the searched range can be in the range now.
  const double change = 2 * dmax + tmax;
  const bool lowok = (mrlo0 <= 0.0) || (rmin - change >= mrlo0);
  const bool highok = (rmax + change <= mrhi0);
  return lowok && highok;
}

// Search the pairs with the bond generator of the structure.

void PairTable::searchPairs(const PeriodicStructureAdapter& stru,
                            const Cell& cell, double rmin, double rmax) {
  // An empty site-offset array never matches a cell (even an empty cell has
  // the end offset). A failed search must be retried before reusing its lists.
  msitebegin0.clear();
  const Lattice& L = stru.getLattice();
  const int cntsites = stru.countSites();
  mrlo0 = rmin - SKIN;
  mrhi0 = rmax + SKIN;
  mpos0 = cell.pos;
  mbase0 = L.base();
  mmno.clear();
  mpairsite0.clear();
  mpairpos.clear();
  mpairtrans.clear();
  // index of each translation by its lattice coordinates
  std::map<std::array<int, 3>, int> transindex;
  auto indexOf = [&](const std::array<int, 3>& mno) {
    const auto [it, added] = transindex.try_emplace(mno, int(mmno.size()));
    if (added) mmno.push_back(mno);
    return it->second;
  };
  auto latticeCoordinates = [&](const R3::Vector& t) {
    const R3::Vector f = L.fractional(t);
    return std::array<int, 3>{int(std::lround(f[0])), int(std::lround(f[1])),
                              int(std::lround(f[2]))};
  };
  // Keep distinct coincident positions in the candidate list. pairAt tests
  // their current distance each time the table is used.
  auto owner = stru.createBondGenerator();
  auto& bnds = static_cast<PeriodicStructureBondGenerator&>(*owner);
  bnds.mkeepcoincident = true;
  bnds.setRmin(mrlo0);
  bnds.setRmax(mrhi0);
  std::vector<std::pair<int, int>> coincident;
  for (int i0 = 0; i0 < cntsites; ++i0) {
    auto append = [&](int position, int translation) {
      mpairsite0.push_back(i0);
      mpairpos.push_back(position);
      mpairtrans.push_back(translation);
    };
    bnds.selectAnchorSite(i0);
    bnds.selectSiteRange(0, i0 + 1);
    for (bnds.rewind(); !bnds.finished(); bnds.next()) {
      const int p = bnds.mcellindex;
      const auto mno = latticeCoordinates(bnds.r1() - cell.pos[p]);
      if (p == cell.sitebegin[i0] && mno == std::array<int, 3>{}) continue;
      const int t = indexOf(mno);
      if (mathutils::eps_eq(bnds.distance(), 0.0))
        coincident.emplace_back(p, t);
      else
        append(p, t);
    }
    // Reserve the trailing slots of each anchor for pairs that may separate.
    // Interleaving them shifts the regular bonds between CUDA reduction rows.
    for (auto [p, t] : coincident) append(p, t);
    coincident.clear();
  }
  msitebegin0 = cell.sitebegin;
  ++msearches;
}

//////////////////////////////////////////////////////////////////////////////
// Built-in peak width models
//////////////////////////////////////////////////////////////////////////////

bool assignPairWidth(PairWidthArgs& rv, const PeakWidthModel& model) {
  rv = {};
  if (typeid(model) == typeid(JeongPeakWidth)) {
    const auto& jeong = static_cast<const JeongPeakWidth&>(model);
    rv = {PairWidthArgs::Jeong, jeong.getDelta1(),           jeong.getDelta2(),
          jeong.getQbroad(),    jeong.getQbroad_seperable(), 0.0};
  } else if (typeid(model) == typeid(DebyeWallerPeakWidth)) {
    rv.kind = PairWidthArgs::DebyeWaller;
  } else if (typeid(model) == typeid(ConstantPeakWidth)) {
    rv.kind = PairWidthArgs::Constant;
    rv.width = static_cast<const ConstantPeakWidth&>(model).getWidth();
  } else {
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// struct PairSumState
//////////////////////////////////////////////////////////////////////////////

PairSumWorkspace* PairSumState::workspaceFor(const CalculationKernel& kernel) {
  if (mfactory != kernel.makeWorkspace) {
    mworkspace.reset(kernel.makeWorkspace());
    mfactory = kernel.makeWorkspace;
  }
  return mworkspace.get();
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy

// End of file
