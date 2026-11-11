/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2013 Brookhaven Science Associates,
 *                   Brookhaven National Laboratory.
 *                   All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 ******************************************************************************
 *
 * class PeriodicStructureAdapter -- universal adapter for structure with
 *     periodic boundary conditions that has no space group symmetry
 *
 * class PeriodicStructureBondGenerator -- bond generator
 *
 *****************************************************************************/

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class PeriodicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Public Methods ------------------------------------------------------------

StructureAdapterPtr PeriodicStructureAdapter::clone() const {
  StructureAdapterPtr rv(new PeriodicStructureAdapter(*this));
  return rv;
}

BaseBondGeneratorPtr PeriodicStructureAdapter::createBondGenerator() const {
  BaseBondGeneratorPtr bnds(
    new PeriodicStructureBondGenerator(shared_from_this()));
  return bnds;
}

double PeriodicStructureAdapter::numberDensity() const {
  const Lattice& L = this->getLattice();
  double rv = this->totalOccupancy() / L.volume();
  return rv;
}

StructureDifference PeriodicStructureAdapter::diff(
  StructureAdapterConstPtr other) const {
  StructureDifference sd = this->StructureAdapter::diff(other);
  if (sd.stru0 == sd.stru1) return sd;
  typedef std::shared_ptr<const class PeriodicStructureAdapter> PPtr;
  PPtr pother = std::dynamic_pointer_cast<PPtr::element_type>(other);
  if (!pother) return sd;
  assert(pother == sd.stru1);
  if (this->getLattice() != pother->getLattice()) return sd;
  sd = this->AtomicStructureAdapter::diff(other);
  return sd;
}

void PeriodicStructureAdapter::setLatPar(double a, double b, double c,
                                         double alphadeg, double betadeg,
                                         double gammadeg) {
  mlattice.setLatPar(a, b, c, alphadeg, betadeg, gammadeg);
}

const Lattice& PeriodicStructureAdapter::getLattice() const { return mlattice; }

void PeriodicStructureAdapter::toCartesian(Atom& a) const {
  const Lattice& L = this->getLattice();
  a.xyz_cartn = L.cartesian(a.xyz_cartn);
  a.uij_cartn = L.cartesianMatrix(a.uij_cartn);
}

void PeriodicStructureAdapter::toFractional(Atom& a) const {
  const Lattice& L = this->getLattice();
  a.xyz_cartn = L.fractional(a.xyz_cartn);
  a.uij_cartn = L.fractionalMatrix(a.uij_cartn);
}

// Comparison functions ------------------------------------------------------

bool operator==(const PeriodicStructureAdapter& stru0,
                const PeriodicStructureAdapter& stru1) {
  const AtomicStructureAdapter& astru0 = stru0;
  const AtomicStructureAdapter& astru1 = stru1;
  bool rv = (astru0 == astru1) && (stru0.getLattice() == stru1.getLattice());
  return rv;
}

bool operator!=(const PeriodicStructureAdapter& stru0,
                const PeriodicStructureAdapter& stru1) {
  return !(stru0 == stru1);
}

//////////////////////////////////////////////////////////////////////////////
// class PeriodicStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Local Constants -----------------------------------------------------------

namespace {

/// Fill candidate buffers with at least this many bonds when possible.
const int CANDIDATE_BATCH = 256;

/// Candidate margin relative to the coordinate scale, well above round-off.
const double CANDIDATE_MARGIN = 1e-9;

/// Major points per kernel call in searchAlongMajor.
const int MAJOR_BLOCK = 512;

/// The kernel runs along the major points when there are fewer minor points
/// than this.  Their candidate flags must fit in mblockbits.
const int MINOR_FEW = 64;

/// Lattice vectors in one chunk of mlattice, which takes 0.8 MB.
const int LATTICE_CHUNK = 1 << 15;

/// Squared limits for BondCandidatesArgs, negative means no bond or limit.
double lowerLimitSquared(double rlo) {
  return (rlo > 0.0) ? (rlo * rlo) : -1.0;
}
double upperLimitSquared(double rhi) {
  return (rhi < 0.0) ? -1.0 : (rhi * rhi);
}

/// Larger of d and dmax, where NaN wins.
double maxOrNaN(double d, double dmax) {
  return (std::isnan(d) || d > dmax) ? d : dmax;
}

/// Narrow the integer range [lo, hi) to [dlo, dhi).  NaN keeps the bound.
void narrowRange(int& lo, int& hi, double dlo, double dhi) {
  if (dlo > lo) lo = (dlo < hi) ? int(dlo) : hi;
  if (dhi < hi) hi = (dhi > lo) ? int(dhi) : lo;
}

}  // namespace

// Constructor ---------------------------------------------------------------

PeriodicStructureBondGenerator::PeriodicStructureBondGenerator(
  StructureAdapterConstPtr adpt)
    : PeriodicStructureBondGenerator(adpt, BondOrder::LatticeFirst) {
  const Lattice& L = mpstructure->getLattice();
  for (int i = 0; i < mpstructure->countSites(); ++i) {
    this->addPosition(i, L.ucvCartesian((*mpstructure)[i].xyz_cartn));
  }
}

PeriodicStructureBondGenerator::PeriodicStructureBondGenerator(
  StructureAdapterConstPtr adpt, BondOrder order)
    : BaseBondGenerator(adpt),
      mpstructure(dynamic_cast<const PeriodicStructureAdapter*>(adpt.get())),
      mbondorder(order) {
  assert(mpstructure);
}

// Public Methods ------------------------------------------------------------

void PeriodicStructureBondGenerator::rewind() {
  // Delay msphere instantiation to here instead of in constructor,
  // so it is possible to use setRmin, setRmax.
  if (!msphere_ready) this->buildSphere();
  // BaseBondGenerator::rewind calls this->rewindSymmetry
  this->BaseBondGenerator::rewind();
}

void PeriodicStructureBondGenerator::selectAnchorSite(int anchor) {
  this->BaseBondGenerator::selectAnchorSite(anchor);
  const int i = msitebegin[anchor];
  mr0 = R3::Vector(mcell.x[i], mcell.y[i], mcell.z[i]);
}

void PeriodicStructureBondGenerator::setRmin(double rmin) {
  // destroy msphere so it will be created on rewind with new rmin
  if (this->getRmin() != rmin) this->releaseSphere();
  this->BaseBondGenerator::setRmin(rmin);
}

void PeriodicStructureBondGenerator::setRmax(double rmax) {
  // destroy msphere so it will be created on rewind with new rmax
  if (this->getRmax() != rmax) this->releaseSphere();
  this->BaseBondGenerator::setRmax(rmax);
}

// Protected Methods ---------------------------------------------------------

// site1 is a lattice vector plus a position of a selected site.  Bonds are
// ordered by major then minor points as in the original nested loops, where
// major are the lattice vectors here and the positions in
// CrystalStructureBondGenerator.  The kernel selects candidates and
// BaseBondGenerator checks their exact distance.
//
// Lattice vectors near the anchor come in chunks of LATTICE_CHUNK points.
// When positions are major, each position goes through all chunks.

void PeriodicStructureBondGenerator::rewindSymmetry() {
  this->selectPositions();
  this->selectAnchorRows();
  this->loadLatticeChunk(true);
  mnextmajor = 0;
  mchunksearched = false;
  mcandcount = 0;
  mcandidx = 0;
  this->seekCandidate();
}

void PeriodicStructureBondGenerator::getNextBond() {
  ++mcandidx;
  this->seekCandidate();
}

void PeriodicStructureBondGenerator::addPosition(int site,
                                                 const R3::Vector& xyz) {
  assert(site + 1 == int(msitebegin.size()) ||
         site + 2 == int(msitebegin.size()));
  msitebegin.resize(site + 2, mcell.size());
  mcell.push_back(xyz);
  mcellsite.push_back(site);
  msitebegin.back() = mcell.size();
  msphere_ready = false;
}

double PeriodicStructureBondGenerator::PointArray::maxDistance(
  const R3::Vector& c) const {
  double d2max = 0.0;
  for (int i = 0; i < this->size(); ++i) {
    const double dx = x[i] - c[0], dy = y[i] - c[1], dz = z[i] - c[2];
    d2max = maxOrNaN(dx * dx + dy * dy + dz * dz, d2max);
  }
  return sqrt(d2max);
}

// Private Methods -----------------------------------------------------------

void PeriodicStructureBondGenerator::buildSphere() {
  const Lattice& L = mpstructure->getLattice();
  const double buffzone = L.ucMaxDiagonalLength();
  PointsInSphere sph(this->getRmin() - buffzone, this->getRmax() + buffzone, L);
  // Join consecutive points along o into rows and store the points of
  // small spheres.  L.cartesian(sph.mno()) gives the bonds exactly the
  // coordinates of the original PointsInSphere loops.
  msphere.clear();
  mspherepoints.clear();
  int npoints = 0;
  for (sph.rewind(); !sph.finished(); sph.next(), ++npoints) {
    if (npoints < LATTICE_CHUNK) {
      mspherepoints.push_back(L.cartesian(sph.mno()));
    }
    const int m = sph.m(), n = sph.n(), o = sph.o();
    if (!msphere.empty()) {
      LatticeRow& last = msphere.back();
      if (last.m == m && last.n == n && last.ohi == o) {
        ++last.ohi;
        continue;
      }
    }
    msphere.push_back({m, n, o, o + 1, npoints});
  }
  if (npoints > LATTICE_CHUNK) {
    mspherepoints = {};
  }
  // The norm along a row is largest at one of its ends.
  double maxlatnorm = 0.0;
  for (const LatticeRow& rw : msphere) {
    for (int o : {rw.olo, rw.ohi - 1}) {
      const double r = R3::norm(L.cartesian(R3::Vector(rw.m, rw.n, o)));
      maxlatnorm = maxOrNaN(r, maxlatnorm);
    }
  }
  // sphere around the unit cell positions
  mcellcenter = R3::zerovector;
  const int n = mcell.size();
  for (int i = 0; i < n; ++i) {
    mcellcenter += R3::Vector(mcell.x[i], mcell.y[i], mcell.z[i]);
  }
  if (n) mcellcenter /= n;
  mcellradius = mcell.maxDistance(mcellcenter);
  // Widen the r-range by the rounding error of the kernel distances.
  const double maxnorm = mcell.maxDistance() + maxlatnorm;
  const double scale =
    2 * maxnorm + fabs(this->getRmin()) + fabs(this->getRmax()) + 1.0;
  mroundmargin = CANDIDATE_MARGIN * scale;
  mcandd2lo = lowerLimitSquared(this->getRmin() - mroundmargin);
  mcandd2hi = upperLimitSquared(this->getRmax() + mroundmargin);
  msphere_ready = true;
}

void PeriodicStructureBondGenerator::releaseSphere() {
  msphere = {};
  mspherepoints = {};
  msphere_ready = false;
}

void PeriodicStructureBondGenerator::selectPositions() {
  const int nsel = msite_last - msite_first;
  assert(nsel > 0);
  const int site0 = *msite_first;
  // consecutive sites have their positions together in mcell
  if (msite_consecutive) {
    mpositionsbegin = msitebegin[site0];
    mpositions = PointSpan(mcell, mpositionsbegin, msitebegin[site0 + nsel]);
    return;
  }
  mgathered.clear();
  mgatheredcell.clear();
  mgatheredsite.clear();
  for (int k = 0; k < nsel; ++k) {
    const int site = msite_first[k];
    for (int i = msitebegin[site]; i < msitebegin[site + 1]; ++i) {
      mgathered.push_back(mcell.x[i], mcell.y[i], mcell.z[i]);
      mgatheredcell.push_back(i);
      mgatheredsite.push_back(k);
    }
  }
  mpositionsbegin = -1;
  mpositions = PointSpan(mgathered);
}

void PeriodicStructureBondGenerator::selectAnchorRows() {
  manchorrows.clear();
  // The sphere around the unit cell positions shifted by the lattice vector
  // must reach the r-range of the anchor.
  const double margin = mcellradius + mroundmargin;
  const double rlo = this->getRmin() - margin;
  const double rhi = this->getRmax() + margin;
  if (rhi < 0.0) return;
  // Row points are v + o * c relative to the anchor, so their squared
  // distance is dmin2 + cc * (o - os)^2.  The o-limits get one more point
  // on each side to cover their rounding.  NaN keeps the whole row.
  const Lattice& L = mpstructure->getLattice();
  const R3::Vector c = L.cartesian(R3::Vector(0, 0, 1));
  const double cc = R3::dot(c, c);
  const R3::Vector center = mcellcenter - mr0;
  for (const LatticeRow& rw : msphere) {
    const R3::Vector v = L.cartesian(R3::Vector(rw.m, rw.n, 0)) + center;
    const double vc = R3::dot(v, c);
    const double os = -vc / cc;
    double dmin2 = R3::dot(v, v) - vc * vc / cc;
    if (dmin2 < 0.0) dmin2 = 0.0;
    if (dmin2 > rhi * rhi) continue;
    // points within rhi
    const double w = sqrt((rhi * rhi - dmin2) / cc);
    int olo = rw.olo;
    int ohi = rw.ohi;
    narrowRange(olo, ohi, ceil(os - w) - 1, floor(os + w) + 2);
    if (olo >= ohi) continue;
    // points well within rlo, they split the row
    if (rlo > 0.0 && dmin2 < rlo * rlo) {
      const double u = sqrt((rlo * rlo - dmin2) / cc);
      int holelo = olo;
      int holehi = ohi;
      narrowRange(holelo, holehi, floor(os - u + 1) + 1, ceil(os + u - 1));
      if (holelo < holehi) {
        if (olo < holelo) addAnchorRow(rw, olo, holelo);
        if (holehi < ohi) addAnchorRow(rw, holehi, ohi);
        continue;
      }
    }
    addAnchorRow(rw, olo, ohi);
  }
}

void PeriodicStructureBondGenerator::addAnchorRow(const LatticeRow& rw, int olo,
                                                  int ohi) {
  manchorrows.push_back({rw.m, rw.n, olo, ohi, rw.point + (olo - rw.olo)});
}

void PeriodicStructureBondGenerator::loadLatticeChunk(bool first) {
  LatticeChunk& ch = mchunk;
  ch.points.clear();
  ch.first = first;
  if (first) {
    ch.nextrow = 0;
    ch.nextpoint = 0;
  }
  const int nrows = manchorrows.size();
  if (!mspherepoints.x.empty()) {
    // copy the stored sphere points, which are all in one chunk
    for (; ch.nextrow < nrows; ++ch.nextrow) {
      const LatticeRow& rw = manchorrows[ch.nextrow];
      ch.points.append(mspherepoints, rw.point, rw.point + (rw.ohi - rw.olo));
    }
  } else {
    // The retained rows already contain the exact integer coordinates.
    const Lattice& L = mpstructure->getLattice();
    while (ch.nextrow < nrows && ch.points.size() < LATTICE_CHUNK) {
      const LatticeRow& rw = manchorrows[ch.nextrow];
      const int count =
        min(rw.ohi - rw.olo - ch.nextpoint, LATTICE_CHUNK - ch.points.size());
      for (int i = 0; i < count; ++i) {
        ch.points.push_back(
          L.cartesian(R3::Vector(rw.m, rw.n, rw.olo + ch.nextpoint + i)));
      }
      ch.nextpoint += count;
      if (ch.nextpoint == rw.ohi - rw.olo) {
        ++ch.nextrow;
        ch.nextpoint = 0;
      }
    }
  }
  ch.last = (ch.nextrow >= nrows);
}

bool PeriodicStructureBondGenerator::refillCandidates() {
  mcandcount = 0;
  mcandidx = 0;
  const PointSpan lattice(mchunk.points);
  if (mbondorder == BondOrder::LatticeFirst) {
    if (mnextmajor >= lattice.n) {
      if (mchunk.last) return false;
      this->loadLatticeChunk(false);
      mnextmajor = 0;
      return true;
    }
    this->searchMajorPoints(lattice, mpositions);
    return true;
  }
  if (mnextmajor >= mpositions.n) return false;
  // All positions see the same lattice vectors when they fit in one chunk.
  if (mchunk.first && mchunk.last) {
    this->searchMajorPoints(mpositions, lattice);
    return true;
  }
  // Otherwise every position goes through all chunks.
  if (mchunksearched) {
    if (!mchunk.last) {
      this->loadLatticeChunk(false);
    } else {
      if (++mnextmajor >= mpositions.n) return false;
      this->loadLatticeChunk(true);
    }
    mchunksearched = false;
    return true;
  }
  this->growCandidates(lattice.n);
  mcandcount =
    this->searchMajorPoint(mnextmajor, mpositions, lattice, mcandminor.data());
  fill_n(mcandmajor.begin(), mcandcount, mnextmajor);
  mchunksearched = true;
  return true;
}

void PeriodicStructureBondGenerator::growCandidates(size_t capacity) {
  if (mcandmajor.size() < capacity) {
    mcandmajor.resize(capacity);
    mcandminor.resize(capacity);
  }
}

int PeriodicStructureBondGenerator::searchMajorPoint(int imajor,
                                                     const PointSpan& major,
                                                     const PointSpan& minor,
                                                     int* idx) const {
  return kernels::activeCalculationKernel().bondCandidates({
    .x = minor.x,
    .y = minor.y,
    .z = minor.z,
    .n = minor.n,
    .ox = major.x[imajor] - mr0[0],
    .oy = major.y[imajor] - mr0[1],
    .oz = major.z[imajor] - mr0[2],
    .d2lo = mcandd2lo,
    .d2hi = mcandd2hi,
    .idx = idx,
  });
}

void PeriodicStructureBondGenerator::searchMajorPoints(const PointSpan& major,
                                                       const PointSpan& minor) {
  if (minor.n == 0) {
    mnextmajor = major.n;
    return;
  }
  // The kernel is faster along the major points when the minor ones are few.
  if (minor.n < MINOR_FEW && minor.n < major.n) {
    this->searchAlongMajor(major, minor);
    return;
  }
  this->growCandidates(CANDIDATE_BATCH + size_t(minor.n));
  for (; mnextmajor < major.n && mcandcount < CANDIDATE_BATCH; ++mnextmajor) {
    int* idx = mcandminor.data() + mcandcount;
    const int m = this->searchMajorPoint(mnextmajor, major, minor, idx);
    fill_n(mcandmajor.begin() + mcandcount, m, mnextmajor);
    mcandcount += m;
  }
}

void PeriodicStructureBondGenerator::searchAlongMajor(const PointSpan& major,
                                                      const PointSpan& minor) {
  this->growCandidates(CANDIDATE_BATCH + size_t(minor.n) * MAJOR_BLOCK);
  // Bit j of mblockbits[i] flags the candidate bond of the major point i
  // in the block and the minor point j.
  mblockidx.resize(MAJOR_BLOCK);
  mblockbits.resize(MAJOR_BLOCK);
  const auto bondCandidates = kernels::activeCalculationKernel().bondCandidates;
  kernels::BondCandidatesArgs args{.d2lo = mcandd2lo, .d2hi = mcandd2hi};
  args.idx = mblockidx.data();
  while (mnextmajor < major.n && mcandcount < CANDIDATE_BATCH) {
    const int nblock = min(MAJOR_BLOCK, major.n - mnextmajor);
    fill_n(mblockbits.begin(), nblock, 0);
    args.x = major.x + mnextmajor;
    args.y = major.y + mnextmajor;
    args.z = major.z + mnextmajor;
    args.n = nblock;
    for (int j = 0; j < minor.n; ++j) {
      args.ox = minor.x[j] - mr0[0];
      args.oy = minor.y[j] - mr0[1];
      args.oz = minor.z[j] - mr0[2];
      const int m = bondCandidates(args);
      for (int k = 0; k < m; ++k) mblockbits[mblockidx[k]] |= uint64_t(1) << j;
    }
    for (int i = 0; i < nblock; ++i) {
      for (uint64_t bits = mblockbits[i]; bits; bits &= bits - 1) {
        mcandmajor[mcandcount] = mnextmajor + i;
        mcandminor[mcandcount] = countr_zero(bits);
        ++mcandcount;
      }
    }
    mnextmajor += nblock;
  }
}

void PeriodicStructureBondGenerator::seekCandidate() {
  while (mcandidx >= mcandcount) {
    if (!this->refillCandidates()) {
      msite_current = msite_last;
      return;
    }
  }
  const int imajor = mcandmajor[mcandidx];
  const int iminor = mcandminor[mcandidx];
  const bool latticefirst = (mbondorder == BondOrder::LatticeFirst);
  const int ipos = latticefirst ? iminor : imajor;
  const int ilat = latticefirst ? imajor : iminor;
  if (mpositionsbegin >= 0) {
    mcellindex = mpositionsbegin + ipos;
    msite_current = msite_first + (mcellsite[mcellindex] - *msite_first);
  } else {
    mcellindex = mgatheredcell[ipos];
    msite_current = msite_first + mgatheredsite[ipos];
  }
  const PointArray& lattice = mchunk.points;
  mr1[0] = lattice.x[ilat] + mpositions.x[ipos];
  mr1[1] = lattice.y[ilat] + mpositions.y[ipos];
  mr1[2] = lattice.z[ilat] + mpositions.z[ipos];
  this->updateDistance();
}

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PeriodicStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::PeriodicStructureAdapter)

// End of file
