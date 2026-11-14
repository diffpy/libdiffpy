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

#ifndef PERIODICSTRUCTUREADAPTER_HPP_INCLUDED
#define PERIODICSTRUCTUREADAPTER_HPP_INCLUDED

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/srreal/Lattice.hpp>

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {

class PointsInSphere;

class DLL_EXPORT PeriodicStructureAdapter : public AtomicStructureAdapter {
 public:
  // methods - overloaded
  virtual StructureAdapterPtr clone() const;
  virtual BaseBondGeneratorPtr createBondGenerator() const;
  virtual double numberDensity() const;
  virtual StructureDifference diff(StructureAdapterConstPtr other) const;

  // methods - own
  void setLatPar(double a, double b, double c, double alphadeg, double betadeg,
                 double gammadeg);
  const Lattice& getLattice() const;
  void toCartesian(Atom&) const;
  void toFractional(Atom&) const;

 private:
  // data
  Lattice mlattice;

  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<AtomicStructureAdapter>(*this);
    ar & mlattice;
  }
};

typedef std::shared_ptr<PeriodicStructureAdapter> PeriodicStructureAdapterPtr;

// Comparison functions

DLL_EXPORT bool operator==(const PeriodicStructureAdapter&,
                           const PeriodicStructureAdapter&);
DLL_EXPORT bool operator!=(const PeriodicStructureAdapter&,
                           const PeriodicStructureAdapter&);

class DLL_EXPORT PeriodicStructureBondGenerator : public BaseBondGenerator {
 public:
  // constructors
  PeriodicStructureBondGenerator(StructureAdapterConstPtr);
  // mpositions points into the generator's own arrays
  PeriodicStructureBondGenerator(const PeriodicStructureBondGenerator&) =
    delete;
  PeriodicStructureBondGenerator& operator=(
    const PeriodicStructureBondGenerator&) = delete;

  // methods
  // loop control
  virtual void rewind();

  // configuration
  virtual void selectAnchorSite(int);
  virtual void setRmin(double);
  virtual void setRmax(double);

 protected:
  /// Outer loop of the bond order, the lattice vectors or the positions of
  /// the selected sites.
  enum class BondOrder { LatticeFirst, PositionFirst };

  /// Cartesian points stored by coordinate for the calculation kernels.
  struct PointArray {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    int size() const { return int(x.size()); }
    void clear() {
      x.clear();
      y.clear();
      z.clear();
    }
    void push_back(double px, double py, double pz) {
      x.push_back(px);
      y.push_back(py);
      z.push_back(pz);
    }
    void push_back(const R3::Vector& v) { this->push_back(v[0], v[1], v[2]); }
    /// Append the points [first, last) of a.
    void append(const PointArray& a, int first, int last) {
      x.insert(x.end(), a.x.begin() + first, a.x.begin() + last);
      y.insert(y.end(), a.y.begin() + first, a.y.begin() + last);
      z.insert(z.end(), a.z.begin() + first, a.z.begin() + last);
    }
    /// Largest distance of the points from c.
    double maxDistance(const R3::Vector& c = R3::zerovector) const;
  };

  /// Range of points in a PointArray.
  struct PointSpan {
    const double* x = nullptr;
    const double* y = nullptr;
    const double* z = nullptr;
    int n = 0;

    PointSpan() = default;
    PointSpan(const PointArray& a, int first, int last)
        : x(a.x.data() + first),
          y(a.y.data() + first),
          z(a.z.data() + first),
          n(last - first) {}
    PointSpan(const PointArray& a) : PointSpan(a, 0, a.size()) {}
  };

  /// Lattice vectors [m, n, o] for o in [olo, ohi).  point is the index
  /// of [m, n, olo] in the sequence of all sphere points.
  struct LatticeRow {
    int m;
    int n;
    int olo;
    int ohi;
    int point;
  };

  /// One chunk of the lattice vectors near the anchor.
  struct LatticeChunk {
    PointArray points;
    /// The chunk is the first or the last one.
    bool first = true;
    bool last = true;
    /// The next chunk starts at this anchor row and offset within it.
    int nextrow = 0;
    int nextpoint = 0;
  };

  /// Create generator without site1 positions, add them with addPosition.
  PeriodicStructureBondGenerator(StructureAdapterConstPtr, BondOrder);

  // data
  const PeriodicStructureAdapter* mpstructure;
  /// Index in the site1 positions of the current bond.
  int mcellindex = 0;

  // methods
  virtual void rewindSymmetry();
  virtual void getNextBond();
  /// Add unit cell position of site1.  Call in the order of the sites.
  void addPosition(int site, const R3::Vector& xyz);

 private:
  // data
  BondOrder mbondorder;
  /// Unit cell positions of site1 for all sites, ordered by site.
  PointArray mcell;
  /// Site of each position in mcell.
  std::vector<int> mcellsite;
  /// Positions of site i are mcell points [msitebegin[i], msitebegin[i + 1]).
  std::vector<int> msitebegin = {0};
  /// Rows of lattice vectors within the r-range plus a unit cell buffer.
  std::vector<LatticeRow> msphere;
  /// All msphere lattice vectors when they fit in one chunk, else empty.
  PointArray mspherepoints;
  bool msphere_ready = false;
  /// Center and radius of a sphere around the mcell points.
  R3::Vector mcellcenter;
  double mcellradius = 0.0;
  /// Widening of the r-range for the rounding error of kernel distances.
  double mroundmargin = 0.0;
  /// Positions of the selected sites.
  PointSpan mpositions;
  /// Index of mpositions[0] in mcell, or -1 when they are in mgathered.
  int mpositionsbegin = 0;
  /// Selected positions gathered when the sites are not consecutive.
  PointArray mgathered;
  std::vector<int> mgatheredcell;
  std::vector<int> mgatheredsite;
  /// Parts of the msphere rows that can make a bond with the anchor.
  std::vector<LatticeRow> manchorrows;
  LatticeChunk mchunk;
  /// Squared distance limits of the candidate bonds.
  double mcandd2lo = -1.0;
  double mcandd2hi = -1.0;
  /// Candidate bonds as indices of the major (outer loop) and minor points.
  std::vector<int> mcandmajor;
  std::vector<int> mcandminor;
  int mcandcount = 0;
  int mcandidx = 0;
  int mnextmajor = 0;
  /// True when the lattice chunk was searched for the major point
  /// mnextmajor, used when positions are major and there are more chunks.
  bool mchunksearched = false;
  /// Kernel output and candidate bits for one block of major points.
  std::vector<int> mblockidx;
  std::vector<std::uint64_t> mblockbits;

  // methods
  void buildSphere();
  void releaseSphere();
  void selectPositions();
  void selectAnchorRows();
  void addAnchorRow(const LatticeRow& rw, int olo, int ohi);
  void loadLatticeChunk(bool first);
  bool refillCandidates();
  void growCandidates(size_t capacity);
  int searchMajorPoint(int imajor, const PointSpan& major,
                       const PointSpan& minor, int* idx) const;
  void searchMajorPoints(const PointSpan& major, const PointSpan& minor);
  void searchAlongMajor(const PointSpan& major, const PointSpan& minor);
  void seekCandidate();
};

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PeriodicStructureAdapter)

#endif  // PERIODICSTRUCTUREADAPTER_HPP_INCLUDED
