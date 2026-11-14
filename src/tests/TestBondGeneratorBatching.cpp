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
 * Compare the periodic and crystal bond generators with a reference loop.
 *
 *****************************************************************************/

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>

#include "test_helpers.hpp"

using namespace diffpy::srreal;

namespace {

struct Bond {
  int site0;
  int site1;
  R3::Vector r1;
  R3::Vector r01;
  double distance;
  R3::Matrix u1;
};

struct Config {
  double rmin;
  double rmax;
};

constexpr std::array<Config, 4> CONFIGS = {
  {{0.0, 9.7}, {5.0, 4.0}, {0.0, 0.0}, {6.1, 6.9}}};

/// Hexagonal, orthorhombic and triclinic cells.
const char* const STRUCTURES[] = {"ZnS_wurtzite.stru", "CaTiO3.stru",
                                  "alpha_K2Bi8Se13.stru"};

/// Site selections exercised for every anchor.
std::vector<SiteIndices> siteSelections(int cntsites, int anchor) {
  std::vector<SiteIndices> rv(4);
  for (int i = 0; i <= anchor; ++i) rv[0].push_back(i);
  for (int i = cntsites - 1; i >= 0; i -= 2) rv[1].push_back(i);
  for (int i = 1; i < cntsites - 1; ++i) rv[2].push_back(i);
  // rv[3] stays empty
  return rv;
}

bool inRange(const Bond& b, const Config& cfg) {
  using diffpy::mathutils::eps_eq;
  const bool outofrange = (b.distance < cfg.rmin) || (b.distance > cfg.rmax);
  return !outofrange && !eps_eq(b.distance, 0.0);
}

/// Reference: msphere points outer loop, selected sites inner loop.
std::vector<Bond> referencePeriodic(const PeriodicStructureAdapter& stru,
                                    int anchor, const SiteIndices& sel,
                                    const Config& cfg) {
  const Lattice& L = stru.getLattice();
  const double buffzone = L.ucMaxDiagonalLength();
  PointsInSphere sph(cfg.rmin - buffzone, cfg.rmax + buffzone, L);
  const R3::Vector r0 = L.ucvCartesian(stru[anchor].xyz_cartn);
  std::vector<Bond> rv;
  for (sph.rewind(); !sph.finished(); sph.next()) {
    const R3::Vector rc = L.cartesian(sph.mno());
    for (int j : sel) {
      Bond b;
      b.site0 = anchor;
      b.site1 = j;
      const R3::Vector xyz = L.ucvCartesian(stru[j].xyz_cartn);
      for (int k = 0; k < 3; ++k) b.r1[k] = rc[k] + xyz[k];
      for (int k = 0; k < 3; ++k) b.r01[k] = b.r1[k] - r0[k];
      b.distance = R3::norm(b.r01);
      b.u1 = stru.siteCartesianUij(j);
      if (inRange(b, cfg)) rv.push_back(b);
    }
  }
  return rv;
}

/// Reference: sites, then symmetry positions, then msphere points.
std::vector<Bond> referenceCrystal(const CrystalStructureAdapter& stru,
                                   int anchor, const SiteIndices& sel,
                                   const Config& cfg) {
  const Lattice& L = stru.getLattice();
  const double buffzone = L.ucMaxDiagonalLength();
  PointsInSphere sph(cfg.rmin - buffzone, cfg.rmax + buffzone, L);
  const R3::Vector r0 = stru.getEquivalentAtoms(anchor)[0].xyz_cartn;
  std::vector<Bond> rv;
  for (int j : sel) {
    for (const Atom& a : stru.getEquivalentAtoms(j)) {
      for (sph.rewind(); !sph.finished(); sph.next()) {
        const R3::Vector rc = L.cartesian(sph.mno());
        Bond b;
        b.site0 = anchor;
        b.site1 = j;
        for (int k = 0; k < 3; ++k) b.r1[k] = rc[k] + a.xyz_cartn[k];
        for (int k = 0; k < 3; ++k) b.r01[k] = b.r1[k] - r0[k];
        b.distance = R3::norm(b.r01);
        b.u1 = a.uij_cartn;
        if (inRange(b, cfg)) rv.push_back(b);
      }
    }
  }
  return rv;
}

std::vector<Bond> generated(BaseBondGenerator& bnds) {
  std::vector<Bond> rv;
  for (bnds.rewind(); !bnds.finished(); bnds.next()) {
    Bond b;
    b.site0 = bnds.site0();
    b.site1 = bnds.site1();
    b.r1 = bnds.r1();
    b.r01 = bnds.r01();
    b.distance = bnds.distance();
    b.u1 = bnds.Ucartesian1();
    rv.push_back(b);
  }
  return rv;
}

/// Equality of the bond vectors.
const diffpy::mathutils::EpsilonEqual allclose(1e-12);

bool sameBond(const Bond& a, const Bond& e) {
  bool rv = (a.site0 == e.site0) && (a.site1 == e.site1) &&
            (std::fabs(a.distance - e.distance) <= 1e-14 * e.distance) &&
            allclose(a.r1, e.r1) && allclose(a.r01, e.r01);
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) rv = rv && (a.u1(k, l) == e.u1(k, l));
  }
  return rv;
}

void checkSameBonds(const std::vector<Bond>& actual,
                    const std::vector<Bond>& expected) {
  REQUIRE(actual.size() == expected.size());
  size_t i = 0;
  while (i < actual.size() && sameBond(actual[i], expected[i])) ++i;
  if (i == actual.size()) return;
  INFO("first different bond " << i);
  const Bond& a = actual[i];
  const Bond& e = expected[i];
  CHECK(a.site0 == e.site0);
  CHECK(a.site1 == e.site1);
  CHECK(allclose(a.r1, e.r1));
  CHECK(allclose(a.r01, e.r01));
  CHECK(std::fabs(a.distance - e.distance) <= 1e-14 * e.distance);
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) CHECK(a.u1(k, l) == e.u1(k, l));
  }
  FAIL("bond sequences differ");
}

/// Crystal structure from a PDFfit file with a few symmetry operations.
CrystalStructureAdapterPtr loadTestCrystalStructure(const std::string& f) {
  StructureAdapterPtr p = loadTestPeriodicStructure(f);
  const PeriodicStructureAdapter& ps =
    dynamic_cast<const PeriodicStructureAdapter&>(*p);
  const Lattice& L = ps.getLattice();
  CrystalStructureAdapterPtr rv(new CrystalStructureAdapter);
  rv->setLatPar(L.a(), L.b(), L.c(), L.alpha(), L.beta(), L.gamma());
  for (int i = 0; i < ps.countSites(); ++i) rv->append(ps[i]);
  R3::Vector t = R3::zerovector;
  rv->addSymOp(R3::identity(), t);
  rv->addSymOp(-1.0 * R3::identity(), t);
  R3::Matrix R = R3::identity();
  R(0, 0) = -1.0;
  R(1, 1) = -1.0;
  t[0] = 0.5;
  t[1] = 0.5;
  rv->addSymOp(R, t);
  return rv;
}

/// Select bonds up to 6 from the site 1, used with CaTiO3.
void nearSite1(BaseBondGenerator& bnds) {
  bnds.setRmax(6.0);
  bnds.selectAnchorSite(1);
}

/// New bond generator for the bonds up to 6 from the site 1.
BaseBondGeneratorPtr nearSite1(const StructureAdapter& stru) {
  BaseBondGeneratorPtr rv = stru.createBondGenerator();
  nearSite1(*rv);
  return rv;
}

template <class Tstru, class Treference>
void compareAllConfigurations(const Tstru& stru, Treference reference) {
  BaseBondGeneratorPtr bnds = stru.createBondGenerator();
  const int n = stru.countSites();
  for (const Config& cfg : CONFIGS) {
    bnds->setRmin(cfg.rmin);
    bnds->setRmax(cfg.rmax);
    for (int anchor = 0; anchor < n; ++anchor) {
      bnds->selectAnchorSite(anchor);
      const std::vector<SiteIndices> sels = siteSelections(n, anchor);
      for (size_t isel = 0; isel < sels.size(); ++isel) {
        INFO("rmin " << cfg.rmin << ", rmax " << cfg.rmax << ", anchor "
                     << anchor << ", selection " << isel);
        if (isel == 0) {
          bnds->selectSiteRange(0, anchor + 1);
        } else {
          bnds->selectSites(sels[isel]);
        }
        checkSameBonds(generated(*bnds),
                       reference(stru, anchor, sels[isel], cfg));
      }
    }
  }
}

}  // namespace

TEST_CASE("TestBondGeneratorBatching::test_periodic",
          "[BondGeneratorBatching]") {
  testAllKernels([] {
    for (const char* f : STRUCTURES) {
      INFO("structure " << f);
      StructureAdapterPtr p = loadTestPeriodicStructure(f);
      const PeriodicStructureAdapter& stru =
        dynamic_cast<const PeriodicStructureAdapter&>(*p);
      compareAllConfigurations(stru, referencePeriodic);
    }
  });
}

TEST_CASE("TestBondGeneratorBatching::test_crystal",
          "[BondGeneratorBatching]") {
  testAllKernels([] {
    for (const char* f : STRUCTURES) {
      INFO("structure " << f);
      CrystalStructureAdapterPtr stru = loadTestCrystalStructure(f);
      // the symmetry operations must give extra equivalent positions
      int cnteq = 0;
      for (int i = 0; i < stru->countSites(); ++i) {
        cnteq += stru->getEquivalentAtoms(i).size();
      }
      CHECK(cnteq > stru->countSites());
      compareAllConfigurations(*stru, referenceCrystal);
    }
  });
}

TEST_CASE("TestBondGeneratorBatching::test_crystal_changes",
          "[BondGeneratorBatching]") {
  // A new bond generator must follow every change, also in place.
  testAllKernels([] {
    CrystalStructureAdapterPtr stru = loadTestCrystalStructure("CaTiO3.stru");
    auto bonds = [](const CrystalStructureAdapter& s) {
      BaseBondGeneratorPtr bnds = s.createBondGenerator();
      bnds->setRmax(6.0);
      bnds->selectAnchorSite(0);
      return generated(*bnds);
    };
    // the same structure rebuilt without cached symmetry positions
    auto rebuilt = [](const CrystalStructureAdapter& s) {
      const Lattice& L = s.getLattice();
      CrystalStructureAdapterPtr rv(new CrystalStructureAdapter);
      rv->setLatPar(L.a(), L.b(), L.c(), L.alpha(), L.beta(), L.gamma());
      rv->assign(s.begin(), s.end());
      rv->setSymmetryPrecision(s.getSymmetryPrecision());
      for (int i = 0; i < s.countSymOps(); ++i) rv->addSymOp(s.getSymOp(i));
      return rv;
    };
    auto same = [](const std::vector<Bond>& a, const std::vector<Bond>& b) {
      return a.size() == b.size() &&
             std::equal(a.begin(), a.end(), b.begin(), sameBond);
    };
    std::vector<Bond> last = bonds(*stru);
    REQUIRE(!last.empty());
    checkSameBonds(last, bonds(*rebuilt(*stru)));
    R3::Matrix R = R3::identity();
    R(2, 2) = -1.0;
    const std::pair<const char*, std::function<void()>> changes[] = {
      {"atom position", [&] { stru->at(1).xyz_cartn[0] += 0.3; }},
      {"atom displacement", [&] { stru->at(2).uij_cartn(0, 0) += 0.01; }},
      {"lattice", [&] { stru->setLatPar(5.5, 7.7, 5.4, 90.0, 95.0, 90.0); }},
      {"new atom", [&] { stru->append(stru->at(0)); }},
      {"removed atom", [&] { stru->erase(0); }},
      {"symmetry operation", [&] { stru->addSymOp(R, R3::zerovector); }},
      {"symmetry precision", [&] { stru->setSymmetryPrecision(1.0); }},
    };
    for (const auto& [name, change] : changes) {
      INFO("changed " << name);
      change();
      const std::vector<Bond> now = bonds(*stru);
      CHECK(!same(last, now));
      checkSameBonds(now, bonds(*rebuilt(*stru)));
      last = now;
    }
  });
}

TEST_CASE("TestBondGeneratorBatching::test_crystal_refresh",
          "[BondGeneratorBatching]") {
  // A generator keeps its bonds when the structure expands symmetry again.
  testAllKernels([] {
    CrystalStructureAdapterPtr stru = loadTestCrystalStructure("CaTiO3.stru");
    BaseBondGeneratorPtr bnds = nearSite1(*stru);
    const std::vector<Bond> expected = generated(*bnds);
    REQUIRE(!expected.empty());
    // An unchanged structure and its copy keep their symmetry positions when
    // they create another generator.
    const Atom* symatoms = stru->getEquivalentAtoms(0).data();
    BaseBondGeneratorPtr bnds1 = stru->createBondGenerator();
    CHECK(symatoms == stru->getEquivalentAtoms(0).data());
    auto copy = std::make_shared<CrystalStructureAdapter>(*stru);
    const Atom* copyatoms = copy->getEquivalentAtoms(0).data();
    BaseBondGeneratorPtr bnds2 = copy->createBondGenerator();
    CHECK(copyatoms == copy->getEquivalentAtoms(0).data());
    // changes of the structure do not reach the first generator
    stru->updateSymmetryPositions();
    checkSameBonds(generated(*bnds), expected);
    stru->at(2).uij_cartn(0, 0) += 0.01;
    stru->at(3).xyz_cartn[1] += 0.2;
    bnds1 = stru->createBondGenerator();
    checkSameBonds(generated(*bnds), expected);
    stru->addSymOp(R3::identity(), R3::Vector(0.5, 0.5, 0.5));
    CHECK(stru->siteMultiplicity(0) > 0);
    checkSameBonds(generated(*bnds), expected);
  });
}

TEST_CASE("TestBondGeneratorBatching::test_crystal_constructor",
          "[BondGeneratorBatching]") {
  // The public constructor finds the symmetry positions when needed.
  testAllKernels([] {
    CrystalStructureAdapterPtr stru = loadTestCrystalStructure("CaTiO3.stru");
    CrystalStructureBondGenerator direct(stru);
    nearSite1(direct);
    const std::vector<Bond> expected = generated(*nearSite1(*stru));
    REQUIRE(!expected.empty());
    checkSameBonds(generated(direct), expected);
  });
  // a structure without symmetry operations is rejected
  StructureAdapterPtr periodic = loadTestPeriodicStructure("CaTiO3.stru");
  CHECK_THROWS_AS(CrystalStructureBondGenerator(periodic),
                  std::invalid_argument);
  // The generators point into their own arrays and cannot be copied.
  static_assert(!std::is_copy_constructible_v<PeriodicStructureBondGenerator>);
  static_assert(!std::is_copy_assignable_v<PeriodicStructureBondGenerator>);
  static_assert(!std::is_copy_constructible_v<CrystalStructureBondGenerator>);
}

TEST_CASE("TestBondGeneratorBatching::test_crystal_multiplicity",
          "[BondGeneratorBatching]") {
  // Site multiplicities follow atoms moved in place to and from a special
  // position without any bond generator.
  auto stru = std::make_shared<CrystalStructureAdapter>();
  stru->setLatPar(4.0, 5.0, 6.0, 90.0, 90.0, 90.0);
  stru->addSymOp(R3::identity(), R3::zerovector);
  stru->addSymOp(-1.0 * R3::identity(), R3::zerovector);
  Atom a;
  a.atomtype = "C";
  stru->append(a);
  stru->append(a);
  stru->at(1).xyz_cartn = R3::Vector(1.0, 1.0, 1.0);
  CHECK(1 == stru->siteMultiplicity(0));
  CHECK(2 == stru->siteMultiplicity(1));
  CHECK(3.0 == stru->totalOccupancy());
  stru->at(0).xyz_cartn = R3::Vector(0.4, 0.5, 0.6);
  CHECK(2 == stru->siteMultiplicity(0));
  CHECK(2u == stru->getEquivalentAtoms(0).size());
  CHECK(4.0 == stru->totalOccupancy());
  // fractional (1/2, 1/2, 1/2) is a special position as well
  stru->at(1).xyz_cartn = R3::Vector(2.0, 2.5, 3.0);
  CHECK(1 == stru->siteMultiplicity(1));
  CHECK(3.0 == stru->totalOccupancy());
  BaseBondGeneratorPtr bnds = stru->createBondGenerator();
  for (int anchor : {0, 1}) {
    bnds->selectAnchorSite(anchor);
    CHECK(stru->siteMultiplicity(anchor) == bnds->multiplicity());
  }
}

TEST_CASE("TestBondGeneratorBatching::test_empty_first_rewind",
          "[BondGeneratorBatching]") {
  // Empty selection at the first rewind.
  testAllKernels([] {
    StructureAdapterPtr p = loadTestPeriodicStructure("CaTiO3.stru");
    const PeriodicStructureAdapter& stru =
      dynamic_cast<const PeriodicStructureAdapter&>(*p);
    CrystalStructureAdapterPtr cstru = loadTestCrystalStructure("CaTiO3.stru");
    BaseBondGeneratorPtr bnds = stru.createBondGenerator();
    BaseBondGeneratorPtr cbnds = cstru->createBondGenerator();
    for (BaseBondGeneratorPtr b : {bnds, cbnds}) {
      b->setRmax(300.0);
      b->selectAnchorSite(1);
      b->selectSiteRange(2, 2);
      b->rewind();
      CHECK(b->finished());
      b->setRmax(6.0);
    }
    const Config cfg = {0.0, 6.0};
    SiteIndices sel = {0, 3, 7};
    bnds->selectSites(sel);
    cbnds->selectSites(sel);
    checkSameBonds(generated(*bnds), referencePeriodic(stru, 1, sel, cfg));
    checkSameBonds(generated(*cbnds), referenceCrystal(*cstru, 1, sel, cfg));
  });
}

TEST_CASE("TestBondGeneratorBatching::test_range_bounds",
          "[BondGeneratorBatching]") {
  // Bonds exactly at rmin or rmax are included, one ulp beyond are not.
  StructureAdapterPtr p = loadTestPeriodicStructure("ZnS_wurtzite.stru");
  CrystalStructureAdapterPtr c = loadTestCrystalStructure("ZnS_wurtzite.stru");
  testAllKernels([&] {
    for (StructureAdapterPtr stru : {p, StructureAdapterPtr(c)}) {
      BaseBondGeneratorPtr bnds = stru->createBondGenerator();
      bnds->setRmax(5.0);
      bnds->selectAnchorSite(0);
      bnds->selectSiteRange(0, stru->countSites());
      const std::vector<Bond> all = generated(*bnds);
      REQUIRE(all.size() > 10);
      for (size_t ib : {size_t(0), all.size() / 2, all.size() - 1}) {
        const double d = all[ib].distance;
        INFO("distance " << d);
        auto count = [&](double rmin, double rmax) {
          bnds->setRmin(rmin);
          bnds->setRmax(rmax);
          int rv = 0;
          for (const Bond& b : generated(*bnds)) {
            rv += (b.site1 == all[ib].site1 && b.r01[0] == all[ib].r01[0] &&
                   b.r01[1] == all[ib].r01[1] && b.r01[2] == all[ib].r01[2]);
          }
          return rv;
        };
        const double below = std::nextafter(d, 0.0);
        const double above = std::nextafter(d, 2 * d);
        CHECK(1 == count(0.0, d));
        CHECK(0 == count(0.0, below));
        CHECK(1 == count(d, 5.0));
        CHECK(0 == count(above, 5.0));
        CHECK(1 == count(d, d));
        CHECK(1 == count(below, above));
      }
    }
  });
}

TEST_CASE("TestBondGeneratorBatching::test_large_sphere",
          "[BondGeneratorBatching]") {
  // Spheres with many points and few candidates per point.
  StructureAdapterPtr p = loadTestPeriodicStructure("Ni.stru");
  const PeriodicStructureAdapter& stru =
    dynamic_cast<const PeriodicStructureAdapter&>(*p);
  CrystalStructureAdapterPtr cstru = loadTestCrystalStructure("Ni.stru");
  const Config cfg = {100.0, 110.0};
  const Lattice& L = stru.getLattice();
  const double buffzone = L.ucMaxDiagonalLength();
  PointsInSphere sph(cfg.rmin - buffzone, cfg.rmax + buffzone, L);
  int npts = 0;
  for (sph.rewind(); !sph.finished(); sph.next()) ++npts;
  REQUIRE(npts > (1 << 16));
  SiteIndices sel = {0, 1, 2, 3};
  testAllKernels([&] {
    BaseBondGeneratorPtr bnds = stru.createBondGenerator();
    BaseBondGeneratorPtr cbnds = cstru->createBondGenerator();
    for (BaseBondGeneratorPtr b : {bnds, cbnds}) {
      b->setRmin(cfg.rmin);
      b->setRmax(cfg.rmax);
      b->selectSites(sel);
    }
    for (int anchor : {0, 3}) {
      INFO("anchor " << anchor);
      bnds->selectAnchorSite(anchor);
      cbnds->selectAnchorSite(anchor);
      const std::vector<Bond> expected =
        referencePeriodic(stru, anchor, sel, cfg);
      REQUIRE(!expected.empty());
      checkSameBonds(generated(*bnds), expected);
      checkSameBonds(generated(*bnds), expected);
      checkSameBonds(generated(*cbnds),
                     referenceCrystal(*cstru, anchor, sel, cfg));
    }
  });
}

namespace {

/// Set isotropic displacement parameters of an atom.
void setUiso(Atom& a, double uiso) {
  a.anisotropy = false;
  a.uij_cartn = R3::zeromatrix();
  for (int k = 0; k < 3; ++k) a.uij_cartn(k, k) = uiso;
}

/// Add two atoms with Uiso = 0.01, which are 1 apart along the x-axis.
template <class Tstru>
void appendDimer(Tstru& stru) {
  Atom a;
  setUiso(a, 0.01);
  stru.append(a);
  a.xyz_cartn[0] = 1.0;
  stru.append(a);
}

/// msd() at the only bond of the dimer after a rewind.
double dimerMSD(BaseBondGenerator& bnds) {
  bnds.rewind();
  REQUIRE(!bnds.finished());
  REQUIRE(1 == bnds.site1());
  return bnds.msd();
}

}  // namespace

TEST_CASE("TestBondGeneratorBatching::test_msd", "[BondGeneratorBatching]") {
  // msd() of a bond along (1, 2, 2) for isotropic, mixed and anisotropic
  // sites in every kind of bond generator.  The anisotropic tensor gives
  // s'Us / s's = 0.18 / 9 = 0.02 along the bond.
  testAllKernels([] {
    R3::Matrix U = R3::zeromatrix();
    U(0, 0) = 0.01;
    U(1, 1) = 0.03;
    U(2, 2) = 0.01;
    U(0, 1) = U(1, 0) = 0.0025;
    auto atomic = std::make_shared<AtomicStructureAdapter>();
    auto periodic = std::make_shared<PeriodicStructureAdapter>();
    periodic->setLatPar(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    auto crystal = std::make_shared<CrystalStructureAdapter>();
    crystal->setLatPar(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    crystal->addSymOp(R3::identity(), R3::zerovector);
    const std::pair<const char*, std::shared_ptr<AtomicStructureAdapter>>
      adapters[] = {
        {"atomic", atomic}, {"periodic", periodic}, {"crystal", crystal}};
    const double eps = 1e-15;
    for (const auto& [name, stru] : adapters) {
      INFO("adapter " << name);
      appendDimer(*stru);
      stru->at(1).xyz_cartn = R3::Vector(1.0, 2.0, 2.0);
      // a new generator, the crystal one keeps earlier symmetry positions
      auto msd = [&stru] {
        BaseBondGeneratorPtr bnds = stru->createBondGenerator();
        bnds->setRmax(3.5);
        bnds->selectAnchorSite(0);
        return dimerMSD(*bnds);
      };
      CHECK(std::fabs(msd() - 0.02) < eps);
      stru->at(1).anisotropy = true;
      stru->at(1).uij_cartn = U;
      CHECK(std::fabs(msd() - 0.03) < eps);
      stru->at(0).anisotropy = true;
      stru->at(0).uij_cartn = U;
      CHECK(std::fabs(msd() - 0.04) < eps);
    }
  });
}

// End of file
