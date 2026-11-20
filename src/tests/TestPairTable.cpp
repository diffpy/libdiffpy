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
 * TestPairTable -- calculators that add all pairs of a structure at once
 *
 * The pairs must be those of the bond generators and the results those of
 * the pairs added one by one, within rounding errors.
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/GaussianProfile.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/NoMetaStructureAdapter.hpp>
#include <diffpy/srreal/NoSymmetryStructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#if __has_include(<diffpy/srreal/kernels/PairTable.hpp>)
#define LIBDIFFPY_TEST_PAIR_TABLE
#include <diffpy/srreal/kernels/PairTable.hpp>
#include <diffpy/srreal/kernels/PairGeometry.hpp>
#endif

#include "pair_test_helpers.hpp"

using namespace diffpy::srreal;
using namespace pairtests;
#ifdef LIBDIFFPY_TEST_PAIR_TABLE
using kernels::PairGroups;
using kernels::PairTable;
#endif

namespace {

const double TOLERANCE = 1e-12;

// Pairs ---------------------------------------------------------------------

#ifdef LIBDIFFPY_TEST_PAIR_TABLE

struct Pair {
  int site0;
  int site1;
  double dist;
  double msd;
  double weight;
};

/// Pairs from the bond generator, as the evaluator adds them one by one.
std::vector<Pair> generatedPairs(const StructureAdapter& stru, double rmin,
                                 double rmax) {
  std::vector<Pair> rv;
  BaseBondGeneratorPtr bnds = stru.createBondGenerator();
  bnds->setRmin(rmin);
  bnds->setRmax(rmax);
  for (int i0 = 0; i0 < stru.countSites(); ++i0) {
    bnds->selectAnchorSite(i0);
    bnds->selectSiteRange(0, i0 + 1);
    for (bnds->rewind(); !bnds->finished(); bnds->next()) {
      const int scale = (i0 == bnds->site1()) ? 1 : 2;
      rv.push_back({i0, bnds->site1(), bnds->distance(), bnds->msd(),
                    double(scale * bnds->multiplicity())});
    }
  }
  return rv;
}

std::vector<Pair> tablePairs(PairTable& table, const StructureAdapter& stru,
                             double rmin, double rmax) {
  std::vector<Pair> rv;
  const auto g = table.geometryArgs(stru, rmin, rmax);
  for (long long i = 0; i < kernels::pairCount(g); ++i) {
    const auto p = kernels::pairAt(g, i);
    if (p.ok) rv.push_back({p.site0, p.site1, p.dist, p.msd, double(p.weight)});
  }
  return rv;
}

/// The pairs are the same within rounding errors.  Distances that differ
/// in their last bits may sort differently, so the distances and the
/// displacements of every pair of sites are compared in their own orders.
void checkSamePairs(const std::vector<Pair>& expected,
                    const std::vector<Pair>& pairs) {
  REQUIRE(expected.size() == pairs.size());
  auto sorted = [](std::vector<Pair> v, double Pair::* member) {
    std::sort(v.begin(), v.end(), [=](const Pair& a, const Pair& b) {
      return std::tie(a.site0, a.site1, a.weight, a.*member) <
             std::tie(b.site0, b.site1, b.weight, b.*member);
    });
    return v;
  };
  int mismatches = 0;
  double maxerror = 0.0;
  for (double Pair::* member : {&Pair::dist, &Pair::msd}) {
    const std::vector<Pair> e = sorted(expected, member);
    const std::vector<Pair> p = sorted(pairs, member);
    for (size_t i = 0; i < p.size(); ++i) {
      mismatches += (e[i].site0 != p[i].site0) || (e[i].site1 != p[i].site1) ||
                    (e[i].weight != p[i].weight);
      const double error = std::fabs(e[i].*member - p[i].*member);
      maxerror = std::max(maxerror, error / std::fabs(e[i].*member));
    }
  }
  CHECK(0 == mismatches);
  CHECK(maxerror <= TOLERANCE);
}

#endif

// Calculators ---------------------------------------------------------------

/// Value and results of the calculator after its setup.
template <class C>
std::vector<QuantityType> resultsOf(C&& calc, const Setup& setup,
                                    StructureAdapterPtr s) {
  calc.setEvaluatorType(BASIC);
  setup(calc);
  calc.eval(s);
  return {calc.value(), calc.getPDF(), calc.getF()};
}

/// Check that a calculator of the class T adds all pairs at once and gives
/// what it does with the pairs one by one.
template <class T>
void checkAllPairs(const Setup& setup, StructureAdapterPtr s) {
  const auto expected = resultsOf(PairByPair<T>(), setup, s);
  AllPairs<T> calc;
  const auto results = resultsOf(calc, setup, s);
  CHECK(calc.addedAllPairs());
  REQUIRE(expected.size() == results.size());
  // the structure has some pairs
  REQUIRE(*std::ranges::max_element(expected[0]) > 0.0);
  for (size_t i = 0; i < results.size(); ++i) {
    CHECK(maxRelativeError(results[i], expected[i]) <= TOLERANCE);
  }
  const auto builtin = resultsOf(T(), setup, s);
  if constexpr (std::is_same_v<T, DebyePDFCalculator>) {
    CHECK(results[0] == builtin[0]);
    // Normalization may round the occupied scattering factors differently.
    for (size_t i = 1; i < results.size(); ++i) {
      CHECK(maxRelativeError(results[i], builtin[i]) <= TOLERANCE);
    }
  } else {
    CHECK(results == builtin);
  }
}

/// Check that a calculator of the class T adds the pairs one by one.
template <class T>
void checkPairByPair(const Setup& setup, StructureAdapterPtr s) {
  const auto expected = resultsOf(PairByPair<T>(), setup, s);
  REQUIRE(*std::ranges::max_element(expected[0]) > 0.0);
  AllPairs<T> calc;
  CHECK(expected == resultsOf(calc, setup, s));
  CHECK(!calc.addedAllPairs());
  CHECK(expected == resultsOf(T(), setup, s));
}

double maxAbsolute(const QuantityType& v) {
  double rv = 0.0;
  for (double x : v) rv = std::max(rv, std::fabs(x));
  return rv;
}

}  // namespace

// Tests ---------------------------------------------------------------------

#ifdef LIBDIFFPY_TEST_PAIR_TABLE
TEST_CASE("TestPairTable::test_pair_groups", "[PairTable]") {
  PairGroups groups;
  CHECK(groups.groups().empty());
  const double d = 2.5;
  const double dnext = std::nextafter(d, 3.0);
  groups.add(d, 0.1, 0, 1.0);
  groups.add(d, 0.1, 0, 2.0);
  // Adjacent representable inputs can fall on opposite sides of a cutoff.
  groups.add(dnext, 0.1, 0, 4.0);
  groups.add(d, std::nextafter(0.1, 1.0), 0, 8.0);
  REQUIRE(3u == groups.groups().size());
  CHECK(d == groups.groups()[0].dist);
  CHECK(0.1 == groups.groups()[0].width);
  CHECK(3.0 == groups.groups()[0].scale);
  CHECK(4.0 == groups.groups()[1].scale);
  CHECK(8.0 == groups.groups()[2].scale);
  groups.add(d * (1 + 1e-12), 0.1, 0, 1.0);
  groups.add(d, 0.1 * (1 + 1e-12), 0, 1.0);
  groups.add(d, 0.1, 1, 1.0);
  groups.add(-d, 0.1, 0, 1.0);
  groups.add(d, 0.0, 0, 1.0);
  CHECK(8u == groups.groups().size());
  CHECK(3.0 == groups.groups()[0].scale);
  // many groups keep their sums
  groups.clear();
  CHECK(groups.groups().empty());
  const int n = 5000;
  for (int repeat = 0; repeat < 3; ++repeat) {
    for (int i = 0; i < n; ++i) groups.add(1.0 + i, 0.1, i % 7, i);
  }
  REQUIRE(size_t(n) == groups.groups().size());
  for (int i = 0; i < n; ++i) {
    const PairGroups::Group& g = groups.groups()[i];
    CHECK((1.0 + i == g.dist && i % 7 == g.tag && 3.0 * i == g.scale));
  }
}

#endif

TEST_CASE("TestPairTable::test_debye_cutoff_of_adjacent_distances",
          "[PairTable]") {
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "H";
  atom.occupancy = std::nextafter(0.001, 1.0);
  stru->append(atom);
  atom.xyz_cartn[0] = 2.0;
  stru->append(atom);
  atom.xyz_cartn = R3::Vector(0, 30, 0);
  stru->append(atom);
  atom.xyz_cartn[0] = 2.0 + 4 * (std::nextafter(2.0, 3.0) - 2.0);
  stru->append(atom);
  DebyePDFCalculator grouped;
  PairByPair<DebyePDFCalculator> individual;
  for (DebyePDFCalculator* calc :
       {&grouped, static_cast<DebyePDFCalculator*>(&individual)}) {
    calc->setRmax(5);
    calc->setMaxExtension(0);
    calc->setQstep(0.1);
    calc->setScatteringFactorTableByType("neutron");
    calc->getScatteringFactorTable()->setCustomAs("H", "H", 1.0);
  }
  testAllKernels([&] {
    const QuantityType reference = individual.eval(stru);
    REQUIRE(maxAbsolute(reference) > 0);
    CHECK(maxRelativeError(grouped.eval(stru), reference) <= TOLERANCE);
  });
}

TEST_CASE("TestPairTable::test_ordered_scattering_at_cutoff", "[PairTable]") {
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "Ni";
  atom.occupancy = 0.0002;
  stru->append(atom);
  atom.xyz_cartn[0] = 3.0;
  atom.occupancy = 7.069469318503158e-05;
  stru->append(atom);
  DebyePDFCalculator batch;
  PairByPair<DebyePDFCalculator> individual;
  for (DebyePDFCalculator* c :
       {&batch, static_cast<DebyePDFCalculator*>(&individual)}) {
    c->setEvaluatorType(BASIC);
    c->setScatteringFactorTableByType("neutron");
    c->setRmax(8);
  }
  // Reversing the site factors rounds the first term above the cutoff.
  testAllKernels([&] {
    CHECK(batch.eval(stru) == individual.eval(stru));
    CHECK(maxAbsolute(batch.value()) == 0.0);
  });
}

TEST_CASE("TestPairTable::test_debye_width_at_cutoff", "[PairTable]") {
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  stru->append(atom);
  atom.xyz_cartn[0] = 2;
  stru->append(atom);
  BaseDebyeSum batch;
  PairByPair<BaseDebyeSum> individual;
  for (BaseDebyeSum* c : {&batch, static_cast<BaseDebyeSum*>(&individual)}) {
    c->setEvaluatorType(BASIC);
    c->setRmax(3);
    c->setQmin(18);
    c->setQmax(19);
    c->setQstep(1);
    c->setPeakWidthModelByType("constant");
    c->getPeakWidthModel()->setDoubleAttr("width", 0.63);
    // Midway between the cutoffs from multiplying and dividing by the
    // Gaussian conversion constants.
    c->setDebyePrecision(9.20959530020065e-6);
  }
  KernelGuard guard;
  for (const auto& name : availableCalculationKernels()) {
    INFO(name);
    setCalculationKernel(name);
    CHECK(individual.eval(stru).at(18) == 0.0);
    CHECK(batch.eval(stru).at(18) == 0.0);
  }
}

TEST_CASE("TestPairTable::test_wrapped_cropped_peak", "[PairTable]") {
  auto stru = std::make_shared<PeriodicStructureAdapter>();
  stru->setLatPar(8, 20, 20, 90, 90, 90);
  Atom atom;
  atom.atomtype = "Ni";
  atom.xyz_cartn[0] = 7.99;
  stru->append(atom);
  atom.xyz_cartn[0] = 2.02;
  stru->append(atom);
  PDFCalculator cached;
  cached.setEvaluatorType(BASIC);
  cached.setRmax(4);
  cached.setMaxExtension(0);
  cached.setPeakWidthModelByType("constant");
  cached.setDoubleAttr("width", 0.2);
  cached.setPeakProfileByType("croppedgaussian");
  cached.setDoubleAttr("peakprecision", 0.5);
  PDFCalculator fresh = cached;
  cached.eval(stru);
  (*stru)[0].xyz_cartn[0] = 0.02;
  CHECK(cached.eval(stru) == fresh.eval(stru));
}

TEST_CASE("TestPairTable::test_incremental_update_after_reordering",
          "[PairTable]") {
  const auto original = asCluster(*loadPeriodic("CaTiO3.stru"));
  auto reversed = std::make_shared<AtomicStructureAdapter>();
  reversed->assign(original->rbegin(), original->rend());
  auto moved = std::make_shared<AtomicStructureAdapter>(*reversed);
  (*moved)[2].xyz_cartn[2] += 0.09;
  auto check = [&](auto& calc) {
    calc.setRmax(8);
    calc.eval(original);
    calc.eval(reversed);
    REQUIRE(calc.getEvaluatorTypeUsed() == OPTIMIZED);
    const QuantityType incremental = calc.eval(moved);
    REQUIRE(calc.getEvaluatorTypeUsed() == OPTIMIZED);
    calc.setEvaluatorType(BASIC);
    CHECK(maxRelativeError(incremental, calc.eval(moved)) <= TOLERANCE);
  };
  PDFCalculator pdf;
  DebyePDFCalculator debye;
  check(pdf);
  check(debye);
}

#ifdef LIBDIFFPY_TEST_PAIR_TABLE
TEST_CASE("TestPairTable::test_pairs_of_bond_generator", "[PairTable]") {
  const std::pair<double, double> ranges[] = {
    {0.0, 9.7}, {5.0, 4.0}, {0.0, 0.0}, {6.1, 6.9}, {2.0, 12.0}};
  for (const auto& [name, stru] : allStructures()) {
    INFO(name);
    CHECK(PairTable::supports(*stru));
    PairTable table;
    for (const auto& [rmin, rmax] : ranges) {
      INFO("range " << rmin << " to " << rmax);
      const std::vector<Pair> expected = generatedPairs(*stru, rmin, rmax);
      checkSamePairs(expected, tablePairs(table, *stru, rmin, rmax));
      // and again from the remembered pairs
      checkSamePairs(expected, tablePairs(table, *stru, rmin, rmax));
    }
  }
}

TEST_CASE("TestPairTable::test_unsupported_structures", "[PairTable]") {
  struct DerivedPeriodic : PeriodicStructureAdapter {};
  struct DerivedAtomic : AtomicStructureAdapter {};
  const auto periodic = loadPeriodic("Ni.stru");
  auto derived = std::make_shared<DerivedPeriodic>();
  derived->setLatPar(3.52, 3.52, 3.52, 90, 90, 90);
  derived->assign(periodic->begin(), periodic->end());
  CHECK(!PairTable::supports(*derived));
  CHECK(!PairTable::supports(*nometa(derived)));
  CHECK(!PairTable::supports(*nosymmetry(derived)));
  auto cluster = std::make_shared<DerivedAtomic>();
  cluster->assign(periodic->begin(), periodic->end());
  CHECK(!PairTable::supports(*cluster));
  // structures without sites have no pairs
  CHECK(PairTable::supports(*emptyStructureAdapter()));
  CHECK(PairTable::supports(*std::make_shared<DerivedAtomic>()));
  PairTable table;
  CHECK(tablePairs(table, *emptyStructureAdapter(), 0.0, 10.0).empty());
  // The calculators add the pairs of such structures one by one.
  checkPairByPair<PDFCalculator>(setupDefault, derived);
  checkPairByPair<DebyePDFCalculator>(setupDefault, cluster);
}

TEST_CASE("TestPairTable::test_remembered_pairs", "[PairTable]") {
  // Small changes of a periodic structure keep its pairs, large ones or
  // a longer range search them again.
  for (const char* kind : {"periodic", "crystal"}) {
    INFO(kind);
    const auto periodic = loadPeriodic("CaTiO3.stru");
    std::shared_ptr<PeriodicStructureAdapter> stru = periodic;
    // The first sites of the crystal are at special positions, where a
    // site that moves gets more symmetry positions and other pairs.
    int i0 = 0;
    if (std::string(kind) == "crystal") {
      stru = asCrystal(*periodic);
      i0 = 4;
    }
    PairTable table;
    const double rmin = 0.0, rmax = 8.0;
    auto check = [&](int added = -1) {
      const int previous = table.countSearches();
      checkSamePairs(generatedPairs(*stru, rmin, rmax),
                     tablePairs(table, *stru, rmin, rmax));
      if (added >= 0) CHECK(previous + added == table.countSearches());
    };
    check(1);
    check(0);
    // changes within the margin of the search
    stru->at(i0 + 1).xyz_cartn[0] += 0.05;
    // An equivalent site of the crystal crosses the unit cell here.
    check(stru == periodic ? 0 : 1);
    stru->at(i0).xyz_cartn[2] += 0.03;
    stru->at(i0 + 1).uij_cartn *= 1.3;
    stru->at(i0 + 2).occupancy = 0.5;
    check(0);
    const Lattice L = stru->getLattice();
    for (Atom& a : *stru) stru->toFractional(a);
    stru->setLatPar(1.0005 * L.a(), 0.9995 * L.b(), L.c(), L.alpha(), 90.02,
                    L.gamma());
    for (Atom& a : *stru) stru->toCartesian(a);
    check();
    // an atom that crosses the unit cell appears at its other side
    const R3::Vector xyz0 = stru->at(i0).xyz_cartn;
    stru->at(i0).xyz_cartn =
      xyz0 - stru->getLattice().cartesian(R3::Vector(0.0, 0.0, 1e-5));
    check();
    stru->at(i0).xyz_cartn -= R3::Vector(0.0, 0.0, 0.06);
    check();
    stru->at(i0).xyz_cartn = xyz0;
    check();
    // a shorter range needs no search, a longer one does
    const int beforeRangeChange = table.countSearches();
    checkSamePairs(generatedPairs(*stru, 1.0, 6.5),
                   tablePairs(table, *stru, 1.0, 6.5));
    CHECK(beforeRangeChange == table.countSearches());
    checkSamePairs(generatedPairs(*stru, 0.0, 9.0),
                   tablePairs(table, *stru, 0.0, 9.0));
    CHECK(beforeRangeChange + 1 == table.countSearches());
    check(0);
    // changes beyond the margin
    stru->at(i0 + 1).xyz_cartn[1] += 0.4;
    check(1);
    for (Atom& a : *stru) stru->toFractional(a);
    stru->setLatPar(1.03 * L.a(), L.b(), L.c(), L.alpha(), L.beta(), L.gamma());
    for (Atom& a : *stru) stru->toCartesian(a);
    check(1);
    // A site that moves by less than the margin can still bring pairs
    // from beyond the margin into the range, with the site that moves
    // the other way.
    stru->at(i0 + 1).xyz_cartn[2] += 0.2;
    check(1);
    stru->at(i0 + 2).xyz_cartn[2] -= 0.1;
    check(0);
    stru->at(i0 + 2).xyz_cartn[2] -= 0.1;
    check(1);
    // other sites
    stru->append(stru->at(i0));
    stru->at(stru->countSites() - 1).xyz_cartn[1] += 1.0;
    check(1);
    stru->erase(1);
    check(1);
    // The pairs below the range of a search are not remembered.
    PairTable other;
    checkSamePairs(generatedPairs(*stru, 4.0, 7.0),
                   tablePairs(other, *stru, 4.0, 7.0));
    checkSamePairs(generatedPairs(*stru, 3.9, 7.1),
                   tablePairs(other, *stru, 3.9, 7.1));
    CHECK(1 == other.countSearches());
    checkSamePairs(generatedPairs(*stru, 3.0, 7.0),
                   tablePairs(other, *stru, 3.0, 7.0));
    CHECK(2 == other.countSearches());
    stru->at(i0).xyz_cartn[0] += 0.07;
    checkSamePairs(generatedPairs(*stru, 2.9, 7.0),
                   tablePairs(other, *stru, 2.9, 7.0));
    CHECK(2 == other.countSearches());
    checkSamePairs(generatedPairs(*stru, 2.85, 7.0),
                   tablePairs(other, *stru, 2.85, 7.0));
    CHECK(3 == other.countSearches());
  }
}

TEST_CASE("TestPairTable::test_strained_lattice", "[PairTable]") {
  // The pairs are searched again when the changes of the positions and of
  // the translations together exceed the margin, even if neither does.
  auto stru = std::make_shared<PeriodicStructureAdapter>();
  stru->setLatPar(3.0, 3.0, 3.0, 90.0, 90.0, 90.0);
  Atom a;
  a.atomtype = "Ni";
  a.uij_cartn = 0.005 * R3::identity();
  stru->append(a);
  a.xyz_cartn = R3::Vector(1.5, 0.0, 0.0);
  stru->append(a);
  PairTable table;
  auto check = [&](double alat, int searches) {
    INFO("lattice parameter " << alat);
    for (Atom& ai : *stru) stru->toFractional(ai);
    stru->setLatPar(alat, 3.0, 3.0, 90.0, 90.0, 90.0);
    for (Atom& ai : *stru) stru->toCartesian(ai);
    checkSamePairs(generatedPairs(*stru, 0.0, 8.0),
                   tablePairs(table, *stru, 0.0, 8.0));
    CHECK(searches == table.countSearches());
  };
  check(3.0, 1);
  // The second atom is up to three cells from the first one within the
  // range: the translation changes by 0.15 and the position by 0.025.
  check(3.05, 1);
  // by 0.21 and by 0.035, which is 0.28 for the distance of two atoms
  check(3.07, 2);
  check(3.05, 2);
}

TEST_CASE("TestPairTable::test_new_neighbors_after_contraction",
          "[PairTable]") {
  auto stru = std::make_shared<PeriodicStructureAdapter>();
  stru->append(Atom());
  // The first search finds either no translations or only translations
  // along b. Shortening a brings a new pair of neighbors into range.
  for (double b : {6.0, 2.0}) {
    PairTable table;
    stru->setLatPar(6, b, 6, 90, 90, 90);
    const auto before = tablePairs(table, *stru, 0, 5);
    stru->setLatPar(4.9, b, 6, 90, 90, 90);
    const auto expected = generatedPairs(*stru, 0, 5);
    REQUIRE(expected.size() > before.size());
    checkSamePairs(expected, tablePairs(table, *stru, 0, 5));
  }
}

TEST_CASE("TestPairTable::test_retry_interrupted_search", "[PairTable]") {
  // A stack adapter cannot supply the shared owner needed by its generator.
  PeriodicStructureAdapter stack;
  stack.setLatPar(3, 3, 3, 90, 90, 90);
  stack.append(Atom{});
  PairTable table;
  REQUIRE_THROWS_AS(table.geometryArgs(stack, 0, 3.1), std::bad_weak_ptr);
  auto owned = std::make_shared<PeriodicStructureAdapter>(stack);
  const auto expected = generatedPairs(*owned, 0, 3.1);
  REQUIRE(expected.size() == 6);
  checkSamePairs(expected, tablePairs(table, *owned, 0, 3.1));
}

TEST_CASE("TestPairTable::test_near_coincidence_on_a_short_axis",
          "[PairTable]") {
  auto stru = std::make_shared<PeriodicStructureAdapter>();
  stru->setLatPar(0.005, 1, 1, 90, 90, 90);
  Atom site;
  site.xyz_cartn = R3::Vector(0.00125, 0.2, 0.2);
  stru->append(site);
  site.xyz_cartn[0] += 1.3e-8;
  stru->append(site);
  PairTable table;
  REQUIRE(generatedPairs(*stru, 0, 0.001).empty());
  REQUIRE(tablePairs(table, *stru, 0, 0.001).empty());
  // The dormant pair follows the regular neighbors of the last anchor,
  // preserving their assignment to reduction rows.
  const auto g = table.geometryArgs(*stru, 0, 0.001);
  REQUIRE(g.cntpairs > 0);
  const auto last = kernels::pairAt(g, g.cntpairs - 1);
  CHECK(last.site0 == 1);
  CHECK(last.site1 == 0);
  CHECK(last.dist <= g.zeroeps);
  stru->at(1).xyz_cartn[0] += 0.0001;
  const auto expected = generatedPairs(*stru, 0, 0.001);
  REQUIRE(expected.size() == 1);
  checkSamePairs(expected, tablePairs(table, *stru, 0, 0.001));
  CHECK(table.countSearches() == 1);
}

TEST_CASE("TestPairTable::test_coinciding_sites", "[PairTable]") {
  // Sites at the same position make no bond until they move apart.
  for (const char* kind : {"periodic", "crystal", "cluster"}) {
    INFO(kind);
    const auto periodic = loadPeriodic("CaTiO3.stru");
    // a second atom on the first site and one a lattice vector away
    periodic->append(periodic->at(0));
    periodic->at(periodic->countSites() - 1).atomtype = "Sr";
    Atom away = periodic->at(3);
    away.xyz_cartn += periodic->getLattice().cartesian(R3::Vector(1, -2, 0));
    periodic->append(away);
    std::shared_ptr<AtomicStructureAdapter> stru = periodic;
    if (std::string(kind) == "crystal") {
      auto crystal = asCrystal(*periodic);
      crystal->append(crystal->at(0));
      stru = crystal;
    }
    if (std::string(kind) == "cluster") stru = asCluster(*periodic);
    PairTable table;
    auto check = [&] {
      checkSamePairs(generatedPairs(*stru, 0.0, 7.0),
                     tablePairs(table, *stru, 0.0, 7.0));
    };
    check();
    const size_t cntpairs = tablePairs(table, *stru, 0.0, 7.0).size();
    const int last = stru->countSites() - 1;
    stru->at(last).xyz_cartn[0] += 0.01;
    stru->at(last - 1).xyz_cartn[1] += 0.02;
    check();
    CHECK(cntpairs < tablePairs(table, *stru, 0.0, 7.0).size());
    // The sites of the crystal get other symmetry positions and their
    // pairs are searched again.
    if (std::string(kind) != "crystal") CHECK(1 >= table.countSearches());
    stru->at(last).xyz_cartn[0] -= 0.01;
    stru->at(last - 1).xyz_cartn[1] -= 0.02;
    check();
    CHECK(cntpairs == tablePairs(table, *stru, 0.0, 7.0).size());
  }
}

TEST_CASE("TestPairTable::test_special_positions", "[PairTable]") {
  // A site that moves to a special position has fewer positions.
  auto stru = std::make_shared<CrystalStructureAdapter>();
  stru->setLatPar(4.0, 4.5, 5.0, 90.0, 90.0, 90.0);
  stru->addSymOp(R3::identity(), R3::zerovector);
  stru->addSymOp(-1.0 * R3::identity(), R3::zerovector);
  Atom a;
  a.atomtype = "Ni";
  a.uij_cartn = 0.005 * R3::identity();
  for (const R3::Vector& xyz :
       {R3::Vector(0.3, 0.0, 0.0), R3::Vector(0.9, 0.6, 1.3),
        R3::Vector(1.7, 1.2, 0.5)}) {
    a.xyz_cartn = xyz;
    stru->append(a);
  }
  PairTable table;
  auto check = [&] {
    checkSamePairs(generatedPairs(*stru, 0.0, 7.0),
                   tablePairs(table, *stru, 0.0, 7.0));
  };
  check();
  CHECK(2 == stru->siteMultiplicity(0));
  stru->at(0).xyz_cartn[0] = 0.0;
  check();
  CHECK(1 == stru->siteMultiplicity(0));
  stru->at(0).xyz_cartn[0] = 0.01;
  check();
  stru->addSymOp(R3::identity(), R3::Vector(0.5, 0.5, 0.5));
  check();
  stru->setSymmetryPrecision(0.5);
  check();
}

#endif

TEST_CASE("TestPairTable::test_results_of_calculators", "[PairTable]") {
  testAllKernels([] {
    forEachCalculatorCase(
      [](const Setup& setup, StructureAdapterPtr structure) {
        checkAllPairs<PDFCalculator>(setup, structure);
        checkAllPairs<DebyePDFCalculator>(setup, structure);
      });
  });
}

TEST_CASE("TestPairTable::test_peak_profiles", "[PairTable]") {
  // Any peak profile is linear in the scale of its peaks.
  struct WideGaussian : GaussianProfile {
    double operator()(double x, double fwhm) const override {
      return this->GaussianProfile::operator()(x, 1.5 * fwhm);
    }
  };
  const auto periodic = loadPeriodic("CaTiO3.stru");
  const Setup cropped = [](PairQuantity& c) {
    setupDefault(c);
    dynamic_cast<PDFCalculator&>(c).setPeakProfileByType("croppedgaussian");
  };
  const Setup wide = [](PairQuantity& c) {
    setupDefault(c);
    PeakProfilePtr pkf(new WideGaussian);
    pkf->setPrecision(1e-5);
    dynamic_cast<PDFCalculator&>(c).setPeakProfile(pkf);
  };
  testAllKernels([&] {
    for (StructureAdapterPtr stru :
         {StructureAdapterPtr(periodic),
          StructureAdapterPtr(asCrystal(*periodic)),
          StructureAdapterPtr(asCluster(*periodic))}) {
      checkAllPairs<PDFCalculator>(cropped, stru);
      checkPairByPair<PDFCalculator>(wide, stru);
    }
  });
}

TEST_CASE("TestPairTable::test_pairs_added_one_by_one", "[PairTable]") {
  // A derived class may add its pairs in another way, as these that add
  // every pair twice.
  struct TwicePDFCalculator : PDFCalculator {
    void addPairContribution(const BaseBondGenerator& bnds,
                             int summationscale) override {
      this->PDFCalculator::addPairContribution(bnds, 2 * summationscale);
    }
  };
  struct TwiceDebyePDFCalculator : DebyePDFCalculator {
    void addPairContribution(const BaseBondGenerator& bnds,
                             int summationscale) override {
      this->DebyePDFCalculator::addPairContribution(bnds, 2 * summationscale);
    }
  };
  const StructureAdapterPtr stru = loadPeriodic("CaTiO3.stru");
  auto twice = [](QuantityType v) {
    for (double& x : v) x *= 2;
    return v;
  };
  const auto once = resultsOf(PDFCalculator(), setupDefault, stru);
  const auto doubled = resultsOf(TwicePDFCalculator(), setupDefault, stru);
  CHECK(maxRelativeError(doubled[0], twice(once[0])) <= TOLERANCE);
  const auto donce = resultsOf(DebyePDFCalculator(), setupDefault, stru);
  const auto ddoubled =
    resultsOf(TwiceDebyePDFCalculator(), setupDefault, stru);
  CHECK(maxRelativeError(ddoubled[0], twice(donce[0])) <= TOLERANCE);
  // The other peak width models need the bond of each pair.
  struct DerivedWidth : JeongPeakWidth {};
  const Setup derivedwidth = [](PairQuantity& c) {
    setupDefault(c);
    dynamic_cast<PeakWidthModelOwner&>(c).setPeakWidthModel(
      PeakWidthModelPtr(new DerivedWidth));
  };
  checkPairByPair<PDFCalculator>(derivedwidth, stru);
  checkPairByPair<DebyePDFCalculator>(derivedwidth, stru);
  // One part of a parallel run has a share of the pairs.
  const Setup parallel = [](PairQuantity& c) {
    setupDefault(c);
    c.setupParallelRun(0, 2);
  };
  checkPairByPair<PDFCalculator>(parallel, stru);
}

TEST_CASE("TestPairTable::test_masks", "[PairTable]") {
  // The masks leave out the same pairs as for the pairs added one by one.
  const std::pair<const char*, Setup> setups[] = {
    {"without the pairs of two types",
     [](PairQuantity& c) {
       setupDefault(c);
       c.setTypeMask("Ti4+", "O2-", false);
     }},
    {"only the pairs of one type with any other",
     [](PairQuantity& c) {
       setupDefault(c);
       c.maskAllPairs(false);
       c.setTypeMask("Ca2+", PairQuantity::ALLATOMSSTR, true);
     }},
    {"without some pairs of sites",
     [](PairQuantity& c) {
       setupDefault(c);
       c.setPairMask(0, 5, false);
       c.setPairMask(2, 2, false);
       c.setPairMask(9, 1, false);
     }},
    {"without one site",
     [](PairQuantity& c) {
       setupDefault(c);
       c.setPairMask(4, PairQuantity::ALLATOMSINT, false);
     }},
    {"only one site with itself",
     [](PairQuantity& c) {
       setupDefault(c);
       c.maskAllPairs(false);
       c.setPairMask(5, 5, true);
     }},
  };
  const auto periodic = loadPeriodic("CaTiO3.stru");
  const std::pair<const char*, StructureAdapterPtr> structures[] = {
    {"periodic", periodic}, {"crystal", asCrystal(*periodic)}};
  for (const auto& [kind, stru] : structures) {
    AllPairs<PDFCalculator> pdfc;
    AllPairs<DebyePDFCalculator> dpdfc;
    setupDefault(pdfc);
    setupDefault(dpdfc);
    const QuantityType vall = pdfc.eval(stru);
    const QuantityType dvall = dpdfc.eval(stru);
    for (const auto& [name, setup] : setups) {
      INFO(kind << " " << name);
      checkAllPairs<PDFCalculator>(setup, stru);
      checkAllPairs<DebyePDFCalculator>(setup, stru);
      // the mask leaves out some pairs
      setup(pdfc);
      setup(dpdfc);
      CHECK(maxRelativeError(pdfc.eval(stru), vall) > 1e-3);
      CHECK(maxRelativeError(dpdfc.eval(stru), dvall) > 1e-3);
      pdfc.maskAllPairs(true);
      dpdfc.maskAllPairs(true);
    }
    // A mask of all pairs changes nothing in the sum.
    const int cntsites = stru->countSites();
    for (PairQuantity* c : {(PairQuantity*)&pdfc, (PairQuantity*)&dpdfc}) {
      c->maskAllPairs(false);
      CHECK(maxAbsolute(c->eval(stru)) == 0.0);
      for (int i = 0; i < cntsites; ++i) {
        for (int j = 0; j <= i; ++j) c->setPairMask(i, j, true);
      }
    }
    CHECK(vall == pdfc.eval(stru));
    CHECK(dvall == dpdfc.eval(stru));
    CHECK(pdfc.addedAllPairs());
    CHECK(dpdfc.addedAllPairs());
  }
}

TEST_CASE("TestPairTable::test_evaluators", "[PairTable]") {
  const StructureAdapterPtr stru = loadPeriodic("CaTiO3.stru");
  const QuantityType expected =
    resultsOf(PairByPair<PDFCalculator>(), setupDefault, stru)[0];
  for (PQEvaluatorType tp : {BASIC, OPTIMIZED, PQEvaluatorType::CHECK}) {
    INFO("evaluator " << tp);
    AllPairs<PDFCalculator> pdfc;
    setupDefault(pdfc);
    pdfc.setEvaluatorType(tp);
    CHECK(!pdfc.addedAllPairs());
    pdfc.eval(stru);
    CHECK(pdfc.addedAllPairs());
    CHECK(maxRelativeError(pdfc.value(), expected) <= TOLERANCE);
    const QuantityType v = pdfc.value();
    pdfc.eval(stru);
    if (tp != BASIC) CHECK(tp == pdfc.getEvaluatorTypeUsed());
    CHECK(v == pdfc.value());
  }
}

TEST_CASE("TestPairTable::test_fast_updates", "[PairTable]") {
  // The OPTIMIZED evaluator updates the sum of all pairs for the sites
  // that change, pair by pair.
  for (const auto& [name, stru0] : allStructures()) {
    INFO(name);
    PairByPair<PDFCalculator> basic;
    PairByPair<DebyePDFCalculator> dbasic;
    AllPairs<PDFCalculator> optimized;
    AllPairs<DebyePDFCalculator> doptimized;
    for (PairQuantity* c : {(PairQuantity*)&basic, (PairQuantity*)&dbasic}) {
      c->setEvaluatorType(BASIC);
    }
    for (PairQuantity* c :
         {(PairQuantity*)&basic, (PairQuantity*)&optimized,
          (PairQuantity*)&dbasic, (PairQuantity*)&doptimized}) {
      setupDefault(*c);
    }
    StructureAdapterPtr stru = stru0->clone();
    for (int step = 0; step < 4; ++step) {
      INFO("step " << step);
      if (step > 0) {
        // one changed site of the source structure
        StructureAdapter* s = stru.get();
        while (!dynamic_cast<AtomicStructureAdapter*>(s)) {
          if (auto* nm = dynamic_cast<NoMetaStructureAdapter*>(s)) {
            s = nm->getSourceStructure().get();
          } else {
            s = dynamic_cast<NoSymmetryStructureAdapter&>(*s)
                  .getSourceStructure()
                  .get();
          }
        }
        auto& atoms = dynamic_cast<AtomicStructureAdapter&>(*s);
        atoms.at(step % atoms.countSites()).xyz_cartn[step % 3] += 0.07;
      }
      basic.eval(stru);
      optimized.eval(stru);
      dbasic.eval(stru);
      doptimized.eval(stru);
      CHECK(maxRelativeError(optimized.getPDF(), basic.getPDF()) <= 1e-10);
      CHECK(maxRelativeError(doptimized.getPDF(), dbasic.getPDF()) <= 1e-10);
      // the first calculation is of all pairs
      if (!step) CHECK(optimized.addedAllPairs());
      if (!step) CHECK(doptimized.addedAllPairs());
    }
  }
}

TEST_CASE("TestPairTable::test_fast_update_or_all_pairs", "[PairTable]") {
  // The OPTIMIZED evaluator adds all pairs again when that takes less work
  // than the pairs of the changed sites one by one.
  struct Case {
    const char* structure;
    /// sites that move
    std::vector<int> sites;
    /// evaluator for these sites when the calculator can add all pairs
    PQEvaluatorType used;
  };
  const Case cases[] = {
    // Four sites, the pairs merge to few peaks.
    {"Ni.stru", {1}, BASIC},
    // Twenty sites: one of them takes less work by itself, the four of
    // one element do not.
    {"CaTiO3.stru", {5}, OPTIMIZED},
    {"CaTiO3.stru", {4, 5, 6, 7}, BASIC},
    // No pairs merge without symmetry.
    {"alpha_K2Bi8Se13.stru", {3}, OPTIMIZED},
    {"alpha_K2Bi8Se13.stru", {2, 3, 4, 5}, OPTIMIZED},
  };
  for (const Case& c : cases) {
    INFO(c.structure << " with " << c.sites.size() << " moved sites");
    const auto stru = loadPeriodic(c.structure);
    // calculators that add all pairs at once, then ones that do not
    auto check = [&](PairQuantity& basic, PairQuantity& optimized,
                     PairQuantity& dbasic, PairQuantity& doptimized,
                     bool pairbypair) {
      INFO((pairbypair ? "pair by pair" : "all pairs"));
      for (PairQuantity* calc : {&basic, &optimized, &dbasic, &doptimized}) {
        setupDefault(*calc);
      }
      basic.setEvaluatorType(BASIC);
      dbasic.setEvaluatorType(BASIC);
      optimized.eval(stru);
      doptimized.eval(stru);
      const PQEvaluatorType used = pairbypair ? OPTIMIZED : c.used;
      for (int step = 0; step < 3; ++step) {
        for (int i : c.sites) stru->at(i).xyz_cartn[step] += 0.02;
        StructureAdapterPtr copy = stru->clone();
        CHECK(maxRelativeError(optimized.eval(copy), basic.eval(copy)) <=
              1e-10);
        CHECK(maxRelativeError(doptimized.eval(copy), dbasic.eval(copy)) <=
              1e-10);
        CHECK(used == optimized.getEvaluatorTypeUsed());
        CHECK(used == doptimized.getEvaluatorTypeUsed());
      }
      // a structure of the same content takes no work in either way
      optimized.eval(stru->clone());
      CHECK(OPTIMIZED == optimized.getEvaluatorTypeUsed());
      // a removed site has half the work of a moved one
      stru->erase(c.sites[0]);
      CHECK(maxRelativeError(optimized.eval(stru), basic.eval(stru)) <= 1e-10);
      CHECK(OPTIMIZED == optimized.getEvaluatorTypeUsed());
      stru->insert(c.sites[0], Atom(stru->at(0)));
    };
    PDFCalculator basic, optimized;
    DebyePDFCalculator dbasic, doptimized;
    check(basic, optimized, dbasic, doptimized, false);
    PairByPair<PDFCalculator> pbasic, poptimized;
    PairByPair<DebyePDFCalculator> pdbasic, pdoptimized;
    check(pbasic, poptimized, pdbasic, pdoptimized, true);
  }
}

TEST_CASE("TestPairTable::test_copies_and_archives", "[PairTable]") {
  const StructureAdapterPtr stru = loadPeriodic("CaTiO3.stru");
  AllPairs<PDFCalculator> pdfc;
  setupDefault(pdfc);
  pdfc.eval(stru);
  const QuantityType v0 = pdfc.value();
  // a copy has pairs of its own
  AllPairs<PDFCalculator> copy = pdfc;
  copy.setEvaluatorType(BASIC);
  copy.eval(nosymmetry(stru));
  CHECK(copy.addedAllPairs());
  CHECK(maxRelativeError(copy.value(), v0) > 0.01);
  pdfc.setEvaluatorType(BASIC);
  pdfc.eval(stru);
  CHECK(maxRelativeError(pdfc.value(), v0) <= TOLERANCE);
  AllPairs<PDFCalculator> assigned;
  assigned = pdfc;
  assigned.eval(stru);
  CHECK(assigned.addedAllPairs());
  CHECK(maxRelativeError(assigned.value(), v0) <= TOLERANCE);
  // and so has a calculator from an archive
  auto saved = std::make_shared<PDFCalculator>();
  setupDefault(*saved);
  saved->eval(stru);
  std::stringstream storage(std::ios::in | std::ios::out | std::ios::binary);
  diffpy::serialization::oarchive oa(storage, std::ios::binary);
  oa << saved;
  diffpy::serialization::iarchive ia(storage, std::ios::binary);
  std::shared_ptr<PDFCalculator> loaded;
  ia >> loaded;
  CHECK(loaded->value() == v0);
  loaded->setEvaluatorType(BASIC);
  loaded->eval(stru);
  CHECK(loaded->value() == v0);
}

TEST_CASE("TestPairTable::test_many_evaluations", "[PairTable]") {
  // A refinement of a crystal: positions, displacements and the lattice
  // change a little at every evaluation.
  const auto periodic = loadPeriodic("CaTiO3.stru");
  for (std::shared_ptr<PeriodicStructureAdapter> stru :
       {periodic,
        std::shared_ptr<PeriodicStructureAdapter>(asCrystal(*periodic))}) {
    PairByPair<PDFCalculator> onebyone;
    PairByPair<DebyePDFCalculator> donebyone;
    AllPairs<PDFCalculator> atonce;
    AllPairs<DebyePDFCalculator> datonce;
    for (PairQuantity* c :
         {(PairQuantity*)&onebyone, (PairQuantity*)&atonce,
          (PairQuantity*)&donebyone, (PairQuantity*)&datonce}) {
      setupDefault(*c);
      c->setEvaluatorType(BASIC);
    }
    const Lattice L = stru->getLattice();
    for (int k = 0; k < 40; ++k) {
      INFO("evaluation " << k);
      const int n = stru->countSites();
      stru->at(k % n).xyz_cartn[k % 3] += 0.004 * ((k % 5) - 2);
      stru->at((3 * k) % n).uij_cartn *= 1.0 + 0.01 * ((k % 3) - 1);
      for (Atom& a : *stru) stru->toFractional(a);
      const double f = 1.0 + 1e-4 * ((k % 7) - 3);
      stru->setLatPar(f * L.a(), L.b(), L.c() / f, L.alpha(), L.beta(),
                      L.gamma());
      for (Atom& a : *stru) stru->toCartesian(a);
      onebyone.setDoubleAttr("delta2", 0.5 + 0.05 * (k % 4));
      atonce.setDoubleAttr("delta2", 0.5 + 0.05 * (k % 4));
      StructureAdapterPtr copy = stru->clone();
      onebyone.eval(copy);
      atonce.eval(copy);
      donebyone.eval(copy);
      datonce.eval(copy);
      CHECK(atonce.addedAllPairs());
      CHECK(datonce.addedAllPairs());
      CHECK(maxRelativeError(atonce.value(), onebyone.value()) <= TOLERANCE);
      CHECK(maxRelativeError(datonce.value(), donebyone.value()) <= TOLERANCE);
    }
  }
}
