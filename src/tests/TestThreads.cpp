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
 * TestThreads -- calculators that evaluate in several threads at once
 *
 * Catch2 assertions are not thread-safe, the threads only collect results
 * and the main thread checks them.
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <latch>
#include <memory>
#include <thread>
#include <vector>

#include <diffpy/EventTicker.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>

#include "test_helpers.hpp"

using namespace diffpy::srreal;

namespace {

const int NTHREADS = 4;

/// Results of a short refinement-like sequence of evaluations.  The job
/// number selects the calculator and makes the sequences differ.
std::vector<QuantityType> evaluateSequence(int job, StructureAdapterPtr stru0) {
  PDFCalculator pdfc;
  DebyePDFCalculator dpdfc;
  const bool usepdf = (job % 2 == 0);
  PairQuantity& calc = usepdf ? static_cast<PairQuantity&>(pdfc)
                              : static_cast<PairQuantity&>(dpdfc);
  calc.setDoubleAttr("rmax", 8.0 + job);
  calc.setDoubleAttr("qmax", 20.0);
  std::vector<QuantityType> rv;
  for (int k = 0; k < 6; ++k) {
    // The first evaluation reads stru0 itself, which other threads may
    // read too.  The later ones get a new structure with a moved atom and
    // changed calculator attributes.
    StructureAdapterPtr stru = stru0;
    if (k > 0) {
      stru = stru0->clone();
      auto& atoms = dynamic_cast<AtomicStructureAdapter&>(*stru);
      atoms[0].xyz_cartn[0] += 0.01 * k * (job + 1);
    }
    if (k == 3) calc.setDoubleAttr("delta2", 1.0 + 0.1 * job);
    if (k == 4) calc.setDoubleAttr("scale", 1.5);
    if (k == 5) calc.setDoubleAttr("qmax", 17.0 + job);
    calc.eval(stru);
    rv.push_back(usepdf ? pdfc.getPDF() : dpdfc.getPDF());
    rv.push_back(usepdf ? pdfc.getF() : dpdfc.getF());
  }
  return rv;
}

/// Cubic crystal with a face-centered Ni site and a body-centered Cu site.
std::shared_ptr<CrystalStructureAdapter> makeCrystal() {
  auto stru = std::make_shared<CrystalStructureAdapter>();
  stru->setLatPar(3.52, 3.52, 3.52, 90, 90, 90);
  const double centering[4][3] = {
    {0.0, 0.0, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
  for (const auto& t : centering) {
    stru->addSymOp(R3::identity(), R3::Vector(t[0], t[1], t[2]));
  }
  Atom a;
  a.atomtype = "Ni";
  a.uij_cartn(0, 0) = a.uij_cartn(1, 1) = a.uij_cartn(2, 2) = 0.004;
  stru->append(a);
  a.atomtype = "Cu";
  a.occupancy = 0.5;
  a.xyz_cartn = R3::Vector(0.5, 0.5, 0.5);
  stru->toCartesian(a);
  stru->append(a);
  return stru;
}

/// Call fn(job) in NTHREADS threads that start together.
template <class Fn>
void runInThreads(Fn fn) {
  std::latch start(NTHREADS);
  std::vector<std::jthread> threads;
  try {
    for (int job = 0; job < NTHREADS; ++job) {
      threads.emplace_back([&, job] {
        start.arrive_and_wait();
        fn(job);
      });
    }
  } catch (...) {
    start.count_down(NTHREADS - threads.size());
    throw;
  }
}

}  // namespace

TEST_CASE("TestThreads::test_separate_calculators", "[Threads]") {
  const StructureAdapterPtr stru = loadTestPeriodicStructure("CaTiO3.stru");
  std::vector<std::vector<QuantityType>> serial(NTHREADS);
  for (int job = 0; job < NTHREADS; ++job) {
    serial[job] = evaluateSequence(job, stru->clone());
  }
  // Every thread has its own calculator and its own structure.
  for (int repeat = 0; repeat < 5; ++repeat) {
    std::vector<std::vector<QuantityType>> threaded(NTHREADS);
    std::vector<StructureAdapterPtr> strus;
    for (int job = 0; job < NTHREADS; ++job) strus.push_back(stru->clone());
    runInThreads(
      [&](int job) { threaded[job] = evaluateSequence(job, strus[job]); });
    CHECK(serial == threaded);
  }
}

TEST_CASE("TestThreads::test_shared_structure", "[Threads]") {
  std::vector<std::vector<QuantityType>> serial(NTHREADS);
  for (int job = 0; job < NTHREADS; ++job) {
    serial[job] = evaluateSequence(job, makeCrystal());
  }
  // The threads read one crystal, whose symmetry positions are expanded on
  // the first use by any of them.
  for (int repeat = 0; repeat < 5; ++repeat) {
    std::vector<std::vector<QuantityType>> threaded(NTHREADS);
    const StructureAdapterPtr shared = makeCrystal();
    runInThreads(
      [&](int job) { threaded[job] = evaluateSequence(job, shared); });
    CHECK(serial == threaded);
  }
}

TEST_CASE("TestThreads::test_event_ticker", "[Threads]") {
  using diffpy::eventticker::EventTicker;
  const int nclicks = 20000;
  EventTicker before, after;
  std::vector<int> unordered(NTHREADS, 0);
  before.click();
  runInThreads([&](int job) {
    EventTicker last, current;
    for (int i = 0; i < nclicks; ++i) {
      current.click();
      if (!(last < current)) ++unordered[job];
      last = current;
    }
  });
  after.click();
  // No click was lost and every click gave a newer value.
  CHECK(before.value().first == after.value().first);
  CHECK(long(NTHREADS) * nclicks + 1 ==
        after.value().second - before.value().second);
  CHECK(std::vector<int>(NTHREADS, 0) == unordered);
}
