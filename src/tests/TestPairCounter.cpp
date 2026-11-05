/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2009 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestPairCounter -- unit tests for PairCounter class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PairCounter.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestPairCounterFixture {
 protected:
  AtomicStructureAdapterPtr mstru;
  AtomicStructureAdapterPtr mline100;

 public:
  TestPairCounterFixture() {
    mstru.reset(new AtomicStructureAdapter);
    mline100.reset(new AtomicStructureAdapter);
    Atom a;
    for (int i = 0; i < 100; ++i) {
      a.xyz_cartn = R3::Vector(1.0 * i, 0.0, 0.0);
      mline100->append(a);
    }
  }
};  // class TestPairCounterFixture

TEST_CASE_METHOD(TestPairCounterFixture, "TestPairCounter::test_call",
                 "[PairCounter]") {
  PairCounter pcount;
  Atom a;
  for (int i = 0; i < 100; ++i) {
    int npairs = i * (i - 1) / 2;
    CHECK(npairs == pcount(mstru));
    a.xyz_cartn = R3::Vector(1.0 * i, 0.0, 0.0);
    mstru->append(a);
  }
  mstru->clear();
  CHECK(0 == pcount(mstru));
}

TEST_CASE_METHOD(TestPairCounterFixture, "TestPairCounter::test_setRmin",
                 "[PairCounter]") {
  PairCounter pcount;
  CHECK(100 * 99 / 2 == pcount(mline100));
  pcount.setRmin(100);
  CHECK(0 == pcount(mline100));
  pcount.setRmin(99.0);
  CHECK(1 == pcount(mline100));
  pcount.setRmin(1.1);
  CHECK(100 * 99 / 2 - 99 == pcount(mline100));
}

TEST_CASE_METHOD(TestPairCounterFixture, "TestPairCounter::test_setRmax",
                 "[PairCounter]") {
  PairCounter pcount;
  pcount.setRmax(0.9);
  CHECK(0 == pcount(mline100));
  pcount.setRmax(1.1);
  CHECK(99 == pcount(mline100));
  pcount.setRmax(98.5);
  CHECK(100 * 99 / 2 - 1 == pcount(mline100));
}

TEST_CASE_METHOD(TestPairCounterFixture, "TestPairCounter::test_maskTicker",
                 "[PairCounter]") {
  using diffpy::eventticker::EventTicker;
  PairCounter pcount;
  EventTicker et0 = pcount.ticker();
  pcount.maskAllPairs(true);
  pcount.setPairMask(pcount.ALLATOMSINT, pcount.ALLATOMSINT, true);
  pcount.setTypeMask("all", "all", true);
  CHECK(et0 == pcount.ticker());
  pcount.setTypeMask("A", "B", true);
  EventTicker et1 = pcount.ticker();
  CHECK(et0 < et1);
  pcount.setTypeMask("A", "B", true);
  CHECK(et1 == pcount.ticker());
  pcount.setTypeMask("A", "B", false);
  EventTicker et2 = pcount.ticker();
  CHECK(et1 < et2);
  pcount.setTypeMask("A", "all", false);
  EventTicker et3 = pcount.ticker();
  CHECK(et2 < et3);
  pcount.setTypeMask("A", "all", false);
  CHECK(et3 == pcount.ticker());
  pcount.setTypeMask("A", "all", true);
  EventTicker et4 = pcount.ticker();
  CHECK(et3 < et4);
  pcount.setPairMask(0, 0, false);
  EventTicker et5 = pcount.ticker();
  CHECK(et4 < et5);
  pcount.setPairMask(0, 0, false);
  CHECK(et5 == pcount.ticker());
  pcount.setPairMask(0, 5, true);
  CHECK(et5 == pcount.ticker());
  pcount.setPairMask(0, 5, false);
  CHECK(et5 < pcount.ticker());
}

TEST_CASE_METHOD(TestPairCounterFixture, "TestPairCounter::test_parallel",
                 "[PairCounter]") {
  const int ncpu = 7;
  PairCounter pmaster;
  PairCounter pslave[ncpu];
  bool allequal = true;
  for (int cpuindex = 0; cpuindex < ncpu; ++cpuindex) {
    pslave[cpuindex].setupParallelRun(cpuindex, ncpu);
    pslave[cpuindex].eval(mline100);
    allequal =
      allequal && (pslave[cpuindex].value()[0] == pslave[0].value()[0]);
  }
  CHECK(!allequal);
  pmaster.setStructure(mline100);
  CHECK(1u == pmaster.value().size());
  CHECK(0.0 == pmaster.value()[0]);
  for (PairCounter* p = pslave; p != pslave + ncpu; ++p) {
    pmaster.mergeParallelData(p->getParallelData(), ncpu);
  }
  CHECK(100 * 99 / 2 == pmaster.value()[0]);
}

// End of file
