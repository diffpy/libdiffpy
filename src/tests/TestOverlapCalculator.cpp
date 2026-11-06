/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2011 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestOverlapCalculator -- unit tests for the OverlapCalculator class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <functional>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/OverlapCalculator.hpp>
#include <diffpy/serialization.ipp>
#include "test_helpers.hpp"
#include "serialization_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;

class TestOverlapCalculatorFixture {
 protected:
  std::shared_ptr<OverlapCalculator> molc;
  StructureAdapterPtr mnacl;
  double meps;

 public:
  TestOverlapCalculatorFixture() {
    meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
    molc.reset(new OverlapCalculator);
    molc->getAtomRadiiTable()->setCustom("", 1.0);
    molc->getAtomRadiiTable()->setCustom("Na1+", 1.5);
    molc->getAtomRadiiTable()->setCustom("Cl1-", 1.8);
    mnacl = loadTestPeriodicStructure("NaCl.stru");
  }
};  // class TestOverlapCalculatorFixture

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_rmin", "[OverlapCalculator]") {
  molc->eval(mnacl);
  CHECK(molc->distances().size());
  double tot = molc->totalSquareOverlap();
  CHECK(0.0 < tot);
  molc->setRmin(10);
  molc->eval();
  CHECK(molc->distances().empty());
  CHECK(0.0 == molc->totalSquareOverlap());
  molc->setRmin(2);
  molc->eval();
  CHECK(tot == molc->totalSquareOverlap());
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_rmax", "[OverlapCalculator]") {
  molc->setRmax(2.5);
  molc->eval(mnacl);
  CHECK(molc->distances().empty());
  CHECK(0.0 == molc->totalSquareOverlap());
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_lineNoTouch",
                 "[OverlapCalculator]") {
  AtomicStructureAdapterPtr stru(new AtomicStructureAdapter);
  Atom a0, a1;
  a0.xyz_cartn = R3::Vector(0.0, 0.0, 0.0);
  a1.xyz_cartn = R3::Vector(0.0, 0.0, 2.0);
  stru->append(a0);
  stru->append(a1);
  molc->eval(stru);
  CHECK(0.0 == molc->totalSquareOverlap());
  QuantityType sqolps = molc->siteSquareOverlaps();
  REQUIRE(2u == sqolps.size());
  CHECK(0.0 == sqolps[0]);
  CHECK(0.0 == sqolps[1]);
  CHECK(0.0 == molc->flipDiffTotal(0, 0));
  CHECK(0.0 == molc->flipDiffTotal(0, 1));
  std::vector<R3::Vector> g = molc->gradients();
  REQUIRE(2u == g.size());
  CHECK(0.0 == R3::norm(g[0]));
  CHECK(0.0 == R3::norm(g[1]));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_lineTouch",
                 "[OverlapCalculator]") {
  AtomicStructureAdapterPtr stru(new AtomicStructureAdapter);
  Atom a0, a1;
  a0.xyz_cartn = R3::Vector(0.0, 0.0, 0.0);
  a1.xyz_cartn = R3::Vector(0.0, 0.0, 1.5);
  stru->append(a0);
  stru->append(a1);
  molc->eval(stru);
  CHECK(0.25 == molc->totalSquareOverlap());
  QuantityType sqolps = molc->siteSquareOverlaps();
  REQUIRE(2u == sqolps.size());
  CHECK(0.125 == sqolps[0]);
  CHECK(0.125 == sqolps[1]);
  CHECK(0.0 == molc->flipDiffTotal(0, 1));
  std::vector<R3::Vector> g = molc->gradients();
  REQUIRE(2u == g.size());
  R3::Vector g0(0.0, 0.0, +1);
  CHECK(0.0 == R3::distance(g0, g[0]));
  R3::Vector g1(0.0, 0.0, -1);
  CHECK(0.0 == R3::distance(g1, g[1]));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_gradients",
                 "[OverlapCalculator]") {
  AtomicStructureAdapterPtr stru(new AtomicStructureAdapter);
  Atom a0, a1;
  a0.xyz_cartn = R3::Vector(0.0, 0.0, 0.0);
  a1.xyz_cartn = R3::Vector(0.17, 0.31, 0.37);
  stru->append(a0);
  stru->append(a1);
  molc->eval(stru);
  std::vector<R3::Vector> g = molc->gradients();
  REQUIRE(2u == g.size());
  double dx = 1e-6;
  double tsqo0 = molc->totalSquareOverlap();
  R3::Vector r1 = stru->siteCartesianPosition(1);
  R3::Vector gn1;
  for (int k = 0; k != R3::Ndim; ++k) {
    R3::Vector& xyz1 = stru->at(1).xyz_cartn;
    xyz1 = r1;
    xyz1[k] += dx;
    molc->eval(stru);
    gn1[k] = (molc->totalSquareOverlap() - tsqo0) / dx;
  }
  CHECK(R3::norm(gn1) > 1.0);
  CHECK_THAT(0.0, WithinAbs(R3::norm(gn1 - g[1]), 1e-5));
  CHECK_THAT(0.0, WithinAbs(R3::norm(g[0] + g[1]), 1e-5));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_bccTouch",
                 "[OverlapCalculator]") {
  PeriodicStructureAdapterPtr bcc(new PeriodicStructureAdapter);
  bcc->setLatPar(2.0, 2.0, 2.0, 90, 90, 90);
  Atom a;
  a.atomtype = "C";
  bcc->append(a);
  a.xyz_cartn = R3::Vector(0.5, 0.5, 0.5);
  bcc->toCartesian(a);
  bcc->append(a);
  molc->getAtomRadiiTable()->setCustom("C", 1.0);
  molc->eval(bcc);
  CHECK(2 == molc->getStructure()->countSites());
  CHECK(8 * pow(2.0 - sqrt(3.0), 2) == molc->totalSquareOverlap());
  // gradient is zero due to bcc symmetry
  std::vector<R3::Vector> g = molc->gradients();
  REQUIRE(2u == g.size());
  CHECK_THAT(0.0, WithinAbs(R3::norm(g[0]), meps));
  CHECK_THAT(0.0, WithinAbs(R3::norm(g[1]), meps));
  // move the second atom so it touches only one ball
  bcc->at(1).xyz_cartn = R3::Vector(-0.25, 0.0, 0.0);
  bcc->toCartesian(bcc->at(1));
  molc->eval(bcc);
  double tsqo = 1.5 * 1.5 + 0.5 * 0.5;
  CHECK_THAT(tsqo, WithinAbs(molc->totalSquareOverlap(), meps));
  g = molc->gradients();
  REQUIRE(2u == g.size());
  double gx = 1.5 + 0.5;
  CHECK_THAT(-gx, WithinAbs(g[0][0], meps));
  CHECK_THAT(0.0, WithinAbs(g[0][1], meps));
  CHECK_THAT(0.0, WithinAbs(g[0][2], meps));
  CHECK_THAT(+gx, WithinAbs(g[1][0], meps));
  CHECK_THAT(0.0, WithinAbs(g[1][1], meps));
  CHECK_THAT(0.0, WithinAbs(g[1][2], meps));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_NaCl_overlap",
                 "[OverlapCalculator]") {
  double olp = pow(3.3 - 5.62 / 2, 2) * 6 / 2;
  molc->eval(mnacl);
  CHECK_THAT(olp, WithinAbs(molc->meanSquareOverlap(), meps));
  QuantityType sqolps = molc->siteSquareOverlaps();
  REQUIRE(8u == sqolps.size());
  for (int i = 0; i < 8; ++i) {
    CHECK_THAT(olp, WithinAbs(sqolps[i], meps));
  }
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_NaCl_flips",
                 "[OverlapCalculator]") {
  molc->eval(mnacl);
  // flipping the same type should not change the cost
  CHECK(0.0 == molc->flipDiffTotal(0, 0));
  CHECK(0.0 == molc->flipDiffTotal(0, 1));
  CHECK(0.0 == molc->flipDiffTotal(0, 2));
  CHECK(0.0 == molc->flipDiffTotal(0, 3));
  // flipping with the second Cl neighbor
  CHECK_THAT(1.08, WithinAbs(molc->flipDiffTotal(0, 4), meps));
  CHECK_THAT(1.08 / 8, WithinAbs(molc->flipDiffMean(0, 4), meps));
  // flipping with the nearest Cl neighbor
  CHECK_THAT(0.72, WithinAbs(molc->flipDiffTotal(0, 5), meps));
  CHECK_THAT(0.72, WithinAbs(molc->flipDiffTotal(0, 6), meps));
  CHECK_THAT(0.72, WithinAbs(molc->flipDiffTotal(0, 7), meps));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_NaCl_gradient",
                 "[OverlapCalculator]") {
  molc->eval(mnacl);
  // default gradients are all zero
  std::vector<R3::Vector> g = molc->gradients();
  REQUIRE(8u == g.size());
  for (int i = 0; i < 8; ++i) {
    CHECK_THAT(0.0, WithinAbs(R3::norm(g[i]), meps));
  }
  StructureAdapterPtr nacl1 = mnacl->clone();
  PeriodicStructureAdapter& nacl1ref =
    static_cast<PeriodicStructureAdapter&>(*nacl1);
  R3::Vector& xyzc0 = nacl1ref[0].xyz_cartn;
  xyzc0 = R3::Vector(0.02, 0.03, 0.07);
  molc->eval(nacl1);
  double c0 = molc->totalSquareOverlap();
  REQUIRE(8u == molc->gradients().size());
  R3::Vector g0 = molc->gradients()[0];
  CHECK(R3::norm(g0) > meps);
  R3::Vector g0n;
  const double dx = 1e-8;
  for (int i = 0; i < R3::Ndim; ++i) {
    xyzc0[i] += dx;
    molc->eval(nacl1);
    g0n[i] = (molc->totalSquareOverlap() - c0) / dx;
    xyzc0[i] -= dx;
  }
  CHECK_THAT(0.0, WithinAbs(R3::norm(g0 - g0n), 1e-6));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_getNeighborSites",
                 "[OverlapCalculator]") {
  molc->eval(mnacl);
  auto nbsites = molc->getNeighborSites(0);
  REQUIRE(3u == nbsites.size());
  for (int j = 5; j < 8; ++j) {
    CHECK(nbsites.count(j));
  }
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_coordinations",
                 "[OverlapCalculator]") {
  using std::placeholders::_1;
  molc->eval(mnacl);
  auto cnums = molc->coordinations();
  REQUIRE(8u == cnums.size());
  auto eq6 = bind(equal_to<double>(), _1, 6.0);
  CHECK(all_of(cnums.begin(), cnums.end(), eq6));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_coordinationByTypes",
                 "[OverlapCalculator]") {
  molc->eval(mnacl);
  auto cbtb0 = molc->coordinationByTypes(0);
  auto cbtb3 = molc->coordinationByTypes(3);
  auto cbtb5 = molc->coordinationByTypes(5);
  REQUIRE(1u == cbtb0.size());
  REQUIRE(1u == cbtb5.size());
  CHECK(cbtb0 == cbtb3);
  CHECK(6 == cbtb0.at("Cl1-"));
  CHECK(6 == cbtb5.at("Na1+"));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_neighborhoods",
                 "[OverlapCalculator]") {
  auto tb = molc->getAtomRadiiTable();
  molc->eval(mnacl);
  auto nbhood = molc->neighborhoods();
  REQUIRE(1u == nbhood.size());
  REQUIRE(8u == nbhood[0].size());
  tb->resetAll();
  molc->eval(mnacl);
  auto nbsep = molc->neighborhoods();
  REQUIRE(8u == nbsep.size());
  REQUIRE_FALSE(nbsep[5].empty());
  int nb5site = *(nbsep[5].begin());
  CHECK(5 == nb5site);
  molc->maskAllPairs(false);
  molc->eval(mnacl);
  auto nbdark = molc->neighborhoods();
  CHECK(nbdark.empty());
  molc->setTypeMask("Na1+", "Na1+", true);
  molc->eval(mnacl);
  auto nbnanasep = molc->neighborhoods();
  REQUIRE(4u == nbnanasep.size());
  REQUIRE_FALSE(nbnanasep[3].empty());
  int nb3site = *(nbnanasep[3].begin());
  CHECK(3 == nb3site);
  tb->setCustom("Na1+", 5);
  tb->setCustom("Cl1-", 5);
  molc->eval(mnacl);
  auto nbnana = molc->neighborhoods();
  REQUIRE(1u == nbnana.size());
  REQUIRE(4u == nbnana[0].size());
  CHECK(nbnana[0].count(0));
  CHECK(!nbnana[0].count(4));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_NaCl_mixed_overlap",
                 "[OverlapCalculator]") {
  StructureAdapterPtr nacl_mixed;
  nacl_mixed = loadTestPeriodicStructure("NaCl_mixed.stru");
  molc->eval(nacl_mixed);
  double olp = pow(3.3 - 5.62 / 2, 2) * 6 / 2;
  CHECK_THAT(olp, WithinAbs(molc->meanSquareOverlap(), meps));
}

TEST_CASE_METHOD(TestOverlapCalculatorFixture,
                 "TestOverlapCalculator::test_serialization",
                 "[OverlapCalculator]") {
  // build customized PDFCalculator
  AtomicStructureAdapterPtr stru(new AtomicStructureAdapter);
  Atom a0, a1;
  a0.xyz_cartn = R3::Vector(0.0, 0.0, 0.0);
  a1.xyz_cartn = R3::Vector(0.0, 0.0, 1.5);
  stru->append(a0);
  stru->append(a1);
  molc->eval(stru);
  std::shared_ptr<OverlapCalculator> olc1;
  olc1 = dumpandload(molc);
  CHECK(molc.get() != olc1.get());
  CHECK(2 == olc1->getStructure()->countSites());
  CHECK(0.25 == olc1->totalSquareOverlap());
  QuantityType sqolps = molc->siteSquareOverlaps();
  REQUIRE(2u == olc1->siteSquareOverlaps().size());
  CHECK(0.125 == olc1->siteSquareOverlaps()[0]);
  CHECK(0.125 == olc1->siteSquareOverlaps()[1]);
}

// End of file
