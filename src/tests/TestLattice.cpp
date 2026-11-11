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
 * TestLattice -- unit tests for the Lattice class
 *
 *****************************************************************************/

#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/Lattice.hpp>
#include <diffpy/srreal/R3linalg.hpp>

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;
using diffpy::mathutils::EpsilonEqual;

class TestLatticeFixture {
 protected:
  Lattice* lattice;
  double precision;
  EpsilonEqual allclose;

 public:
  TestLatticeFixture() {
    lattice = new Lattice();
    precision = 1.0e-8;
    allclose = EpsilonEqual(precision);
  }

  ~TestLatticeFixture() { delete lattice; }
};  // class TestLatticeFixture

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_Lattice", "[Lattice]") {
  // default lattice should be cartesian
  CHECK(allclose(lattice->base(), R3::identity()));
  // lattice parameters constructor
  Lattice lattice1(1.0, 2.0, 3.0, 90, 90, 120);
  R3::Vector va = lattice1.va();
  R3::Vector vb = lattice1.vb();
  R3::Vector vc = lattice1.vc();
  double adotb = R3::dot(va, vb);
  double adotc = R3::dot(va, vc);
  double bdotc = R3::dot(vb, vc);
  CHECK_THAT(1.0, WithinAbs(R3::norm(va), precision));
  CHECK_THAT(2.0, WithinAbs(R3::norm(vb), precision));
  CHECK_THAT(3.0, WithinAbs(R3::norm(vc), precision));
  CHECK_THAT(-0.5 * 1.0 * 2.0, WithinAbs(adotb, precision));
  CHECK_THAT(0.0, WithinAbs(adotc, precision));
  CHECK_THAT(0.0, WithinAbs(bdotc, precision));
  CHECK_THAT(sqrt(4.0 / 3), WithinAbs(lattice1.ar(), precision));
  CHECK_THAT(sqrt(1.0 / 3), WithinAbs(lattice1.br(), precision));
  CHECK_THAT(1.0 / 3, WithinAbs(lattice1.cr(), precision));
  CHECK_THAT(90.0, WithinAbs(lattice1.alphar(), precision));
  CHECK_THAT(90.0, WithinAbs(lattice1.betar(), precision));
  CHECK_THAT(60.0, WithinAbs(lattice1.gammar(), precision));
  // lattice vectors constructor
  va = R3::Vector(1.0, 1.0, 0.0);
  vb = R3::Vector(0.0, 1.0, 1.0);
  vc = R3::Vector(1.0, 0.0, 1.0);
  Lattice lattice2(va, vb, vc);
  CHECK_THAT(sqrt(2.0), WithinAbs(lattice2.a(), precision));
  CHECK_THAT(60.0, WithinAbs(lattice2.alpha(), precision));
}

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_setLatPar",
                 "[Lattice]") {
  lattice->setLatPar(1.0, 2.0, 3.0, 90, 90, 120);
  R3::Matrix base_check(sqrt(0.75), -0.5, 0.0, 0.0, +2.0, 0.0, 0.0, +0.0, 3.0);
  CHECK(allclose(base_check, lattice->base()));
  R3::Matrix recbase_check(sqrt(4.0 / 3), sqrt(1.0 / 12), 0.0, 0.0, 0.5, 0.0,
                           0.0, 0.0, 1.0 / 3.0);
  CHECK(allclose(recbase_check, lattice->recbase()));
}

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_setLatBase",
                 "[Lattice]") {
  R3::Vector va(1.0, 1.0, 0.0);
  R3::Vector vb(0.0, 1.0, 1.0);
  R3::Vector vc(1.0, 0.0, 1.0);
  lattice->setLatBase(va, vb, vc);
  CHECK_THAT(sqrt(2.0), WithinAbs(lattice->a(), precision));
  CHECK_THAT(sqrt(2.0), WithinAbs(lattice->b(), precision));
  CHECK_THAT(sqrt(2.0), WithinAbs(lattice->c(), precision));
  CHECK_THAT(60.0, WithinAbs(lattice->alpha(), precision));
  CHECK_THAT(60.0, WithinAbs(lattice->beta(), precision));
  CHECK_THAT(60.0, WithinAbs(lattice->gamma(), precision));
  // check determinant of rotation matrix
  double detRot0 = R3::determinant(lattice->baserot());
  CHECK_THAT(1.0, WithinAbs(detRot0, precision));
  // check if rotation matrix works
  R3::Matrix base_check(va[0], va[1], va[2], vb[0], vb[1], vb[2], vc[0], vc[1],
                        vc[2]);
  CHECK(allclose(base_check, lattice->base()));
  R3::Matrix recbase_check(0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5);
  CHECK(allclose(recbase_check, lattice->recbase()));
  lattice->setLatPar(lattice->a(), lattice->b(), lattice->c(), 44.0, 66.0,
                     88.0);
  CHECK(!allclose(base_check, lattice->base()));
  CHECK(!allclose(recbase_check, lattice->recbase()));
  lattice->setLatPar(lattice->a(), lattice->b(), lattice->c(), 60.0, 60.0,
                     60.0);
  CHECK(allclose(base_check, lattice->base()));
  CHECK(allclose(recbase_check, lattice->recbase()));
}

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_invalid_unchanged",
                 "[Lattice]") {
  // A degenerate lattice throws and leaves the lattice as it was.
  lattice->setLatPar(1.0, 2.0, 3.0, 80, 90, 100);
  const Lattice saved = *lattice;
  const R3::Vector v(0.3, -1.7, 2.9);
  auto unchanged = [&]() {
    return lattice->a() == saved.a() && lattice->b() == saved.b() &&
           lattice->c() == saved.c() && lattice->alpha() == saved.alpha() &&
           lattice->beta() == saved.beta() &&
           lattice->gamma() == saved.gamma() && lattice->ar() == saved.ar() &&
           lattice->base() == saved.base() &&
           lattice->recbase() == saved.recbase() &&
           lattice->baserot() == saved.baserot() &&
           lattice->metrics() == saved.metrics() &&
           lattice->fractional(v) == saved.fractional(v) &&
           lattice->cartesian(v) == saved.cartesian(v);
  };
  REQUIRE(unchanged());
  CHECK_THROWS_AS(lattice->setLatPar(0.0, 4.0, 4.0, 90, 90, 90),
                  std::invalid_argument);
  CHECK(unchanged());
  const R3::Vector va(1.0, 0.0, 0.0);
  const R3::Vector vc(0.0, 0.0, 1.0);
  CHECK_THROWS_AS(lattice->setLatBase(va, va, vc), std::invalid_argument);
  CHECK(unchanged());
  CHECK_THROWS_AS(lattice->setLatBase(va, vc, R3::Vector(0.0, 1.0, 0.0)),
                  std::invalid_argument);
  CHECK(unchanged());
}

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_dist", "[Lattice]") {
  R3::Vector va(1.0, 2.0, 2.0);
  R3::Vector vb(0.0, 0.0, 0.0);
  CHECK_THAT(3.0, WithinAbs(lattice->distance(va, vb), precision));
  lattice->setLatPar(2.0, 2.0, 2.0, 90.0, 90.0, 90.0);
  CHECK_THAT(6.0, WithinAbs(lattice->distance(va, vb), precision));
}

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_angle", "[Lattice]") {
  R3::Vector va(1.0, 0.0, 0.0);
  R3::Vector vb(0.0, 1.0, 0.0);
  CHECK_THAT(90.0, WithinAbs(lattice->angledeg(va, vb), precision));
  lattice->setLatPar(2.0, 2.0, 2.0, 90.0, 90.0, 120.0);
  CHECK_THAT(120.0, WithinAbs(lattice->angledeg(va, vb), precision));
}

TEST_CASE_METHOD(TestLatticeFixture, "TestLattice::test_ucvCartesian",
                 "[Lattice]") {
  R3::Vector ucv;
  ucv = lattice->ucvCartesian(R3::Vector(0.1, 0.2, 0.3));
  R3::Vector ucv_check(0.1, 0.2, 0.3);
  CHECK(allclose(ucv_check, ucv));
  ucv = lattice->ucvCartesian(R3::Vector(1.1, 13.2, -0.7));
  CHECK(allclose(ucv_check, ucv));
  ucv = lattice->ucvCartesian(R3::Vector(0.5, 0.5, 0.5));
  ucv_check = R3::Vector(0.5, 0.5, 0.5);
  CHECK(allclose(ucv_check, ucv));
  lattice->setLatPar(13, 17, 19, 37, 41, 47);
  ucv = lattice->ucvCartesian(R3::Vector(100.0, 100.0, 100.0));
  ucv_check = R3::Vector(8.09338442077, 9.55056747225, 30.0389043325);
  CHECK(allclose(ucv_check, ucv));
}

// End of file
