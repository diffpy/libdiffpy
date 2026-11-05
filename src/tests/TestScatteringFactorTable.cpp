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
 * TestScatteringFactorTable -- unit tests for implementations
 *     of the ScatteringFactorTable class
 *
 *****************************************************************************/

#include <typeinfo>
#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/mathutils.hpp>
#include "serialization_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;

class TestScatteringFactorTableFixture {
 protected:
  double mtol;
  double meps;
  ScatteringFactorTablePtr msftb;

 public:
  TestScatteringFactorTableFixture() {
    mtol = 1.0e-4;
    meps = 10 * diffpy::mathutils::DOUBLE_EPS;
  }
};  // class TestScatteringFactorTableFixture

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_factory",
                 "[ScatteringFactorTable]") {
  ScatteringFactorTablePtr sfx0, sfx1, sfn0, sfn1;
  CHECK_THROWS_AS(ScatteringFactorTable::createByType("invalid"),
                  invalid_argument);
  sfx0 = ScatteringFactorTable::createByType("xray");
  sfx1 = ScatteringFactorTable::createByType("X");
  CHECK(sfx0.get());
  CHECK(sfx1.get());
  CHECK(sfx0->type() == sfx1->type());
  ScatteringFactorTable& r_sfx0 = *sfx0;
  ScatteringFactorTable& r_sfx1 = *sfx1;
  CHECK(typeid(r_sfx0) == typeid(r_sfx1));
  sfn0 = ScatteringFactorTable::createByType("neutron");
  sfn1 = ScatteringFactorTable::createByType("N");
  CHECK(sfn0.get());
  CHECK(sfn1.get());
  CHECK(sfn0->type() == sfn1->type());
  ScatteringFactorTable& r_sfn0 = *sfn0;
  ScatteringFactorTable& r_sfn1 = *sfn1;
  CHECK(typeid(r_sfn0) == typeid(r_sfn1));
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_setCustomAs",
                 "[ScatteringFactorTable]") {
  msftb = ScatteringFactorTable::createByType("X");
  CHECK_THROWS_AS(msftb->lookup(""), invalid_argument);
  CHECK_THAT(6.0, WithinAbs(msftb->lookup("C"), 0.01));
  msftb->setCustomAs("C", "C", 6.3);
  CHECK_THAT(6.3, WithinAbs(msftb->lookup("C"), meps));
  msftb->setCustomAs("C", "C", 6.4);
  CHECK_THAT(6.4, WithinAbs(msftb->lookup("C"), meps));
  CHECK_THROWS_AS(msftb->lookup("Ccustom"), invalid_argument);
  msftb->setCustomAs("Ccustom", "C", 6.5);
  CHECK_THAT(6.5, WithinAbs(msftb->lookup("Ccustom"), meps));
  msftb->resetCustom("C");
  CHECK_THAT(6.5, WithinAbs(msftb->lookup("Ccustom"), meps));
  CHECK_THAT(6.0, WithinAbs(msftb->lookup("C"), 0.01));
  msftb->resetAll();
  CHECK_THROWS_AS(msftb->lookup("Ccustom"), invalid_argument);
  CHECK_THAT(6.0, WithinAbs(msftb->lookup("C"), 0.01));
  msftb->setCustomAs("Calias", "C");
  CHECK(msftb->lookup("C", 0) == msftb->lookup("Calias", 0));
  CHECK_THAT(msftb->lookup("C", 2.5),
             WithinAbs(msftb->lookup("Calias", 2.5), meps));
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_getCustomSymbols",
                 "[ScatteringFactorTable]") {
  msftb = ScatteringFactorTable::createByType("X");
  CHECK(msftb->getCustomSymbols().empty());
  msftb->setCustomAs("C", "C", 6.1);
  CHECK(1u == msftb->getCustomSymbols().size());
  CHECK(1u == msftb->getCustomSymbols().count("C"));
  msftb->setCustomAs("C", "C", 6.3);
  CHECK(1u == msftb->getCustomSymbols().size());
  CHECK(1u == msftb->getCustomSymbols().count("C"));
  CHECK_THAT(6.3, WithinAbs(msftb->lookup("C"), meps));
  ScatteringFactorTablePtr sftb1 = msftb->clone();
  msftb->resetCustom("C");
  CHECK(msftb->getCustomSymbols().empty());
  CHECK(1u == sftb1->getCustomSymbols().size());
  CHECK(1u == sftb1->getCustomSymbols().count("C"));
  CHECK_THAT(6.3, WithinAbs(sftb1->lookup("C"), meps));
  sftb1->resetAll();
  CHECK(msftb->getCustomSymbols().empty());
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_ticker",
                 "[ScatteringFactorTable]") {
  using diffpy::eventticker::EventTicker;
  msftb = ScatteringFactorTable::createByType("X");
  EventTicker e0 = msftb->ticker();
  CHECK(e0 == msftb->clone()->ticker());
  CHECK(e0 == dumpandload(msftb)->ticker());
  msftb->setCustomAs("C", "C", 6.1);
  CHECK(e0 < msftb->ticker());
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_SFTXray",
                 "[ScatteringFactorTable]") {
  msftb = ScatteringFactorTable::createByType("X");
  CHECK_THAT(1.0, WithinAbs(msftb->lookup("H"), 0.01));
  CHECK_THAT(8.0, WithinAbs(msftb->lookup("O"), 0.01));
  CHECK_THAT(10.0, WithinAbs(msftb->lookup("O2-"), 0.01));
  CHECK_THAT(11.0, WithinAbs(msftb->lookup("Na"), 0.01));
  CHECK_THAT(10.0, WithinAbs(msftb->lookup("Na+"), 0.01));
  CHECK(msftb->lookup("Na+") == msftb->lookup("Na1+"));
  CHECK_THAT(74.0, WithinAbs(msftb->lookup("W"), 0.04));
  CHECK_THAT(88.0, WithinAbs(msftb->lookup("Ra"), 0.04));
  CHECK(msftb->lookup("Si") == msftb->lookup("Si0+"));
  CHECK("xray" == msftb->clone()->type());
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_SFTElectron",
                 "[ScatteringFactorTable]") {
  using diffpy::mathutils::DOUBLE_MAX;
  msftb = ScatteringFactorTable::createByType("E");
  CHECK(DOUBLE_MAX == msftb->lookup("H"));
  CHECK(DOUBLE_MAX == msftb->lookup("Ra"));
  CHECK_THAT(3.42104, WithinAbs(msftb->lookup("Na", 1), 1e-5));
  CHECK_THAT(1.34868, WithinAbs(msftb->lookup("Na", 3), 1e-5));
  CHECK_THAT(0.832158, WithinAbs(msftb->lookup("Na", 5), 1e-5));
  CHECK_THROWS_AS(msftb->lookup("H4+"), invalid_argument);
  CHECK_THROWS_AS(msftb->lookup("H4+", 3), invalid_argument);
  CHECK("electron" == msftb->clone()->type());
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_SFTNeutron",
                 "[ScatteringFactorTable]") {
  msftb = ScatteringFactorTable::createByType("N");
  CHECK_THAT(3.63, WithinAbs(msftb->lookup("Na"), mtol));
  CHECK_THAT(-3.37, WithinAbs(msftb->lookup("Ti"), mtol));
  CHECK_THAT(5.805, WithinAbs(msftb->lookup("O"), mtol));
  CHECK_THAT(6.6484, WithinAbs(msftb->lookup("C"), mtol));
  CHECK(msftb->lookup("Na") == msftb->lookup("Na1+"));
  CHECK(msftb->lookup("Ge") == msftb->lookup("Ge0+"));
  CHECK("neutron" == msftb->clone()->type());
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_ElectronNumber",
                 "[ScatteringFactorTable]") {
  msftb = ScatteringFactorTable::createByType("electronnumber");
  CHECK(8.0 == msftb->lookup("O"));
  CHECK(10.0 == msftb->lookup("O2-"));
  CHECK(18.0 == msftb->lookup("K+"));
  CHECK(18.0 == msftb->lookup("K1+"));
  CHECK(68.0 == msftb->lookup("W6+"));
  CHECK_THROWS_AS(msftb->lookup("H4+"), invalid_argument);
  CHECK_THROWS_AS(msftb->lookup("O3+"), invalid_argument);
  CHECK("electronnumber" == msftb->clone()->type());
}

TEST_CASE_METHOD(TestScatteringFactorTableFixture,
                 "TestScatteringFactorTable::test_serialization",
                 "[ScatteringFactorTable]") {
  ScatteringFactorTablePtr sftb1;
  msftb = ScatteringFactorTable::createByType("electronnumber");
  msftb->setCustomAs("H", "H", 1.23);
  sftb1 = dumpandload(msftb);
  CHECK(string("electronnumber") == sftb1->type());
  CHECK(1.23 == sftb1->lookup("H"));
  CHECK(1u == sftb1->getCustomSymbols().size());
  sftb1 = dumpandload(ScatteringFactorTable::createByType("N"));
  CHECK(string("neutron") == sftb1->type());
  sftb1 = dumpandload(ScatteringFactorTable::createByType("X"));
  CHECK(string("xray") == sftb1->type());
  sftb1 = dumpandload(ScatteringFactorTable::createByType("E"));
  CHECK(string("electron") == sftb1->type());
}

// End of file
