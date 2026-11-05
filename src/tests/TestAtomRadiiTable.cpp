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
 * TestAtomRadiiTable -- unit tests for the AtomRadiiTable class
 *
 *****************************************************************************/

#include <stdexcept>
#include <catch2/catch_test_macros.hpp>

#include "serialization_helpers.hpp"
#include <diffpy/srreal/AtomRadiiTable.hpp>
#include <diffpy/srreal/ConstantRadiiTable.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestAtomRadiiTableFixture {
 protected:
  AtomRadiiTablePtr mrtb;

 public:
  TestAtomRadiiTableFixture() {
    mrtb = AtomRadiiTable::createByType("constant");
  }
};  // class TestAtomRadiiTableFixture

TEST_CASE_METHOD(TestAtomRadiiTableFixture, "TestAtomRadiiTable::test_lookup",
                 "[AtomRadiiTable]") {
  CHECK(0.0 == mrtb->lookup("everything"));
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture,
                 "TestAtomRadiiTable::test_setCustom", "[AtomRadiiTable]") {
  mrtb->setCustom("C", 1.23);
  CHECK(1.23 == mrtb->lookup("C"));
  mrtb->setCustom("C", 1.5);
  CHECK(1.5 == mrtb->lookup("C"));
  mrtb->resetAll();
  CHECK(0.0 == mrtb->lookup("C"));
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture,
                 "TestAtomRadiiTable::test_fromString", "[AtomRadiiTable]") {
  CHECK(0u == mrtb->getAllCustom().size());
  mrtb->fromString("A:1.1, B:1.2,    C:  1.3  \t");
  CHECK(3u == mrtb->getAllCustom().size());
  CHECK(1.1 == mrtb->lookup("A"));
  CHECK(1.2 == mrtb->lookup("B"));
  CHECK(1.3 == mrtb->lookup("C"));
  CHECK_THROWS_AS(mrtb->fromString("X12"), invalid_argument);
  CHECK(3u == mrtb->getAllCustom().size());
  mrtb->fromString("D:4.3");
  CHECK(4.3 == mrtb->lookup("D"));
  CHECK(4u == mrtb->getAllCustom().size());
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture,
                 "TestAtomRadiiTable::test_resetCustom", "[AtomRadiiTable]") {
  mrtb->fromString("A:1.1, B:1.2,    C:  1.3");
  CHECK(3u == mrtb->getAllCustom().size());
  mrtb->resetCustom("B");
  CHECK(1.1 == mrtb->lookup("A"));
  CHECK(1.3 == mrtb->lookup("C"));
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture, "TestAtomRadiiTable::test_resetAll",
                 "[AtomRadiiTable]") {
  mrtb->fromString("A:1.1, B:1.2,    C:  1.3");
  CHECK(3u == mrtb->getAllCustom().size());
  mrtb->resetAll();
  CHECK(mrtb->getAllCustom().empty());
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture, "TestAtomRadiiTable::test_toString",
                 "[AtomRadiiTable]") {
  CHECK(mrtb->toString().empty());
  mrtb->fromString("A:1.1,    C:  1.3, B:1.2");
  CHECK(string("A:1.1,B:1.2,C:1.3") == mrtb->toString());
  CHECK(string("A:1.1, B:1.2, C:1.3") == mrtb->toString(", "));
  mrtb->resetAll();
  CHECK(mrtb->toString().empty());
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture,
                 "TestAtomRadiiTable::test_setDefault", "[AtomRadiiTable]") {
  ConstantRadiiTable* crtb = dynamic_cast<ConstantRadiiTable*>(mrtb.get());
  crtb->setDefault(7.1);
  CHECK(7.1 == mrtb->lookup("everything"));
  crtb->setDefault(3.1);
  CHECK(3.1 == mrtb->lookup("Na"));
}

TEST_CASE_METHOD(TestAtomRadiiTableFixture,
                 "TestAtomRadiiTable::test_serialization", "[AtomRadiiTable]") {
  AtomRadiiTablePtr rtb1;
  mrtb->setCustom("H", 1.23);
  ConstantRadiiTable* crtb = dynamic_cast<ConstantRadiiTable*>(mrtb.get());
  crtb->setDefault(0.5);
  rtb1 = dumpandload(mrtb);
  CHECK(1.23 == rtb1->lookup("H"));
  CHECK(1u == rtb1->getAllCustom().size());
  CHECK(mrtb.get() != rtb1.get());
  ConstantRadiiTable* crtb1 = dynamic_cast<ConstantRadiiTable*>(rtb1.get());
  CHECK(0.5 == crtb1->getDefault());
}

// End of file
