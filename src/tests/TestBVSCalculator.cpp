/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2010 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestBVSCalculator -- unit tests suite
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>
#include "test_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using diffpy::srreal::BVParametersTablePtr;
using diffpy::srreal::BVSCalculator;
using diffpy::srreal::PeriodicStructureAdapter;
using diffpy::srreal::PeriodicStructureAdapterPtr;
using diffpy::srreal::StructureAdapterPtr;

//////////////////////////////////////////////////////////////////////////////
// TestBVSCalculator
//////////////////////////////////////////////////////////////////////////////

class TestBVSCalculatorFixture {
 protected:
  StructureAdapterPtr mnacl;
  std::shared_ptr<BVSCalculator> mbvc;

 public:
  TestBVSCalculatorFixture() {
    mnacl = loadTestPeriodicStructure("NaCl.stru");
    mbvc.reset(new BVSCalculator);
  }
};  // class TestBVSCalculatorFixture

TEST_CASE_METHOD(TestBVSCalculatorFixture, "TestBVSCalculator::test_NaCl",
                 "[BVSCalculator]") {
  const double eps = 1e-4;
  mbvc->eval(mnacl);
  CHECK(8u == mbvc->value().size());
  CHECK_THAT(+1.01352, WithinAbs(mbvc->value()[0], eps));
  CHECK_THAT(+1.01352, WithinAbs(mbvc->value()[1], eps));
  CHECK_THAT(+1.01352, WithinAbs(mbvc->value()[2], eps));
  CHECK_THAT(+1.01352, WithinAbs(mbvc->value()[3], eps));
  CHECK_THAT(-1.01352, WithinAbs(mbvc->value()[4], eps));
  CHECK_THAT(-1.01352, WithinAbs(mbvc->value()[5], eps));
  CHECK_THAT(-1.01352, WithinAbs(mbvc->value()[6], eps));
  CHECK_THAT(-1.01352, WithinAbs(mbvc->value()[7], eps));
}

TEST_CASE_METHOD(TestBVSCalculatorFixture, "TestBVSCalculator::test_NaCl_mixed",
                 "[BVSCalculator]") {
  StructureAdapterPtr nacl_mixed = loadTestPeriodicStructure("NaCl_mixed.stru");
  mbvc->eval(mnacl);
  BVSCalculator bvc;
  bvc.eval(nacl_mixed);
  CHECK_THAT(mbvc->bvrmsdiff(), WithinAbs(bvc.bvrmsdiff(), 1e-12));
}

TEST_CASE_METHOD(TestBVSCalculatorFixture,
                 "TestBVSCalculator::test_customAtomValences",
                 "[BVSCalculator]") {
  using diffpy::srreal::Atom;
  const double eps = 1e-4;
  mbvc->eval(mnacl);
  CHECK_THAT(0.01352, WithinAbs(mbvc->bvrmsdiff(), eps));
  BVParametersTablePtr bvtb = mbvc->getBVParamTable();
  bvtb->setAtomValence("Cl1-", 0);
  mbvc->eval(mnacl);
  CHECK(0.0 == mbvc->value()[0]);
  CHECK(0.0 == mbvc->value()[4]);
  // create structure with bare atom symbols "Na", "Cl".
  PeriodicStructureAdapterPtr naclbare =
    std::dynamic_pointer_cast<PeriodicStructureAdapter>(mnacl->clone());
  for (int i = 0; i < naclbare->countSites(); ++i) {
    Atom& a = naclbare->at(i);
    a.atomtype = a.atomtype.substr(0, 2);
  }
  CHECK(string("Na") == naclbare->siteAtomType(0));
  // verify valence sums are zero for a standard setup.
  mbvc->eval(naclbare);
  CHECK(0.0 == mbvc->value()[0]);
  CHECK(0.0 == mbvc->value()[4]);
  // verify valence sums with custom atom valences.
  bvtb->setAtomValence("Na", +1);
  bvtb->setAtomValence("Cl", -1);
  mbvc->eval(naclbare);
  CHECK_THAT(0.01352, WithinAbs(mbvc->bvrmsdiff(), eps));
}

TEST_CASE_METHOD(TestBVSCalculatorFixture,
                 "TestBVSCalculator::test_setValencePrecision",
                 "[BVSCalculator]") {
  CHECK_THROWS_AS(mbvc->setValencePrecision(0), invalid_argument);
  CHECK_THROWS_AS(mbvc->setValencePrecision(-1), invalid_argument);
  CHECK_THROWS_AS(mbvc->setDoubleAttr("valenceprecision", 0), invalid_argument);
}

TEST_CASE_METHOD(TestBVSCalculatorFixture,
                 "TestBVSCalculator::test_getRmaxUsed", "[BVSCalculator]") {
  mbvc->eval(mnacl);
  double rmaxused = mbvc->getDoubleAttr("rmaxused");
  CHECK(rmaxused == mbvc->getRmaxUsed());
  CHECK(rmaxused < mbvc->getRmax());
  mbvc->setDoubleAttr("valenceprecision",
                      mbvc->getDoubleAttr("valenceprecision") / 10.0);
  CHECK(rmaxused < mbvc->getDoubleAttr("rmaxused"));
  mbvc->setValencePrecision(1e-7);
  mbvc->setDoubleAttr("rmax", 5.0);
  CHECK(5.0 == mbvc->getDoubleAttr("rmaxused"));
  // check if value updates with changes in the table.
  BVParametersTablePtr bvtb = mbvc->getBVParamTable();
  bvtb->setCustom("Na", 1, "Cl", -1, 0.0, 0);
  CHECK(0.0 == mbvc->getDoubleAttr("rmaxused"));
}

TEST_CASE_METHOD(TestBVSCalculatorFixture,
                 "TestBVSCalculator::test_serialization", "[BVSCalculator]") {
  mbvc->eval(mnacl);
  stringstream storage(ios::in | ios::out | ios::binary);
  diffpy::serialization::oarchive oa(storage, ios::binary);
  oa << mbvc;
  diffpy::serialization::iarchive ia(storage, ios::binary);
  std::shared_ptr<BVSCalculator> bvc1;
  CHECK(!bvc1.get());
  ia >> bvc1;
  CHECK(mbvc.get() != bvc1.get());
  const double eps = 1e-4;
  CHECK(8u == bvc1->value().size());
  CHECK_THAT(+1.01352, WithinAbs(bvc1->value()[0], eps));
  CHECK_THAT(+1.01352, WithinAbs(bvc1->value()[1], eps));
  CHECK_THAT(+1.01352, WithinAbs(bvc1->value()[2], eps));
  CHECK_THAT(+1.01352, WithinAbs(bvc1->value()[3], eps));
  CHECK_THAT(-1.01352, WithinAbs(bvc1->value()[4], eps));
  CHECK_THAT(-1.01352, WithinAbs(bvc1->value()[5], eps));
  CHECK_THAT(-1.01352, WithinAbs(bvc1->value()[6], eps));
  CHECK_THAT(-1.01352, WithinAbs(bvc1->value()[7], eps));
}

// End of file
