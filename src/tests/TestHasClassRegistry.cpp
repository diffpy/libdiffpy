/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2016 Brookhaven Science Associates,
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
 * TestHasClassRegistry -- test registration machinery for named classes.
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>

using diffpy::srreal::ScatteringFactorTable;
using diffpy::srreal::ScatteringFactorTablePtr;

//////////////////////////////////////////////////////////////////////////////
// TestHasClassRegistry
//////////////////////////////////////////////////////////////////////////////

class TestHasClassRegistryFixture {
 protected:
  typedef diffpy::HasClassRegistry<ScatteringFactorTable> HCRBase;
  ScatteringFactorTablePtr msftb;
  HCRBase* mpreg;

 public:
  TestHasClassRegistryFixture() {
    msftb = ScatteringFactorTable::createByType("xray");
    mpreg = msftb.get();
  }

  ~TestHasClassRegistryFixture() {
    // restore SFTXray registration if removed in some test.
    msftb->deregisterType(msftb->type());
    msftb->registerThisType();
    CHECK(ScatteringFactorTable::isRegisteredType("X"));
    CHECK(ScatteringFactorTable::isRegisteredType("xray"));
  }
};  // class TestHasClassRegistryFixture

TEST_CASE_METHOD(TestHasClassRegistryFixture,
                 "TestHasClassRegistry::test_deregisterType",
                 "[HasClassRegistry]") {
  CHECK(0 == ScatteringFactorTable::deregisterType("invalid"));
  CHECK(2 == ScatteringFactorTable::deregisterType("X"));
  CHECK(!(mpreg->isRegisteredType("xray")));
  CHECK(!(mpreg->isRegisteredType("X")));
}

TEST_CASE_METHOD(TestHasClassRegistryFixture,
                 "TestHasClassRegistry::test_isRegisteredType",
                 "[HasClassRegistry]") {
  CHECK(msftb->isRegisteredType("xray"));
  CHECK(ScatteringFactorTable::isRegisteredType("X"));
  CHECK(false == mpreg->isRegisteredType("invalid"));
}

TEST_CASE_METHOD(TestHasClassRegistryFixture,
                 "TestHasClassRegistry::test_getAliasedTypes",
                 "[HasClassRegistry]") {
  std::map<std::string, std::string> atps;
  atps = ScatteringFactorTable::getAliasedTypes();
  CHECK(4 == atps.size());
  CHECK(atps.count("X"));
  CHECK("xray" == atps["X"]);
  CHECK("electronnumber" == atps["EN"]);
}

// End of file
