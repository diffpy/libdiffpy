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
 * TestEmptyStructureAdapter -- test adapter that is always empty
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// TestEmptyStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestEmptyStructureAdapterFixture {
 protected:
  StructureAdapterPtr mstru;

 public:
  TestEmptyStructureAdapterFixture() { mstru = emptyStructureAdapter(); }
};  // class TestEmptyStructureAdapterFixture

TEST_CASE_METHOD(TestEmptyStructureAdapterFixture,
                 "TestEmptyStructureAdapter::test_empty_instance",
                 "[EmptyStructureAdapter]") {
  CHECK(0 == mstru->countSites());
  CHECK(mstru.get() == mstru->clone().get());
  CHECK(mstru == emptyStructureAdapter());
  CHECK_THROWS_AS(mstru->siteAtomType(0), std::out_of_range);
  CHECK_THROWS_AS(mstru->siteCartesianPosition(0), std::out_of_range);
  CHECK_THROWS_AS(mstru->siteMultiplicity(0), std::out_of_range);
  CHECK_THROWS_AS(mstru->siteOccupancy(0), std::out_of_range);
  CHECK_THROWS_AS(mstru->siteAnisotropy(0), std::out_of_range);
  CHECK_THROWS_AS(mstru->siteCartesianPosition(0), std::out_of_range);
}

TEST_CASE_METHOD(TestEmptyStructureAdapterFixture,
                 "TestEmptyStructureAdapter::test_serialization",
                 "[EmptyStructureAdapter]") {
  StructureAdapterPtr stru1;
  stru1 = dumpandload(mstru);
  AtomicStructureAdapterPtr astru1 =
    std::dynamic_pointer_cast<AtomicStructureAdapter>(stru1);
  CHECK(!astru1);
  CHECK(0 == stru1->countSites());
  CHECK(stru1 == stru1->clone());
  CHECK_THROWS_AS(mstru->siteAtomType(0), std::out_of_range);
}

}  // namespace srreal
}  // namespace diffpy

// End of file
