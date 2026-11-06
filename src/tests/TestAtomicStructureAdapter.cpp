/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2013 Brookhaven Science Associates,
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
 * TestAtomicStructureAdapter -- unit tests for an adapter that
 *     stores data in a series of Atom objects
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <memory>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// TestAtomicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestAtomicStructureAdapterFixture {
 protected:
  StructureAdapterPtr mstru;
  AtomicStructureAdapterPtr mpstru;

 public:
  TestAtomicStructureAdapterFixture() {
    mstru = StructureAdapterPtr(new AtomicStructureAdapter);
    mpstru = std::dynamic_pointer_cast<AtomicStructureAdapter>(mstru);
  }
};  // class TestAtomicStructureAdapterFixture

TEST_CASE_METHOD(TestAtomicStructureAdapterFixture,
                 "TestAtomicStructureAdapter::test_diff",
                 "[AtomicStructureAdapter]") {
  typedef StructureDifference::Method DM;
  const DM::Type& NONE = DM::NONE;
  const DM::Type& SIDEBYSIDE = DM::SIDEBYSIDE;
  const DM::Type& SORTED = DM::SORTED;
  StructureDifference sd;
  sd = mstru->diff(emptyStructureAdapter());
  CHECK(sd.add1.empty());
  CHECK(!sd.allowsfastupdate());
  sd = mstru->diff(StructureAdapterPtr());
  CHECK(NONE == sd.diffmethod);
  CHECK(!sd.allowsfastupdate());
  Atom ai;
  ai.atomtype = "C";
  const int SZ = 10;
  for (int i = 0; i < SZ; ++i) {
    ai.xyz_cartn[0] = i;
    mpstru->append(ai);
  }
  AtomicStructureAdapterPtr cpstru =
    std::make_shared<AtomicStructureAdapter>(*mpstru);
  sd = mstru->diff(mstru);
  CHECK(sd.allowsfastupdate());
  sd = mstru->diff(cpstru);
  CHECK(SIDEBYSIDE == sd.diffmethod);
  CHECK(sd.allowsfastupdate());
  CHECK(sd.pop0.empty());
  CHECK(sd.add1.empty());
  (*cpstru)[0].atomtype = "N";
  sd = mstru->diff(cpstru);
  CHECK(SIDEBYSIDE == sd.diffmethod);
  CHECK(sd.allowsfastupdate());
  CHECK(1u == sd.pop0.size());
  CHECK(1u == sd.add1.size());
  for (int i = 1; i < (1 - sqrt(0.5)) * SZ; ++i) {
    cpstru->erase(0);
    sd = mstru->diff(cpstru);
    CHECK(SORTED == sd.diffmethod);
    CHECK(sd.allowsfastupdate());
    CHECK(i == int(sd.pop0.size()));
    CHECK(sd.add1.empty());
  }
  cpstru->erase(0);
  sd = mstru->diff(cpstru);
  CHECK(!sd.allowsfastupdate());
  Atom a2;
  a2.atomtype = "N";
  cpstru->append(a2);
  sd = mstru->diff(cpstru);
  CHECK(!sd.allowsfastupdate());
  CHECK(1u == sd.add1.size());
}

TEST_CASE_METHOD(TestAtomicStructureAdapterFixture,
                 "TestAtomicStructureAdapter::test_serialization",
                 "[AtomicStructureAdapter]") {
  Atom ai;
  ai.atomtype = "C";
  ai.xyz_cartn = R3::Vector(1, 2, 3);
  ai.anisotropy = true;
  mpstru->append(ai);
  ai.atomtype = "H";
  ai.xyz_cartn = R3::Vector(4, 5, 6);
  ai.anisotropy = false;
  mpstru->append(ai);
  StructureAdapterPtr stru1;
  stru1 = dumpandload(mstru);
  AtomicStructureAdapterPtr astru1 =
    std::dynamic_pointer_cast<AtomicStructureAdapter>(stru1);
  CHECK(2 == astru1->countSites());
  CHECK((*mpstru)[0] == (*astru1)[0]);
  CHECK((*mpstru)[1] == (*astru1)[1]);
}

TEST_CASE_METHOD(TestAtomicStructureAdapterFixture,
                 "TestAtomicStructureAdapter::test_comparison",
                 "[AtomicStructureAdapter]") {
  Atom ai;
  ai.atomtype = "C";
  const int SZ = 10;
  for (int i = 0; i < SZ; ++i) {
    ai.xyz_cartn[0] = i;
    mpstru->append(ai);
  }
  AtomicStructureAdapterPtr cpstru =
    std::make_shared<AtomicStructureAdapter>(*mpstru);
  CHECK(*mpstru == *cpstru);
  CHECK(!(*mpstru != *cpstru));
  cpstru->at(0).atomtype = "H";
  CHECK(*mpstru != *cpstru);
  CHECK(!(*mpstru == *cpstru));
}

}  // namespace srreal
}  // namespace diffpy

// End of file
