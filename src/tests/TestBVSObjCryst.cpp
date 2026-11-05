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
 * TestBVSObjCryst -- unit tests for BVS calculation for an ObjCryst
 *   crystal structure
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>

#include "objcryst_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;

//////////////////////////////////////////////////////////////////////////////
// class Test
//////////////////////////////////////////////////////////////////////////////

class TestBVSObjCrystFixture {
 protected:
  const ObjCryst::Crystal& mnacl = cachedTestCrystal("NaCl.cif");
  BVSCalculator mbvc;
};  // class TestBVSObjCrystFixture

TEST_CASE_METHOD(TestBVSObjCrystFixture, "TestBVSObjCryst::test_NaCl",
                 "[BVSObjCryst]") {
  const double eps = 1e-4;
  mbvc.eval(mnacl);
  CHECK(2u == mbvc.value().size());
  CHECK_THAT(+1.01352, WithinAbs(mbvc.value()[0], eps));
  CHECK_THAT(-1.01352, WithinAbs(mbvc.value()[1], eps));
}

TEST_CASE_METHOD(TestBVSObjCrystFixture, "TestBVSObjCryst::test_NaCl_mixed",
                 "[BVSObjCryst]") {
  unique_ptr<ObjCryst::Crystal> nacl_mixed;
  nacl_mixed.reset(loadTestCrystal("NaCl_mixed.cif"));
  mbvc.eval(mnacl);
  BVSCalculator bvc;
  bvc.eval(*nacl_mixed);
  CHECK_THAT(mbvc.bvrmsdiff(), WithinAbs(bvc.bvrmsdiff(), 1e-12));
}

// End of file
