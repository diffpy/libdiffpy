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
 * TestSphericalShapeEnvelope -- unit tests for SphericalShapeEnvelope
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/PDFEnvelope.hpp>
#include "serialization_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;

class TestSphericalShapeEnvelopeFixture {
 protected:
  PDFEnvelopePtr menvelope;

 public:
  TestSphericalShapeEnvelopeFixture() {
    menvelope = PDFEnvelope::createByType("sphericalshape");
  }
};  // class TestSphericalShapeEnvelopeFixture

TEST_CASE_METHOD(TestSphericalShapeEnvelopeFixture,
                 "TestSphericalShapeEnvelope::test_create",
                 "[SphericalShapeEnvelope]") {
  CHECK(0.0 == menvelope->getDoubleAttr("spdiameter"));
  menvelope->setDoubleAttr("spdiameter", 13.0);
  CHECK(13.0 == menvelope->getDoubleAttr("spdiameter"));
  PDFEnvelopePtr e1 = menvelope->create();
  CHECK(0.0 == e1->getDoubleAttr("spdiameter"));
}

TEST_CASE_METHOD(TestSphericalShapeEnvelopeFixture,
                 "TestSphericalShapeEnvelope::test_copy",
                 "[SphericalShapeEnvelope]") {
  menvelope->setDoubleAttr("spdiameter", 13.0);
  CHECK(13.0 == menvelope->getDoubleAttr("spdiameter"));
  PDFEnvelopePtr e1 = menvelope->clone();
  CHECK(13.0 == e1->getDoubleAttr("spdiameter"));
}

TEST_CASE_METHOD(TestSphericalShapeEnvelopeFixture,
                 "TestSphericalShapeEnvelope::test_type",
                 "[SphericalShapeEnvelope]") {
  CHECK("sphericalshape" == menvelope->type());
}

TEST_CASE_METHOD(TestSphericalShapeEnvelopeFixture,
                 "TestSphericalShapeEnvelope::test_parentheses_operator",
                 "[SphericalShapeEnvelope]") {
  const PDFEnvelope& fne = *menvelope;
  CHECK(1.0 == fne(0.0));
  CHECK(1.0 == fne(100.0));
  menvelope->setDoubleAttr("spdiameter", 10.0);
  CHECK(1.0 == fne(0.0));
  CHECK(0.0 == fne(10.0));
  CHECK(0.0 == fne(100));
  CHECK(fne(9.99) > 0.0);
  CHECK_THAT(0.3125, WithinAbs(fne(5), 1e-8));
}

TEST_CASE_METHOD(TestSphericalShapeEnvelopeFixture,
                 "TestSphericalShapeEnvelope::test_serialization",
                 "[SphericalShapeEnvelope]") {
  menvelope->setDoubleAttr("spdiameter", 13.1);
  PDFEnvelopePtr e1 = dumpandload(menvelope);
  CHECK(string("sphericalshape") == e1->type());
  CHECK(13.1 == e1->getDoubleAttr("spdiameter"));
}

// End of file
