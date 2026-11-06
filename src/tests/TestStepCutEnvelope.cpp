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
 * TestStepCutEnvelope -- unit tests for StepCutEnvelope class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/PDFEnvelope.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

class TestStepCutEnvelopeFixture {
 protected:
  PDFEnvelopePtr menvelope;

 public:
  TestStepCutEnvelopeFixture() {
    menvelope = PDFEnvelope::createByType("stepcut");
  }
};  // class TestStepCutEnvelopeFixture

TEST_CASE_METHOD(TestStepCutEnvelopeFixture, "TestStepCutEnvelope::test_create",
                 "[StepCutEnvelope]") {
  CHECK(0.0 == menvelope->getDoubleAttr("stepcut"));
  menvelope->setDoubleAttr("stepcut", 13.0);
  CHECK(13.0 == menvelope->getDoubleAttr("stepcut"));
  PDFEnvelopePtr e1(menvelope->create());
  CHECK(0.0 == e1->getDoubleAttr("stepcut"));
}

TEST_CASE_METHOD(TestStepCutEnvelopeFixture, "TestStepCutEnvelope::test_copy",
                 "[StepCutEnvelope]") {
  menvelope->setDoubleAttr("stepcut", 13.0);
  CHECK(13.0 == menvelope->getDoubleAttr("stepcut"));
  PDFEnvelopePtr e1(menvelope->clone());
  CHECK(13.0 == e1->getDoubleAttr("stepcut"));
}

TEST_CASE_METHOD(TestStepCutEnvelopeFixture, "TestStepCutEnvelope::test_type",
                 "[StepCutEnvelope]") {
  CHECK("stepcut" == menvelope->type());
}

TEST_CASE_METHOD(TestStepCutEnvelopeFixture,
                 "TestStepCutEnvelope::test_parentheses_operator",
                 "[StepCutEnvelope]") {
  const PDFEnvelope& fne = *menvelope;
  CHECK(1.0 == fne(-1.0));
  CHECK(1.0 == fne(0.0));
  CHECK(1.0 == fne(+1.0));
  menvelope->setDoubleAttr("stepcut", 1.0);
  CHECK(1.0 == fne(-1.0));
  CHECK(1.0 == fne(0.0));
  CHECK(1.0 == fne(+1.0));
  CHECK(0.0 == fne(+1.0001));
  CHECK(0.0 == fne(+2.0));
}

TEST_CASE_METHOD(TestStepCutEnvelopeFixture,
                 "TestStepCutEnvelope::test_serialization",
                 "[StepCutEnvelope]") {
  menvelope->setDoubleAttr("stepcut", 13.1);
  PDFEnvelopePtr e1 = dumpandload(menvelope);
  CHECK(string("stepcut") == e1->type());
  CHECK(13.1 == e1->getDoubleAttr("stepcut"));
}

// End of file
