/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2019 Brookhaven Science Associates,
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
 * TestPDFBaseline -- test PDFBaseline and derived classes
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/PDFBaseline.hpp>
#include <diffpy/srreal/LinearBaseline.hpp>
#include <diffpy/srreal/ZeroBaseline.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// TestPDFBaseline
//////////////////////////////////////////////////////////////////////////////

class TestPDFBaselineFixture {
 protected:
  // data
  PDFBaselinePtr ml;
  LinearBaseline* mcl;
  PDFBaselinePtr mz;
  ZeroBaseline* mcz;

 public:
  TestPDFBaselineFixture() {
    ml = PDFBaseline::createByType("linear");
    mcl = static_cast<LinearBaseline*>(ml.get());
    mz = PDFBaseline::createByType("zero");
    mcz = static_cast<ZeroBaseline*>(mz.get());
  }
};  // class TestPDFBaselineFixture

TEST_CASE_METHOD(TestPDFBaselineFixture, "TestPDFBaseline::test_clone",
                 "[PDFBaseline]") {
  CHECK("linear" == ml->clone()->type());
  CHECK("zero" == mz->clone()->type());
}

TEST_CASE_METHOD(TestPDFBaselineFixture, "TestPDFBaseline::test_attributes",
                 "[PDFBaseline]") {
  CHECK(1 == ml->namesOfDoubleAttributes().size());
  CHECK(0 == mz->namesOfDoubleAttributes().size());
  mcl->setSlope(0.1);
  CHECK(0.1 == ml->getDoubleAttr("slope"));
  ml->setDoubleAttr("slope", 0.2);
  CHECK(0.2 == mcl->getSlope());
}

TEST_CASE_METHOD(TestPDFBaselineFixture, "TestPDFBaseline::test_calculate",
                 "[PDFBaseline]") {
  CHECK(0.0 == (*ml)(2.0));
  CHECK(0.0 == (*mz)(2.0));
  mcl->setSlope(-2);
  CHECK(-4.0 == (*ml)(2.0));
}

TEST_CASE_METHOD(TestPDFBaselineFixture, "TestPDFBaseline::test_serialization",
                 "[PDFBaseline]") {
  mcl->setSlope(0.3);
  PDFBaselinePtr bl1 = dumpandload(ml);
  LinearBaseline* cbl1 = dynamic_cast<LinearBaseline*>(bl1.get());
  CHECK(cbl1);
  CHECK(0.3 == cbl1->getSlope());
  PDFBaselinePtr bz1 = dumpandload(mz);
  ZeroBaseline* cbz1 = dynamic_cast<ZeroBaseline*>(bz1.get());
  CHECK(cbz1);
}

}  // namespace srreal
}  // namespace diffpy

// End of file
