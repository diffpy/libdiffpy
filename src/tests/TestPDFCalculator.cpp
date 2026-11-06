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
 * TestPDFCalculator -- unit tests for PDFCalculator class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>
#include <diffpy/serialization.hpp>
#include "test_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;

class TestPDFCalculatorFixture {
 protected:
  std::shared_ptr<PDFCalculator> mpdfc;
  StructureAdapterPtr memptystru;
  double meps;
  double mepsdb;

 public:
  TestPDFCalculatorFixture() {
    memptystru = emptyStructureAdapter();
    meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
    mepsdb = 10 * diffpy::mathutils::DOUBLE_EPS;
    mpdfc.reset(new PDFCalculator);
  }
};  // class TestPDFCalculatorFixture

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_setPeakWidthModel",
                 "[PDFCalculator]") {
  const PeakWidthModelPtr& jpw0 = mpdfc->getPeakWidthModel();
  CHECK(0.0 == jpw0->getDoubleAttr("delta1"));
  CHECK(0.0 == jpw0->getDoubleAttr("delta2"));
  CHECK(0.0 == jpw0->getDoubleAttr("qbroad"));
  CHECK(0.0 == jpw0->getDoubleAttr("qbroad_seperable"));
  JeongPeakWidth jpw;
  jpw.setDelta1(1.0);
  jpw.setDelta2(2.0);
  jpw.setQbroad(3.0);
  jpw.setQbroad_seperable(4.0);
  mpdfc->setPeakWidthModel(jpw.clone());
  const PeakWidthModelPtr& jpw1 = mpdfc->getPeakWidthModel();
  CHECK(1.0 == jpw1->getDoubleAttr("delta1"));
  CHECK(2.0 == jpw1->getDoubleAttr("delta2"));
  CHECK(3.0 == jpw1->getDoubleAttr("qbroad"));
  CHECK(4.0 == jpw1->getDoubleAttr("qbroad_seperable"));
}

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_getPeakWidthModel",
                 "[PDFCalculator]") {
  string tp = "jeong";
  CHECK(tp == mpdfc->getPeakWidthModel()->type());
  PeakWidthModelPtr pwm(PeakWidthModel::createByType("debye-waller"));
  mpdfc->setPeakWidthModel(pwm);
  tp = "debye-waller";
  CHECK(tp == mpdfc->getPeakWidthModel()->type());
  mpdfc->setPeakWidthModel(ConstantPeakWidth().clone());
  tp = "constant";
  CHECK(tp == mpdfc->getPeakWidthModel()->type());
}

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_access_PeakProfile",
                 "[PDFCalculator]") {
  PeakProfilePtr pkf(mpdfc->getPeakProfile()->clone());
  pkf->setPrecision(1.1);
  mpdfc->setPeakProfile(pkf);
  CHECK(1.1 == mpdfc->getPeakProfile()->getPrecision());
  CHECK(1.1 == mpdfc->getDoubleAttr("peakprecision"));
  mpdfc->setDoubleAttr("peakprecision", 0.2);
  CHECK(0.2 == mpdfc->getDoubleAttr("peakprecision"));
  CHECK(pkf->type() == mpdfc->getPeakProfile()->type());
  mpdfc->setPeakProfileByType("gaussian");
  CHECK(0.2 == mpdfc->getDoubleAttr("peakprecision"));
  CHECK_THROWS_AS(mpdfc->setPeakProfileByType("invalid"), logic_error);
}

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_access_Envelopes",
                 "[PDFCalculator]") {
  CHECK(2u == mpdfc->usedEnvelopeTypes().size());
  CHECK(1.0 == mpdfc->getDoubleAttr("scale"));
  CHECK(0.0 == mpdfc->getDoubleAttr("qdamp"));
  mpdfc->setDoubleAttr("scale", 3.0);
  CHECK(3.0 == mpdfc->getDoubleAttr("scale"));
  mpdfc->addEnvelopeByType("scale");
  CHECK(1.0 == mpdfc->getDoubleAttr("scale"));
  QResolutionEnvelope qdamp4;
  qdamp4.setQdamp(4);
  mpdfc->addEnvelope(qdamp4.clone());
  CHECK(4.0 == mpdfc->getDoubleAttr("qdamp"));
  CHECK_THROWS_AS(mpdfc->addEnvelopeByType("invalid"), logic_error);
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getPDF",
                 "[PDFCalculator]") {
  QuantityType pdf;
  pdf = mpdfc->getPDF();
  CHECK(1000u == pdf.size());
  mpdfc->setRmin(2.0);
  mpdfc->setRmax(0.0);
  mpdfc->eval(memptystru);
  pdf = mpdfc->getPDF();
  CHECK(pdf.empty());
  mpdfc->setRmax(2.0);
  mpdfc->eval(memptystru);
  pdf = mpdfc->getPDF();
  CHECK(pdf.empty());
  mpdfc->setRmax(2.01001);
  mpdfc->eval(memptystru);
  pdf = mpdfc->getPDF();
  CHECK(2u == pdf.size());
  CHECK(0.0 == pdf[0]);
  CHECK(0.0 == pdf[1]);
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getRDF",
                 "[PDFCalculator]") {
  QuantityType rdf = mpdfc->getRDF();
  CHECK(1000u == rdf.size());
  CHECK(0.0 == *min_element(rdf.begin(), rdf.end()));
  CHECK(0.0 == *max_element(rdf.begin(), rdf.end()));
  mpdfc->eval(memptystru);
  rdf = mpdfc->getRDF();
  CHECK(1000u == rdf.size());
  CHECK(0.0 == *min_element(rdf.begin(), rdf.end()));
  CHECK(0.0 == *max_element(rdf.begin(), rdf.end()));
  mpdfc->setRmin(2.0);
  mpdfc->setRmax(0.0);
  mpdfc->eval(memptystru);
  rdf = mpdfc->getRDF();
  CHECK(rdf.empty());
  mpdfc->setRmax(2.01001);
  mpdfc->eval(memptystru);
  rdf = mpdfc->getRDF();
  CHECK(2u == rdf.size());
  CHECK(0.0 == rdf[0]);
  CHECK(0.0 == rdf[1]);
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getF",
                 "[PDFCalculator]") {
  QuantityType fq = mpdfc->getF();
  CHECK(1024u == fq.size());
  CHECK(0.0 == *min_element(fq.begin(), fq.end()));
  CHECK(0.0 == *max_element(fq.begin(), fq.end()));
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getQgrid",
                 "[PDFCalculator]") {
  CHECK(1024u == mpdfc->getQgrid().size());
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getQmin",
                 "[PDFCalculator]") {
  CHECK(0.0 == mpdfc->getQmin());
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getQmax",
                 "[PDFCalculator]") {
  CHECK_THAT(100 * M_PI, WithinAbs(mpdfc->getQmax(), meps));
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getQstep",
                 "[PDFCalculator]") {
  const double qstep0 = 100 * M_PI / 1024;
  const double qstep1 = 100 * M_PI / 2048;
  const double qstep2 = 100 * M_PI / 4096;
  const double qstep3 = M_PI / (1024 * 0.03);
  CHECK_THAT(qstep0, WithinAbs(mpdfc->getQstep(), meps));
  mpdfc->setRmax(20);
  CHECK_THAT(qstep1, WithinAbs(mpdfc->getQstep(), meps));
  mpdfc->setQmax(10);
  CHECK_THAT(qstep2, WithinAbs(mpdfc->getQstep(), meps));
  mpdfc->setRstep(0.03);
  CHECK_THAT(qstep3, WithinAbs(mpdfc->getQstep(), meps));
  // test qstep after evaluation of some structure.
  StructureAdapterPtr catio3;
  catio3 = loadTestPeriodicStructure("CaTiO3.stru");
  mpdfc->eval(catio3);
  CHECK_THAT(qstep3, WithinAbs(mpdfc->getQstep(), meps));
  // now qstep is not updated because of non-trivial structure.
  mpdfc->setRstep(0.01);
  CHECK(qstep2 != mpdfc->getQstep());
  // qstep is correct after calculation.
  mpdfc->eval();
  CHECK_THAT(qstep2, WithinAbs(mpdfc->getQstep(), meps));
}

TEST_CASE_METHOD(TestPDFCalculatorFixture, "TestPDFCalculator::test_getRgrid",
                 "[PDFCalculator]") {
  QuantityType rgrid0 = mpdfc->getRgrid();
  CHECK(rgrid0 == mpdfc->getRgrid());
  mpdfc->setRmin(5);
  mpdfc->setRmax(4);
  CHECK(mpdfc->getRgrid().empty());
}

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_getExtendedRmin", "[PDFCalculator]") {
  // empty structure should not extend the grid at all.
  CHECK(0.0 == mpdfc->getExtendedRmin());
  mpdfc->setRmin(5);
  mpdfc->eval(memptystru);
  CHECK_THAT(5.0, WithinAbs(mpdfc->getExtendedRmin(), mepsdb));
}

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_getExtendedRmax", "[PDFCalculator]") {
  // empty structure should not extend the grid at all.
  CHECK_THAT(10.0, WithinAbs(mpdfc->getExtendedRmax(), mepsdb));
  mpdfc->setRmax(7);
  mpdfc->eval(memptystru);
  CHECK_THAT(7.0, WithinAbs(mpdfc->getExtendedRmax(), mepsdb));
}

TEST_CASE_METHOD(TestPDFCalculatorFixture,
                 "TestPDFCalculator::test_serialization", "[PDFCalculator]") {
  // build customized PDFCalculator
  mpdfc->setPeakWidthModelByType("constant");
  mpdfc->setDoubleAttr("width", 0.123);
  mpdfc->setPeakProfileByType("gaussian");
  mpdfc->setDoubleAttr("peakprecision", 0.011);
  mpdfc->setScatteringFactorTableByType("electronnumber");
  mpdfc->getScatteringFactorTable()->setCustomAs("H", "H", 1.1);
  // dump it to string
  stringstream storage(ios::in | ios::out | ios::binary);
  diffpy::serialization::oarchive oa(storage, ios::binary);
  oa << mpdfc;
  diffpy::serialization::iarchive ia(storage, ios::binary);
  std::shared_ptr<PDFCalculator> pdfc1;
  ia >> pdfc1;
  CHECK(pdfc1.get() != mpdfc.get());
  CHECK(string("constant") == pdfc1->getPeakWidthModel()->type());
  CHECK(0.123 == pdfc1->getDoubleAttr("width"));
  CHECK(0.011 == pdfc1->getDoubleAttr("peakprecision"));
  CHECK(string("electronnumber") == pdfc1->getScatteringFactorTable()->type());
  CHECK(1.1 == pdfc1->getScatteringFactorTable()->lookup("H"));
}

// End of file
