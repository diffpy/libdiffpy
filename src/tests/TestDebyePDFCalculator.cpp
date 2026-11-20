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
 * TestDebyePDFCalculator -- unit tests for DebyePDFCalculator class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <functional>
#include <memory>
#include <sstream>
#include <vector>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>
#include <diffpy/serialization.hpp>

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;

class TestDebyePDFCalculatorFixture {
 protected:
  std::shared_ptr<DebyePDFCalculator> mpdfc;
  AtomicStructureAdapterPtr memptystru;
  AtomicStructureAdapterPtr mstru10;
  AtomicStructureAdapterPtr mstru10d1;
  AtomicStructureAdapterPtr mstru10r;
  AtomicStructureAdapterPtr mstru9;
  diffpy::mathutils::EpsilonEqual allclose;
  double meps;

 public:
  TestDebyePDFCalculatorFixture() {
    const int SZ = 10;
    meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
    mpdfc.reset(new DebyePDFCalculator);
    memptystru = std::make_shared<AtomicStructureAdapter>();
    mstru10 = std::make_shared<AtomicStructureAdapter>();
    Atom ai;
    ai.atomtype = "C";
    ai.uij_cartn = R3::identity();
    ai.uij_cartn(0, 0) = ai.uij_cartn(1, 1) = ai.uij_cartn(2, 2) = 0.004;
    for (int i = 0; i < SZ; ++i) {
      ai.xyz_cartn[0] = i;
      mstru10->append(ai);
    }
    mstru10d1 = std::make_shared<AtomicStructureAdapter>(*mstru10);
    (*mstru10d1)[0].atomtype = "Au";
    mstru10r = std::make_shared<AtomicStructureAdapter>();
    mstru10r->assign(mstru10->rbegin(), mstru10->rend());
    mstru9 = std::make_shared<AtomicStructureAdapter>(*mstru10);
    mstru9->erase(9);
  }
};  // class TestDebyePDFCalculatorFixture

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_setQmax",
                 "[DebyePDFCalculator]") {
  const double dq0 = mpdfc->getQstep();
  const double qmax0 = mpdfc->getQmax();
  const double qmax1 = 2 * dq0;
  mpdfc->setQmax(qmax1);
  CHECK(qmax1 / 4.0 == mpdfc->getQstep());
  mpdfc->setQmax(qmax0);
  CHECK(dq0 == mpdfc->getQstep());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_setQstep",
                 "[DebyePDFCalculator]") {
  mpdfc->setRstep(2.0);
  mpdfc->setQstep(0.02);
  CHECK(0.02 == mpdfc->getQstep());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_access_Envelopes",
                 "[DebyePDFCalculator]") {
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

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getPDF",
                 "[DebyePDFCalculator]") {
  QuantityType pdf;
  CHECK(1000u == mpdfc->getPDF().size());
  mpdfc->setRmin(2.0);
  mpdfc->setRmax(0.0);
  mpdfc->eval(memptystru);
  pdf = mpdfc->getPDF();
  CHECK(mpdfc->getPDF().empty());
  mpdfc->setRmax(2.0);
  mpdfc->eval(memptystru);
  pdf = mpdfc->getPDF();
  CHECK(mpdfc->getPDF().empty());
  mpdfc->setRmax(2.01001);
  mpdfc->eval(memptystru);
  pdf = mpdfc->getPDF();
  CHECK(2u == pdf.size());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getRDF",
                 "[DebyePDFCalculator]") {
  // getRDF is not yet implemented and just returns empty array.
  QuantityType rdf = mpdfc->getRDF();
  CHECK(1000u == rdf.size());
  CHECK(0.0 == *std::min_element(rdf.begin(), rdf.end()));
  CHECK(0.0 == *std::max_element(rdf.begin(), rdf.end()));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getF", "[DebyePDFCalculator]") {
  QuantityType fq = mpdfc->getF();
  CHECK(92u == fq.size());
  CHECK(0.0 == *min_element(fq.begin(), fq.end()));
  CHECK(0.0 == *max_element(fq.begin(), fq.end()));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getQgrid",
                 "[DebyePDFCalculator]") {
  CHECK(92u == mpdfc->getQgrid().size());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getQmin",
                 "[DebyePDFCalculator]") {
  CHECK(0.0 == mpdfc->getQmin());
  mpdfc->setQmin(1.0);
  const DebyePDFCalculator& dpc = *mpdfc;
  CHECK(1.0 == dpc.getQmin());
  CHECK(0.0 == dpc.BaseDebyeSum::getQmin());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getQmax",
                 "[DebyePDFCalculator]") {
  CHECK_THAT(25.0, WithinAbs(mpdfc->getQmax(), meps));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getQstep",
                 "[DebyePDFCalculator]") {
  const double rcalchi = ceil((10.0 + 12 * M_PI / 25.0) / 0.01) * 0.01;
  const double qstep = M_PI / rcalchi;
  CHECK_THAT(qstep, WithinAbs(mpdfc->getQstep(), meps));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_getRgrid",
                 "[DebyePDFCalculator]") {
  QuantityType rgrid0 = mpdfc->getRgrid();
  CHECK(rgrid0 == mpdfc->getRgrid());
  mpdfc->setRmin(5);
  mpdfc->setRmax(4);
  CHECK(mpdfc->getRgrid().empty());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_serialization",
                 "[DebyePDFCalculator]") {
  // build customized DebyePDFCalculator
  mpdfc->setPeakWidthModelByType("constant");
  mpdfc->setDoubleAttr("width", 0.123);
  mpdfc->setDoubleAttr("debyeprecision", 0.00011);
  mpdfc->setScatteringFactorTableByType("electronnumber");
  mpdfc->getScatteringFactorTable()->setCustomAs("H", "H", 1.1);
  // dump it to string
  stringstream storage(ios::in | ios::out | ios::binary);
  diffpy::serialization::oarchive oa(storage, ios::binary);
  oa << mpdfc;
  diffpy::serialization::iarchive ia(storage, ios::binary);
  std::shared_ptr<DebyePDFCalculator> pdfc1;
  ia >> pdfc1;
  CHECK(pdfc1.get() != mpdfc.get());
  CHECK(string("constant") == pdfc1->getPeakWidthModel()->type());
  CHECK(0.123 == pdfc1->getDoubleAttr("width"));
  CHECK(0.00011 == pdfc1->getDoubleAttr("debyeprecision"));
  CHECK(string("electronnumber") == pdfc1->getScatteringFactorTable()->type());
  CHECK(1.1 == pdfc1->getScatteringFactorTable()->lookup("H"));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_change_atom",
                 "[DebyePDFCalculator]") {
  using std::placeholders::_1;
  mpdfc->setQmin(1.0);
  DebyePDFCalculator pdfcb = *mpdfc;
  DebyePDFCalculator pdfco = *mpdfc;
  CHECK(1.0 == pdfcb.getQmin());
  pdfcb.setEvaluatorType(BASIC);
  pdfco.setEvaluatorType(OPTIMIZED);
  CHECK(NONE == pdfcb.getEvaluatorTypeUsed());
  CHECK(NONE == pdfco.getEvaluatorTypeUsed());
  pdfcb.eval(mstru10);
  pdfco.eval(mstru10);
  QuantityType gb = pdfcb.getPDF();
  QuantityType go = pdfco.getPDF();
  CHECK(!gb.empty());
  CHECK(gb == go);
  int cnonzero =
    count_if(gb.begin(), gb.end(), bind(not_equal_to<double>(), _1, 0.0));
  CHECK(cnonzero);
  CHECK(BASIC == pdfcb.getEvaluatorType());
  CHECK(OPTIMIZED == pdfco.getEvaluatorType());
  // first call of pdfco should use the BASIC evaluation
  CHECK(BASIC == pdfcb.getEvaluatorTypeUsed());
  CHECK(BASIC == pdfco.getEvaluatorTypeUsed());
  // test second call on the same structure
  pdfco.eval(mstru10);
  go = pdfco.getPDF();
  CHECK(gb == go);
  CHECK(OPTIMIZED == pdfco.getEvaluatorTypeUsed());
  // test structure with one different atom
  pdfcb.eval(mstru10d1);
  pdfco.eval(mstru10d1);
  QuantityType gb1 = pdfcb.getPDF();
  QuantityType go1 = pdfco.getPDF();
  CHECK(!allclose(gb, gb1));
  CHECK(OPTIMIZED == pdfco.getEvaluatorTypeUsed());
  CHECK(allclose(gb1, go1));
  // change position of 1 atom
  mstru10d1->at(0).xyz_cartn[1] = 0.5;
  pdfcb.eval(mstru10d1);
  pdfco.eval(mstru10d1);
  QuantityType gb2 = pdfcb.getPDF();
  QuantityType go2 = pdfco.getPDF();
  CHECK(!allclose(gb1, gb2));
  CHECK(OPTIMIZED == pdfco.getEvaluatorTypeUsed());
  CHECK(allclose(gb2, go2));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_reverse_atoms",
                 "[DebyePDFCalculator]") {
  mpdfc->setQmin(1.1);
  DebyePDFCalculator pdfcb = *mpdfc;
  DebyePDFCalculator pdfco = *mpdfc;
  pdfcb.setEvaluatorType(BASIC);
  pdfco.setEvaluatorType(OPTIMIZED);
  pdfcb.eval(mstru10);
  pdfco.eval(mstru10);
  pdfco.eval(mstru10r);
  QuantityType gb = pdfcb.getPDF();
  QuantityType go = pdfco.getPDF();
  CHECK(OPTIMIZED == pdfco.getEvaluatorTypeUsed());
  CHECK(allclose(gb, go));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_remove_atom",
                 "[DebyePDFCalculator]") {
  mpdfc->setQmin(1.2);
  DebyePDFCalculator pdfcb = *mpdfc;
  DebyePDFCalculator pdfco = *mpdfc;
  pdfcb.setEvaluatorType(BASIC);
  pdfco.setEvaluatorType(OPTIMIZED);
  pdfcb.eval(mstru9);
  pdfco.eval(mstru10);
  pdfco.eval(mstru9);
  QuantityType gb = pdfcb.getPDF();
  QuantityType go = pdfco.getPDF();
  CHECK(OPTIMIZED == pdfco.getEvaluatorTypeUsed());
  CHECK(allclose(gb, go));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_qmin_click",
                 "[DebyePDFCalculator]") {
  mpdfc->setEvaluatorType(OPTIMIZED);
  CHECK(NONE == mpdfc->getEvaluatorTypeUsed());
  mpdfc->eval(mstru10);
  mpdfc->eval(mstru10);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  mpdfc->setQmin(4);
  mpdfc->eval(mstru10);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  mpdfc->setQmin(0);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_qmax_click",
                 "[DebyePDFCalculator]") {
  mpdfc->setEvaluatorType(OPTIMIZED);
  mpdfc->eval(mstru10);
  mpdfc->eval(mstru10);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  mpdfc->setQmax(mpdfc->getQmax());
  mpdfc->eval(mstru10);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  mpdfc->setQmax(mpdfc->getQmax() - 1.0);
  mpdfc->eval(mstru10);
  CHECK(BASIC == mpdfc->getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_delta12_click",
                 "[DebyePDFCalculator]") {
  mpdfc->setEvaluatorType(OPTIMIZED);
  mpdfc->eval(mstru10);
  mpdfc->setDoubleAttr("delta1", mpdfc->getDoubleAttr("delta1"));
  mpdfc->setDoubleAttr("delta2", mpdfc->getDoubleAttr("delta2"));
  mpdfc->eval(mstru10);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  mpdfc->setDoubleAttr("delta1", 0.5 + mpdfc->getDoubleAttr("delta1"));
  mpdfc->eval(mstru10);
  CHECK(BASIC == mpdfc->getEvaluatorTypeUsed());
  mpdfc->eval(mstru10);
  mpdfc->setDoubleAttr("delta2", 0.5 + mpdfc->getDoubleAttr("delta2"));
  mpdfc->eval(mstru10);
  CHECK(BASIC == mpdfc->getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_SFTB_click",
                 "[DebyePDFCalculator]") {
  mpdfc->setEvaluatorType(OPTIMIZED);
  mpdfc->eval(mstru10);
  mpdfc->setScatteringFactorTable(mpdfc->getScatteringFactorTable());
  mpdfc->eval(mstru10);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  mpdfc->setScatteringFactorTableByType("neutron");
  mpdfc->eval(mstru10);
  CHECK(BASIC == mpdfc->getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_site_occupancies",
                 "[DebyePDFCalculator]") {
  // Sites of one element with different occupancies scatter as sites of
  // different atom types with those occupancies.
  AtomicStructureAdapterPtr mixed =
    std::make_shared<AtomicStructureAdapter>(*mstru10);
  AtomicStructureAdapterPtr named =
    std::make_shared<AtomicStructureAdapter>(*mstru10);
  for (int i : {0, 3, 4, 8}) {
    mixed->at(i).occupancy = 0.4;
    named->at(i).occupancy = 0.4;
    named->at(i).atomtype = "C2";
  }
  mpdfc->getScatteringFactorTable()->setCustomAs("C2", "C");
  for (PQEvaluatorType tp : {BASIC, OPTIMIZED}) {
    INFO("evaluator " << tp);
    mpdfc->setEvaluatorType(tp);
    mpdfc->eval(named);
    const QuantityType gnamed = mpdfc->getPDF();
    const QuantityType fnamed = mpdfc->getF();
    mpdfc->eval(mstru10);
    CHECK(!allclose(gnamed, mpdfc->getPDF()));
    mpdfc->eval(mixed);
    CHECK(allclose(gnamed, mpdfc->getPDF()));
    CHECK(allclose(fnamed, mpdfc->getF()));
  }
  // A new occupancy of the first site of an element is a change of that
  // site only, which the fast update must get right.
  mpdfc->setEvaluatorType(OPTIMIZED);
  mpdfc->eval(mixed);
  mixed->at(0).occupancy = 0.7;
  named->at(0).occupancy = 0.7;
  mpdfc->eval(mixed);
  CHECK(OPTIMIZED == mpdfc->getEvaluatorTypeUsed());
  const QuantityType gfast = mpdfc->getPDF();
  mpdfc->setEvaluatorType(BASIC);
  mpdfc->eval(named);
  CHECK(allclose(gfast, mpdfc->getPDF()));
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_retry_unknown_atom",
                 "[DebyePDFCalculator]") {
  const QuantityType expected = mpdfc->eval(mstru10);
  mstru10->at(1).atomtype = "Bogus";
  CHECK_THROWS_AS(mpdfc->eval(mstru10), invalid_argument);
  CHECK_THROWS_AS(mpdfc->eval(mstru10), invalid_argument);
  mstru10->at(1).atomtype = "C";
  CHECK(mpdfc->eval(mstru10) == expected);
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_replace_table_by_reference",
                 "[DebyePDFCalculator]") {
  DebyePDFCalculator reference;
  reference.setScatteringFactorTableByType("neutron");
  const auto neutron = reference.eval(mstru10);
  for (const auto evaluator : {BASIC, OPTIMIZED}) {
    DebyePDFCalculator calc;
    calc.setEvaluatorType(evaluator);
    auto& table = calc.getScatteringFactorTable();
    REQUIRE(calc.eval(mstru10) != neutron);
    table = reference.getScatteringFactorTable();
    CHECK(calc.eval(mstru10) == neutron);
    CHECK(calc.getEvaluatorTypeUsed() == BASIC);
  }
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_DBPDF_kept_scattering_factors",
                 "[DebyePDFCalculator]") {
  // The kept scattering factors must follow every change of the calculator
  // and of the occupancies.
  typedef std::function<void(DebyePDFCalculator&)> Config;
  std::vector<Config> history;
  auto configure = [&](Config cfg) {
    cfg(*mpdfc);
    history.push_back(cfg);
  };
  // the calculator that is checked
  DebyePDFCalculator* calc = mpdfc.get();
  auto check = [&](StructureAdapterPtr stru) {
    DebyePDFCalculator fresh;
    for (const Config& cfg : history) cfg(fresh);
    fresh.setEvaluatorType(BASIC);
    fresh.eval(stru);
    calc->eval(stru);
    CHECK(allclose(fresh.getPDF(), calc->getPDF()));
    CHECK(allclose(fresh.getF(), calc->getF()));
  };
  check(mstru10);
  check(mstru10d1);
  configure(
    [](DebyePDFCalculator& c) { c.setScatteringFactorTableByType("neutron"); });
  check(mstru10d1);
  configure([](DebyePDFCalculator& c) {
    c.getScatteringFactorTable()->setCustomAs("C", "Au");
  });
  check(mstru10d1);
  configure([](DebyePDFCalculator& c) { c.setQmax(15.0); });
  check(mstru10d1);
  configure([](DebyePDFCalculator& c) { c.setRstep(0.02); });
  check(mstru10d1);
  // Many occupancies in turn, which leave arrays that no site uses.
  AtomicStructureAdapterPtr stru =
    std::make_shared<AtomicStructureAdapter>(*mstru10d1);
  // a calculator of a derived class adds its pairs one by one
  struct DerivedDebyePDFCalculator : DebyePDFCalculator {};
  DerivedDebyePDFCalculator onebyone;
  for (const Config& cfg : history) cfg(onebyone);
  for (DebyePDFCalculator* c : {(DebyePDFCalculator*)&onebyone, mpdfc.get()}) {
    const bool pairbypair = (c == &onebyone);
    calc = c;
    check(stru);
    for (int k = 0; k < 40; ++k) {
      INFO((pairbypair ? "pair by pair" : "all pairs") << " step " << k);
      stru->at(3).occupancy = 0.5 + 0.01 * k;
      stru->at(7).occupancy = 0.9 - 0.01 * k;
      check(stru);
      if (pairbypair) CHECK(OPTIMIZED == calc->getEvaluatorTypeUsed());
    }
  }
  calc = mpdfc.get();
  // a restored calculator builds its scattering factors again
  stringstream storage(ios::in | ios::out | ios::binary);
  diffpy::serialization::oarchive oa(storage, ios::binary);
  oa << mpdfc;
  diffpy::serialization::iarchive ia(storage, ios::binary);
  std::shared_ptr<DebyePDFCalculator> pdfc1;
  ia >> pdfc1;
  const QuantityType g = mpdfc->getPDF();
  stru->at(3).occupancy = 0.3;
  pdfc1->eval(stru);
  CHECK(!allclose(g, pdfc1->getPDF()));
  mpdfc.swap(pdfc1);
  check(stru);
}

TEST_CASE_METHOD(TestDebyePDFCalculatorFixture,
                 "TestDebyePDFCalculator::test_structure_dependent_scattering",
                 "[DebyePDFCalculator]") {
  struct SiteFactor : DebyePDFCalculator {
    double sfSiteAtQ(int site, const double&) const override {
      // Both site-specific and structure-wide dependence.
      return mstructure->siteCartesianUij(site)(0, 0) +
             mstructure->siteCartesianUij(0)(1, 1);
    }
  };
  for (auto evaluator : {BASIC, OPTIMIZED}) {
    SiteFactor calc;
    calc.setEvaluatorType(evaluator);
    auto stru = std::make_shared<AtomicStructureAdapter>(*mstru10);
    for (double u : {0.004, 0.05, 0.02, 0.004}) {
      (*stru)[0].uij_cartn = u * R3::identity();
      // Give the reference sites unique keys. It must always evaluate each
      // site's curve, including sites whose own parameters stayed unchanged.
      auto named = std::make_shared<AtomicStructureAdapter>(*stru);
      for (int i = 0; i < named->countSites(); ++i)
        (*named)[i].atomtype = std::to_string(i);
      SiteFactor fresh;
      fresh.setEvaluatorType(BASIC);
      CHECK(allclose(calc.eval(stru), fresh.eval(named)));
    }
  }
}

// End of file
