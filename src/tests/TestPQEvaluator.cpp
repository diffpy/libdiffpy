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
 * TestPQEvaluator -- unit tests for the PairQuantity evaluator class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <functional>
#include <memory>

#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/PairCounter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/OverlapCalculator.hpp>
#include "test_helpers.hpp"

namespace diffpy {
namespace srreal {

using namespace std;
using diffpy::mathutils::EpsilonEqual;

// calculator with invalid support for OPTIMIZED evaluation

class BadPairCounter : public PairCounter {
 protected:
  // intentionally faulty support for PQEvaluatorOptimized
  virtual void stashPartialValue() {}
  virtual void restorePartialValue() {}
};

// calculator that fails after a given number of pair contributions

class FailingPDFCalculator : public PDFCalculator {
 public:
  /// Pair contributions left before the failure, negative for no failure.
  int countdown = -1;
  /// Pair contributions added so far.
  int calls = 0;

 protected:
  virtual void addPairContribution(const BaseBondGenerator& bnds,
                                   int summationscale) {
    if (countdown == 0) throw runtime_error("failed pair contribution");
    if (countdown > 0) --countdown;
    ++calls;
    this->PDFCalculator::addPairContribution(bnds, summationscale);
  }
};

// Crystal with an inversion center at the origin.  The first atom sits
// on it at x = 0 and has two symmetry positions elsewhere.

CrystalStructureAdapterPtr inversionCrystal(double x) {
  CrystalStructureAdapterPtr rv(new CrystalStructureAdapter);
  rv->setLatPar(4.0, 4.5, 5.0, 90.0, 90.0, 90.0);
  rv->addSymOp(R3::identity(), R3::zerovector);
  rv->addSymOp(-1.0 * R3::identity(), R3::zerovector);
  const double positions[][3] = {{x, 0.0, 0.0},
                                 {0.9, 0.6, 1.3},
                                 {1.7, 1.2, 0.5},
                                 {0.4, 1.9, 2.1},
                                 {1.3, 0.3, 1.9}};
  Atom a;
  a.atomtype = "Ni";
  a.uij_cartn = 0.005 * R3::identity();
  for (const auto& xyz : positions) {
    a.xyz_cartn = R3::Vector(xyz[0], xyz[1], xyz[2]);
    rv->append(a);
  }
  return rv;
}

//////////////////////////////////////////////////////////////////////////////
// TestPQEvaluator
//////////////////////////////////////////////////////////////////////////////

class TestPQEvaluatorFixture {
 protected:
  // data
  EpsilonEqual allclose;
  PDFCalculator mpdfcb;
  PDFCalculator mpdfco;
  QuantityType mzeros;
  AtomicStructureAdapterPtr mstru10;
  AtomicStructureAdapterPtr mstru10d1;
  AtomicStructureAdapterPtr mstru10r;
  AtomicStructureAdapterPtr mstru9;

  // methods
  QuantityType pdfcdiff(StructureAdapterPtr stru) {
    mpdfcb.eval(stru);
    mpdfco.eval(stru);
    QuantityType gb = mpdfcb.getPDF();
    QuantityType go = mpdfco.getPDF();
    assert(mpdfcb.getRgrid() == mpdfco.getRgrid());
    QuantityType rv(gb.size());
    transform(gb.begin(), gb.end(), go.begin(), rv.begin(), minus<double>());
    return rv;
  }

 public:
  TestPQEvaluatorFixture() {
    // setup PDFCalculator instances and the zero vector for comparison
    mpdfcb.setEvaluatorType(BASIC);
    mpdfco.setEvaluatorType(OPTIMIZED);
    mzeros.assign(mpdfcb.getRgrid().size(), 0.0);
    // setup structure instances
    const int SZ = 10;
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
};  // class TestPQEvaluatorFixture

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_initial_evaluator_type_used",
                 "[PQEvaluator]") {
  PDFCalculator pdfc;
  CHECK(NONE == pdfc.getEvaluatorTypeUsed());
  pdfc.setEvaluatorType(OPTIMIZED);
  CHECK(NONE == pdfc.getEvaluatorTypeUsed());
  pdfc.setEvaluatorType(PQEvaluatorType::CHECK);
  CHECK(NONE == pdfc.getEvaluatorTypeUsed());
  pdfc.setEvaluatorType(BASIC);
  CHECK(NONE == pdfc.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_PDF_change_atom", "[PQEvaluator]") {
  using std::placeholders::_1;
  CHECK(BASIC == mpdfcb.getEvaluatorType());
  CHECK(OPTIMIZED == mpdfco.getEvaluatorType());
  CHECK(mzeros == this->pdfcdiff(mstru10));
  // first call of mpdfco should use the BASIC evaluation
  CHECK(BASIC == mpdfcb.getEvaluatorTypeUsed());
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
  // verify there are some non-zero values in the PDF.
  QuantityType gb = mpdfcb.getPDF();
  CHECK(!gb.empty());
  int cnonzero =
    count_if(gb.begin(), gb.end(), bind(not_equal_to<double>(), _1, 0.0));
  CHECK(cnonzero);
  // test second call on the same structure
  mpdfco.eval(mstru10);
  CHECK(gb == mpdfco.getPDF());
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  // test structure with one different atom
  CHECK(allclose(mzeros, this->pdfcdiff(mstru10d1)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb1 = mpdfcb.getPDF();
  CHECK(!allclose(gb, gb1));
  // change position of 1 atom
  mstru10d1->at(0).xyz_cartn[1] = 0.5;
  CHECK(allclose(mzeros, this->pdfcdiff(mstru10d1)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb2 = mpdfcb.getPDF();
  CHECK(!allclose(gb1, gb2));
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_PDF_reverse_atoms", "[PQEvaluator]") {
  mpdfco.eval(mstru10);
  CHECK(allclose(mzeros, this->pdfcdiff(mstru10r)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_PDF_remove_atom", "[PQEvaluator]") {
  mpdfco.eval(mstru10);
  CHECK(allclose(mzeros, this->pdfcdiff(mstru9)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture, "TestPQEvaluator::test_PDF_type_mask",
                 "[PQEvaluator]") {
  mpdfcb.setTypeMask("O2-", "all", false);
  mpdfco.setTypeMask("O2-", "all", false);
  PeriodicStructureAdapterPtr litao =
    std::dynamic_pointer_cast<PeriodicStructureAdapter>(
      loadTestPeriodicStructure("LiTaO3.stru"));
  CHECK(mzeros == this->pdfcdiff(litao));
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb0 = mpdfcb.getPDF();
  // change masked-away oxygen
  Atom& o29 = litao->at(29);
  CHECK("O2-" == o29.atomtype);
  o29.xyz_cartn = R3::Vector(0.1, 0.2, 0.3);
  CHECK(mzeros == this->pdfcdiff(litao));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb1 = mpdfcb.getPDF();
  CHECK(gb0 == gb1);
  // change active lithium
  Atom& li0 = litao->at(0);
  li0.occupancy = 0.1;
  CHECK(allclose(mzeros, this->pdfcdiff(litao)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb2 = mpdfcb.getPDF();
  CHECK(!allclose(gb0, gb2));
}

TEST_CASE_METHOD(TestPQEvaluatorFixture, "TestPQEvaluator::test_PDF_pair_mask",
                 "[PQEvaluator]") {
  mpdfcb.setPairMask(0, 3, false);
  mpdfco.setPairMask(0, 3, false);
  CHECK(mzeros == this->pdfcdiff(mstru10));
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb0 = mpdfcb.getPDF();
  // check structure with one changed atom
  CHECK(allclose(mzeros, this->pdfcdiff(mstru10d1)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb1 = mpdfcb.getPDF();
  CHECK(!allclose(gb0, gb1));
  // check structure with one removed atom
  CHECK(allclose(mzeros, this->pdfcdiff(mstru9)));
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb2 = mpdfcb.getPDF();
  CHECK(!allclose(gb1, gb2));
  // check full evaluation for structure with reversed atoms
  mpdfco.eval(mstru10);
  CHECK(mzeros == this->pdfcdiff(mstru10r));
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
  QuantityType gb3 = mpdfcb.getPDF();
  CHECK(allclose(gb0, gb3));
  // check effect of pair mask updates
  mpdfco.eval(mstru10r);
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  mpdfco.setPairMask(0, 3, false);
  mpdfco.eval(mstru10r);
  CHECK(OPTIMIZED == mpdfco.getEvaluatorTypeUsed());
  mpdfco.setPairMask(0, 3, true);
  mpdfco.eval(mstru10r);
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_optimized_supported", "[PQEvaluator]") {
  mpdfcb.eval(mstru10);
  CHECK(BASIC == mpdfcb.getEvaluatorTypeUsed());
  mpdfcb.setEvaluatorType(OPTIMIZED);
  CHECK(BASIC == mpdfcb.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_optimized_unsupported",
                 "[PQEvaluator]") {
  OverlapCalculator olc;
  CHECK(BASIC == olc.getEvaluatorType());
  CHECK_THROWS_AS(olc.setEvaluatorType(OPTIMIZED), invalid_argument);
  CHECK(BASIC == olc.getEvaluatorType());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_checkevaluator_unsupported",
                 "[PQEvaluator]") {
  OverlapCalculator olc;
  CHECK(BASIC == olc.getEvaluatorType());
  CHECK_THROWS_AS(olc.setEvaluatorType(PQEvaluatorType::CHECK),
                  invalid_argument);
  CHECK(BASIC == olc.getEvaluatorType());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_checkevaluator_passes",
                 "[PQEvaluator]") {
  mpdfco.setEvaluatorType(PQEvaluatorType::CHECK);
  CHECK(PQEvaluatorType::CHECK == mpdfco.getEvaluatorType());
  mpdfco.eval(mstru10);
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
  mpdfco.eval(mstru10);
  CHECK(PQEvaluatorType::CHECK == mpdfco.getEvaluatorTypeUsed());
  mpdfco.eval(mstru9);
  CHECK(PQEvaluatorType::CHECK == mpdfco.getEvaluatorTypeUsed());
  mpdfco.eval(emptyStructureAdapter());
  CHECK(BASIC == mpdfco.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_checkevaluator_fails",
                 "[PQEvaluator]") {
  BadPairCounter badcounter;
  badcounter.setEvaluatorType(PQEvaluatorType::CHECK);
  CHECK(45 == badcounter(mstru10));
  CHECK(BASIC == badcounter.getEvaluatorTypeUsed());
  CHECK_THROWS_AS(badcounter(mstru10), logic_error);
  CHECK(PQEvaluatorType::CHECK == badcounter.getEvaluatorTypeUsed());
}

TEST_CASE_METHOD(TestPQEvaluatorFixture, "TestPQEvaluator::test_failed_update",
                 "[PQEvaluator]") {
  // A value left by a failed update must not pass for a valid one.
  mpdfcb.eval(mstru10);
  const QuantityType v10 = mpdfcb.value();
  mpdfcb.eval(mstru10d1);
  const QuantityType v10d1 = mpdfcb.value();
  // count the pair contributions in a fast update
  FailingPDFCalculator counter;
  counter.eval(mstru10);
  counter.calls = 0;
  counter.eval(mstru10d1);
  REQUIRE(OPTIMIZED == counter.getEvaluatorTypeUsed());
  const int nfast = counter.calls;
  REQUIRE(nfast > 4);
  FailingPDFCalculator pdfc;
  SECTION("full calculation") {
    // all atoms are different, which needs a full calculation
    AtomicStructureAdapterPtr shifted =
      std::make_shared<AtomicStructureAdapter>(*mstru10);
    for (Atom& a : *shifted) a.xyz_cartn[1] += 0.5;
    pdfc.eval(mstru10);
    pdfc.countdown = 3;
    CHECK_THROWS_AS(pdfc.eval(shifted), runtime_error);
    pdfc.countdown = -1;
    pdfc.eval(mstru10);
    CHECK(BASIC == pdfc.getEvaluatorTypeUsed());
    CHECK(allclose(v10, pdfc.value()));
  }
  SECTION("fast update") {
    pdfc.eval(mstru10);
    pdfc.countdown = nfast / 2;
    CHECK_THROWS_AS(pdfc.eval(mstru10d1), runtime_error);
    pdfc.countdown = -1;
    pdfc.eval(mstru10);
    CHECK(BASIC == pdfc.getEvaluatorTypeUsed());
    CHECK(allclose(v10, pdfc.value()));
  }
  SECTION("check of a fast update") {
    // fail in the full calculation that follows the fast update
    pdfc.setEvaluatorType(PQEvaluatorType::CHECK);
    pdfc.eval(mstru10);
    pdfc.countdown = nfast + 2;
    CHECK_THROWS_AS(pdfc.eval(mstru10d1), runtime_error);
    pdfc.countdown = -1;
    pdfc.eval(mstru10d1);
    CHECK(allclose(v10d1, pdfc.value()));
  }
}

TEST_CASE_METHOD(TestPQEvaluatorFixture,
                 "TestPQEvaluator::test_structure_changed_in_place",
                 "[PQEvaluator]") {
  // The first evaluation after an atom moves off a special position must
  // include its new symmetry positions.
  mpdfcb.setRmax(7.0);
  mpdfcb.eval(inversionCrystal(0.3));
  const QuantityType expected = mpdfcb.getPDF();
  for (PQEvaluatorType tp : {BASIC, OPTIMIZED, PQEvaluatorType::CHECK}) {
    INFO("evaluator " << tp);
    CrystalStructureAdapterPtr stru = inversionCrystal(0.0);
    PDFCalculator pdfc;
    pdfc.setEvaluatorType(tp);
    pdfc.setRmax(7.0);
    pdfc.eval(stru);
    REQUIRE(1 == stru->siteMultiplicity(0));
    REQUIRE(!allclose(expected, pdfc.getPDF()));
    stru->at(0).xyz_cartn[0] = 0.3;
    for (int i = 0; i < 2; ++i) {
      pdfc.eval(stru);
      // the first evaluation is a fast update unless the type is BASIC
      if (i == 0) CHECK(tp == pdfc.getEvaluatorTypeUsed());
      CHECK(2 == stru->siteMultiplicity(0));
      CHECK(allclose(expected, pdfc.getPDF()));
    }
  }
}

}  // namespace srreal
}  // namespace diffpy

// End of file
