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
 * TestPeakWidthModel -- test PeakWidthModel and derived classes
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/DebyeWallerPeakWidth.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include "serialization_helpers.hpp"

using Catch::Matchers::WithinAbs;

namespace diffpy {
namespace srreal {

const double tuiso = 0.004;
const double UtoB = 8 * M_PI * M_PI;

//////////////////////////////////////////////////////////////////////////////
// TestPeakWidthModel
//////////////////////////////////////////////////////////////////////////////

class TestPeakWidthModelFixture {
 protected:
  // data
  PeakWidthModelPtr mcnpw;
  ConstantPeakWidth* mccnpw;
  PeakWidthModelPtr mdwpw;
  PeakWidthModelPtr mjepw;
  JeongPeakWidth* mcjepw;
  double meps;

  // methods
  static AtomicStructureAdapterPtr dimer() {
    AtomicStructureAdapterPtr adpt(new AtomicStructureAdapter);
    Atom a;
    a.atomtype = "C";
    a.uij_cartn = R3::identity() * tuiso;
    adpt->append(a);
    a.xyz_cartn[2] = 1;
    adpt->append(a);
    return adpt;
  }

 public:
  TestPeakWidthModelFixture() {
    mcnpw = PeakWidthModel::createByType("constant");
    mccnpw = static_cast<ConstantPeakWidth*>(mcnpw.get());
    mdwpw = PeakWidthModel::createByType("debye-waller");
    mjepw = PeakWidthModel::createByType("jeong");
    mcjepw = static_cast<JeongPeakWidth*>(mjepw.get());
    meps = 10 * diffpy::mathutils::DOUBLE_EPS;
  }
};  // class TestPeakWidthModelFixture

TEST_CASE_METHOD(TestPeakWidthModelFixture, "TestPeakWidthModel::test_clone",
                 "[PeakWidthModel]") {
  CHECK("constant" == mcnpw->clone()->type());
  CHECK("debye-waller" == mdwpw->clone()->type());
  CHECK("jeong" == mjepw->clone()->type());
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_attributes", "[PeakWidthModel]") {
  CHECK(0 == mdwpw->namesOfDoubleAttributes().size());
  CHECK(4 == mjepw->namesOfDoubleAttributes().size());
  mcjepw->setDelta1(0.1);
  mcjepw->setDelta2(0.2);
  mcjepw->setQbroad(0.3);
  mcjepw->setQbroad_seperable(0.4);
  CHECK(0.1 == mjepw->getDoubleAttr("delta1"));
  CHECK(0.2 == mjepw->getDoubleAttr("delta2"));
  CHECK(0.3 == mjepw->getDoubleAttr("qbroad"));
  CHECK(0.4 == mjepw->getDoubleAttr("qbroad_seperable"));
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_calculate", "[PeakWidthModel]") {
  using diffpy::mathutils::GAUSS_SIGMA_TO_FWHM;
  AtomicStructureAdapterPtr adpt = dimer();
  BaseBondGeneratorPtr bnds = adpt->createBondGenerator();
  bnds->setRmax(2.0);
  bnds->rewind();
  const double w = GAUSS_SIGMA_TO_FWHM * sqrt(2 * tuiso);
  CHECK_THAT(w, WithinAbs(mdwpw->calculate(*bnds), meps));
  CHECK_THAT(w, WithinAbs(mjepw->calculate(*bnds), meps));
  mcjepw->setDelta2(0.75);
  CHECK_THAT(0.5 * w, WithinAbs(mjepw->calculate(*bnds), meps));
  mcjepw->setDelta2(0);
  mcjepw->setQbroad(sqrt(3));
  mcjepw->setQbroad_seperable(0);
  CHECK_THAT(2 * w, WithinAbs(mjepw->calculate(*bnds), meps));
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_serialization", "[PeakWidthModel]") {
  PeakWidthModelPtr pwm1 = dumpandload(mdwpw);
  DebyeWallerPeakWidth* cpwm1 = dynamic_cast<DebyeWallerPeakWidth*>(pwm1.get());
  CHECK(cpwm1);
  mcjepw->setDelta1(1.1);
  mcjepw->setDoubleAttr("delta2", 2.2);
  mcjepw->setDoubleAttr("qbroad", 0.03);
  mcjepw->setDoubleAttr("qbroad_seperable", 0.003);
  PeakWidthModelPtr pwm2 = dumpandload(mjepw);
  JeongPeakWidth* cpwm2 = dynamic_cast<JeongPeakWidth*>(pwm2.get());
  CHECK(cpwm2);
  CHECK(1.1 == pwm2->getDoubleAttr("delta1"));
  CHECK(2.2 == pwm2->getDoubleAttr("delta2"));
  CHECK(0.03 == pwm2->getDoubleAttr("qbroad"));
  CHECK(0.003 == pwm2->getDoubleAttr("qbroad_seperable"));
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_CPWM_attributes",
                 "[PeakWidthModel]") {
  CHECK("constant" == mcnpw->type());
  CHECK(mcnpw->hasDoubleAttr("width"));
  CHECK(mcnpw->hasDoubleAttr("bisowidth"));
  CHECK(mcnpw->hasDoubleAttr("uisowidth"));
  CHECK(3 == mcnpw->namesOfDoubleAttributes().size());
  mccnpw->setWidth(0.1);
  CHECK(0.1 == mcnpw->getDoubleAttr("width"));
  CHECK(0.1 == mcnpw->maxWidth(dimer(), 0, 4));
  const double bw0 = mcnpw->getDoubleAttr("bisowidth");
  const double uw0 = mcnpw->getDoubleAttr("uisowidth");
  mccnpw->setWidth(2 * mccnpw->getWidth());
  CHECK_THAT(bw0, WithinAbs(UtoB * uw0, meps));
  CHECK_THAT(4 * bw0, WithinAbs(mcnpw->getDoubleAttr("bisowidth"), meps));
  CHECK_THAT(4 * uw0, WithinAbs(mcnpw->getDoubleAttr("uisowidth"), meps));
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_CPWM_bisowidth",
                 "[PeakWidthModel]") {
  const double biso = UtoB * tuiso;
  mcnpw->setDoubleAttr("bisowidth", biso);
  CHECK_THAT(tuiso, WithinAbs(mcnpw->getDoubleAttr("uisowidth"), meps));
  mcnpw->setDoubleAttr("uisowidth", -tuiso);
  CHECK_THAT(-biso, WithinAbs(mcnpw->getDoubleAttr("bisowidth"), meps));
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_CPWM_uisowidth",
                 "[PeakWidthModel]") {
  StructureAdapterPtr adpt = dimer();
  BaseBondGeneratorPtr bnds = adpt->createBondGenerator();
  PeakWidthModelPtr dwpw(new DebyeWallerPeakWidth);
  bnds->setRmax(2.0);
  bnds->rewind();
  const double w0 = dwpw->calculate(*bnds);
  CHECK(0.0 == mcnpw->calculate(*bnds));
  mcnpw->setDoubleAttr("uisowidth", tuiso);
  CHECK_THAT(w0, WithinAbs(mcnpw->calculate(*bnds), meps));
  mcnpw->setDoubleAttr("uisowidth", -tuiso);
  CHECK_THAT(-tuiso, WithinAbs(mcnpw->getDoubleAttr("uisowidth"), meps));
  CHECK_THAT(-w0, WithinAbs(mcnpw->getDoubleAttr("width"), meps));
}

TEST_CASE_METHOD(TestPeakWidthModelFixture,
                 "TestPeakWidthModel::test_CPWM_serialization",
                 "[PeakWidthModel]") {
  PeakWidthModelPtr pwm1 = dumpandload(mcnpw);
  CHECK("constant" == pwm1->type());
  CHECK(0.0 == pwm1->getDoubleAttr("width"));
  mcnpw->setDoubleAttr("width", 0.1);
  double uw0 = mcnpw->getDoubleAttr("uisowidth");
  PeakWidthModelPtr pwm2 = dumpandload(mcnpw);
  CHECK(0.1 == pwm2->getDoubleAttr("width"));
  CHECK_THAT(uw0, WithinAbs(pwm2->getDoubleAttr("uisowidth"), meps));
}

}  // namespace srreal
}  // namespace diffpy

// End of file
