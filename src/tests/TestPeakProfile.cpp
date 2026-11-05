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
 * TestPeakProfile -- unit tests for various PeakProfile classes
 *
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/PeakProfile.hpp>
#include <diffpy/srreal/GaussianProfile.hpp>
#include <diffpy/srreal/CroppedGaussianProfile.hpp>
#include "serialization_helpers.hpp"

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;
using diffpy::mathutils::eps_eq;

class TestPeakProfileFixture {
 protected:
  PeakProfilePtr mpkgauss;
  PeakProfilePtr mpkgcrop;
  double meps;

 public:
  TestPeakProfileFixture() {
    mpkgauss = PeakProfile::createByType("gaussian");
    mpkgcrop = PeakProfile::createByType("croppedgaussian");
    meps = 10 * diffpy::mathutils::DOUBLE_EPS;
  }
};  // class TestPeakProfileFixture

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_factory",
                 "[PeakProfile]") {
  CHECK_THROWS_AS(PeakProfile::createByType("invalid"), invalid_argument);
  CHECK(string("gaussian") == mpkgauss->type());
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_clone",
                 "[PeakProfile]") {
  CHECK("gaussian" == mpkgauss->clone()->type());
  CHECK("croppedgaussian" == mpkgcrop->clone()->type());
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_attributes",
                 "[PeakProfile]") {
  CHECK(1 == mpkgauss->namesOfDoubleAttributes().size());
  CHECK(1 == mpkgcrop->namesOfDoubleAttributes().size());
  mpkgauss->setDoubleAttr("peakprecision", 0.05);
  const double lo0 = mpkgauss->xboundlo(0.2);
  GaussianProfile g;
  CHECK(!eps_eq(lo0, g.xboundlo(0.2)));
  g.setPrecision(0.05);
  CHECK(lo0 == g.xboundlo(0.2));
  mpkgcrop->setDoubleAttr("peakprecision", 0.05);
  CroppedGaussianProfile fcrop;
  CHECK(!eps_eq(fcrop(0.1, 0.2), (*mpkgcrop)(0.1, 0.2)));
  fcrop.setPrecision(0.05);
  CHECK(fcrop(0.1, 0.2) == (*mpkgcrop)(0.1, 0.2));
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_y",
                 "[PeakProfile]") {
  double Afwhm1 = 2 * sqrt(M_LN2 / M_PI);
  const PeakProfile& pkgauss = *mpkgauss;
  CHECK_THAT(Afwhm1, WithinAbs(pkgauss(0, 1), meps));
  CHECK_THAT(Afwhm1 / 2, WithinAbs(pkgauss(+0.5, 1), meps));
  CHECK_THAT(Afwhm1 / 2, WithinAbs(pkgauss(-0.5, 1), meps));
  const PeakProfile& pkgcrop = *mpkgcrop;
  mpkgauss->setDoubleAttr("peakprecision", 1e-6);
  mpkgcrop->setDoubleAttr("peakprecision", 1e-6);
  CHECK(pkgcrop(0, 1) > pkgauss(0, 1) + meps);
  CHECK(0 == pkgcrop(1, 0.1));
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_xboundlo",
                 "[PeakProfile]") {
  const double epsy = 1e-8;
  mpkgauss->setPrecision(epsy);
  const double Afwhm1 = 2 * sqrt(M_LN2 / M_PI);
  const double xblo1 = mpkgauss->xboundlo(1);
  const double Afwhm3 = Afwhm1 / 3;
  const double xblo3 = mpkgauss->xboundlo(3);
  CHECK_THAT(epsy, WithinAbs((*mpkgauss)(xblo1, 1) / Afwhm1, meps));
  CHECK_THAT(epsy, WithinAbs((*mpkgauss)(xblo3, 3) / Afwhm3, meps));
  CHECK(xblo1 < 0);
  CHECK(xblo3 < 0);
  mpkgauss->setPrecision(10);
  CHECK(0.0 == mpkgauss->xboundlo(1));
  mpkgauss->setPrecision(0.1);
  CHECK(0.0 == mpkgauss->xboundlo(0));
  CHECK(0.0 == mpkgauss->xboundlo(-1));
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_xboundhi",
                 "[PeakProfile]") {
  const double epsy = 1e-8;
  mpkgauss->setPrecision(epsy);
  const double Afwhm1 = 2 * sqrt(M_LN2 / M_PI);
  const double xbhi1 = mpkgauss->xboundhi(1);
  const double Afwhm3 = Afwhm1 / 3;
  const double xbhi3 = mpkgauss->xboundhi(3);
  CHECK_THAT(epsy, WithinAbs((*mpkgauss)(xbhi1, 1) / Afwhm1, meps));
  CHECK_THAT(epsy, WithinAbs((*mpkgauss)(xbhi3, 3) / Afwhm3, meps));
  CHECK(xbhi1 > 0);
  CHECK(xbhi3 > 0);
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_setPrecision",
                 "[PeakProfile]") {
  using diffpy::mathutils::DOUBLE_EPS;
  double epsy = 1e-7;
  CHECK_THAT(0.0, WithinAbs(mpkgauss->getPrecision(), 10 * DOUBLE_EPS));
  mpkgauss->setPrecision(epsy);
  CHECK(epsy == mpkgauss->getPrecision());
  double xbhi1 = mpkgauss->xboundhi(1);
  mpkgauss->setPrecision(1e-4);
  CHECK(xbhi1 != mpkgauss->xboundhi(1));
}

TEST_CASE_METHOD(TestPeakProfileFixture, "TestPeakProfile::test_serialization",
                 "[PeakProfile]") {
  mpkgauss->setPrecision(0.0123);
  PeakProfilePtr pk1 = dumpandload(mpkgauss);
  CHECK(pk1.get());
  CHECK(mpkgauss.get() != pk1.get());
  CHECK(string("gaussian") == pk1->type());
  CHECK(0.0123 == pk1->getPrecision());
  mpkgcrop->setPrecision(0.004);
  PeakProfilePtr pk2 = dumpandload(mpkgcrop);
  CHECK("croppedgaussian" == pk2->type());
  CHECK(0.004 == pk2->getPrecision());
  CHECK((*mpkgcrop)(0.1, 0.2) == (*pk2)(0.1, 0.2));
}

// End of file
