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
 * TestAttributes -- test the helper Attributes class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/Attributes.hpp>

// Local example class -------------------------------------------------------

namespace {

class Example : public diffpy::Attributes {
 public:
  Example() : ma(0.0) {
    this->registerDoubleAttribute("a", this, &Example::geta, &Example::seta);
    this->registerDoubleAttribute("b", this, &Example::geta);
  }

  double geta() const { return ma; }
  void seta(double a) { ma = a; }

 private:
  double ma;
};

}  // namespace

// ---------------------------------------------------------------------------

class TestAttributesFixture {
 protected:
  // data
  std::unique_ptr<Example> mobj;

 public:
  TestAttributesFixture() { mobj.reset(new Example); }
};  // class TestAttributesFixture

TEST_CASE_METHOD(TestAttributesFixture, "TestAttributes::test_names",
                 "[Attributes]") {
  CHECK(2 == mobj->namesOfDoubleAttributes().size());
  CHECK(1 == mobj->namesOfWritableDoubleAttributes().size());
}

TEST_CASE_METHOD(TestAttributesFixture, "TestAttributes::test_getsetdoubleattr",
                 "[Attributes]") {
  using diffpy::attributes::DoubleAttributeError;
  CHECK(0.0 == mobj->getDoubleAttr("a"));
  mobj->setDoubleAttr("a", 1.3);
  CHECK(1.3 == mobj->getDoubleAttr("a"));
  CHECK_THROWS_AS(mobj->setDoubleAttr("b", 3), DoubleAttributeError);
  CHECK(1.3 == mobj->getDoubleAttr("b"));
  const Example ex1;
  CHECK(0.0 == ex1.getDoubleAttr("a"));
  CHECK_THROWS_AS(ex1.getDoubleAttr("bad"), DoubleAttributeError);
}

// End of file
