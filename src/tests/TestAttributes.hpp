/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2019 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestAttributes -- test the helper Attributes class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/Attributes.hpp>

// Local example class -------------------------------------------------------

namespace {

class Example : public diffpy::Attributes
{
    public:

        Example() : ma(0.0)
        {
            this->registerDoubleAttribute("a", this,
                    &Example::geta, &Example::seta);
            this->registerDoubleAttribute("b", this,
                    &Example::geta);
        }

        double geta() const  { return ma; }
        void seta(double a)  { ma = a; }

    private:

        double ma;

};

}   // namespace

// ---------------------------------------------------------------------------

class TestAttributes : public CxxTest::TestSuite
{
    private:

        // data
        std::unique_ptr<Example> mobj;

    public:

        void setUp()
        {
            mobj.reset(new Example);
        }


        void test_names()
        {
            TS_ASSERT_EQUALS(2, mobj->namesOfDoubleAttributes().size());
            TS_ASSERT_EQUALS(1, mobj->namesOfWritableDoubleAttributes().size());
        }


        void test_getsetdoubleattr()
        {
            using diffpy::attributes::DoubleAttributeError;
            TS_ASSERT_EQUALS(0.0, mobj->getDoubleAttr("a"));
            mobj->setDoubleAttr("a", 1.3);
            TS_ASSERT_EQUALS(1.3, mobj->getDoubleAttr("a"));
            TS_ASSERT_THROWS(mobj->setDoubleAttr("b", 3), DoubleAttributeError);
            TS_ASSERT_EQUALS(1.3, mobj->getDoubleAttr("b"));
            const Example ex1;
            TS_ASSERT_EQUALS(0.0, ex1.getDoubleAttr("a"));
            TS_ASSERT_THROWS(ex1.getDoubleAttr("bad"), DoubleAttributeError);
        }

};  // class TestAttributes

// End of file
