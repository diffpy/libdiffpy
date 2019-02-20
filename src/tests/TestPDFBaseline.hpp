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
* class TestPDFBaseline -- test PDFBaseline and derived classes
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/PDFBaseline.hpp>
#include <diffpy/srreal/LinearBaseline.hpp>
#include <diffpy/srreal/ZeroBaseline.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class TestPDFBaseline
//////////////////////////////////////////////////////////////////////////////

class TestPDFBaseline : public CxxTest::TestSuite
{
    private:

        // data
        PDFBaselinePtr ml;
        LinearBaseline* mcl;
        PDFBaselinePtr mz;
        ZeroBaseline* mcz;

    public:

        void setUp()
        {
            ml = PDFBaseline::createByType("linear");
            mcl = static_cast<LinearBaseline*>(ml.get());
            mz = PDFBaseline::createByType("zero");
            mcz = static_cast<ZeroBaseline*>(mz.get());
        }


        void test_clone()
        {
            TS_ASSERT_EQUALS("linear", ml->clone()->type());
            TS_ASSERT_EQUALS("zero", mz->clone()->type());
        }


        void test_attributes()
        {
            TS_ASSERT_EQUALS(1, ml->namesOfDoubleAttributes().size());
            TS_ASSERT_EQUALS(0, mz->namesOfDoubleAttributes().size());
            mcl->setSlope(0.1);
            TS_ASSERT_EQUALS(0.1, ml->getDoubleAttr("slope"));
            ml->setDoubleAttr("slope", 0.2);
            TS_ASSERT_EQUALS(0.2, mcl->getSlope());
        }


        void test_calculate()
        {
            TS_ASSERT_EQUALS(0.0, (*ml)(2.0));
            TS_ASSERT_EQUALS(0.0, (*mz)(2.0));
            mcl->setSlope(-2);
            TS_ASSERT_EQUALS(-4.0, (*ml)(2.0));
        }


        void test_serialization()
        {
            mcl->setSlope(0.3);
            PDFBaselinePtr bl1 = dumpandload(ml);
            LinearBaseline* cbl1 = dynamic_cast<LinearBaseline*>(bl1.get());
            TS_ASSERT(cbl1);
            TS_ASSERT_EQUALS(0.3, cbl1->getSlope());
            PDFBaselinePtr bz1 = dumpandload(mz);
            ZeroBaseline* cbz1 = dynamic_cast<ZeroBaseline*>(bz1.get());
            TS_ASSERT(cbz1);
        }

};  // class TestPDFBaseline

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestPDFBaseline;

// End of file
