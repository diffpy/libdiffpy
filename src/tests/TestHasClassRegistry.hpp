/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2016 Brookhaven Science Associates,
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
* class TestHasClassRegistry -- test registration machinery for named classes.
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/ScatteringFactorTable.hpp>

using diffpy::srreal::ScatteringFactorTable;
using diffpy::srreal::ScatteringFactorTablePtr;

//////////////////////////////////////////////////////////////////////////////
// class TestHasClassRegistry
//////////////////////////////////////////////////////////////////////////////

class TestHasClassRegistry : public CxxTest::TestSuite
{
    private:

        typedef diffpy::HasClassRegistry<ScatteringFactorTable> HCRBase;
        ScatteringFactorTablePtr msftb;
        HCRBase* mpreg;

    public:

        void setUp()
        {
            msftb = ScatteringFactorTable::createByType("xray");
            mpreg = msftb.get();
        }


        void test_isRegisteredType()
        {
            TS_ASSERT(msftb->isRegisteredType("xray"));
            TS_ASSERT(ScatteringFactorTable::isRegisteredType("X"));
            TS_ASSERT_EQUALS(false, mpreg->isRegisteredType("invalid"));
        }

};  // class TestHasClassRegistry

// End of file
