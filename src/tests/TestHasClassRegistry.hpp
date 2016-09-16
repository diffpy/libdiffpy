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

        void tearDown()
        {
            // restore SFTXray registration if removed in some test.
            msftb->deregisterType(msftb->type());
            msftb->registerThisType();
            TS_ASSERT(ScatteringFactorTable::isRegisteredType("X"));
            TS_ASSERT(ScatteringFactorTable::isRegisteredType("xray"));
        }


        void test_deregisterType()
        {
            TS_ASSERT_EQUALS(0,
                    ScatteringFactorTable::deregisterType("invalid"));
            TS_ASSERT_EQUALS(2,
                    ScatteringFactorTable::deregisterType("X"));
            TS_ASSERT(!(mpreg->isRegisteredType("xray")));
            TS_ASSERT(!(mpreg->isRegisteredType("X")));
        }


        void test_isRegisteredType()
        {
            TS_ASSERT(msftb->isRegisteredType("xray"));
            TS_ASSERT(ScatteringFactorTable::isRegisteredType("X"));
            TS_ASSERT_EQUALS(false, mpreg->isRegisteredType("invalid"));
        }


        void test_getAliasedTypes()
        {
            std::map<std::string, std::string> atps;
            atps = ScatteringFactorTable::getAliasedTypes();
            TS_ASSERT_EQUALS(4, atps.size());
            TS_ASSERT(atps.count("X"));
            TS_ASSERT_EQUALS("xray", atps["X"]);
            TS_ASSERT_EQUALS("electronnumber", atps["EN"]);
        }

};  // class TestHasClassRegistry

// End of file
