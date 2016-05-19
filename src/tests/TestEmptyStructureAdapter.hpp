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
* class TestEmptyStructureAdapter -- test adapter that is always empty
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// class TestEmptyStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestEmptyStructureAdapter : public CxxTest::TestSuite
{
    private:

        StructureAdapterPtr mstru;

    public:

        void setUp()
        {
            mstru = emptyStructureAdapter();
        }


        void test_empty_instance()
        {
            TS_ASSERT_EQUALS(0, mstru->countSites());
            TS_ASSERT_EQUALS(mstru.get(), mstru->clone().get());
            TS_ASSERT_EQUALS(mstru, emptyStructureAdapter());
            TS_ASSERT_THROWS(mstru->siteAtomType(0), std::out_of_range);
            TS_ASSERT_THROWS(mstru->siteCartesianPosition(0),
                    std::out_of_range);
            TS_ASSERT_THROWS(mstru->siteMultiplicity(0), std::out_of_range);
            TS_ASSERT_THROWS(mstru->siteOccupancy(0), std::out_of_range);
            TS_ASSERT_THROWS(mstru->siteAnisotropy(0), std::out_of_range);
            TS_ASSERT_THROWS(mstru->siteCartesianPosition(0),
                    std::out_of_range);
        }


        void test_serialization()
        {
            StructureAdapterPtr stru1;
            stru1 = dumpandload(mstru);
            AtomicStructureAdapterPtr astru1 =
                boost::dynamic_pointer_cast<AtomicStructureAdapter>(stru1);
            TS_ASSERT(!astru1);
            TS_ASSERT_EQUALS(0, stru1->countSites());
            TS_ASSERT_EQUALS(stru1, stru1->clone());
            TS_ASSERT_THROWS(mstru->siteAtomType(0), std::out_of_range);
        }

};  // class TestEmptyStructureAdapter

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestEmptyStructureAdapter;

// End of file
