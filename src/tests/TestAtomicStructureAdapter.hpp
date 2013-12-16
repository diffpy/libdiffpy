/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestAtomicStructureAdapter -- unit tests for an adapter that
*     stores data in a series of Atom objects
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <boost/make_shared.hpp>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// class TestAtomicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestAtomicStructureAdapter : public CxxTest::TestSuite
{
    private:

        StructureAdapterPtr mstru;
        AtomicStructureAdapterPtr mpstru;

    public:

        void setUp()
        {
            mstru = StructureAdapterPtr(new AtomicStructureAdapter);
            mpstru =
                boost::dynamic_pointer_cast<AtomicStructureAdapter>(mstru);
        }


        void test_diff()
        {
            typedef StructureDifference::Method DM;
            const DM::Type& NONE = DM::NONE;
            const DM::Type& SIDEBYSIDE = DM::SIDEBYSIDE;
            const DM::Type& SORTED = DM::SORTED;
            StructureDifference sd;
            sd = mstru->diff(emptyStructureAdapter());
            TS_ASSERT(sd.add1.empty());
            TS_ASSERT(!sd.allowsfastupdate());
            sd = mstru->diff(StructureAdapterPtr());
            TS_ASSERT_EQUALS(NONE, sd.diffmethod);
            TS_ASSERT(!sd.allowsfastupdate());
            Atom ai;
            ai.atomtype = "C";
            const int SZ = 10;
            for (int i = 0; i < SZ; ++i)
            {
                ai.cartesianposition[0] = i;
                mpstru->append(ai);
            }
            AtomicStructureAdapterPtr cpstru =
                boost::make_shared<AtomicStructureAdapter>(*mpstru);
            sd = mstru->diff(mstru);
            TS_ASSERT(sd.allowsfastupdate())
            sd = mstru->diff(cpstru);
            TS_ASSERT_EQUALS(SIDEBYSIDE, sd.diffmethod);
            TS_ASSERT(sd.allowsfastupdate())
            TS_ASSERT(sd.pop0.empty());
            TS_ASSERT(sd.add1.empty());
            (*cpstru)[0].atomtype = "N";
            sd = mstru->diff(cpstru);
            TS_ASSERT_EQUALS(SIDEBYSIDE, sd.diffmethod);
            TS_ASSERT(sd.allowsfastupdate())
            TS_ASSERT_EQUALS(1u, sd.pop0.size());
            TS_ASSERT_EQUALS(1u, sd.add1.size());
            for (int i = 1; i < (1 - sqrt(0.5)) * SZ; ++i)
            {
                cpstru->erase(0);
                sd = mstru->diff(cpstru);
                TS_ASSERT_EQUALS(SORTED, sd.diffmethod);
                TS_ASSERT(sd.allowsfastupdate());
                TS_ASSERT_EQUALS(i, int(sd.pop0.size()));
                TS_ASSERT(sd.add1.empty());
            }
            cpstru->erase(0);
            sd = mstru->diff(cpstru);
            TS_ASSERT(!sd.allowsfastupdate());
            Atom a2;
            a2.atomtype = "N";
            cpstru->append(a2);
            sd = mstru->diff(cpstru);
            TS_ASSERT(!sd.allowsfastupdate());
            TS_ASSERT_EQUALS(1u, sd.add1.size());
        }


        void test_serialization()
        {
            Atom ai;
            ai.atomtype = "C";
            ai.cartesianposition = R3::Vector(1, 2, 3);
            ai.anisotropy = true;
            mpstru->append(ai);
            ai.atomtype = "H";
            ai.cartesianposition = R3::Vector(4, 5, 6);
            ai.anisotropy = false;
            mpstru->append(ai);
            StructureAdapterPtr stru1;
            stru1 = dumpandload(mstru);
            AtomicStructureAdapterPtr astru1 =
                boost::dynamic_pointer_cast<AtomicStructureAdapter>(stru1);
            TS_ASSERT_EQUALS(2, astru1->countSites());
            TS_ASSERT_EQUALS((*mpstru)[0], (*astru1)[0]);
            TS_ASSERT_EQUALS((*mpstru)[1], (*astru1)[1]);
        }


        void test_comparison()
        {
            Atom ai;
            ai.atomtype = "C";
            const int SZ = 10;
            for (int i = 0; i < SZ; ++i)
            {
                ai.cartesianposition[0] = i;
                mpstru->append(ai);
            }
            AtomicStructureAdapterPtr cpstru =
                boost::make_shared<AtomicStructureAdapter>(*mpstru);
            TS_ASSERT_EQUALS(*mpstru, *cpstru);
            TS_ASSERT(!(*mpstru != *cpstru));
            cpstru->at(0).atomtype = "H";
            TS_ASSERT_DIFFERS(*mpstru, *cpstru);
            TS_ASSERT(!(*mpstru == *cpstru));
        }

};  // class TestAtomicStructureAdapter

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestAtomicStructureAdapter;

// End of file
