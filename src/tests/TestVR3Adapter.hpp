/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestVR3Adapter -- unit tests for an adapter to the VR3Structure class
*
*****************************************************************************/

#include <typeinfo>
#include <sstream>
#include <cxxtest/TestSuite.h>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/VR3Structure.hpp>

using namespace std;
using namespace boost;
using namespace diffpy::srreal;

//////////////////////////////////////////////////////////////////////////////
// class TestVR3Adapter
//////////////////////////////////////////////////////////////////////////////

class TestVR3Adapter : public CxxTest::TestSuite
{
    private:

        StructureAdapterPtr mempty;
        StructureAdapterPtr mthree;

    public:

        void setUp()
        {
            if (!mempty.get())
            {
                mempty.reset(new VR3Adapter);
            }
            if (!mthree.get())
            {
                VR3Structure struthree;
                struthree.push_back(R3::Vector(0, 0, 0));
                struthree.push_back(R3::Vector(1, 2, 3));
                struthree.push_back(R3::Vector(4, 5, 6));
                mthree = createStructureAdapter(struthree);
            }
        }


        void test_typeid()
        {
            TS_ASSERT(typeid(VR3Adapter) == typeid(*mempty));
        }


        void test_countSites()
        {
            TS_ASSERT_EQUALS(0, mempty->countSites());
            TS_ASSERT_EQUALS(3, mthree->countSites());
        }


        void test_totalOccupancy()
        {
            TS_ASSERT_EQUALS(0.0, mempty->totalOccupancy());
            TS_ASSERT_EQUALS(3.0, mthree->totalOccupancy());
        }


        void test_numberDensity()
        {
            TS_ASSERT_EQUALS(0.0, mempty->numberDensity());
            TS_ASSERT_EQUALS(0.0, mthree->numberDensity());
        }


        void test_siteCartesianPosition()
        {
            TS_ASSERT_EQUALS(4.0, mthree->siteCartesianPosition(2)[0]);
            TS_ASSERT_EQUALS(5.0, mthree->siteCartesianPosition(2)[1]);
            TS_ASSERT_EQUALS(6.0, mthree->siteCartesianPosition(2)[2]);
        }


        void test_siteAnisotropy()
        {
            for (int i = 0; i < mthree->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(false, mthree->siteAnisotropy(i));
            }
        }


        void test_siteCartesianUij()
        {
            // Uij is always zero
            TS_ASSERT_EQUALS(R3::zeromatrix(), mthree->siteCartesianUij(0));
        }


        void test_siteAtomType()
        {
            for (int i = 0; i < mthree->countSites(); ++i)
            {
                TS_ASSERT(mthree->siteAtomType(0).empty());
            }
        }


        void test_serialization()
        {
            stringstream storage(ios::in | ios::out | ios::binary);
            diffpy::serialization::oarchive oa(storage, ios::binary);
            oa << mthree;
            diffpy::serialization::iarchive ia(storage, ios::binary);
            StructureAdapterPtr three1;
            TS_ASSERT(!three1.get());
            ia >> three1;
            TS_ASSERT_DIFFERS(mthree.get(), three1.get());
            TS_ASSERT_EQUALS(3, three1->countSites());
            TS_ASSERT_EQUALS(3.0, three1->totalOccupancy());
            TS_ASSERT_EQUALS(0.0, three1->numberDensity());
            for (int i = 0; i < mthree->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(0.0, R3::distance(
                            mthree->siteCartesianPosition(i),
                            three1->siteCartesianPosition(i)));
            }
        }

};  // class TestVR3Adapter

// End of file
