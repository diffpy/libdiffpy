/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestPairCounter -- unit tests for PairCounter class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PairCounter.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestPairCounter : public CxxTest::TestSuite
{

private:

    AtomicStructureAdapterPtr mstru;
    AtomicStructureAdapterPtr mline100;

public:

    void setUp()
    {
        mstru.reset(new AtomicStructureAdapter);
        if (!mline100)
        {
            mline100.reset(new AtomicStructureAdapter);
            Atom a;
            for (int i = 0; i < 100; ++i)
            {
                a.xyz_cartn = R3::Vector(1.0*i, 0.0, 0.0);
                mline100->append(a);
            }
        }
    }


    void test_call()
    {
        PairCounter pcount;
        Atom a;
        for (int i = 0; i < 100; ++i)
        {
            int npairs = i * (i - 1) / 2;
            TS_ASSERT_EQUALS(npairs, pcount(mstru));
            a.xyz_cartn = R3::Vector(1.0*i, 0.0, 0.0);
            mstru->append(a);
        }
        mstru->clear();
        TS_ASSERT_EQUALS(0, pcount(mstru));
    }


    void test_setRmin()
    {
        PairCounter pcount;
        TS_ASSERT_EQUALS(100*99/2, pcount(mline100));
        pcount.setRmin(100);
        TS_ASSERT_EQUALS(0, pcount(mline100));
        pcount.setRmin(99.0);
        TS_ASSERT_EQUALS(1, pcount(mline100));
        pcount.setRmin(1.1);
        TS_ASSERT_EQUALS(100*99/2 - 99, pcount(mline100));
    }


    void test_setRmax()
    {
        PairCounter pcount;
        pcount.setRmax(0.9);
        TS_ASSERT_EQUALS(0, pcount(mline100));
        pcount.setRmax(1.1);
        TS_ASSERT_EQUALS(99, pcount(mline100));
        pcount.setRmax(98.5);
        TS_ASSERT_EQUALS(100*99/2 - 1, pcount(mline100));
    }


    void test_parallel()
    {
        const int ncpu = 7;
        PairCounter pmaster;
        PairCounter pslave[ncpu];
        bool allequal = true;
        for (int cpuindex = 0; cpuindex < ncpu; ++cpuindex)
        {
            pslave[cpuindex].setupParallelRun(cpuindex, ncpu);
            pslave[cpuindex].eval(mline100);
            allequal = allequal &&
                (pslave[cpuindex].value()[0] == pslave[0].value()[0]);
        }
        TS_ASSERT(!allequal);
        pmaster.setStructure(mline100);
        TS_ASSERT_EQUALS(1u, pmaster.value().size());
        TS_ASSERT_EQUALS(0.0, pmaster.value()[0]);
        for (PairCounter* p = pslave; p != pslave + ncpu; ++p)
        {
            pmaster.mergeParallelData(p->getParallelData(), ncpu);
        }
        TS_ASSERT_EQUALS(100 * 99 / 2, pmaster.value()[0]);
    }

};  // class TestPairCounter

// End of file
