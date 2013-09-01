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
* class TestPQEvaluator -- unit tests for the PairQuantity evaluator class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <functional>

#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>

namespace diffpy {
namespace srreal {

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// class TestPQEvaluator
//////////////////////////////////////////////////////////////////////////////

class TestPQEvaluator : public CxxTest::TestSuite
{
    private:

        typedef boost::shared_ptr<AtomicStructureAdapter> AtomicAdapterPtr;
        AtomicAdapterPtr mstru10;
        AtomicAdapterPtr mstru10d1;
        AtomicAdapterPtr mstru10r;
        AtomicAdapterPtr mstru9;

     public:

        void setUp()
        {
            const int SZ = 10;
            mstru10.reset(new AtomicStructureAdapter);
            Atom ai;
            ai.atomtype = "C";
            ai.cartesianuij = R3::identity();
            ai.cartesianuij(0, 0) = ai.cartesianuij(1, 1) =
                ai.cartesianuij(2, 2) = 0.004;
            for (int i = 0; i < SZ; ++i)
            {
                ai.cartesianposition[0] = i;
                mstru10->append(ai);
            }
            mstru10d1.reset(new AtomicStructureAdapter(*mstru10));
            (*mstru10d1)[0].atomtype = "Au";
            mstru10r.reset(new AtomicStructureAdapter);
            mstru10r->assign(mstru10->rbegin(), mstru10->rend());
            mstru9.reset(new AtomicStructureAdapter(*mstru10));
            mstru9->remove(9);
        }

        void test_PDFCalculator()
        {
            PDFCalculator pdfcb;
            PDFCalculator pdfco;
            pdfcb.setEvaluator(BASIC);
            pdfco.setEvaluator(OPTIMIZED);
            pdfcb.eval(mstru10);
            pdfco.eval(mstru10);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT(!gb.empty());
            TS_ASSERT_EQUALS(gb, go);
            int cnonzero = count_if(gb.begin(), gb.end(),
                    bind1st(not_equal_to<double>(), 0.0));
            TS_ASSERT(cnonzero);
            TS_ASSERT_EQUALS(BASIC, pdfcb.mevaluator->typeint());
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->typeint());
            // first call of pdfco should use the BASIC evaluation
            TS_ASSERT_EQUALS(BASIC, pdfcb.mevaluator->mtypeused);
            TS_ASSERT_EQUALS(BASIC, pdfco.mevaluator->mtypeused);
            // test second call on the same structure
            pdfco.eval(mstru10);
            go = pdfco.getPDF();
            TS_ASSERT_EQUALS(gb, go);
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->mtypeused);
        }

};  // class TestPQEvaluator

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestPQEvaluator;

// End of file
