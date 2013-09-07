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
#include <boost/make_shared.hpp>

#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>

namespace diffpy {
namespace srreal {

using namespace std;
using diffpy::mathutils::EpsilonEqual;

//////////////////////////////////////////////////////////////////////////////
// class TestPQEvaluator
//////////////////////////////////////////////////////////////////////////////

class TestPQEvaluator : public CxxTest::TestSuite
{
    private:

        EpsilonEqual allclose;
        AtomicStructureAdapterPtr mstru10;
        AtomicStructureAdapterPtr mstru10d1;
        AtomicStructureAdapterPtr mstru10r;
        AtomicStructureAdapterPtr mstru9;

     public:

        void setUp()
        {
            const int SZ = 10;
            mstru10 = boost::make_shared<AtomicStructureAdapter>();
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
            mstru10d1 = boost::make_shared<AtomicStructureAdapter>(*mstru10);
            (*mstru10d1)[0].atomtype = "Au";
            mstru10r = boost::make_shared<AtomicStructureAdapter>();
            mstru10r->assign(mstru10->rbegin(), mstru10->rend());
            mstru9 = boost::make_shared<AtomicStructureAdapter>(*mstru10);
            mstru9->remove(9);
        }

        void test_PDF_change_atom()
        {
            PDFCalculator pdfcb;
            PDFCalculator pdfco;
            pdfcb.setEvaluatorType(BASIC);
            pdfco.setEvaluatorType(OPTIMIZED);
            pdfcb.eval(mstru10);
            pdfco.eval(mstru10);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT(!gb.empty());
            TS_ASSERT_EQUALS(gb, go);
            int cnonzero = count_if(gb.begin(), gb.end(),
                    bind1st(not_equal_to<double>(), 0.0));
            TS_ASSERT(cnonzero);
            TS_ASSERT_EQUALS(BASIC, pdfcb.getEvaluatorType());
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorType());
            // first call of pdfco should use the BASIC evaluation
            TS_ASSERT_EQUALS(BASIC, pdfcb.mevaluator->mtypeused);
            TS_ASSERT_EQUALS(BASIC, pdfco.mevaluator->mtypeused);
            // test second call on the same structure
            pdfco.eval(mstru10);
            go = pdfco.getPDF();
            TS_ASSERT_EQUALS(gb, go);
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->mtypeused);
            // test structure with one different atom
            pdfcb.eval(mstru10d1);
            pdfco.eval(mstru10d1);
            QuantityType gb1 = pdfcb.getPDF();
            QuantityType go1 = pdfco.getPDF();
            TS_ASSERT(!allclose(gb, gb1));
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->mtypeused);
            TS_ASSERT(allclose(gb1, go1));
            // change position of 1 atom
            AtomicStructureAdapterPtr stru10d1s =
                boost::make_shared<AtomicStructureAdapter>(*mstru10d1);
            (*stru10d1s)[0].cartesianposition[1] = 0.5;
            pdfcb.eval(stru10d1s);
            pdfco.eval(stru10d1s);
            QuantityType gb2 = pdfcb.getPDF();
            QuantityType go2 = pdfco.getPDF();
            TS_ASSERT(!allclose(gb1, gb2));
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->mtypeused);
            TS_ASSERT(allclose(gb2, go2));
        }


        void test_PDF_reverse_atoms()
        {
            PDFCalculator pdfcb;
            PDFCalculator pdfco;
            pdfcb.setEvaluatorType(BASIC);
            pdfco.setEvaluatorType(OPTIMIZED);
            pdfcb.eval(mstru10);
            pdfco.eval(mstru10);
            pdfco.eval(mstru10r);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->mtypeused);
            TS_ASSERT(allclose(gb, go));
        }


        void test_PDF_remove_atom()
        {
            PDFCalculator pdfcb;
            PDFCalculator pdfco;
            pdfcb.setEvaluatorType(BASIC);
            pdfco.setEvaluatorType(OPTIMIZED);
            pdfcb.eval(mstru9);
            pdfco.eval(mstru10);
            pdfco.eval(mstru9);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.mevaluator->mtypeused);
            TS_ASSERT(allclose(gb, go));
        }

};  // class TestPQEvaluator

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestPQEvaluator;

// End of file
